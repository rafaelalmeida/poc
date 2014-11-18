#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <opencv2/core/core.hpp>
#include <opencv2/gpu/gpu.hpp>

#include "classification.h"
#include "common.h"
#include "config.h"
#include "description_lwir.h"
#include "description_vis.h"
#include "ensemble.h"
#include "gdal_driver.h"
#include "logging.h"
#include "models.h"
#include "segmentation.h"
#include "statistics.h"
#include "utils.h"

using namespace std;
using namespace cv;

using namespace segmentation;
using namespace classification;

bool verbose;
Logger *logger = NULL;

// Signatures
void rescale(Mat& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	InterpolationMode interpolationMode);

void setupClassifiers(Ensemble& e, Mat vis, LWIRImage& lwir);

void printClassHistogram(ThematicMap& M);

// Main function
int main(int argc, char **argv) {
	Configuration conf;
	config::parse(argv, argc, conf);
	verbose = conf.verbose;

	// Setup logger
	logger = new Logger(conf.logPath);
	logger->saveArguments(argc, argv);

	// Load images
	log("loading VIS image...");
	Mat vis = gdal_driver::loadVIS(conf.pathVIS);
	log("loading LWIR image...");
	LWIRImage lwir = gdal_driver::loadLWIR(conf.pathLWIR);
	log("loading training data...");
	Mat training = gdal_driver::loadTrainingData(conf.pathTraining);

	// Ignore unclassed pixels to speed up classification, if enabled
	if (conf.ignoreUnclassedPixels) {
		Mat trainingMask = training > 0;
		Mat trainingMask3C;
		cvtColor(trainingMask, trainingMask3C, CV_GRAY2BGR);

		vis = vis & trainingMask3C;
	}

	// Rescale images if necessary
	if (conf.scaleVIS != 1.0 || conf.scaleLWIR != 1.0) {
		log("rescaling images...");
		rescale(vis, lwir, conf.scaleVIS, conf.scaleLWIR, 
			conf.interpolationMode);
	}

	log("creating training thematic maps...");

	// Create training map objects
	ThematicMap trainingMapVIS(training);
	ThematicMap trainingMapLWIR(training);

	// Scale the training map
	trainingMapVIS.resize(vis.size());
	trainingMapLWIR.resize(lwir.size());

	// Reduce dimensionality of LWIR image
	log("reducing LWIR dimensionality...");
	lwir.reduceDimensionality(LWIR_BANDS_TO_KEEP_ON_PCA);

	// Create k-fold splits
	log("creating k-fold splits...");
	vector<ThematicMap> splits = trainingMapVIS.split(K_FOLDS);
	int c = 0;
	for (auto s : splits) {
		string n = "training-split-";
		n += (c + '0');

		logger->saveImage(n.c_str(), s.coloredMap());
		c++;
	}

	// Save training map for debugging
	logger->saveImage("training", blend(vis, trainingMapVIS.coloredMap()));

	// Show classes counts in the whole training map
	cerr << "Amount of pixels per class in the whole map: " << endl;
	printClassHistogram(trainingMapVIS);

	// Set ROI if there is one
	Rect roi;
	if (conf.roiWidth > 0 && conf.roiHeight > 0) {
		roi = Rect(conf.roiX, conf.roiY, conf.roiWidth, conf.roiHeight);
		vis = vis(roi);
		training = training(roi);

		log("applying ROI to LWIR...");
		lwir.setRoi(roi);
	}

	log("segmenting image...");
	Segmentation segmentationVIS;
	if (conf.segmentationMode == GRID) {
		segmentationVIS = segmentation::segmentVISGrid(vis, conf.gridTileSize);
	}
	else if (conf.segmentationMode == SLIC) {
		// Try to guess rescaled parameters for SLIC, if desired
		if (conf.slicAutoScaleParameters) {
			float s = conf.scaleVIS;

			conf.slicRegionSize *= s;
			conf.slicMinRegionSize *= s;
		}

		segmentationVIS = segmentation::segmentVIS_SLIC(vis, 
			conf.slicRegionSize, conf.slicMinRegionSize, 
			conf.slicRegularization);
	}
	else {
		assert(false && "Unsupported segmentation mode");
	}

	Segmentation segmentationLWIR = segmentation::segmentLWIRPixelated(
		lwir, vis);

	// Save segmentation representation
	logger->saveImage("segmentation", segmentationVIS.representation());

	// Open result file
	ofstream *results = logger->makeFile("results.txt");
	*results << "MAP AGREEMENT KAPPA" << endl;

	// Run k-fold cross validation
	float bestKappa = FLT_MIN;
	int bestFold = -1; // Sentinel
	vector<pair<string, Mat> > bestClassifications;
	vector<ThematicMap> foldValidationMaps;
	ThematicMap bestConsensus;

	cerr << "running k-fold cross-validation (k = " << K_FOLDS << ")" << endl;
	for (int fold = 0; fold < K_FOLDS; fold++) {
		// Report progress
		cerr << "running fold " << (fold+1) << " of " << K_FOLDS << endl;

		// Build training map (made by all folds except this one)
		ThematicMap T(trainingMapVIS.size());
		for (int i = 0; i < K_FOLDS; i++) {
			if (i != fold) {
				T.combine(splits[i]);
			}
		}

		// Resize this fold's training map to the correct sizes
		ThematicMap T_VIS = T.clone();
		ThematicMap T_LWIR = T.clone();
		T_VIS.resize(vis.size());
		T_LWIR.resize(lwir.size());

		// Show classes counts in this fold
		cerr << "Amount of pixels per class in fold: " << endl;
		printClassHistogram(T_VIS);

		// Build the ensemble
		Ensemble E(MAJORITY_VOTING, segmentationVIS, segmentationLWIR,
		T_VIS, T_LWIR);
		E.setParallel(conf.parallel);
		setupClassifiers(E, vis, lwir);

		// Train the ensemble
		E.train();

		// Run the classification
		ThematicMap C = E.classify();

		// Use the current fold for validation
		Mat V = splits[fold].asMat();
		foldValidationMaps.push_back(V);

		// Calculate performance metrics
		Mat X = C.asMat();
		float a = statistics::agreement(V, X);
		float k = statistics::kappa(V, X);

		// Show statistics
		cerr << "agreement = " << a << ", kappa = " << k << endl;

		// Log results
		string imageName("fold_");
		imageName += (fold+'0');
		logger->saveImage(imageName.c_str(), blend(vis, C.coloredMap()));

		// See if there is improvement
		if (k > bestKappa) {
			bestKappa = k;
			bestFold = fold;
			bestClassifications = E.individualClassifications();
			bestConsensus = C;
		}
	}

	cerr << "best kappa = " << bestKappa << " on fold " << (bestFold) << 
		endl;

	// Process best individual classifications of best performing ensemble
	log("calculating individual statistics...");
	Mat G = foldValidationMaps[bestFold].asMat();
	for (auto c : bestClassifications) {
		// Save results
		Mat X = c.second;

		float agreement = statistics::agreement(G, X);
		float kappa = statistics::kappa(G, X);

		*results << c.first << " " << agreement << " " << kappa 
			<< endl;

		// Save images
		string imageName("classifier-");
		imageName += c.first;
		logger->saveImage(imageName.c_str(), 
			blend(vis, ThematicMap(c.second).coloredMap()));
	}

	// Calculate consensus kappa
	log("calculating consensus statistics...");

	Mat C = bestConsensus.asMat();
	float agreement = statistics::agreement(G, C);
	float kappa = statistics::kappa(G, C);

	cerr << agreement << " " << kappa << endl;
	*results << "MAJORITY " << agreement << " " << kappa << endl;

	// Close result file
	results->close();
	delete results;
	
	// Destroy logger
	delete logger;

	return 0;
}

void rescale(Mat& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	InterpolationMode interpolationMode) {

	// Scale LWIR
	lwir.rescale(scaleLWIR, interpolationMode);

	// Scale VIS
	Mat visR;
	resize(vis, visR, Size(), scaleVIS, scaleVIS, 
		translateInterpolationMode(interpolationMode));
	vis = visR;
}

void setupClassifiers(Ensemble& ensemble, Mat vis, LWIRImage& lwir) {
	// Register ensemble classifiers - SVM
	/*ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new GCHDescriptor("GCH")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new ACCDescriptor("ACC")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new BICDescriptor("BIC")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new LCHDescriptor("LCH")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new SIGDescriptor("SIG")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new REDUCEDSIGDescriptor("RSIG")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new MOMENTSDescriptor("MMT")));*/

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new GCHDescriptor("GCH")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::RTREES, 
		vis, new GCHDescriptor("GCH")));
}

void printClassHistogram(ThematicMap& M) {
	auto counts = M.getClassesCounts();
	for (int i = 0; i < counts.size(); i++) {
		cerr << counts[(unsigned char) i] << " ";
	}
	cerr << endl;
}
