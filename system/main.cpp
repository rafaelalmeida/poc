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
void rescale(VISImage& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	InterpolationMode interpolationMode);

void setupClassifiers(Ensemble& e, VISImage& vis, LWIRImage& lwir, 
	vector<Descriptor*> descriptors);

void printClassHistogram(ThematicMap& M);

// Main function
int main(int argc, char **argv) {
	// Start timer
	Stopwatch swatchMain;
	swatchMain.start();

	// Parse configuration
	Configuration conf;
	config::parse(argv, argc, conf);
	verbose = conf.verbose;

	// Setup logger
	logger = new Logger(conf.logPath);
	logger->saveArguments(argc, argv);

	// Load images
	log("loading VIS image...");
	VISImage vis = gdal_driver::loadVIS(conf.pathVIS);
	log("loading LWIR image...");
	LWIRImage lwir = gdal_driver::loadLWIR(conf.pathLWIR);
	log("loading training data...");
	Mat training = gdal_driver::loadTrainingData(conf.pathTraining);

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
	logger->saveImage("training", 
	                  blend(vis.asMat(), trainingMapVIS.coloredMap()));

	// Show classes counts in the whole training map
	cerr << "Amount of pixels per class in the whole map: " << endl;
	printClassHistogram(trainingMapVIS);

	// Set ROI if there is one
	Rect roi;
	if (conf.roiWidth > 0 && conf.roiHeight > 0) {
		roi = Rect(conf.roiX, conf.roiY, conf.roiWidth, conf.roiHeight);
		vis.setRoi(roi);
		training = training(roi);

		log("applying ROIs to LWIR...");
		lwir.setRoi(roi);
	}

	// Create segmentation - VIS
	log("creating segmentation...");
	Segmentation segmentationVIS;
	if (conf.segmentationMode == GRID) {
		segmentationVIS = segmentation::segmentVISGrid(vis.asMat(), 
			conf.gridTileSize);
	}
	else if (conf.segmentationMode == SLIC) {
		// Try to guess rescaled parameters for SLIC, if desired
		if (conf.slicAutoScaleParameters) {
			float s = conf.scaleVIS;

			conf.slicRegionSize *= s;
			conf.slicMinRegionSize *= s;
		}

		segmentationVIS = segmentation::segmentVIS_SLIC(vis.asMat(), 
			conf.slicRegionSize, conf.slicMinRegionSize, 
			conf.slicRegularization);
	}
	else {
		assert(false && "Unsupported segmentation mode");
	}

	// Create segmentation - LWIR
	Segmentation segmentationLWIR = segmentation::segmentLWIRPixelated(
		lwir, vis);

	// Save segmentation representation
	logger->saveImage("segmentation", segmentationVIS.representation());

	// Assign images to segmentations
	segmentationVIS.setImage(&vis);
	segmentationLWIR.setImage(&lwir);

	// Create descriptors
	vector<Descriptor*> descriptors;
	descriptors.push_back(new GCHDescriptor("GCH"));
	descriptors.push_back(new BICDescriptor("BIC"));
	descriptors.push_back(new LCHDescriptor("LCH"));
	descriptors.push_back(new SIGDescriptor("SIG"));
	descriptors.push_back(new REDUCEDSIGDescriptor("RSIG"));
	descriptors.push_back(new MOMENTSDescriptor("MMT"));

	// Describe segmentation objects
	for (auto d : descriptors) {
		if (d->getType() == VIS) {
			segmentationVIS.describe(d);
		}
		else {
			segmentationLWIR.describe(d);
		}
	}

	// Open result file
	ofstream *results = logger->makeFile("results.txt");
	*results << "MAP AGREEMENT KAPPA" << endl;

	// Initialize time loggers
	double descriptionTime = 0;
	double trainingTime = 0;
	double classificationTime = 0;

	// Run k-fold cross validation
	float bestKappa = -numeric_limits<float>::max();
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
		setupClassifiers(E, vis, lwir, descriptors);

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
		logger->saveImage(imageName.c_str(), 
		                  blend(vis.asMat(), C.coloredMap()));

		// Register time taken
		descriptionTime += E.getTotalDescriptionTime();
		trainingTime += E.getTotalTrainingTime();
		classificationTime += E.getTotalClassificationTime();

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
			blend(vis.asMat(), ThematicMap(c.second).coloredMap()));
	}

	// Calculate consensus kappa
	log("calculating consensus statistics...");

	Mat C = bestConsensus.asMat();
	float agreement = statistics::agreement(G, C);
	float kappa = statistics::kappa(G, C);

	cerr << agreement << " " << kappa << endl;
	*results << "MAJORITY " << agreement << " " << kappa << endl;

	// Stop timer and log total execution time
	swatchMain.stop();
	double totalTime = swatchMain.read();
	*results << endl << "Total wall time: " << totalTime << endl;

	// Log detailed execution times
	*results << endl << 
		"Description time: " << descriptionTime << endl <<
		"Training time: " << trainingTime << endl <<
		"Classification time: " << classificationTime << endl;

	// Close result file
	results->close();
	delete results;
	
	// Destroy logger
	delete logger;

	// Destroy descriptors
	for (auto d : descriptors) {
		delete d;
	}

	return 0;
}

void rescale(VISImage& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	InterpolationMode interpolationMode) {

	// Scale LWIR
	lwir.rescale(scaleLWIR, interpolationMode);

	// Scale VIS
	vis.rescale(scaleVIS, interpolationMode);
}

void setupClassifiers(Ensemble& ensemble, VISImage& vis, LWIRImage& lwir,
	vector<Descriptor*> descriptors) {
	// List classifier engines
	vector<ClassifierEngine> engines = {
		ClassifierEngine::SVM, 
		ClassifierEngine::KNN, 
		ClassifierEngine::DTREE, 
		ClassifierEngine::MLP
	};

	// Create the classifier <-> descriptor pairs and add them to the ensemble
	for (auto engine : engines) {
		for (auto descriptor : descriptors) {
			Classifier *classifier;
			if (descriptor->getType() == ImageType::VIS) {
				classifier = new Classifier(engine, &vis, descriptor);
			}
			else {
				classifier = new Classifier(engine, &lwir, descriptor);
			}

			cerr << "Adding classifier " << classifier->getID() << endl;
			ensemble.addClassifier(classifier);
		}
	}
}

void printClassHistogram(ThematicMap& M) {
	auto counts = M.getClassesCounts();
	for (int i = 0; i < counts.size(); i++) {
		cerr << counts[(unsigned char) i] << " ";
	}
	cerr << endl;
}
