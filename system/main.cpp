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

void rescale(Mat& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	ResamplingMethod resamplingMethod);

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

	// Rescale images if necessary
	if (conf.scaleVIS != 1.0 || conf.scaleLWIR != 1.0) {
		log("rescaling images...");
		rescale(vis, lwir, conf.scaleVIS, conf.scaleLWIR, 
			conf.resamplingMethod);
	}
	
	// Scale training map to the correct sizes (one for VIS and one for LWIR)
	Mat trainingVIS, trainingLWIR;
	resize(training, trainingVIS, vis.size(), 0, 0, 
		TRAINING_INTERPOLATION_MODE);
	resize(training, trainingLWIR, lwir.size(), 0, 0, 
		TRAINING_INTERPOLATION_MODE);

	// Reduce dimensionality of LWIR image
	log("reducing LWIR dimensionality...");
	lwir.reduceDimensionality(LWIR_BANDS_TO_KEEP_ON_PCA);

	// Create training thematic maps
	log("creating training thematic maps...");
	ThematicMap trainingMapVIS(trainingVIS);
	ThematicMap trainingMapLWIR(trainingLWIR);

	// DEBUG: test k-fold
	vector<ThematicMap> splits = trainingMapVIS.split(5);
	int c = 0;
	for (auto s : splits) {
		string n = "training-split-";
		n += (c + '0');

		logger->saveImage(n.c_str(), s.coloredMap());
		c++;
	}

	// Save training map for debugging
	logger->saveImage("training", blend(vis, trainingMapVIS.coloredMap()));

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

	// Setup classifier ensemble
	Ensemble ensemble(MAJORITY_VOTING, segmentationVIS, segmentationLWIR,
		trainingMapVIS, trainingMapLWIR);
	ensemble.setLogger(logger);
	ensemble.setParallel(conf.parallel);

	// Register ensemble classifiers
	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new GCHDescriptor("GCH")));

	/*ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new ACCDescriptor("ACC")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new BICDescriptor("BIC")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new LCHDescriptor("LCH")));

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		vis, new UnserDescriptor("UNS")));*/

	/*ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new SIGDescriptor("SIG")));*/

	/*ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new REDUCEDSIGDescriptor("RSIG")));*/

	/*ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new MOMENTSDescriptor("MMT")));*/

	// Train ensemble
	log("starting ensemble training...");
	ensemble.train();

	// Classify ensemble
	log("classifying image...");
	ThematicMap classification = ensemble.classify();
	log("classifying image... done     ");

	// Grab individual classifications
	auto allClassifications = ensemble.individualClassifications();

	// Open result file
	ofstream results = logger->makeFile("results.txt");
	results << "MAP AGREEMENT KAPPA" << endl;

	// Log classification maps
	for (auto c : allClassifications) {
		char path[16];
		sprintf(path, "classifier_%s", c.first.c_str());

		logger->saveImage(path, blend(vis, 
			ThematicMap(c.second).coloredMap()));
	}

	Mat coloredMap = classification.coloredMap();
	logger->saveImage("consensus", blend(vis, coloredMap));

	// Calculate kappa for all classifications
	log("calculating individual statistics...");
	Mat G = trainingMapVIS.asMat();
	for (auto c : allClassifications) {
		Mat X = c.second;

		float agreement = statistics::agreement(G, X);
		float kappa = statistics::kappa(G, X);

		results << c.first << " " << agreement << " " << kappa 
			<< endl;
	}

	// Calculate consensus kappa
	log("calculating consensus statistics...");

	Mat C = classification.asMat();
	float agreement = statistics::agreement(G, C);
	float kappa = statistics::kappa(G, C);

	cerr << agreement << " " << kappa << endl;
	results << "MAJORITY " << agreement << " " << kappa << endl;

	// Close result file
	results.close();
	
	// Destroy logger
	delete logger;

	return 0;
}

void rescale(Mat& vis, LWIRImage& lwir, float scaleVIS, float scaleLWIR, 
	ResamplingMethod resamplingMethod) {

	// Scale LWIR
	lwir.rescale(scaleLWIR, resamplingMethod);

	// Scale VIS
	Mat visR;
	resize(vis, visR, Size(), scaleVIS, scaleVIS, 
		translateInterpolationMode(resamplingMethod));
	vis = visR;
}
