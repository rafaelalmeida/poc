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

void rescale(Mat& vis, LWIRImage& lwir, Mat trainingOriginal, Mat& trainingVIS, 
	Mat& trainingLWIR, float scaleVIS, float scaleLWIR, 
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
	Mat visFull = gdal_driver::loadVIS(conf.pathVIS);
	log("loading LWIR image...");
	LWIRImage lwir = gdal_driver::loadLWIR(conf.pathLWIR);
	log("loading training data...");
	Mat trainingFull = gdal_driver::loadTrainingData(conf.pathTraining);

	// Matrixes with default settings
	Mat vis = visFull, training = trainingFull, trainingVIS, trainingLWIR;

	// Rescale images if necessary
	if (conf.scaleVIS != 1.0 || conf.scaleLWIR != 1.0) {
		log("rescaling images...");
		rescale(vis, lwir, training, trainingVIS, trainingLWIR, 
			conf.scaleVIS, conf.scaleLWIR, conf.resamplingMethod);
	}

	// Reduce dimensionality of LWIR image
	lwir.reduceDimensionality(LWIR_BANDS_TO_KEEP_ON_PCA);

	// Create training thematic maps
	ThematicMap trainingMapVIS(trainingVIS);
	ThematicMap trainingMapLWIR(trainingLWIR);

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
		segmentationVIS = segmentation::segmentVIS_SLIC(vis);
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

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, 
		&lwir, new MOMENTSDescriptor("MMT")));

	// Train ensemble
	log("training ensemble...");
	ensemble.train();

	// Classify ensemble
	log("classifying image...");
	ThematicMap classification = ensemble.classify();
	log("classifying image... done     ");

	// Grab individual classifications
	auto allClassifications = ensemble.individualClassifications();

	// Open result file
	ofstream results = logger->makeFile("results.txt");

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
	log("calculating individual kappas...");
	for (auto c : allClassifications) {
		Mat map = c.second;
		float kappa = statistics::kappa(trainingMapVIS.asMat(), map);
		results << "kappa " << c.first << " " << kappa << endl;
	}

	// Calculate consensus kappa
	log("calculating consensus kappa...");
	float k = statistics::kappa(trainingMapVIS.asMat(), 
		classification.asMat());
	cerr << k << endl;
	results << "kappa consensus " << k << endl;

	// Close result file
	results.close();
	
	// Destroy logger
	delete logger;

	return 0;
}

void rescale(Mat& vis, LWIRImage& lwir, Mat trainingOriginal, Mat& trainingVIS, 
	Mat& trainingLWIR, float scaleVIS, float scaleLWIR, 
	ResamplingMethod resamplingMethod) {

	// Scale LWIR
	if (scaleLWIR != 1.0) {
		// Rescale the LWIR image
		lwir.rescale(scaleLWIR, resamplingMethod);

		// Rescale LWIR training map
		resize(trainingOriginal, trainingLWIR, lwir.size(), 0, 0,
			TRAINING_INTERPOLATION_MODE);
	}
	else {
		trainingLWIR = trainingOriginal;
	}

	// Scale VIS
	if (scaleVIS != 1.0) {
		// Rescale the VIS image
		Mat visR;
		resize(vis, visR, Size(), scaleVIS, scaleVIS, 
			translateInterpolationMode(resamplingMethod));
		vis = visR;

		// Rescale VIS training map
		resize(trainingOriginal, trainingVIS, Size(), scaleVIS, scaleVIS,
			TRAINING_INTERPOLATION_MODE);
	}
	else {
		trainingVIS = trainingOriginal;
	}
}
