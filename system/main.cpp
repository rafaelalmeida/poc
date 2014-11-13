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
	if (conf.logEnabled) {
		logger = new Logger(conf.logPath);
	}

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
	if (logger) {
		logger->saveImage("training", 
			blend(vis, trainingMapVIS.coloredMap()));
	}

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
	Segmentation segmentation;
	if (conf.segmentationMode == GRID) {
		segmentation = segmentation::segmentVISGrid(vis, conf.gridTileSize);
	}
	else {
		assert(false && "Unsupported segmentation mode");
	}

	// Save segmentation representation
	if (logger) {
		logger->saveImage("segmentation", segmentation.representation());
	}

	// Setup classifier ensemble
	Ensemble ensemble(MAJORITY_VOTING, segmentation, trainingMapVIS, 
		trainingMapLWIR);
	ensemble.setLogger(logger);
	ensemble.setParallel(conf.parallel);

	// Register ensemble classifiers
	ensemble.addClassifier(new Classifier("SVM-GCH", ClassifierEngine::SVM, 
		vis, new GCHDescriptor()));

	/*ensemble.addClassifier(new Classifier("SVM-ACC", ClassifierEngine::SVM, 
		vis, new ACCDescriptor()));

	ensemble.addClassifier(new Classifier("SVM-BIC", ClassifierEngine::SVM, 
		vis, new BICDescriptor()));

	ensemble.addClassifier(new Classifier("SVM-LCH", ClassifierEngine::SVM, 
		vis, new LCHDescriptor()));

	ensemble.addClassifier(new Classifier("SVM-Unser", ClassifierEngine::SVM, 
		vis, new UnserDescriptor()));*/

	/*ensemble.addClassifier(new Classifier("SVM-SIG", ClassifierEngine::SVM, 
		&lwir, new SIGDescriptor()));*/

	/*ensemble.addClassifier(new Classifier("SVM-SIG", ClassifierEngine::SVM, 
		&lwir, new REDUCEDSIGDescriptor()));*/

	ensemble.addClassifier(new Classifier("SVM-MOMENTS", ClassifierEngine::SVM, 
		&lwir, new MOMENTSDescriptor()));

	// Train ensemble
	log("training ensemble...");
	ensemble.train();

	// Classify ensemble
	log("classifying image...");
	ThematicMap classification = ensemble.classify();
	log("classifying image... done     ");

	// Log classification results
	if (logger) {
		vector<pair<string, Mat> > allClassifications = 
			ensemble.individualClassifications();

		for (auto c : allClassifications) {
			char path[16];
			sprintf(path, "classifier_%s", c.first.c_str());

			logger->saveImage(path, blend(vis, 
				ThematicMap(c.second).coloredMap()));
		}

		Mat coloredMap = classification.coloredMap();
		logger->saveImage("consensus", blend(vis, coloredMap));
	}
	else {
		Mat coloredMap = classification.coloredMap();
		showImage(blend(vis, coloredMap));	
	}

	// Calculate statistics
	log("calculating kappa...");
	float k = statistics::kappa(training, classification.asMat());
	cout << k << endl;
	
	// Destroy logger
	if (logger) {
		delete logger;
	}

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
