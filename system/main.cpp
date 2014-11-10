#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <opencv2/core/core.hpp>
#include <opencv2/gpu/gpu.hpp>

#include "classification.h"
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

int main(int argc, char **argv) {
	Configuration conf;
	config::parse(argv, argc, conf);
	verbose = conf.verbose;

	// Setup logger
	Logger *logger = NULL;
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
	Mat vis = visFull, training = trainingFull;

	// Set sampling mode
	if (conf.samplingMode == UPSAMPLE_LWIR) {
		log("upscaling LWIR...");
		lwir.upscale(visFull.size());
	}
	else {
		log("downsampling VIS...");
		Mat resizedVis, resizedTraining;
		resize(visFull, resizedVis, lwir.size());
		resize(trainingFull, resizedTraining, lwir.size());

		vis = resizedVis;
		training = resizedTraining;

		if (logger) {
			logger->saveImage("vis", vis);
		}
	}

	// Save training map for debugging
	if (logger) {
		logger->saveImage("training", blend(vis, CoverMap(training).coloredMap()));
	}

	// Set ROI if exists
	Rect roi;
	if (conf.roiWidth > 0 && conf.roiHeight > 0) {
		roi = Rect(conf.roiX, conf.roiY, conf.roiWidth, conf.roiHeight);
		vis = vis(roi);
		training = training(roi);

		log("applying ROI to LWIR...");
		lwir.setRoi(roi);
	}

	// Log cropped images
	if (logger) {
		logger->saveImage("training-ROI", blend(vis, CoverMap(training).coloredMap()));
		logger->saveImage("vis-ROI", vis);
	}

	log("segmenting image...");
	Segmentation segmentation;
	if (conf.segmentationMode == GRID) {
		segmentation = segmentation::segmentVISGrid(vis);
		if (logger) {
			logger->saveImage("segmentation", segmentation.representation());
		}
	}
	else {
		assert(false && "Unsupported segmentation mode");
	}

	// Setup classifier ensemble
	CoverMap tMap(training);
	Ensemble ensemble(MAJORITY_VOTING, segmentation, tMap);
	ensemble.setLogger(logger);
	ensemble.setParallel(conf.parallel);

	ensemble.addClassifier(new Classifier("SVM-GCH", ClassifierEngine::SVM, vis, new GCHDescriptor()));
	ensemble.addClassifier(new Classifier("SVM-ACC", ClassifierEngine::SVM, vis, new ACCDescriptor()));
	//ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, &lwir, new SIGDescriptor()));

	log("training ensemble...");
	ensemble.train();

	log("classifying image...");
	CoverMap classification = ensemble.classify();
	log("classifying image... done     ");

	if (logger) {
		vector<pair<string, Mat> > allClassifications = ensemble.individualClassifications();

		for (auto c : allClassifications) {
			char path[16];
			sprintf(path, "classifier_%s", c.first.c_str());

			logger->saveImage(path, blend(vis, CoverMap(c.second).coloredMap()));
		}

		Mat coloredMap = classification.coloredMap();
		logger->saveImage("consensus", blend(vis, coloredMap));
	}
	else {
		Mat coloredMap = classification.coloredMap();
		showImage(blend(vis, coloredMap));	
	}

	log("calculating kappa...");
	float k = statistics::kappa(training, classification.asMat());
	cout << k << endl;
	
	// Destroy logger
	if (logger) {
		delete logger;
	}

	return 0;
}
