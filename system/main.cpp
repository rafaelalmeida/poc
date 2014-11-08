#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <opencv2/core/core.hpp>

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
	}

	int lowVal;
	int highVal;

	namedWindow("win");
	createTrackbar("low threshold", "win", &lowVal, 1000*10);
	createTrackbar("high threshold", "win", &highVal, 1000*10);

	while (true) {
		double low = lowVal / 10.0;
		double high = highVal / 10.0;

		cout << 
			"low threshold = " << low << endl <<
			"high threshold = " << high << endl;


		Mat gray;
		cvtColor(vis, gray, CV_BGR2GRAY);

		Mat edges(gray.size(), gray.type());
		Canny(gray, edges, low, high);

		imshow("win", edges);

		// Wait for enter
		int key = waitKey(0);
		if (key == 'q') break;
	}

	return 0;
}

int main3(int argc, char **argv) {
	Configuration conf;
	config::parse(argv, argc, conf);
	verbose = conf.verbose;

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
	}

	int spVal;
	int srVal;

	namedWindow("win");
	createTrackbar("spatial radius", "win", &spVal, 255*10);
	createTrackbar("color radius", "win", &srVal, 255*10);

	while (true) {
		double spatial = spVal / 10.0;
		double color = srVal / 10.0;

		cout << 
			"spatial radius = " << spatial << endl <<
			"color radius = " << color << endl;


		Mat dst(vis.size(), vis.type());
		pyrMeanShiftFiltering(vis, dst, spatial, color);

		//Segmentation seg(getColorBlobs(dst));

		imshow("win", dst);

		// Wait for enter
		int key = waitKey(0);
		if (key == 'q') break;
	}

	return 0;
}

int main2(int argc, char **argv) {
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

	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, vis, new GCHDescriptor()));
	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, vis, new ACCDescriptor()));
	//ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, &lwir, new SIGDescriptor()));

	log("training classifier...");
	ensemble.train();

	log("classifying image...");
	CoverMap classification = ensemble.classify();

	if (logger) {
		vector<Mat> allClassifications = ensemble.individualClassifications();

		int i = 0;
		for (auto c : allClassifications) {
			char path[16];
			sprintf(path, "classifier_%d", i);

			logger->saveImage(path, blend(vis, CoverMap(c).coloredMap()));
			i++;
		}

		Mat coloredMap = classification.coloredMap();
		logger->saveImage("consensus", blend(vis, coloredMap));
	}
	else {
		Mat coloredMap = classification.coloredMap();
		showImage(blend(vis, coloredMap));	
	}

	// Destroy logger
	if (logger) {
		delete logger;
	}

	return 0;

	/*log("calculating kappa...");
	float k = statistics::kappa(training, map);
	cout << k << endl;*/
}
