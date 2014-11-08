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

extern bool verbose;

int main0(int argc, char **argv) {
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
	//createTrackbar("low threshold", "win", &lowVal, 1000*10);
	//createTrackbar("high threshold", "win", &highVal, 1000*10);

	while (true) {
		Mat gray;
		cvtColor(vis, gray, CV_BGR2GRAY);

		Mat gradient(gray.size(), gray.type());
		Sobel(gray, gradient, -1, 1, 1);

		gradient = 255 - gradient;

		Mat bin(gray.size(), gray.type());
		threshold(gradient, bin, 0, 255, THRESH_BINARY + THRESH_OTSU);

		Mat dist(gray.size(), gray.type());
		distanceTransform(bin, dist, CV_DIST_L2, CV_DIST_MASK_PRECISE);

		imshow("win", dist);

		// Wait for enter
		int key = waitKey(0);
		if (key == 'q') break;
	}

	return 0;
}

int main5(int argc, char **argv) {
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
	//createTrackbar("low threshold", "win", &lowVal, 1000*10);
	//createTrackbar("high threshold", "win", &highVal, 1000*10);

	while (true) {
		Mat gray;
		cvtColor(vis, gray, CV_BGR2GRAY);

		Mat bin(gray.size(), gray.type());
		threshold(gray, bin, 0, 255, THRESH_BINARY + THRESH_OTSU);

		imshow("win", bin);

		// Wait for enter
		int key = waitKey(0);
		if (key == 'q') break;
	}

	return 0;
}

int main4(int argc, char **argv) {
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