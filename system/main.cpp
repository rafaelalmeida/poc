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
#include "gdal_driver.h"
#include "models.h"
#include "segmentation.h"
#include "statistics.h"
#include "utils.h"

using namespace std;
using namespace cv;

using namespace segmentation;

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

	// Upscale LWIR image
	log("upscaling LWIR...");
	lwir.upscale(visFull.size());

	// Set ROI if exists
	Rect roi;
	Mat vis = visFull, training = trainingFull;
	if (conf.roiX > 0 || conf.roiY > 0 || conf.roiWidth > 0 || conf.roiHeight > 0) {
		roi = Rect(conf.roiX, conf.roiY, conf.roiWidth, conf.roiHeight);
		vis = visFull(roi);
		training = trainingFull(roi);

		log("applying ROI to LWIR...");
		lwir = lwir(roi);
	}

	log("training classifier...");
	CvSVM *svm = classification::trainSVM(vis, training, description_vis::GCH);

	log("segmenting image...");
	Segmentation segments = segmentation::segmentVISGrid(vis);

	log("classifying image...");
	Mat map = classification::predict(vis, segments, svm);

	log("calculating kappa...");
	float k = statistics::kappa(training, map);
	cout << k << endl;

	return 0;
}
