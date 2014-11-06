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
		lwir.setRoi(roi);
	}

	log("segmenting image...");
	Segmentation segmentation = segmentation::segmentVISGrid(vis);

	// Setup classifier ensemble
	CoverMap tMap(training);
	Ensemble ensemble(MAJORITY_VOTING, segmentation, tMap);
	ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, vis, new GCHDescriptor()));
	//ensemble.addClassifier(new Classifier(ClassifierEngine::SVM, &lwir, new SIGDescriptor()));

	log("training classifier...");
	ensemble.train();

	log("classifying image...");
	CoverMap classification = ensemble.classify();

	Mat coloredMap = classification.coloredMap();
	showImage(blend(vis, coloredMap));

	return 0;

	/*log("calculating kappa...");
	float k = statistics::kappa(training, map);
	cout << k << endl;*/
}
