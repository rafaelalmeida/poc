#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/ml/ml.hpp>

#include "gdal_driver.h"
#include "utils.h"
#include "statistics.h"
#include "segmentation.h"
#include "description_lwir.h"
#include "description_vis.h"

using namespace std;
using namespace cv;

int main() {
	cerr << "loading VIS image..." << endl;
	Mat vis = gdal_driver::loadVIS("data/subset/TelopsDatasetCityVisible_20cm_Subset.img");
	cerr << "loading training data..." << endl;
	Mat training = gdal_driver::loadTrainingData("data/subset/TrainingMap_ENVI_RAW_format.raw");
	
	cerr << "processing valid regions..." << endl;
	list<Mat> masks = segmentation::makeSegmentMasksFromPosterizedImage(training);

	list<float> labels;
	list<Mat> validSegments;

	cerr << "recovering region labels..." << endl;
	for (list<Mat>::iterator it = masks.begin(); it != masks.end(); ++it) {
		float label = segmentation::getSegmentLabel(training, *it);
		if (label != 0) {
			labels.push_back(label);
			validSegments.push_back(*it);
		}
	}

	Mat labelsMat(labels.size(), 1, CV_32FC1);
	float c = 0;
	for (list<float>::iterator it = labels.begin(); it != labels.end(); ++it) {
		labelsMat.at<float>(c) = *it;
		c++;
	}

	cerr << "describing regions using GCH..." << endl;
	Mat features = description_vis::GCH(vis, validSegments);

	cerr << "training SVM..." << endl;

	CvSVM SVM;

	CvSVMParams params;
    params.svm_type    = CvSVM::C_SVC;
    params.kernel_type = CvSVM::LINEAR;
    params.term_crit   = cvTermCriteria(CV_TERMCRIT_ITER, 100, 1e-6);

    SVM.train(features, labelsMat, Mat(), Mat(), params);

	return 0;
}

int main2() {
	Mat training = gdal_driver::loadTrainingData("data/subset/TrainingMap_ENVI_RAW_format.raw");
	Mat vis = gdal_driver::loadVIS("data/subset/TelopsDatasetCityVisible_20cm_Subset.img");
	vector<Mat> matx = gdal_driver::loadLWIR("data/subset/TelopsDatasetCityLWIR_Subset.img");

	Mat M = matx[0];
	Mat avg = Mat::zeros(M.rows, M.cols, M.type());

	for (int c = 0; c < matx.size(); c++) {
		avg += matx[c];
	}
	avg = avg / matx.size();

	list<Mat> segments = segmentation::segmentLWIRCanny(avg);

	Mat repr = floatImageTo8UC3Image(avg);

	// Plot contours
	for (list<Mat>::iterator i = segments.begin(); i != segments.end(); ++i)
	{
		vector<vector<Point> > contours;
  		vector<Vec4i> hierarchy;

  		findContours(*i, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
  		drawContours(repr, contours, -1, Scalar(0, 0, 255), 1);
	}

	//showImage("win", repr);

	// Calculate and show stats
	vector<float> sizes(segments.size());
	int c = 0;
	for (list<Mat>::iterator i = segments.begin(); i != segments.end(); ++i)
	{
		sizes[c] = countNonZero(*i);
		c++;
	}

	return 0;
}
