#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/ml/ml.hpp>

#include "classification.h"
#include "description_lwir.h"
#include "description_vis.h"
#include "gdal_driver.h"
#include "segmentation.h"
#include "statistics.h"
#include "utils.h"

using namespace std;
using namespace cv;

int main() {
	cerr << "loading VIS image..." << endl;
	Mat visOrig = gdal_driver::loadVIS("data/subset/TelopsDatasetCityVisible_20cm_Subset.img");
	cerr << "loading training data..." << endl;
	Mat trainingOrig = gdal_driver::loadTrainingData("data/subset/TrainingMap_ENVI_RAW_format.raw");
	
	Rect roi;
	roi.x = 940;
	roi.y = 2619;
	roi.width = 549;
	roi.height = 675;

	Mat vis = visOrig(roi);
	Mat training = trainingOrig(roi);

	list<Mat> segments = segmentation::segmentVISGrid(vis);
	Mat repr = segmentation::representSegmentation(segments);
	showImage("win", repr);

	return 0;

	classification::trainSVM(vis, training, description_vis::GCH);

    cerr << "drawing land cover map..." << endl;
    Mat coverMap = makeLandCoverMap(training);
    showImage("win", coverMap);
    showImage("win", vis);

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
