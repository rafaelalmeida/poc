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
	cerr << "loading LWIR image..." << endl;
	vector<Mat> matx = gdal_driver::loadLWIR("data/subset/TelopsDatasetCityLWIR_Subset.img");
	cerr << "loading training data..." << endl;
	Mat trainingOrig = gdal_driver::loadTrainingData("data/subset/TrainingMap_ENVI_RAW_format.raw");
	
	Rect roi;
	roi.x = 940;
	roi.y = 2619;
	roi.width = 549;
	roi.height = 675;

	Mat vis = visOrig(roi);
	Mat training = trainingOrig(roi);

	CvSVM *svm = classification::trainSVM(vis, training, description_vis::GCH);

	list<Mat> segments = segmentation::segmentVISGrid(vis);
	Mat map = classification::predict(vis, segments, svm);

	return 0;
}
