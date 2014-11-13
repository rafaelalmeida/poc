#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "common.h"
#include "models.h"
#include "utils.h"

extern "C" {
	#include "include/vl/generic.h"
	#include "include/vl/slic.h"
}

using namespace cv;
using namespace std;

// Forward declarations
class LWIRImage;

namespace segmentation {
	class Segmentation {
		list<SparseMat> _masks;

		public:
			int segmentCount();
			list<SparseMat> getSegments();
			Mat representation();
			Segmentation pixelize();
			Segmentation() {};
			Segmentation(list<SparseMat> masks);
			Size getMapSize();
	};

	Segmentation segmentVISGrid(Mat M, int tileSize=10);
	Segmentation segmentVISWatershed(Mat M, int tileSize=100);
	Segmentation segmentVIS_SLIC(Mat M, int regionSize, int minRegionSize, 
		float regularization);

	Segmentation segmentLWIRPixelated(LWIRImage& lwir, Mat vis);

	Mat makeNonMissingDataMask(Mat vis);

	float getSegmentLabel(Mat classificationMap, Mat mask);
	list<SparseMat> getColorBlobs(Mat posterized);
	vector<SparseMat> pixelSegmentation(Mat image);
}

#endif