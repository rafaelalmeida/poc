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

#include "models.h"
#include "utils.h"

using namespace cv;
using namespace std;

enum SegmentationMode {
	GRID
};

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

	Segmentation segmentLWIRMeanShift(Mat M);
	Segmentation segmentVISMeanShift(Mat M);
	Segmentation segmentLWIRCanny(Mat M);
	Segmentation segmentVISCanny(Mat M);
	Segmentation segmentVISGrid(Mat M, int tileSize=10);

	float getSegmentLabel(Mat classificationMap, Mat mask);
	list<SparseMat> getColorBlobs(Mat posterized);
	vector<SparseMat> pixelSegmentation(Mat image);
}

#endif