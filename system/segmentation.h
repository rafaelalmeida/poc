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
#include "description.h"
#include "models.h"
#include "utils.h"

extern "C" {
	#include "include/vl/generic.h"
	#include "include/vl/slic.h"
}

using namespace cv;
using namespace std;

// Forward declarations
class Descriptor;
class LWIRImage;
class Region;
class Segmentation;

// Class to represent a region
class Region {
	// The sparse matrix that contains a mask in which non-zero values 
	// represent the pixels of the image that belong to this region.
	SparseMat _mask;

	// A pointer to the parent segmentation, because many useful data (like
	// feature vectors) are there.
	Segmentation *_parentSegmentation;

	// The index of this region in the parent segmentation. This avoids having
	// to look up in where the feature vectors of this region are in the
	// segmentation feature matrix.
	int _parentIdx;

	public:
		Region(SparseMat mask, int parentIdx);
		Mat getDescription(string descriptorID);
		SparseMat getMask();
};

// Class to represent a Segmentation (collection of regions)
class Segmentation {
	// Regions that make up this segmentation
	list<Region> _regions;

	// Descriptions that describe the regions. Each collection item is a 
	// feature matrix: a collection of feature vectors, made by different 
	// descriptors. Items in this collection are identified by their descriptor 
	// ID, and represented by a 32FC1 OpenCV dense matrix, containing as many 
	// rows as the number of regions in this segmentation, and where each row 
	// is a feature vector created by that descriptor. The n-th line in any 
	// feature matrix is the feature vector of the n-th region of the region
	// list.
	map<string, Mat> _descriptions;

	public:
		Segmentation() {};
		Segmentation(list<SparseMat> masks);

		int regionCount();
		list<Region> getRegions();
		Mat representation();
		Segmentation pixelize();
		Size getMapSize();

		Mat getDescription(string descriptorID);
		void describe(Descriptor *descriptor);
};

namespace segmentation {
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