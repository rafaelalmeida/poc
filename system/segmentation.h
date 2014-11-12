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

enum SegmentationMode {
	GRID
};

namespace segmentation {
	class Segmentation {
		std::list<cv::SparseMat> _masks;

		public:
			Segmentation() {};
			Segmentation(std::list<cv::SparseMat> masks);
			std::list<cv::SparseMat> getSegments();
			cv::Mat representation();
			cv::Size getMapSize();
	};

	Segmentation segmentLWIRMeanShift(cv::Mat M);
	Segmentation segmentVISMeanShift(cv::Mat M);
	Segmentation segmentLWIRCanny(cv::Mat M);
	Segmentation segmentVISCanny(cv::Mat M);
	Segmentation segmentVISGrid(cv::Mat M, int tileSize=10);

	float getSegmentLabel(cv::Mat classificationMap, cv::Mat mask);
	std::list<cv::SparseMat> getColorBlobs(cv::Mat posterized);
}

#endif