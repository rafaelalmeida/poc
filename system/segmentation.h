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
		std::list<cv::Mat> _masks;

		public:
			Segmentation(std::list<cv::Mat> masks);
			std::list<cv::Mat> getSegments();
			cv::Mat representation();
	};

	Segmentation segmentLWIRMeanShift(cv::Mat M);
	Segmentation segmentLWIRCanny(cv::Mat M);
	Segmentation segmentVISGrid(cv::Mat M);

	float getSegmentLabel(cv::Mat classificationMap, cv::Mat mask);
	std::list<cv::Mat> getColorBlobs(cv::Mat posterized);
	cv::Mat representSegmentation(std::list<cv::Mat> masks);
}

#endif