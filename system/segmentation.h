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

#include "utils.h"

namespace segmentation {
	std::list<cv::Mat> segmentLWIRMeanShift(cv::Mat M);
	std::list<cv::Mat> segmentLWIRCanny(cv::Mat M);
}

#endif