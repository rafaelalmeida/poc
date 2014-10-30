#ifndef GDAL_DRIVER_H
#define GDAL_DRIVER_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace gdal_driver {
	std::vector<cv::Mat> loadLWIR(const char*);
	cv::Mat loadVIS(const char*);
	cv::Mat loadTrainingData(const char*);
}

#endif