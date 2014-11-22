#ifndef GDAL_DRIVER_H
#define GDAL_DRIVER_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>

#include <gdal_priv.h>
#include <opencv2/opencv.hpp>

#include "models.h"

namespace gdal_driver {
	LWIRImage loadLWIR(const char*);
	VISImage loadVIS(const char*);
	cv::Mat loadTrainingData(const char*);
}

#endif
