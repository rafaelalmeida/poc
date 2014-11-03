#ifndef DESCRIPTION_VIS_H
#define DESCRIPTION_VIS_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "utils.h"

extern "C" {
	#include "descriptors/image.h"
	#include "descriptors/cimage.h"

	#include "descriptors/gch.h"
}

namespace description_vis {
	// TODO DESCRIPTOR LIST
	// ACC
	// BIC
	// CBC
	// EMD
	// EOAC
	// Gabor
	// IRM
	// LCH
	// LAS
	// QCCH
	// Spytec
	// Steerable Pyramid
	// Unser

	cv::Mat GCH(const cv::Mat image, const std::list<cv::Mat> masks);
}

#endif