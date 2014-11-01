#ifndef DESCRIPTION_LWIR_H
#define DESCRIPTION_LWIR_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

namespace description_lwir {
	cv::Mat spectralSignature(std::vector<cv::Mat> lwirChannels, cv::Mat mask);
}

#endif