#ifndef MODELS_H
#define MODELS_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "utils.h"

class LWIRImage {
	public:
		// Members
		std::vector<cv::Mat> bands;
		cv::Rect roi;

		// Constructors
		LWIRImage();
		LWIRImage(std::vector<cv::Mat> bands);

		// Methods
		cv::Mat average();
		cv::Mat equalized();
		int numBands();
		cv::Mat spectralSignature(cv::Mat mask);
		void upscale(cv::Size size);
		void setRoi(cv::Rect roi);
};

class CoverMap {
	cv::Mat _map;

	public:
		CoverMap(cv::Mat training);

		cv::Mat asMat();
		float getRegionClass(cv::Mat mask);
		cv::Mat coloredMap();
};

#endif
