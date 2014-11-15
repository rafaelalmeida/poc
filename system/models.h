#ifndef MODELS_H
#define MODELS_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <stdio.h>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "config.h"
#include "utils.h"

class LWIRImage {
	public:
		// Members
		std::vector<cv::Mat> bands;
		std::vector<cv::Mat> reducedBands;
		cv::Rect roi;

		// Extreme values, for normalization purposes
		float minVal;
		float maxVal;

		float minValReduced;
		float maxValReduced;

		// Constructors
		LWIRImage();
		LWIRImage(std::vector<cv::Mat> bands);

		// Methods
		cv::Mat average();
		cv::Mat equalized();
		void minMaxAcrossBands(std::vector<cv::Mat> bands, float *minVal, 
			float *maxVal);
		cv::Mat normalizedSpectralSignature(cv::Mat mask, bool reduced=false);
		cv::Mat normalizedSpectralSignature(cv::Point point, 
			bool reduced=false);
		cv::Mat spectralSignature(cv::Mat mask, bool reduced=false);
		cv::Size size();
		int numBands();
		int numReducedBands();
		void reduceDimensionality(int keep);
		void setRoi(cv::Rect roi);
		void rescale(float scale, InterpolationMode mode);
};

class ThematicMap {
	cv::Mat _map;

	public:
		ThematicMap();
		ThematicMap(cv::Mat training);
		ThematicMap(cv::Size size);

		cv::Mat asMat();
		ThematicMap clone();
		cv::Mat coloredMap();
		cv::Size size();
		float getRegionClass(cv::Mat mask);
		std::map<unsigned char, int> getClassesCounts();
		std::list<std::pair<cv::SparseMat, int> > enumerateRegions();
		std::vector<ThematicMap> split(int k);
		void resize(cv::Size newSize);
		void combine(ThematicMap T);
};

#endif
