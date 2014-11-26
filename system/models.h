#ifndef MODELS_H
#define MODELS_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <random>
#include <stdio.h>

#include <opencv2/opencv.hpp>

#include "common.h"
#include "config.h"
#include "counter.h"
#include "statistics.h"
#include "utils.h"

// Forward declarations
class KFolder;

// Base class to represent a remote sensing image
class RSImage {
	ImageType _type;

	public:
		RSImage(ImageType type) : _type(type) {};
		ImageType getType() { return _type; };
		virtual cv::Size size() = 0;
};

// Class to represent a VIS image
class VISImage : public RSImage {
	cv::Mat _vis;

	public:
		VISImage();
		VISImage(cv::Mat vis);
		cv::Size size() override;
		cv::Mat asMat();
		void rescale(float scale, InterpolationMode mode);
		void setRoi(cv::Rect roi);
};

// Class to represent a LWIR image
class LWIRImage : public RSImage {
	private:
		// Members
		std::vector<cv::Mat> bands;
		std::vector<cv::Mat> reducedBands;
		cv::Rect roi;

		bool wasReduced = false;

		// Extreme values, for normalization purposes
		float minVal;
		float maxVal;

		float minValReduced;
		float maxValReduced;

		cv::Mat doGetSpectralSignature(cv::Mat *mask=NULL, 
			cv::Point *point=NULL, bool reduced=false);

		void doNormalizeSpectralSignature(cv::Mat signature, bool reduced);

	public:
		// Constructors
		LWIRImage();
		LWIRImage(std::vector<cv::Mat> bands);

		// Methods
		cv::Mat average();
		cv::Mat equalized();
		void minMaxAcrossBands(std::vector<cv::Mat> bands, float *minVal, 
			float *maxVal);

		// Spectral signature calculation
		cv::Mat normalizedSpectralSignature(cv::Mat mask, bool reduced=false);
		cv::Mat normalizedSpectralSignature(cv::Point point, 
			bool reduced=false);
		cv::Mat spectralSignature(cv::Mat mask, bool reduced=false);
		cv::Mat spectralSignature(cv::Point point, bool reduced=false);

		cv::Size size() override;
		int numBands();
		int numReducedBands();
		void reduceDimensionality(int keep);
		void setRoi(cv::Rect roi);
		void rescale(float scale, InterpolationMode mode);
};

// Class to represent a thematic map
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

		std::map<unsigned char, int> getClassesCounts();
		std::list<std::pair<cv::SparseMat, int> > enumerateRegions();
		cv::SparseMat getFullMask();

		// Splits the map using a KFolder object. Return form is a list
		// of k (training, validation) thematic maps.
		std::vector<std::pair<ThematicMap, ThematicMap> > split(KFolder folder);

		// Calculate a vector of class probabilities for this map, 
		// disconsidering unclassed pixels
		std::vector<float> *classProbabilities();

		void resize(cv::Size newSize);
		void combine(ThematicMap T);
};

#endif
