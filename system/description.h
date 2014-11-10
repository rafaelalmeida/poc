#ifndef DESCRIPTION_H
#define DESCRIPTION_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/ml/ml.hpp>

#include "constants.h"
#include "models.h"
#include "segmentation.h"

extern "C" {
	#include "descriptors/image.h"
	#include "descriptors/cimage.h"

	#include "descriptors/acc.h"
	#include "descriptors/bic.h"
	#include "descriptors/gch.h"
	#include "descriptors/lch.h"
	#include "descriptors/unser.h"
}

class Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, cv::Mat mask);
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks);

		virtual cv::Mat describe(LWIRImage image, cv::Mat mask);
		virtual cv::Mat describe(LWIRImage image, std::list<cv::Mat> masks);
};

// CBC
// EMD
// EOAC
// Gabor
// IRM
// LAS
// QCCH
// Spytec
// Steerable Pyramid

class GCHDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks) override;
};

class ACCDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks) override;
};

class BICDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks) override;
};

class LCHDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks) override;
};

class UnserDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(cv::Mat image, std::list<cv::Mat> masks) override;
};

class SIGDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(LWIRImage image, std::list<cv::Mat> masks) override;
};

class ENERGYDescriptor : public Descriptor {
	public:
		virtual cv::Mat describe(LWIRImage image, std::list<cv::Mat> masks) override;
};

#endif
