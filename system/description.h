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
#include "statistics.h"

extern "C" {
	#include "descriptors/image.h"
	#include "descriptors/cimage.h"

	#include "descriptors/acc.h"
	#include "descriptors/bic.h"
	#include "descriptors/cbc.h"
	#include "descriptors/gch.h"
	#include "descriptors/lch.h"
	#include "descriptors/unser.h"
}

// Import namespaces symbols for readability
using namespace cv;
using namespace std;

// Base descriptor class - subclass this and implement the describe() method
// for the VIS or the LWIR image, depending on the type of descritor. The 
// classifier will know which method to call. Implement only the method which
// receives a list of masks, the one that receives a single mask will
// automatically use the other implementation.
class Descriptor {
	public:
		virtual Mat describe(Mat image, SparseMat mask);
		virtual Mat describe(Mat image, list<SparseMat> masks);

		virtual Mat describe(LWIRImage image, SparseMat mask);
		virtual Mat describe(LWIRImage image, list<SparseMat> masks);
};

/* Classifiers left to build wrappers (TODO)
   - EMD
   - EOAC
   - Gabor
   - IRM
   - LAS
   - QCCH
   - Spytec
   - Steerable Pyramid
*/

// Descriptor wrappers
class GCHDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class ACCDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class BICDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class LCHDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class UnserDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class CBCDescriptor : public Descriptor {
	public:
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class SIGDescriptor : public Descriptor {
	public:
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

class REDUCEDSIGDescriptor : public Descriptor {
	public:
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

class MOMENTSDescriptor : public Descriptor {
	public:
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

// Util functions
Mat convertHistogramColor(Mat image, list<SparseMat> masks, 
	int dimensions, Histogram *(*descriptor)(CImage*, Image*));

#endif
