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

// Forward declarations
class VISImage;
class LWIRImage;

// Descriptor types
enum DescriptorType {
	VIS,
	LWIR
};

// Base descriptor class - subclass this and implement the describe() method
// for the VIS or the LWIR image, depending on the type of descritor. The 
// classifier will know which method to call. Implement only the method which
// receives a list of masks, the one that receives a single mask will
// automatically use the other implementation.
class Descriptor {
	// Members
	string _id;
	DescriptorType _type;

	public:
		Descriptor(const char *id, DescriptorType type);

		virtual Mat describe(Mat image, SparseMat mask);
		virtual Mat describe(Mat image, list<SparseMat> masks);

		virtual Mat describe(LWIRImage image, SparseMat mask);
		virtual Mat describe(LWIRImage image, list<SparseMat> masks);

		DescriptorType getType() { return _type; };
		std::string getID();
};

// Wrapper classes for each classifier type, for easier polymorphism
class VISDescriptor : public Descriptor {
	public:
		VISDescriptor(const char *id) : Descriptor(id, DescriptorType::VIS) {};
};

class LWIRDescriptor : public Descriptor {
	public:
		LWIRDescriptor(const char *id) : Descriptor(id, DescriptorType::LWIR) {};
};

// Descriptor wrappers
class GCHDescriptor : public VISDescriptor {
	public:
		GCHDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class ACCDescriptor : public VISDescriptor {
	public:
		ACCDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class BICDescriptor : public VISDescriptor {
	public:
		BICDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class LCHDescriptor : public VISDescriptor {
	public:
		LCHDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class UnserDescriptor : public VISDescriptor {
	public:
		UnserDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class CBCDescriptor : public VISDescriptor {
	public:
		CBCDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(Mat image, list<SparseMat> masks) override;
};

class SIGDescriptor : public LWIRDescriptor {
	public:
		SIGDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

class REDUCEDSIGDescriptor : public LWIRDescriptor {
	public:
		REDUCEDSIGDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

class MOMENTSDescriptor : public LWIRDescriptor {
	public:
		MOMENTSDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, list<SparseMat> masks) override;
};

// Util functions
Mat convertHistogramColor(Mat image, list<SparseMat> masks, 
	int dimensions, Histogram *(*descriptor)(CImage*, Image*));

#endif
