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
class Segmentation;
class Region;

// Base descriptor class - subclass this and implement the describe() method
// for the VIS or the LWIR image, depending on the type of descritor. The 
// classifier will know which method to call.
class Descriptor {
	// Members
	string _id;
	ImageType _type;

	public:
		Descriptor(const char *id, ImageType type);

		virtual Mat describe(VISImage vis, Segmentation *S);
		virtual Mat describe(LWIRImage image, Segmentation *S);

		ImageType getType() { return _type; };
		std::string getID();
};

// Wrapper classes for each classifier type, for easier polymorphism
class VISDescriptor : public Descriptor {
	public:
		VISDescriptor(const char *id) : Descriptor(id, ImageType::VIS) {};
};

class LWIRDescriptor : public Descriptor {
	public:
		LWIRDescriptor(const char *id) : Descriptor(id, ImageType::LWIR) {};
};

// Descriptor wrappers
class GCHDescriptor : public VISDescriptor {
	public:
		GCHDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(VISImage vis, Segmentation *S) override;
};

class ACCDescriptor : public VISDescriptor {
	public:
		ACCDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(VISImage vis, Segmentation *S) override;
};

class BICDescriptor : public VISDescriptor {
	public:
		BICDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(VISImage vis, Segmentation *S) override;
};

class LCHDescriptor : public VISDescriptor {
	public:
		LCHDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(VISImage vis, Segmentation *S) override;
};

class UnserDescriptor : public VISDescriptor {
	public:
		UnserDescriptor(const char* id) : VISDescriptor(id) {};
		virtual Mat describe(VISImage vis, Segmentation *S) override;
};

class SIGDescriptor : public LWIRDescriptor {
	public:
		SIGDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, Segmentation *S) override;
};

class REDUCEDSIGDescriptor : public LWIRDescriptor {
	public:
		REDUCEDSIGDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, Segmentation *S) override;
};

class MOMENTSDescriptor : public LWIRDescriptor {
	public:
		MOMENTSDescriptor(const char* id) : LWIRDescriptor(id) {};
		virtual Mat describe(LWIRImage image, Segmentation *S) override;
};

// Util functions
Mat convertHistogramColor(Descriptor *D, VISImage vis, 
	Segmentation *S, int dimensions, 
	Histogram *(*descriptor)(CImage*, Image*));

#endif
