#include "description.h"

using namespace cv;
using namespace std;

// ==================
// Base class methods
// ==================
Mat Descriptor::describe(VISImage vis, Segmentation S) {
	assert(false && "Not impelemented. Maybe this is a LWIR descriptor?");
	return Mat();
}

Mat Descriptor::describe(LWIRImage image, Segmentation S) {
	assert(false && "Not impelemented. Maybe this is a VIS descriptor?");
	return Mat();
}

// ===============================
// LWIR Descriptors implementation
// ===============================
Mat getSpectralSignature(LWIRImage image, Segmentation S, bool reduced) {
	int cols = reduced ? image.numReducedBands() : image.numBands();
	Mat samples(S.regionCount(), cols, CV_32FC1);

	int i = 0;
	for (auto mask : S.getRegionMasks()) {
		Mat sig = image.normalizedSpectralSignature(densify(mask), reduced);
		sig.row(0).copyTo(samples.row(i++));
	}

	return samples;
}

Mat SIGDescriptor::describe(LWIRImage image, Segmentation S) {
	return getSpectralSignature(image, S.getRegionMasks(), false);
}

Mat REDUCEDSIGDescriptor::describe(LWIRImage image, Segmentation S) {
	return getSpectralSignature(image, S.getRegionMasks(), true);
}

Mat MOMENTSDescriptor::describe(LWIRImage image, Segmentation S) {
	const int MAX_MOMENT_ORDER = 4;

	Mat samples(S.regionCount(), MAX_MOMENT_ORDER, CV_32FC1);

	int i = 0;
	for (auto mask : S.getRegionMasks()) {
		Mat sig = image.spectralSignature(densify(mask));

		const float *p = sig.ptr<float>(0);
		vector<float> values(p, p + sig.cols);

		vector<float> moments = statistics::moments(values, MAX_MOMENT_ORDER);

		int j = 0;
		for (auto m : moments) {
			samples.at<float>(i, j++) = moments[j];
		}

		i++;
	}

	return samples;
}

// =======================
// VIS descriptor wrappers
// =======================
Mat GCHDescriptor::describe(VISImage vis, Segmentation S) {
	return convertHistogramColor(vis.asMat(), S.getRegionMasks(), 
		GCHDimensions(), &GCH);
}

Mat ACCDescriptor::describe(VISImage vis, Segmentation S) {
	return convertHistogramColor(vis.asMat(), S.getRegionMasks(), 
		ACCDimensions(), &ACC);
}

Mat BICDescriptor::describe(VISImage vis, Segmentation S) {
	return convertHistogramColor(vis.asMat(), S.getRegionMasks(), 
		BICDimensions(), &BIC);
}

Mat LCHDescriptor::describe(VISImage vis, Segmentation S) {
	return convertHistogramColor(vis.asMat(), S.getRegionMasks(), 
		LCHDimensions(), &LCH);
}

Mat UnserDescriptor::describe(VISImage vis, Segmentation S) {
	// TODO: refactor to leverage convertHistogramColor code

	int dimensions = UnserDimensions();
	
	// Convert to grayscale
	Mat gray;
	cvtColor(vis.asMat(), gray, CV_BGR2GRAY);

	Image *img = matToRawGray(gray);

	Mat samples(S.regionCount(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : S.getRegionMasks()) {
		Image *cMask = matToRawGray(densify(mask));

		Histogram *hist = Unser(img, cMask);
		assert(hist->n == dimensions);

		for (int i = 0; i < hist->n; i++) {
			samples.at<float>(currentMask, i) = (float) hist->v[i];
		}

		DestroyHistogram(&hist);
		DestroyImage(&cMask);

		currentMask++;
	}

	DestroyImage(&img);

	return samples;
}

Mat convertHistogramColor(VISImage vis, Segmentation S, 
	int dimensions, Histogram *(*descriptor)(CImage*, Image*)) {
	
	CImage *cimg = matToRawColor(vis.asMat());

	Mat samples(S.regionCount(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : S.getRegionMasks()) {
		Image *cMask = matToRawGray(densify(mask));

		Histogram *hist = descriptor(cimg, cMask);
		assert(hist->n == dimensions);

		for (int i = 0; i < hist->n; i++) {
			samples.at<float>(currentMask, i) = (float) hist->v[i];
		}

		DestroyHistogram(&hist);
		DestroyImage(&cMask);

		currentMask++;
	}

	DestroyCImage(&cimg);

	return samples;
}

std::string Descriptor::getID() {
	return _id;
}

Descriptor::Descriptor(const char *id, ImageType type) 
	: _id(string(id)), _type(type) {}
