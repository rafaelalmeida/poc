#include "description.h"

using namespace cv;
using namespace std;

// ==================
// Base class methods
// ==================
Mat Descriptor::describe(VISImage vis, Segmentation *S) {
	FATAL_ERROR("Not impelemented. Maybe this is a LWIR descriptor?");
	return Mat();
}

Mat Descriptor::describe(LWIRImage image, Segmentation *S) {
	FATAL_ERROR("Not impelemented. Maybe this is a VIS descriptor?");
	return Mat();
}

// ===============================
// LWIR Descriptors implementation
// ===============================
Mat getSpectralSignature(Descriptor *D, LWIRImage image, Segmentation *S, 
	bool reduced) {

	int cols = reduced ? image.numReducedBands() : image.numBands();
	Mat samples(S->regionCount(), cols, CV_32FC1);

	int i = 0;
	for (auto region : S->getRegions()) {
		Mat sig;
		if (region.getRepresentationMode() == PIXEL) {
			Point P = region.getPoint();
			sig = image.normalizedSpectralSignature(P, reduced);
		}
		else {
			SparseMat M = region.getMask();
			sig = image.normalizedSpectralSignature(densify(M), reduced);
		}

		sig.row(0).copyTo(samples.row(i++));
		S->upcountDescription(D);
	}

	return samples;
}

Mat SIGDescriptor::describe(LWIRImage image, Segmentation *S) {
	return getSpectralSignature(this, image, S, false);
}

Mat REDUCEDSIGDescriptor::describe(LWIRImage image, Segmentation *S) {
	return getSpectralSignature(this, image, S, true);
}

Mat MOMENTSDescriptor::describe(LWIRImage image, Segmentation *S) {
	const int MAX_MOMENT_ORDER = 4;

	Mat samples(S->regionCount(), MAX_MOMENT_ORDER, CV_32FC1);

	int i = 0;
	for (auto region : S->getRegions()) {
		// Get the original spectral signature
		Mat sig;
		if (region.getRepresentationMode() == PIXEL) {
			Point P = region.getPoint();
			sig = image.spectralSignature(P);
		}
		else {
			SparseMat M = region.getMask();
			sig = image.spectralSignature(densify(M));
		}

		// Convert into a vector
		const float *p = sig.ptr<float>(0);
		vector<float> values(p, p + sig.cols);

		// Calculate statistical moments of this vector
		vector<float> moments = statistics::moments(values, MAX_MOMENT_ORDER);

		// Put the moments in the feature matrix
		int j = 0;
		for (auto m : moments) {
			samples.at<float>(i, j++) = moments[j];
		}

		// Update and report progress
		i++;
		S->upcountDescription(this);
	}

	return samples;
}

// =======================
// VIS descriptor wrappers
// =======================
Mat GCHDescriptor::describe(VISImage vis, Segmentation *S) {
	return convertHistogramColor(this, vis.asMat(), S, GCHDimensions(), &GCH);
}

Mat ACCDescriptor::describe(VISImage vis, Segmentation *S) {
	return convertHistogramColor(this, vis.asMat(), S, ACCDimensions(), &ACC);
}

Mat BICDescriptor::describe(VISImage vis, Segmentation *S) {
	return convertHistogramColor(this, vis.asMat(), S, BICDimensions(), &BIC);
}

Mat LCHDescriptor::describe(VISImage vis, Segmentation *S) {
	return convertHistogramColor(this, vis.asMat(), S, LCHDimensions(), &LCH);
}

Mat UnserDescriptor::describe(VISImage vis, Segmentation *S) {
	// TODO: refactor to leverage convertHistogramColor code

	int dimensions = UnserDimensions();
	
	// Convert to grayscale
	Mat gray;
	cvtColor(vis.asMat(), gray, CV_BGR2GRAY);

	Image *img = matToRawGray(gray);

	Mat samples(S->regionCount(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : S->getRegionMasks()) {
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

Mat convertHistogramColor(Descriptor *D, VISImage vis, 
	Segmentation *S, int dimensions, 
	Histogram *(*descriptor)(CImage*, Image*)) {
	
	CImage *cimg = matToRawColor(vis.asMat());

	Mat samples(S->regionCount(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : S->getRegionMasks()) {
		Image *cMask = matToRawGray(densify(mask));

		Histogram *hist = descriptor(cimg, cMask);
		assert(hist->n == dimensions);

		for (int i = 0; i < hist->n; i++) {
			samples.at<float>(currentMask, i) = (float) hist->v[i];
		}

		DestroyHistogram(&hist);
		DestroyImage(&cMask);

		currentMask++;
		S->upcountDescription(D);
	}

	DestroyCImage(&cimg);

	return samples;
}

std::string Descriptor::getID() {
	return _id;
}

Descriptor::Descriptor(const char *id, ImageType type) 
	: _id(string(id)), _type(type) {}
