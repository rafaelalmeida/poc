#include "description.h"

using namespace cv;
using namespace std;

Mat Descriptor::describe(Mat image, SparseMat mask) {
	list<SparseMat> masks = {mask};
	return this->describe(image, masks);
}

Mat Descriptor::describe(LWIRImage image, SparseMat mask) {
	list<SparseMat> masks = {mask};
	return this->describe(image, masks);
}

Mat Descriptor::describe(Mat image, list<SparseMat> masks) {
	assert(false && "Not impelemented. Maybe this is a LWIR descriptor?");
	return Mat();
}

Mat Descriptor::describe(LWIRImage image, list<SparseMat> masks) {
	assert(false && "Not impelemented. Maybe this is a VIS descriptor?");
	return Mat();
}

// LWIR Descriptors implementation
Mat getSpectralSignature(LWIRImage image, list<SparseMat> masks, bool reduced) {
	int cols = reduced ? image.numReducedBands() : image.numBands();
	Mat samples(masks.size(), cols, CV_32FC1);

	int i = 0;
	for (auto mask : masks) {
		Mat sig = image.normalizedSpectralSignature(densify(mask), reduced);

		sig.row(0).copyTo(samples.row(i));

		i++;
	}

	return samples;
}

Mat SIGDescriptor::describe(LWIRImage image, list<SparseMat> masks) {
	return getSpectralSignature(image, masks, false);
}

Mat REDUCEDSIGDescriptor::describe(LWIRImage image, list<SparseMat> masks) {
	return getSpectralSignature(image, masks, true);
}

Mat MOMENTSDescriptor::describe(LWIRImage image, list<SparseMat> masks) {
	const int MAX_MOMENT_ORDER = 4;

	Mat samples(masks.size(), MAX_MOMENT_ORDER, CV_32FC1);

	int i = 0;
	for (auto mask : masks) {
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

// VIS descriptor wrappers

Mat GCHDescriptor::describe(Mat image, list<SparseMat> masks) {
	return convertHistogramColor(image, masks, GCHDimensions(), &GCH);
}

Mat ACCDescriptor::describe(Mat image, list<SparseMat> masks) {
	return convertHistogramColor(image, masks, ACCDimensions(), &ACC);
}

Mat BICDescriptor::describe(Mat image, list<SparseMat> masks) {
	return convertHistogramColor(image, masks, BICDimensions(), &BIC);
}

Mat LCHDescriptor::describe(Mat image, list<SparseMat> masks) {
	return convertHistogramColor(image, masks, LCHDimensions(), &LCH);
}

Mat UnserDescriptor::describe(Mat image, list<SparseMat> masks) {
	// TODO: refactor to leverage convertHistogramColor code

	int dimensions = UnserDimensions();
	
	// Convert to grayscale
	Mat gray;
	cvtColor(image, gray, CV_BGR2GRAY);

	Image *img = matToRawGray(gray);

	Mat samples(masks.size(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : masks) {
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

// TODO: fix (dimensions reported differently for training and classification)
Mat CBCDescriptor::describe(Mat image, list<SparseMat> masks) {
	assert(masks.size() > 0);

	CImage *cimg = matToRawColor(image);

	list<Ap_FeatureVector1D*> fvs;

	// Extract all the feature vectors
	int nFeatures = 0;

	for (auto mask : masks) {
		Image *cMask = matToRawGray(densify(mask));

		Ap_FeatureVector1D *fv = CBC(cimg, cMask, &nFeatures);
		fvs.push_back(fv);

		DestroyImage(&cMask);
	}

	cerr << "FEATURES = " << nFeatures << endl;

	DestroyCImage(&cimg);

	// Create the samples matrix with the dimensions reported by the descriptor
	Mat samples(masks.size(), nFeatures, CV_32FC1);

	// Converts the feature vectors to the Mat
	int currentRow = 0;
	for (auto fv : fvs) {
		for (int i = 0; i < (*fv)->n; i++) {
			samples.at<float>(currentRow, i) = (float) (*fv)->X[i];
		}

		// TODO: fix segfault
		//DestroyFeatureVector1D(fv);
	}

	return samples;
}

Mat convertHistogramColor(Mat image, list<SparseMat> masks, 
	int dimensions, Histogram *(*descriptor)(CImage*, Image*)) {
	
	CImage *cimg = matToRawColor(image);

	Mat samples(masks.size(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : masks) {
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

Descriptor::Descriptor(const char *id) : _id(string(id)) {}
