#include "description.h"

using namespace cv;
using namespace std;

cv::Mat Descriptor::describe(cv::Mat image, cv::Mat mask) {
	list<Mat> masks = {mask};
	return this->describe(image, masks);
}

cv::Mat Descriptor::describe(LWIRImage image, cv::Mat mask) {
	list<Mat> masks = {mask};
	return this->describe(image, masks);
}

cv::Mat Descriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	assert(false && "Not impelemented. Maybe this is a LWIR descriptor?");
	return Mat();
}

cv::Mat Descriptor::describe(LWIRImage image, std::list<cv::Mat> masks) {
	assert(false && "Not impelemented. Maybe this is a VIS descriptor?");
	return Mat();
}

// LWIR Descriptors implementation

cv::Mat SIGDescriptor::describe(LWIRImage image, std::list<cv::Mat> masks) {
	Mat samples(masks.size(), image.numBands(), CV_32FC1);

	int i = 0;
	for (auto mask : masks) {
		Mat sig = image.normalizedSpectralSignature(mask);

		sig.row(0).copyTo(samples.row(i));

		i++;
	}

	return samples;
}

cv::Mat MOMENTSDescriptor::describe(LWIRImage image, std::list<cv::Mat> masks) {
	const int MAX_MOMENT_ORDER = 4;

	Mat samples(masks.size(), MAX_MOMENT_ORDER, CV_32FC1);

	int i = 0;
	for (auto mask : masks) {
		Mat sig = image.spectralSignature(mask);

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

cv::Mat GCHDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	return convertHistogramColor(image, masks, GCHDimensions(), &GCH);
}

cv::Mat ACCDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	return convertHistogramColor(image, masks, ACCDimensions(), &ACC);
}

cv::Mat BICDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	return convertHistogramColor(image, masks, BICDimensions(), &BIC);
}

cv::Mat LCHDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	return convertHistogramColor(image, masks, LCHDimensions(), &LCH);
}

cv::Mat UnserDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	// TODO: refactor to leverage convertHistogramColor code

	int dimensions = UnserDimensions();
	
	// Convert to grayscale
	Mat gray;
	cvtColor(image, gray, CV_BGR2GRAY);

	Image *img = matToRawGray(gray);

	Mat samples(masks.size(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : masks) {
		Image *cMask = matToRawGray(mask);

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
cv::Mat CBCDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	assert(masks.size() > 0);

	CImage *cimg = matToRawColor(image);

	list<Ap_FeatureVector1D*> fvs;

	// Extract all the feature vectors
	int nFeatures = 0;

	for (auto mask : masks) {
		Image *cMask = matToRawGray(mask);

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

cv::Mat convertHistogramColor(cv::Mat image, std::list<cv::Mat> masks, 
	int dimensions, Histogram *(*descriptor)(CImage*, Image*)) {
	
	CImage *cimg = matToRawColor(image);

	Mat samples(masks.size(), dimensions, CV_32FC1);

	int currentMask = 0;
	for (auto mask : masks) {
		Image *cMask = matToRawGray(mask);

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
