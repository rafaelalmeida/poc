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
		Mat sig = image.spectralSignature(mask);

		Mat thisFeatureVector = sig.t();
		thisFeatureVector.row(0).copyTo(samples.row(i));

		i++;
	}

	return samples;
}

cv::Mat ENERGYDescriptor::describe(LWIRImage image, std::list<cv::Mat> masks) {
	Mat samples(masks.size(), 1, CV_32FC1);

	int i = 0;
	for (auto mask : masks) {
		Mat sig = image.spectralSignature(mask);
		float meanEnergy = mean(sig)[0];

		samples.at<float>(i, 0) = meanEnergy;

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
