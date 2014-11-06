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

cv::Mat GCHDescriptor::describe(cv::Mat image, std::list<cv::Mat> masks) {
	const int HIST_BINS = 64;

	CImage *cimg = matToRawColor(image);

	Mat samples(masks.size(), HIST_BINS, CV_32FC1);

	int currentMask = 0;
	for (auto mask : masks) {
		Image *cMask = matToRawGray(mask);

		Histogram *hist = GCH(cimg, cMask);
		assert(hist->n == HIST_BINS);

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