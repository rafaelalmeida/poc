#include "models.h"

using namespace cv;

LWIRImage::LWIRImage(std::vector<cv::Mat> bands) {
	this->bands = bands;
}

cv::Mat LWIRImage::spectralSignature(cv::Mat mask) {
	int bands = this->bands.size();
	Mat sig(bands, 1, CV_32FC1);

	int c = 1;
	for (auto band : this->bands) {
		Mat roi = band & mask;
		sig.at<float>(c, 1) = mean(roi)[0];
		c++;
	}

	return sig;
}

void LWIRImage::upscale(cv::Size size) {
	for (auto& band : bands) {
		Mat resized;
		resize(band, resized, size);

		band = resized;
	}
}

LWIRImage LWIRImage::operator()(cv::Rect roi) {
	vector<Mat> newBands;
	newBands.reserve(bands.size());

	for (auto band : bands) {
		Mat newBand = band.clone();
		newBands.push_back(newBand(roi));
	}

	return LWIRImage(newBands);
}

cv::Mat LWIRImage::average() {
	Mat M = *(this->bands.begin());
	Mat avg(M.rows, M.cols, M.type());

	for (auto&& band : this->bands) {
		avg += band;
	}

	avg /= this->bands.size();

	return avg;
}

cv::Mat LWIRImage::equalized() {
	Mat src = floatImageTo8UC3Image(this->average());
	Mat channels[3];
	cv::split(src, channels);
	src = channels[0];

	Mat dst(src.rows, src.cols, src.type());

	equalizeHist(src, dst);

	return dst;
}
