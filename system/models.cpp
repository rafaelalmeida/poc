#include "models.h"

using namespace cv;

LWIRImage::LWIRImage(std::vector<cv::Mat> bands) {
	this->bands = bands;

	// Finds the extremes for normalization
	Mat extremes(1, bands.size() * 2, CV_32FC1);
	int i = 0;
	for (auto band : bands) {
		double minValL, maxValL;
		minMaxIdx(band, &minValL, &maxValL);

		extremes.at<float>(0, i) = (float) minValL;
		extremes.at<float>(0, i+1) = (float) maxValL;
		i += 2;
	}

	// Saves the extreme values
	double minValGlobal, maxValGlobal;
	minMaxIdx(extremes, &minValGlobal, &maxValGlobal);
	this->minVal = (float) minValGlobal;
	this->maxVal = (float) maxValGlobal;
}

cv::Mat LWIRImage::spectralSignature(cv::Mat mask) {
	int bands = this->bands.size();
	Mat sig(1, bands, CV_32FC1);

	int c = 0;
	for (auto band : this->bands) {
		Mat cropped;
		if (this->roi.area() > 0) {
			cropped = band(this->roi);
		}
		else {
			cropped = band;
		}
		
		assert((cropped.size() == mask.size()) && 
			"Mask is not the correct size");

		sig.at<float>(0, c) = mean(cropped, mask)[0];

		c++;
	}

	return sig;
}

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Mat mask) {
	Mat sig = this->spectralSignature(mask);
	cerr << sig << endl;

	for (int col = 0; col < sig.cols; col++) {
		float val = sig.at<float>(0, col);
		float nVal = (val - minVal) / (maxVal - minVal);

		sig.at<float>(0, col) = nVal;
	}

	return sig;
}

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Point point) {
	Mat mask = Mat::zeros(bands[0].size(), CV_8UC1);
	mask.at<unsigned char>(point) = 255;

	return this->normalizedSpectralSignature(mask);
}

void LWIRImage::upscale(cv::Size size) {
	for (auto& band : bands) {
		Mat resized;
		resize(band, resized, size);

		band = resized;
	}
}

cv::Mat LWIRImage::average() {
	Mat full = *(this->bands.begin());
	Mat M = full(roi);
	Mat avg(M.rows, M.cols, M.type());

	for (auto&& band : this->bands) {
		avg += band(roi);
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

int LWIRImage::numBands() {
	return this->bands.size();
}

void LWIRImage::setRoi(cv::Rect roi) {
	this->roi = roi;
}

CoverMap::CoverMap(Mat training) {
	_map = training;
}

cv::Mat CoverMap::asMat() {
	return _map;
}

float CoverMap::getRegionClass(cv::Mat mask) {
	assert(_map.type() == CV_8UC1);
	assert(mask.type() == CV_8UC1);
	assert((_map.rows == mask.rows) && (_map.cols = mask.cols));

	typedef unsigned char uchar;

	Counter<uchar> counter;

	for (int row = 0; row < _map.rows; row++) {
		for (int col = 0; col < _map.cols; col++) {
			if (mask.at<uchar>(row, col)) {
				uchar val = _map.at<uchar>(row, col);
				counter.inc(val);
			}
		}
	}

	return (float) counter.top();
}

cv::Mat CoverMap::coloredMap() {
	assert(_map.type() == CV_8UC1);

	Mat map(_map.rows, _map.cols, CV_8UC3);

	list<Mat> regions = segmentation::getColorBlobs(_map);
	for (auto r : regions) {
		vector<Mat> contours;
		findContours(r, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
		float label = segmentation::getSegmentLabel(_map, r);

		Scalar color;
		if (label == 0) { // UNCLASSIFIED
			color = Scalar(0, 0, 0); // BLACK
		}
		else if (label == 1) { // ROAD
			color = Scalar(255, 0, 255); // MAGENTA
		}
		else if (label == 2) { // TREES
			color = Scalar(0, 255, 0); // GREEN
		}
		else if (label == 3) { // RED ROOF
			color = Scalar(0, 0, 255); // RED
		}
		else if (label == 4) { // GREY ROOF
			color = Scalar(255, 255, 0); // CYAN
		}
		else if (label == 5) { // CONCRETE ROOF
			color = Scalar(128, 0, 128); // PURPLE
		}
		else if (label == 6) { // VEGETATION
			color = Scalar(87, 139, 46); // SEA GREEN
		}
		else if (label == 7) { // BARE SOIL
			color = Scalar(0, 255, 255); // YELLOW
		}

		drawContours(map, contours, -1, color, CV_FILLED);
	}

	return map;
}

std::map<unsigned char, int> CoverMap::getClassesCounts() {
	Counter<unsigned char> counter;

	for (int row = 0; row < _map.rows; row++) {
		for (int col = 0; col < _map.cols; col++) {
			counter.inc(_map.at<unsigned char>(row, col));
		}
	}

	return counter.getCounts();
}

cv::Size LWIRImage::size() {
	return bands[0].size();
}
