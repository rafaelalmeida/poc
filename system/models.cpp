#include "models.h"

using namespace cv;

LWIRImage::LWIRImage(std::vector<cv::Mat> bands) {
	// Saves the bands
	this->bands = bands;

	// Finds the extremes for normalization
	this->minMaxAcrossBands(bands, &(this->minVal), &(this->maxVal));

	// Reduce dimensionality via PCA
	reduceDimensionality(PCA_COMPONENTS);
}

cv::Mat LWIRImage::spectralSignature(cv::Mat mask, bool reduced) {
	vector<Mat> &myBands = reduced ? this->reducedBands : this->bands;

	int numBands = myBands.size();
	Mat sig(1, numBands, CV_32FC1);

	int c = 0;
	for (auto band : myBands) {
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

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Mat mask, bool reduced) {
	Mat sig = this->spectralSignature(mask, reduced);

	// Determines the correct extreme values for normalization
	float minVal, maxVal;
	if (reduced) {
		minVal = this->minValReduced;
		maxVal = this->maxValReduced;
	}
	else {
		minVal = this->minVal;
		maxVal = this->maxVal;
	}

	// Normalizes the values
	for (int col = 0; col < sig.cols; col++) {
		float val = sig.at<float>(0, col);
		float nVal = (val - minVal) / (maxVal - minVal);

		sig.at<float>(0, col) = nVal;
	}

	return sig;
}

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Point point, bool reduced) {
	Mat mask = Mat::zeros(this->size(), CV_8UC1);
	mask.at<unsigned char>(point) = 255;

	return this->normalizedSpectralSignature(mask, reduced);
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

int LWIRImage::numReducedBands() {
	return this->reducedBands.size();
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

void LWIRImage::reduceDimensionality(int keep) {
	Mat data = formatImagesForPCA(bands);
	PCA pca(data, Mat(), CV_PCA_DATA_AS_ROW, keep);

	// Converts the eigenvectors back into images
	Size bandSize = this->size();
	reducedBands.resize(keep);

	for (int i = 0; i < keep; i++) {
		reducedBands[i] = pca.eigenvectors.row(i).reshape(1, bandSize.height);
	}

	// Finds extremes for normalization
	this->minMaxAcrossBands(reducedBands, &(this->minValReduced), 
		&(this->maxValReduced));
}

void LWIRImage::minMaxAcrossBands(std::vector<cv::Mat> bands, float *minVal, 
	float *maxVal) {

	// Find the extreme of each band
	Mat extremes(1, bands.size() * 2, CV_32FC1);
	int i = 0;
	for (auto band : bands) {
		double minValL, maxValL;
		minMaxIdx(band, &minValL, &maxValL);

		extremes.at<float>(0, i) = (float) minValL;
		extremes.at<float>(0, i+1) = (float) maxValL;
		i += 2;
	}

	// Find the global extremes
	double minValD, maxValD;
	minMaxIdx(extremes, &minValD, &maxValD);
	*minVal = (float) minValD;
	*maxVal = (float) maxValD;
}