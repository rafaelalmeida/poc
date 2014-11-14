#include "models.h"

using namespace cv;

LWIRImage::LWIRImage(std::vector<cv::Mat> bands) {
	// Saves the bands
	this->bands = bands;

	// Finds the extremes for normalization
	this->minMaxAcrossBands(bands, &(this->minVal), &(this->maxVal));

	// Init the ROI as blank
	this->roi = Rect(0, 0, 0, 0);
}

cv::Mat LWIRImage::spectralSignature(cv::Mat mask, bool reduced) {
	vector<Mat> &myBands = reduced ? this->reducedBands : this->bands;

	int numBands = myBands.size();
	Mat sig(1, numBands, CV_32FC1);

	int c = 0;
	for (auto band : myBands) {
		Mat B = (this->roi.area() > 0) ? band(this->roi) : band;

		assert((B.size() == mask.size()) && 
			"Mask is not the correct size");

		sig.at<float>(0, c) = mean(B, mask)[0];

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

void LWIRImage::rescale(float scale, ResamplingMethod mode) {
	for (auto& band : bands) {
		Mat resized;
		resize(band, resized, Size(), scale, scale, 
			translateInterpolationMode(mode));

		band = resized;
	}
}

cv::Mat LWIRImage::average() {
	Mat full = this->bands.front();

	// Apply ROI, if there is one
	Mat M = full;
	if (roi.area() > 0) {
		M = full(roi);
	}
	
	Mat avg(M.rows, M.cols, M.type());

	for (auto band : this->bands) {
		if (roi.area() > 0) {
			avg += band(roi);
		}
		else {
			avg += band;
		}
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

ThematicMap::ThematicMap(Mat training) {
	_map = training;
}

cv::Mat ThematicMap::asMat() {
	return _map;
}

float ThematicMap::getRegionClass(cv::Mat mask) {
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

cv::Mat ThematicMap::coloredMap() {
	Mat M = _map;
	assert(M.type() == CV_8UC1);
	Mat map = Mat::zeros(M.rows, M.cols, CV_8UC3);

	for (auto region : this->enumerateRegions()) {
		// Determine which region to paint
		Mat denseRegion = densify(region.first);
		vector<Mat> contours;
		findContours(denseRegion, contours, CV_RETR_EXTERNAL, 
			CV_CHAIN_APPROX_SIMPLE);

		// Select the right color to paint
		int label = region.second;
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

		// Paint the region
		drawContours(map, contours, -1, color, CV_FILLED);
	}

	return map;
}

std::map<unsigned char, int> ThematicMap::getClassesCounts() {
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

cv::Size ThematicMap::size() {
	return this->asMat().size();
}

std::list<std::pair<cv::SparseMat, int> > ThematicMap::enumerateRegions() {
	assert(_map.size().area() > 0 && "Error: map seems to be blank");
	assert(_map.type() == CV_8UC1 && "Error: map seems to be the wrong type.");

	// Use the fact that all non-undefined labels are positive to threshold
	// the image and enumerate the regions using contour finding
	Mat bin(_map.size(), CV_8UC1);
	threshold(_map, bin, 1, 255, THRESH_BINARY);
	vector<vector<Point> > contours;
	findContours(bin, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

	// Iterate over the found contours and build the region masks
	list<pair<SparseMat, int> > regions;
	for (int idx = 0; idx < contours.size(); idx++) {
		// Build the mask
		Mat mask = Mat::zeros(_map.size(), CV_8UC1);
		drawContours(mask, contours, idx, Scalar(255), CV_FILLED);

		// Get the region label. We assume, by construction, that all pixels
		// inside a region have the same label. So, for efficiency, we take
		// the first one we find and look at its label.
		int label = 0;
		bool keepLooking = true;

		for (int row = 0; keepLooking && (row < _map.rows); row++) {
			for (int col = 0; keepLooking && (col < _map.cols); col++) {
				if (mask.at<unsigned char>(row, col) == 255) {
					label = (int) _map.at<unsigned char>(row, col);
					keepLooking = false;
				}
			}
		}

		// By construction, we only look at regions whose label is positive,
		// so if label is still the sentinel value, something wrong happened.
		assert(label > 0);

		regions.push_back(make_pair(SparseMat(mask), label));
	}

	return regions;
}

std::vector<ThematicMap> ThematicMap::split(int k) {
	// Grab the list of regions and get it into a vector
	list<pair<SparseMat, int> > regions = this->enumerateRegions();
	vector<pair<SparseMat, int> > regionsV(regions.begin(), regions.end());

	// Create a vector of folds
	vector<int> folds(regions.size());

	// Assign the folds
	int foldToGo = 0;
	for (int i = 0; i < folds.size(); i++) {
		folds[i] = foldToGo++;

		// Circle back to first fold
		if (foldToGo == k) {
			foldToGo = 0;
		}
	}
	
	// Shuffle the assigned folds
	auto engine = default_random_engine(RANDOM_SEED);
	shuffle(begin(folds), end(folds), engine);

	// Group the shuffled segments together
	vector<list<pair<SparseMat, int> > > splitRegionGroups(k);
	for (int i = 0; i < folds.size(); i++) {
		splitRegionGroups[folds[i]].push_back(regionsV[i]);
	}

	// Convert the grouped segments into matrixes with the correct maps
	vector<ThematicMap> splits;
	splits.reserve(k);

	for (auto& g : splitRegionGroups) {
		Mat M = Mat::zeros(_map.size(), _map.type());

		// Paint the regions into the blank map
		for (auto r : g) {
			Mat R = densify(r.first);
			M += r.second * R / 255;
		}

		splits.push_back(ThematicMap(M));
	}

	return splits;
}
