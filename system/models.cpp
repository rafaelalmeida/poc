#include "models.h"

using namespace cv;

LWIRImage::LWIRImage(std::vector<cv::Mat> bands) : RSImage(LWIR) {
	// Saves the bands
	this->bands = bands;

	// Finds the extremes for normalization
	this->minMaxAcrossBands(bands, &(this->minVal), &(this->maxVal));

	// Init the ROI as blank
	this->roi = Rect(0, 0, 0, 0);
}

LWIRImage::LWIRImage() : RSImage(LWIR) {}

cv::Mat LWIRImage::doGetSpectralSignature(cv::Mat *mask, cv::Point *point, 
	bool reduced) {

	// Check arguments
	if (mask == NULL && point == NULL) {
		FATAL_ERROR("Either mask or point needs to be provided to calculate \
			spectral signature");
	}

	vector<Mat> &myBands = reduced ? this->reducedBands : this->bands;

	int numBands = myBands.size();
	Mat sig(1, numBands, CV_32FC1);

	int c = 0;
	for (auto band : myBands) {
		Mat B = (this->roi.area() > 0) ? band(this->roi) : band;

		// Get the signature from mask or point
		if (mask != NULL) {
			assert((B.size() == mask->size()) && "Mask has wrong size");
			sig.at<float>(0, c) = mean(B, *mask)[0];
		}
		else {
			sig.at<float>(0, c) = B.at<float>(point->y, point->x);
		}

		c++;
	}

	return sig;
}

cv::Mat LWIRImage::spectralSignature(cv::Mat mask, bool reduced) {
	return doGetSpectralSignature(&mask, NULL, reduced);
}

cv::Mat LWIRImage::spectralSignature(cv::Point point, bool reduced) {
	return doGetSpectralSignature(NULL, &point, reduced);
}

void LWIRImage::doNormalizeSpectralSignature(cv::Mat signature, bool reduced) {
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
	for (int col = 0; col < signature.cols; col++) {
		float val = signature.at<float>(0, col);
		float nVal = (val - minVal) / (maxVal - minVal);

		signature.at<float>(0, col) = nVal;
	}
}

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Mat mask, bool reduced) {
	Mat sig = this->spectralSignature(mask, reduced);
	doNormalizeSpectralSignature(sig, reduced);

	return sig;
}

cv::Mat LWIRImage::normalizedSpectralSignature(cv::Point point, bool reduced) {
	Mat sig = this->spectralSignature(point, reduced);
	doNormalizeSpectralSignature(sig, reduced);

	return sig;
}

void LWIRImage::rescale(float scale, InterpolationMode mode) {
	// Resize normal bands
	for (auto& band : bands) {
		Mat resized;
		resize(band, resized, Size(), scale, scale, 
			translateInterpolationMode(mode));

		band = resized;
	}

	// Resize reduced bands
	for (auto& band : reducedBands) {
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
	_map = training.clone();
}

ThematicMap::ThematicMap(cv::Size size) {
	_map = Mat::zeros(size, CV_8UC1);
}

ThematicMap::ThematicMap() {
	_map = Mat(); // Dummy map, just so we can create the object
}

cv::Mat ThematicMap::asMat() {
	return _map;
}

cv::Mat ThematicMap::coloredMap() {
	assert(_map.type() == CV_8UC1);
	Mat map = Mat::zeros(_map.size(), CV_8UC3);

	vector<vector<Point> > contours;
	for (int row = 0; row < _map.rows; row++) {
		for (int col = 0; col < _map.cols; col++) {
			int label = (int) _map.at<unsigned char>(row, col);

			// Select the right color to paint
			Vec3b color;
			if (label == 0) { // UNCLASSIFIED
				color = Vec3b(0, 0, 0); // BLACK
			}
			else if (label == 1) { // ROAD
				color = Vec3b(255, 0, 255); // MAGENTA
			}
			else if (label == 2) { // TREES
				color = Vec3b(0, 255, 0); // GREEN
			}
			else if (label == 3) { // RED ROOF
				color = Vec3b(0, 0, 255); // RED
			}
			else if (label == 4) { // GREY ROOF
				color = Vec3b(255, 255, 0); // CYAN
			}
			else if (label == 5) { // CONCRETE ROOF
				color = Vec3b(128, 0, 128); // PURPLE
			}
			else if (label == 6) { // VEGETATION
				color = Vec3b(87, 139, 46); // SEA GREEN
			}
			else if (label == 7) { // BARE SOIL
				color = Vec3b(0, 255, 255); // YELLOW
			}
			else {
				FATAL_ERROR("Error! Unknown label");
			}

			// Paint the pixel
			map.at<Vec3b>(row, col) = color;
		}
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

	// Finds the connected color components representing each region
	list<SparseMat> connectedComponents = segmentation::getColorBlobs(_map);

	// Fill region list
	list<pair<SparseMat, int> > regions;
	for (auto cc : connectedComponents) {
		// Build the mask
		Mat mask = densify(cc);

		// Skip connected components whose value is zero
		if (countNonZero(_map & mask) == 0) {
			continue;
		}

		// Get the region label. We assume, by construction, that all pixels
		// inside a region have the same label. So, for efficiency, we take
		// the first one we find and look at its label.
		int label = -1;
		bool keepLooking = true;

		for (int row = 0; keepLooking && (row < _map.rows); row++) {
			for (int col = 0; keepLooking && (col < _map.cols); col++) {
				if (mask.at<unsigned char>(row, col) == 255) {
					label = (int) _map.at<unsigned char>(row, col);
					keepLooking = false;
				}
			}
		}

		// By construction, we only look at regions whose label is non-negative,
		// so if label is still the sentinel value, something wrong happened.
		assert(label > -1);

		regions.push_back(make_pair(cc, label));
	}

	return regions;
}

std::vector<std::pair<ThematicMap, ThematicMap> > ThematicMap::split(
	KFolder folder) {

	// Grab the list of regions and get it into a vector
	list<pair<SparseMat, int> > regions = this->enumerateRegions();
	vector<pair<SparseMat, int> > regionsV(regions.begin(), regions.end());

	// Create the split pairs (for now in index form)
	list<KFoldIndexes> foldsIdx = folder.makeFolds(regions.size());

	// Create the vector to hold the split pairs
	vector<pair<ThematicMap, ThematicMap> > folds;
	folds.reserve(foldsIdx.size());
	for (auto f : foldsIdx) {
		// The training and validation thematic maps (we will populate them
		// below)
		ThematicMap T(this->size());
		ThematicMap V(this->size());

		// Populate the training map
		for (auto idx : f.first) {
			pair<SparseMat, int> region = regionsV[idx];
			T._map += region.second * (densify(region.first) / 255);
		}

		// Populate the validation map
		for (auto idx : f.second) {
			pair<SparseMat, int> region = regionsV[idx];
			V._map += region.second * (densify(region.first) / 255);
		}

		// Add the pair to the list
		folds.push_back(make_pair(T, V));
	}

	return folds;
}

void ThematicMap::resize(cv::Size newSize) {
	Mat R;
	cv::resize(_map, R, newSize, 0, 0, TRAINING_INTERPOLATION_MODE);
	_map = R;
}

void ThematicMap::combine(ThematicMap T) {
	_map += T.asMat();
}

ThematicMap ThematicMap::clone() {
	return ThematicMap(_map);
}

VISImage::VISImage(cv::Mat vis) 
	: _vis(vis),
	  RSImage(VIS) {}

VISImage::VISImage() : RSImage(VIS) {}

cv::Size VISImage::size() {
	return _vis.size();
}

cv::Mat VISImage::asMat() {
	return _vis;
}

void VISImage::rescale(float scale, InterpolationMode mode) {
	Mat visR;
	resize(_vis, visR, Size(), scale, scale, translateInterpolationMode(mode));
	_vis = visR;
}

void VISImage::setRoi(cv::Rect roi) {
	_vis = _vis(roi);
}

cv::SparseMat ThematicMap::getFullMask() {
	Mat mask(_map.rows, _map.cols, _map.type());
	threshold(_map, mask, 0, 255, THRESH_BINARY);
	return SparseMat(mask);
}

std::vector<float> *ThematicMap::classProbabilities() {
	vector<float> *P = new vector<float>();
	P->reserve(CLASS_COUNT);

	// Count the classes
	auto counts = getClassesCounts();
	for (int i = 1; i < CLASS_COUNT; i++) {
		if (counts.find((unsigned char) i) != counts.end()) {
			P->push_back((float) counts[(unsigned char) i]);
		}
	}

	// Get the total, for normalization
	float sum = 0;
	for (auto p : *P) {
		sum += p;
	}

	// Normalize so all probabilities add up to 1
	for (auto& p : *P) {
		p /= sum;
	}

	return P;
}
