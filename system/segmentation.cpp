#include "segmentation.h"

using namespace cv;
using namespace std;

using namespace segmentation;

Segmentation segmentation::segmentVIS_SLIC(Mat M, int regionSize, 
	int minRegionSize, float regularization) {
	// Adapted from https://github.com/davidstutz/vlfeat-slic-example/blob/
	// 79bc823000ec1bc7fb204b0613578e13d3dff1ae/vlfeat_slic_cli/main.cpp

	assert(M.channels() == 3);

	// Convert image to one-dimensional array.
	float* image = new float[M.rows*M.cols*M.channels()];
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.cols; ++j) {
			image[j + M.cols*i + M.cols*M.rows*0] = M.at<cv::Vec3b>(i, j)[0];
			image[j + M.cols*i + M.cols*M.rows*1] = M.at<cv::Vec3b>(i, j)[1];
			image[j + M.cols*i + M.cols*M.rows*2] = M.at<cv::Vec3b>(i, j)[2];
		}
	}

	// The algorithm will store the final segmentation in a one-dimensional 
	// array.
	vl_uint32* segmentation = new vl_uint32[M.rows*M.cols];
	vl_size height = M.rows;
	vl_size width = M.cols;
	vl_size channels = M.channels();

	// The region size defines the number of superpixels obtained.
	// Regularization describes a trade-off between the color term and the
	// spatial term.
	vl_size region = regionSize;
	vl_size minRegion = minRegionSize;
	float _regularization = regularization;

	// Run the actual segmentation procedure
	vl_slic_segment(segmentation, image, width, height, channels, region, 
		_regularization, minRegion);
            
    // Convert segmentation.
    int** labels = new int*[M.rows];
    for (int i = 0; i < M.rows; ++i) {
        labels[i] = new int[M.cols];
                
        for (int j = 0; j < M.cols; ++j) {
            labels[i][j] = (int) segmentation[j + M.cols*i];
        }
    }

    // Free memory used by float image and segmentation (we already converted
    // it above)
    delete [] image;
    delete [] segmentation;

    // Finds the total number of regions. The elements of the labels
    // matrix are the region labels starting from zero.
    int maxLabel = 0;
    for (int i = 0; i < M.rows; i++) {
    	for (int j = 0; j < M.cols; j++) {
    		int l = labels[i][j];
    		if (l > maxLabel) {
    			maxLabel = l;
    		}
    	}
    }

    int totalRegions = maxLabel + 1;

    // Create the segment masks
    vector<SparseMat> masks;
    for (int i = 0; i < totalRegions; i++) {
    	int sizes[2] = {M.rows, M.cols};
    	SparseMat SM(2, sizes, CV_8UC1);
    	masks.push_back(SM);
    }

    // Fill the segment masks where original data is non missing
    Mat nonMissingPixels = segmentation::makeNonMissingDataMask(M);
    for (int i = 0; i < M.rows; i++) {
    	for (int j = 0; j < M.cols; j++) {
    		if (nonMissingPixels.at<unsigned char>(i, j) == 255) {
    			int label = labels[i][j];
    			*(masks[label].ptr(i, j, true)) = 255;
    		}
    	}
    }

    // Delete the labels array
    for (int i = 0; i < M.rows; i++) {
    	delete [] labels[i];
    }
    delete [] labels;

    // Convert the mask vector into a list, filtering empty segments
    // (segments that were created by SLIC but on which all pixels are missing)
    list<SparseMat> masksList;
    for (auto m : masks) {
    	if (m.nzcount() > 0) {
    		masksList.push_back(m);
    	}
    }

    // Create the segmentation and return it
    return Segmentation(masksList);
}

Segmentation segmentation::segmentVISGrid(cv::Mat M, int tileSize) {
	list<SparseMat> segments;

	int regionsPerLine = (M.cols / tileSize) + (M.cols % tileSize);
	int regionsPerColumn = (M.rows / tileSize) + (M.rows % tileSize);

	Mat nonMissing = makeNonMissingDataMask(M);

	for (int i = 0; i < regionsPerLine; i++) {
		for (int j = 0; j < regionsPerColumn; j++) {
			Mat segment = Mat::zeros(M.size(), CV_8UC1);
			Rect region(i*tileSize, j*tileSize, tileSize, tileSize);

			rectangle(segment, region, Scalar(255), CV_FILLED);

			Mat filteredSegment = segment & nonMissing;
			if (countNonZero(filteredSegment) > 0) {
				segments.push_back(SparseMat(filteredSegment));
			}
		}
	}

	return Segmentation(segments);
}

Segmentation segmentation::segmentLWIRPixelated(LWIRImage& lwir, 
	VISImage& vis) {

	// Make a mask of non-missing pixels using the VIS image. We only consider
	// those pixels to create the LWIR segmentation.
	Mat nonMissing = makeNonMissingDataMask(vis.asMat());
	Mat M;
	resize(nonMissing, M, lwir.size(), 0, 0, INTER_NEAREST);

	// Adds every pixel
	list<Point> pixels;
	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			if (M.at<unsigned char>(row, col) == 255) {
				pixels.push_back(Point(col, row));
			}
		}
	}

	// Create and return the segmentation
	Segmentation S(pixels);
	S.setPixalated();
	return S;
}

// Create a binary mask to show where we have missing data and ignore
// those pixels later
Mat segmentation::makeNonMissingDataMask(Mat vis) {
	Mat gray(vis.size(), CV_8UC1);
	cvtColor(vis, gray, CV_BGR2GRAY);

	Mat bin(gray.size(), CV_8UC1);
	threshold(gray, bin, 0, 255, CV_THRESH_BINARY);

	// TODO? fill holes in mask to account for pixels inside non missing region
	// that have zero value

	return bin;
}

float segmentation::getSegmentLabel(Mat classificationMap, Mat mask) {
	assert(classificationMap.type() == CV_8UC1);
	assert(mask.type() == CV_8UC1);
	assert((classificationMap.rows == mask.rows) && (classificationMap.cols = mask.cols));

	return (float) mean(classificationMap, mask)[0];
}

list<SparseMat> segmentation::getColorBlobs(Mat posterized) {
	list<SparseMat> segments;

	Mat clone = posterized.clone();
	Mat mask(posterized.rows+2, posterized.cols+2, CV_8UC1, Scalar::all(0));

	for (int y = 0; y < clone.rows; y++) {
		for (int x = 0; x < clone.cols; x++) {
			// Find connected components
			if (mask.at<uchar>(y+1, x+1) == 0) {
				Mat previousMask = mask.clone();

				floodFill(clone, mask, Point(x,y), Scalar(0, 0, 0), 0, 
						  Scalar::all(0), Scalar::all(0));

				Mat difference = previousMask ^ mask;
				Mat segmentMask = 255 * difference(Range(1, difference.rows - 1),
												   Range(1, difference.cols - 1));

				segments.push_back(SparseMat(segmentMask));
			}
		}
	}

	return segments;
}

list<SparseMat> segmentation::getGrayBlobs(Mat posterized, int threshold) {
	// Alias and checking
	Mat M = posterized;
	assert(M.type() == CV_8UC1);

	// List to hold the blobs
	list<SparseMat> blobs;

	// Scan all levels from the threshold to find the blobs
	for (int level = threshold; level < 256; level++) {
		Mat mask = (M == level);
		if (countNonZero(mask) > 0) {
			vector<vector<Point> > contours;
			findContours(mask, contours, CV_RETR_EXTERNAL, 
				CV_CHAIN_APPROX_NONE);

			for (int i = 0; i < contours.size(); i++) {
				Mat blob = Mat::zeros(M.rows, M.cols, M.type());
				drawContours(blob, contours, i, Scalar(255), CV_FILLED);
				blobs.push_back(SparseMat(blob));
			}
		}
	}

	return blobs;
}

Segmentation::Segmentation(list<SparseMat> masks) {
	int idx = 0;
	for (auto m : masks) {
		this->_regions.push_back(Region(m, idx++));
	}

	this->setRegionsParent();
}

Segmentation::Segmentation(list<Point> points) {
	int idx = 0;
	for (auto p : points) {
		this->_regions.push_back(Region(p, idx++));
	}

	this->setRegionsParent();
}

Segmentation::Segmentation(ThematicMap M) {
	list<pair<SparseMat, int> > regions1 = M.enumerateRegions();
	int idx = 0;
	for (auto r1 : regions1) {
		_regions.push_back(Region(r1.first, idx++));
	}

	this->setRegionsParent();
}

Segmentation::Segmentation(const Segmentation& other)
	: _regions(other._regions),
	  _descriptions(other._descriptions),
	  _image(other._image),
	  _pixelated(other._pixelated) {

	// Update regions parent reference
	this->setRegionsParent();
}

Segmentation Segmentation::operator=(const Segmentation& rhs) {
	_regions = rhs._regions;
	_descriptions = rhs._descriptions;
	_image = rhs._image;
	_pixelated = rhs._pixelated;

	// Update regions parent reference
	this->setRegionsParent();

	return *this;
}

list<Region> Segmentation::getRegions() {
	return _regions;
}

list<SparseMat> Segmentation::getRegionMasks() {
	list<SparseMat> masks;
	for (auto& r : _regions) {
		masks.push_back(r.getMask());
	}

	return masks;
}

cv::Mat Segmentation::representation() {
	assert(_regions.size() > 0);

	Mat repr = Mat::zeros(getMapSize(), CV_8UC3);
	for (auto region : _regions) {
		Mat clone = densify(region.getMask());
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
}

cv::Size Segmentation::getMapSize() {
	assert(_image != NULL);
	return _image->size();
	return densify(_regions.front().getMask()).size();
}

int Segmentation::regionCount() {
	return _regions.size();
}

Mat Segmentation::getDescription(string descriptorID) {
	return _descriptions[descriptorID];
}

SparseMat Region::getMask() {
	return _mask;
}

Point Region::getPoint() {
	return _point;
}

Region::Region(SparseMat mask, int parentIdx) 
	: _mask(mask),
	  _parentIdx(parentIdx),
	  _representationMode(MASK) {}

Region::Region(Point point, int parentIdx) 
	: _point(point),
	  _parentIdx(parentIdx),
	  _representationMode(PIXEL) {}

Mat Region::getDescription(string descriptorID) {
	// Check if parent has been 
	assert (_parentSegmentation != NULL && 
		"Region construction error: no parent has been set");

	// Get features
	Mat features = _parentSegmentation->getDescription(descriptorID);
	return features.row(_parentIdx);
}

void Segmentation::setImage(RSImage *image) {
	_image = image;
}

RSImage *Segmentation::getImage() {
	if (_image == NULL) {
		FATAL_ERROR("Error: no image defined for this segmentation.");
	}

	return _image;
}

void Segmentation::describe(Descriptor *descriptor) {
	// Extract features
	Mat features;
	if (_image->getType() == VIS) {
		features = descriptor->describe(*((VISImage*) _image), this);
	}
	else {
		features = descriptor->describe(*((LWIRImage*) _image), this);
	}

	// Record features
	_mtx.lock();
	_descriptions[descriptor->getID()] = features;
	_mtx.unlock();
}

void Segmentation::setPixalated(bool v) {
	_pixelated = v;
}

bool Segmentation::isPixelated() {
	return _pixelated;
}

Segmentation Segmentation::cloneWithMask(SparseMat mask) {
	// Check eligibility
	if (!_pixelated) {
		FATAL_ERROR(
			"Error: cloning with mask only works on pixelated segmentations.");

		return Segmentation(); // Dummy
	}

	// Densify the mask
	Mat maskD = densify(mask);

	// Find the regions to clone
	list<int> clonedRegionsIdx;
	list<SparseMat> clonedRegions;
	int idx = 0;
	for (auto& r : _regions) {
		// Check if region is in mask.
		if (r.getRepresentationMode() == PIXEL) {
			Point P = r.getPoint();

			// Check if pixel is in mask
			if (maskD.at<unsigned char>(P) == 255) {
				clonedRegions.push_back(r.getMask());
				clonedRegionsIdx.push_back(idx);	
			}
		}
		else {
			if (countNonZero(densify(r.getMask()) & maskD) > 0) {
				clonedRegions.push_back(r.getMask());
				clonedRegionsIdx.push_back(idx);
			}
		}

		// Increase the index
		idx++;
	}

	// Create the feature matrixes (FMs) that correspond to the cloned regions
	int fmRows = clonedRegions.size();
	map<string, Mat> subFMs;
	for (auto d : _descriptions) {
		string id = d.first;
		Mat FM = d.second;

		Mat subFM(fmRows, FM.cols, FM.type());
		int i = 0;
		for (auto rIdx : clonedRegionsIdx) {
			FM.row(rIdx).copyTo(subFM.row(i++));
		}

		subFMs[id] = subFM;
	}

	// Create the Segmentation object
	Segmentation cloned(clonedRegions);
	cloned.setPixalated();
	cloned._descriptions = subFMs;

	return cloned;
}

Segmentation Segmentation::cloneWithIndexes(vector<int> indexes) {
	// Put the regions into a vector
	vector<Region> regionsV(_regions.begin(), _regions.end());

	// Clone the desired indexes
	list<SparseMat> clonedRegions;
	for (auto index : indexes) {
		clonedRegions.push_back(regionsV[index].getMask());
	}

	// Create the feature matrixes (FMs) that correspond to the cloned regions
	int fmRows = indexes.size();
	map<string, Mat> subFMs;
	for (auto d : _descriptions) {
		string id = d.first;
		Mat FM = d.second;

		Mat subFM(fmRows, FM.cols, FM.type());
		int i = 0;
		for (auto index : indexes) {
			FM.row(index).copyTo(subFM.row(i++));
		}

		subFMs[id] = subFM;
	}

	// Create the Segmentation object
	assert(_image != NULL);
	Segmentation cloned(clonedRegions);
	cloned._descriptions = subFMs;
	cloned.setImage(_image);
	cloned.setRegionsParent();

	return cloned;
}

std::vector<std::pair<Segmentation, Segmentation> > Segmentation::split(
	KFolder folder) {

	// Get the regions into a vector
	vector<Region> regions(_regions.begin(), _regions.end());

	// Create the split pairs (for now in index form)
	list<KFoldIndexes> foldsIdx = folder.makeFolds(regions.size());

	// Create the vector to hold the split pairs
	vector<pair<Segmentation, Segmentation> > folds;
	folds.reserve(foldsIdx.size());
	for (auto f : foldsIdx) {
		Segmentation T = this->cloneWithIndexes(f.first);
		Segmentation V = this->cloneWithIndexes(f.second);

		folds.push_back(make_pair(T, V));
	}

	return folds;
}

void Segmentation::setRegionsParent() {
	for (auto& r : _regions) {
		r.setParent(this);
	}
}

void Region::setParent(Segmentation *parent) {
	_parentSegmentation = parent;
}

void Segmentation::showDescriptionStats() {
	for (auto& p : _descriptions) {
		cerr << p.first << ": " << p.second.size() << endl;
	}
}

void Segmentation::upcountDescription(Descriptor *descriptor) {
	_regionsDescribed.inc(descriptor->getID());
}

Counter<string> Segmentation::getDescriptionCounts() {
	return _regionsDescribed;
}

RegionRepresentationMode Region::getRepresentationMode() {
	return _representationMode;
}
