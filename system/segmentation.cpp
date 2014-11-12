#include "segmentation.h"

using namespace cv;
using namespace std;

using namespace segmentation;

Segmentation segmentation::segmentVISGrid(cv::Mat M, int tileSize) {
	list<SparseMat> segments;

	int regionsPerLine = (M.cols / tileSize) + (M.cols % tileSize);
	int regionsPerColumn = (M.rows / tileSize) + (M.rows % tileSize);

	Mat gray(M.size(), CV_8UC1);
	cvtColor(M, gray, CV_BGR2GRAY);

	Mat bin(gray.size(), CV_8UC1);
	threshold(gray, bin, 1, 255, CV_THRESH_BINARY);

	for (int i = 0; i < regionsPerLine; i++) {
		for (int j = 0; j < regionsPerColumn; j++) {
			Mat segment = Mat::zeros(M.size(), CV_8UC1);
			Rect region(i*tileSize, j*tileSize, tileSize, tileSize);

			rectangle(segment, region, Scalar(255), CV_FILLED);

			Mat filteredSegment = segment & bin;
			if (countNonZero(filteredSegment) > 0) {
				segments.push_back(SparseMat(filteredSegment));
			}
		}
	}

	return Segmentation(segments);
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
			// Report progress
			if ((y*clone.rows + x) % 1000 == 0) {
				int progress = 100.0 * ((y*clone.cols + x)) / 
					(clone.rows * clone.cols);

				cerr << "finding blobs... " << progress << 
					"%...     " << "\r" << flush;
			}

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

	cerr << "finding blobs... done" << endl;

	return segments;
}

Segmentation::Segmentation(list<SparseMat> masks) {
	this->_masks = masks;
}

list<SparseMat> Segmentation::getSegments() {
	return _masks;
}

cv::Mat Segmentation::representation() {
	assert(_masks.size() > 0);

	Mat repr = Mat::zeros(getMapSize(), CV_8UC3);
	for (auto mask : _masks) {
		Mat clone = densify(mask);
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
}

cv::Size Segmentation::getMapSize() {
	return densify(_masks.front()).size();
}

vector<SparseMat> segmentation::pixelSegmentation(Mat image) {
	vector<SparseMat> segments;
	segments.reserve(image.rows * image.cols);

	for (int row = 0; row < image.rows; row++) {
		for (int col = 0; image.cols; col++) {
			Mat mask(image.size(), CV_8UC1);
			mask.at<unsigned char>(row, col) = 255;

			SparseMat sparseMask(mask);

			segments.push_back(sparseMask);
		}
	}

	return segments;
}

Segmentation Segmentation::pixelize() {
	Mat fullMask = Mat::zeros(this->getMapSize(), CV_8UC1);
	for (auto m : _masks) {
		fullMask += densify(m);
	}

	list<SparseMat> pixelatedSegments;

	for (int row = 0; row < fullMask.rows; row++) {
		for (int col = 0; col < fullMask.cols; col++) {
			if (fullMask.at<unsigned char>(row, col) != 0) {
				Mat mask = Mat::zeros(fullMask.size(), CV_8UC1);
				mask.at<unsigned char>(row, col) = 255;

				pixelatedSegments.push_back(SparseMat(mask));
			}
		}
	}

	return Segmentation(pixelatedSegments);
}

int Segmentation::segmentCount() {
	return _masks.size();
}