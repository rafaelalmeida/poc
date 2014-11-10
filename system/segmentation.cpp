#include "segmentation.h"

using namespace cv;
using namespace std;

using namespace segmentation;

Segmentation segmentation::segmentLWIRMeanShift(Mat M) {
	// Convert to 8UC3 image
	Mat img8uc3 = floatImageTo8UC3Image(M);

	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(img8uc3, res, 9, 9, 0);

	return Segmentation(getColorBlobs(res));
}

Segmentation segmentation::segmentVISMeanShift(Mat M) {
	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(M, res, 9, 9, 0);

	return Segmentation(getColorBlobs(res));
}

Segmentation segmentation::segmentLWIRCanny(Mat M) {
	list<Mat> ret;

	Mat conv = floatImageTo8UC3Image(M);

	Mat edges;
	int thres = 10;
	Canny(conv, edges, thres, 3*thres, 5);
	
	return Segmentation(ret);
}

Segmentation segmentation::segmentVISCanny(Mat M) {
	list<Mat> ret;

	Mat gray(M.size(), CV_8UC1);
	cvtColor(M, gray, CV_BGR2GRAY);

	Mat edges;
	int thres = 50;
	Canny(gray, edges, thres, 3*thres, 3);

	vector<vector<Point> > contours;
	findContours(edges, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

	Mat cont = Mat::zeros(edges.size(), edges.type());
	drawContours(cont, contours, -1, Scalar(255), CV_FILLED);
	
	showImage(cont, 0.2);

	return Segmentation(ret);
}

Segmentation segmentation::segmentVISGrid(cv::Mat M, int tileSize) {
	list<Mat> segments;

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
				segments.push_back(filteredSegment);
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

std::list<cv::Mat> segmentation::getColorBlobs(cv::Mat posterized) {
	list<Mat> segments;

	Mat clone = posterized.clone();
	Mat mask(posterized.rows+2, posterized.cols+2, CV_8UC1, Scalar::all(0));

	for (int y = 0; y < clone.rows; y++) {
		for (int x = 0; x < clone.cols; x++) {
			if ((y*clone.rows + x) % 1000 == 0) {
				int progress = 100.0 * ((y*clone.cols + x)) / (clone.rows * clone.cols);
				cerr << "finding blobs... " << progress << "%...     " << "\r" << flush;
			}

			if (mask.at<uchar>(y+1, x+1) == 0) {
				Mat previousMask = mask.clone();

				floodFill(clone, mask, Point(x,y), Scalar(0, 0, 0), 0, 
				          Scalar::all(0), Scalar::all(0));

				Mat difference = previousMask ^ mask;
				Mat segmentMask = 255 * difference(Range(1, difference.rows - 1),
				                                   Range(1, difference.cols - 1));

				segments.push_back(segmentMask);
			}
		}
	}

	cerr << "finding blobs... done" << endl;

	return segments;
}

Segmentation::Segmentation(list<Mat> masks) {
	this->_masks = masks;
}

list<Mat> Segmentation::getSegments() {
	return _masks;
}

cv::Mat Segmentation::representation() {
	assert(_masks.size() > 0);
	Mat M = *_masks.begin();

	Mat repr = Mat::zeros(M.size(), CV_8UC3);
	for (auto mask : _masks) {
		Mat clone = mask.clone();
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
}

cv::Size Segmentation::getMapSize() {
	return _masks.begin()->size();
}
