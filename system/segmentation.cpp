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

	return Segmentation(makeSegmentMasksFromPosterizedImage(res));
}

Segmentation segmentation::segmentLWIRCanny(Mat M) {
	list<Mat> ret;

	Mat conv = floatImageTo8UC3Image(M);

	Mat edges;
	int thres = 10;
	Canny(conv, edges, thres, 3*thres, 3);
	
	return Segmentation(ret);
}

Segmentation segmentation::segmentVISGrid(cv::Mat M) {
	const int GRID_SIZE = 32;

	list<Mat> segments;

	int regionsPerLine = M.cols / GRID_SIZE;
	int regionsPerColumn = M.rows / GRID_SIZE;

	for (int i = 0; i < regionsPerLine; i++) {
		for (int j = 0; j < regionsPerColumn; j++) {
			Mat segment(M.rows, M.cols, CV_8UC1);
			Rect region(i*GRID_SIZE, j*GRID_SIZE, GRID_SIZE, GRID_SIZE);

			rectangle(segment, region, Scalar(255), CV_FILLED);
			segments.push_back(segment);
		}
	}

	return Segmentation(segments);
}

float segmentation::getSegmentLabel(Mat classificationMap, Mat mask) {
	assert(classificationMap.type() == CV_8UC1);
	assert(mask.type() == CV_8UC1);
	assert((classificationMap.rows == mask.rows) && (classificationMap.cols = mask.cols));

	typedef unsigned char uchar;

	Counter<uchar> counter;

	for (int row = 0; row < classificationMap.rows; row++) {
		for (int col = 0; col < classificationMap.cols; col++) {
			if (mask.at<uchar>(row, col)) {
				uchar val = classificationMap.at<uchar>(row, col);
				counter.inc(val);
			}
		}
	}

	return (float) counter.top();
}

std::list<cv::Mat> segmentation::makeSegmentMasksFromPosterizedImage(cv::Mat posterized) {
	list<Mat> segments;

	Mat clone = posterized.clone();
	Mat mask(posterized.rows+2, posterized.cols+2, CV_8UC1, Scalar::all(0));

	for (int y = 0; y < clone.rows; y++)
	{
		for (int x = 0; x < clone.cols; x++)
		{
			if (mask.at<uchar>(y+1, x+1) == 0)
			{
				Mat previousMask = mask.clone();

				floodFill(clone, mask, Point(x,y), Scalar(0, 0, 0), 0, 
				          Scalar::all(1), Scalar::all(1));

				Mat difference = previousMask ^ mask;
				Mat segmentMask = 255 * difference(Range(1, difference.rows - 1),
				                             Range(1, difference.cols - 1));

				segments.push_back(segmentMask);
			}
		}
	}

	return segments;
}

cv::Mat segmentation::representSegmentation(std::list<cv::Mat> masks) {
	assert(masks.size() > 0);
	Mat M = *masks.begin();

	Mat repr(M.rows, M.cols, CV_8UC3);
	for (list<Mat>::iterator it = masks.begin(); it != masks.end(); ++it) {
		Mat clone = it->clone();
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
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

	Mat repr(M.rows, M.cols, CV_8UC3);
	for (list<Mat>::iterator it = _masks.begin(); it != _masks.end(); ++it) {
		Mat clone = it->clone();
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
}