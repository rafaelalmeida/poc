#include "segmentation.h"

using namespace cv;
using namespace std;

list<Mat> segmentation::segmentLWIRMeanShift(Mat M) {
	// Convert to 8UC3 image
	Mat img8uc3 = floatImageTo8UC3Image(M);

	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(img8uc3, res, 9, 9, 0);

	return makeSegmentMasksFromPosterizedImage(res);
}

list<Mat> segmentation::segmentLWIRCanny(Mat M) {
	list<Mat> ret;

	Mat conv = floatImageTo8UC3Image(M);

	Mat edges;
	int thres = 10;
	Canny(conv, edges, thres, 3*thres, 3);
	showImage("img", edges);
	
	return ret;
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
