#include "segmentation.h"

using namespace cv;
using namespace std;

list<Mat> segmentation::segmentLWIRMeanShift(Mat M) {
	// Convert to 8UC3 image
	Mat img8uc3 = floatImageTo8UC3Image(M);

	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(img8uc3, res, 9, 9, 0);

	// Separates the segments
	list<Mat> segments;
	Mat mask(res.rows+2, res.cols+2, CV_8UC1, Scalar::all(0));
	for (int y = 0; y < res.rows; y++)
	{
		for (int x = 0; x < res.cols; x++)
		{
			if (mask.at<uchar>(y+1, x+1) == 0)
			{
				Mat previousMask = mask.clone();

				floodFill(res, mask, Point(x,y), Scalar(0, 0, 0), 0, Scalar::all(1), Scalar::all(1));

				Mat difference = previousMask ^ mask;
				segments.push_back(difference);
			}
		}
	}

	return segments;
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
	typedef unsigned char uchar;

	Mat roi = classificationMap & mask;
	Counter<uchar> counter;

	for (int row = 0; row < roi.rows; row++) {
		for (int col = 0; col < roi.cols; col++) {
			counter.inc(roi.at<uchar>(row, col));
		}
	}

	return (float) counter.top();
}
