#include "utils.h"

using namespace cv;
using namespace std;

Mat floatImageTo8UC3Image(Mat floatImage) {
	Mat M = floatImage;

	double max, min;
	minMaxIdx(M, &min, &max);
	Mat ret(M.rows, M.cols, CV_8UC3);
	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			uchar val = (uchar) (255 * M.at<float>(row, col) / max);
			ret.at<Vec3b>(row, col) = Vec3b(val, val, val);
		}
	}

	return ret;
}

void showImage(const char *winname, Mat img, int delay) {
	namedWindow(winname);
	imshow(winname, img);
	waitKey(delay);
}

Image *matToRawGray(cv::Mat gray) {
	assert(gray.type() == CV_8UC1);

	Image *img = CreateImage(gray.cols, gray.rows);
	for (int row = 0; row < gray.rows; row++) {
		for (int col = 0; col < gray.cols; col++) {
			img->val[row*gray.cols + col] = (int) gray.at<unsigned char>(row, col);
		}
	}

	return img;
}

CImage *matToRawColor(cv::Mat color) {
	assert(color.type() == CV_8UC3);

	Mat channels[3];
	split(color, channels);

	CImage *cimg = CreateCImage(channels[0].cols, channels[0].rows);
	cimg->C[0] = matToRawGray(channels[2]);
	cimg->C[1] = matToRawGray(channels[1]);
	cimg->C[2] = matToRawGray(channels[0]);

	return cimg;
}

cv::Mat onesLike(cv::Mat M) {
	return 255 * Mat::ones(M.rows, M.cols, CV_8UC1);
}

cv::Mat makeLandCoverMap(cv::Mat labels) {
	assert(labels.type() == CV_8UC1);

	Mat labelsClone = labels.clone();
	Mat map(labels.rows, labels.cols, CV_8UC3);

	list<Mat> regions = segmentation::makeSegmentMasksFromPosterizedImage(labels);
	for (list<Mat>::iterator it = regions.begin(); it != regions.end(); ++it) {
		vector<Mat> contours;
		findContours(*it, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
		float label = segmentation::getSegmentLabel(labels, *it);

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

cv::Mat blend(cv::Mat M1, cv::Mat M2) {
	assert(M1.rows == M2.rows && M1.cols == M2.cols);
	assert(M1.type() == CV_8UC3 && M2.type() == CV_8UC3);

	cv::Mat dst(M1.rows, M1.cols, CV_8UC3);
	addWeighted(M1, 0.5, M2, 0.5, 0.0, dst);

	return dst;
}

void log(const char *msg) {
	if (verbose) {
		cerr << msg << endl;
	}
}
