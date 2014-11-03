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
