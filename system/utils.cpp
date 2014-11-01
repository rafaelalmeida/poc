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


