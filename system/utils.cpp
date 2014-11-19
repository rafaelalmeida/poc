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

void showImage(Mat img, float scale, const char *winname, int delay) {
	if (winname == NULL) {
		winname = "win";
	}

	Mat display = img;
	if (scale != 1.0) {
		Mat aux;
		resize(display, aux, aux.size(), scale, scale, INTER_NEAREST);
		display = aux;
	}

	namedWindow(winname);
	imshow(winname, display);
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

cv::Mat mergeVISandLWIR(cv::Mat vis, cv::Mat lwirAvg) {
	Mat visGray(vis.rows, vis.cols, CV_8UC1);
	cvtColor(vis, visGray, CV_BGR2GRAY);

	return visGray;
}

void colorReduce(cv::Mat& image, int div) {    
    int nl = image.rows;                    // number of lines
    int nc = image.cols * image.channels(); // number of elements per line

    for (int j = 0; j < nl; j++)
    {
        // get the address of row j
        uchar* data = image.ptr<uchar>(j);

        for (int i = 0; i < nc; i++)
        {
            // process each pixel
            data[i] = data[i] / div * div + div / 2;
        }
    }
}

Mat formatImagesForPCA(const vector<Mat> &data)
{
    Mat dst(static_cast<int>(data.size()), data[0].rows*data[0].cols, CV_32F);
    for (unsigned int i = 0; i < data.size(); i++) {
        Mat image_row = data[i].clone().reshape(1,1);
        Mat row_i = dst.row(i);
        image_row.convertTo(row_i, CV_32F);
    }

    return dst;
}

cv::Mat densify(cv::SparseMat sm) {
	Mat dense = Mat::zeros(sm.size(0), sm.size(1), sm.type());
	sm.convertTo(dense, sm.type());
	
	return dense;
}

cv::Rect scaleROI(cv::Size original, cv::Size dst, cv::Rect roi) {
	double fx = (float) dst.width / original.width;
	double fy = (float) dst.height / original.height;

	Rect newRoi(fx * roi.x, fy * roi.y, fx * roi.width, fy * roi.height);

	return newRoi;
}

int translateInterpolationMode(InterpolationMode method) {
	if (method == NEAREST_NEIGHBOR) {
		return INTER_NEAREST;
	}
	else if (method == LINEAR) {
		return INTER_LINEAR;
	}
	else if (method == CUBIC) {
		return INTER_CUBIC;
	}
	else {
		// Should not get here
		assert(false && "Invalid interpolation mode.");
		return -1;
	}
}

double Stopwatch::read() { 
	return (double) ((_finish - _start).count()) / 1000000;
}

void Stopwatch::start() { 
	_start = chrono::system_clock::now(); 
}

void Stopwatch::stop() { 
	_finish = chrono::system_clock::now(); 
}
