#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <cstdint>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>
#include <time.h>

#include <opencv2/opencv.hpp>

#include "common.h"
#include "segmentation.h"

extern "C" {
	#include "descriptors/image.h"
	#include "descriptors/cimage.h"
}

using namespace std;

cv::Mat blend(cv::Mat M1, cv::Mat M2);
cv::Mat floatImageTo8UC3Image(cv::Mat floatImage);
cv::Mat formatImagesForPCA(const std::vector<cv::Mat> &data);
cv::Mat makeLandThematicMap(cv::Mat labels);
cv::Mat onesLike(cv::Mat M);
void colorReduce(cv::Mat& image, int div=64);
void showImage(cv::Mat img, float scale=1.0, const char *winname=NULL, int delay=0);
cv::Mat densify(cv::SparseMat sm);
cv::Rect scaleROI(cv::Size original, cv::Size dst, cv::Rect roi);
int translateInterpolationMode(InterpolationMode method);

cv::Mat mergeVISandLWIR(cv::Mat vis, cv::Mat lwirAvg);

Image *matToRawGray(cv::Mat gray);
CImage *matToRawColor(cv::Mat color);

inline int streq(const char *a, const char *b) {
	return strcmp(a, b) == 0;
}

// Helper class to measure time
class Stopwatch {
	private:
		std::chrono::time_point<std::chrono::system_clock> _start;
		std::chrono::time_point<std::chrono::system_clock> _finish;

	public:
		double read();
		void start();
		void stop();
};

#endif