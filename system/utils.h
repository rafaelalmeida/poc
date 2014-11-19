#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <cstdint>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>
#include <time.h>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "common.h"
#include "segmentation.h"

extern "C" {
	#include "descriptors/image.h"
	#include "descriptors/cimage.h"
}

using namespace std;

extern bool verbose;

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

void log(const char *msg);

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

// Helper class to count things
template <typename T>
class Counter {
	std::map<T, int> _map;

	public:
		void inc(T item);
		T top();
		std::map<T, int> getCounts();
		int getCount(T item);
};

template <typename T>
void Counter<T>::inc(T item) {
	if (_map.count(item) == 0) {
		_map[item] = 0;
	}

	_map[item]++;
}

template <typename T> 
T Counter<T>::top() {
	int maxVal = 0;
	T cMax;

	for (typename map<T, int>::iterator i = _map.begin(); i != _map.end(); ++i) {
		if (i->second > maxVal) {
			maxVal = i->second;
			cMax = i->first;
		}
	}

	return cMax;
}

template <typename T>
std::map<T, int> Counter<T>::getCounts() {
	return _map;
}

template <typename T>
int Counter<T>::getCount(T item) {
	auto val = _map.find(item);
	if (val == _map.end()) {
		return 0;
	}
	else {
		return val->second;
	}
}

#endif