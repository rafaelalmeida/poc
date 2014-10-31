#ifndef STAT_H
#define STAT_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#define MAX_CATEGORIES 256;

namespace statistics {
	std::vector<float> summary(std::vector<float> samples);
	double kappa(cv::Mat A, cv::Mat B);
}

#endif