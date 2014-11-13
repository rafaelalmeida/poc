#ifndef STAT_H
#define STAT_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "utils.h"

#define MAX_CATEGORIES 256;

namespace statistics {
	double agreement(cv::Mat A, cv::Mat B);
	double kappa(cv::Mat A, cv::Mat B);
	
	std::vector<float> moments(std::vector<float> samples, int maxOrder);
	
	std::vector<float> summary(std::vector<float> samples);
}

#endif