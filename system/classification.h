#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/ml/ml.hpp>

#include "segmentation.h"

namespace classification {
	CvSVM *trainSVM(cv::Mat image, cv::Mat trainingMap, cv::Mat (*descriptor)(cv::Mat, const std::list<cv::Mat>));
	
}

#endif