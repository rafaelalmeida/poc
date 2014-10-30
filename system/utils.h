#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

void showImage(const char *winname, cv::Mat img, int delay=0);
cv::Mat floatImageTo8UC3Image(cv::Mat floatImage);

#endif