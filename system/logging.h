#ifndef LOGGING_H
#define LOGGING_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

class Logger {
	std::string _fullPath;

	// Constructors
	public:
		Logger(const char *basePath);
		void saveImage(const char *name, cv::Mat img);
		std::ofstream* makeFile(const char *name);
		void saveArguments(int argc, char **argv);
};

#endif
