#ifndef LOGGING_H
#define LOGGING_H

#include <cstdint>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <stdio.h>
#include <sys/stat.h>

#include <opencv2/opencv.hpp>

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
