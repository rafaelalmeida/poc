#include "logging.h"

using namespace std;

Logger::Logger(const char *basePath) {
	time_t rawtime;
	struct tm * timeinfo;

	time (&rawtime);
	timeinfo = localtime (&rawtime);

	char buffer[256];
	strftime(buffer, 256, "%F %T", timeinfo);

	char fullPath[256];
	strcpy(fullPath, basePath);
	strcat(fullPath, "/");
	strcat(fullPath, buffer);

	mkdir(fullPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	_fullPath = string(fullPath);
}

void Logger::saveImage(const char *name, cv::Mat img) {
	string s = _fullPath + "/" + name + ".png";
	imwrite(s.c_str(), img);
}
