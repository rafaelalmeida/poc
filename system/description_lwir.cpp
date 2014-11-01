#include "description_lwir.h"

using namespace std;
using namespace cv;

Mat description_lwir::spectralSignature(vector<Mat> lwirChannels, Mat mask) {
	int bands = lwirChannels.size();
	Mat sig(bands, 1, CV_32FC1);

	int band = 1;
	for (vector<Mat>::iterator i = lwirChannels.begin(); i != lwirChannels.end(); ++i) {
		Mat roi = *i & mask;
		sig.at<float>(band, 1) = mean(roi)[0];
		band++;
	}

	return sig;
}