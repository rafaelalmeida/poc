#include "description_vis.h"

using namespace cv;

cv::Mat description_vis::GCH(const cv::Mat image, const cv::Mat mask) {
	const int HIST_BINS = 64;

	CImage *cimg = matToRawColor(image);
	Image *maskRaw = matToRawGray(mask);

	Mat samples(1, HIST_BINS, CV_32FC1);

	Histogram *hist = GCH(cimg, maskRaw);
	assert(hist->n == HIST_BINS);

	for (int i = 0; i < hist->n; i++) {
		samples.at<float>(0, i) = (float) hist->v[i];
	}

	DestroyHistogram(&hist);

	DestroyCImage(&cimg);
	DestroyImage(&maskRaw);

	return samples;
}
