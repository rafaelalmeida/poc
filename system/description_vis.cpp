#include "description_vis.h"

using namespace cv;

cv::Mat description_vis::GCH(const cv::Mat image, const std::list<cv::Mat> masks) {
	const int HIST_BINS = 64;

	CImage *cimg = matToRawColor(image);

	Mat samples(masks.size(), HIST_BINS, CV_32FC1);

	int currentMask = 0;
	for (list<Mat>::const_iterator it = masks.begin(); it != masks.end(); it++) {
		Image *mask = matToRawGray(*it);

		Histogram *hist = GCH(cimg, mask);
		assert(hist->n == HIST_BINS);

		for (int i = 0; i < hist->n; i++) {
			samples.at<float>(currentMask, i) = (float) hist->v[i];
		}

		DestroyHistogram(&hist);
		DestroyImage(&mask);

		currentMask++;
	}

	DestroyCImage(&cimg);

	return samples;
}

cv::Mat description_vis::GCH(const cv::Mat image, const cv::Mat mask) {
	list<Mat> masks;
	masks.push_back(mask);
	return GCH(image, masks);
}