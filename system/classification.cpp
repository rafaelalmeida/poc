#include "classification.h"

using namespace cv;
using namespace std;

using namespace segmentation;

CvSVM *classification::trainSVM(cv::Mat image, cv::Mat trainingMap, cv::Mat (*descriptor)(cv::Mat, const std::list<cv::Mat>)) {
	// Extract regions from training map
	list<Mat> masks = segmentation::makeSegmentMasksFromPosterizedImage(trainingMap);

	// Recover region labels
	list<float> labels;
	list<Mat> validSegments;
	for (list<Mat>::iterator it = masks.begin(); it != masks.end(); ++it) {
		float label = segmentation::getSegmentLabel(trainingMap, *it);
		if (label != 0) { // Disconsider unclassified regions
			labels.push_back(label);
			validSegments.push_back(*it);
		}
	}

	// Create labels training matrix
	Mat labelsMat(labels.size(), 1, CV_32FC1);
	float c = 0;
	for (list<float>::iterator it = labels.begin(); it != labels.end(); ++it) {
		labelsMat.at<float>(c) = *it;
		c++;
	}

	// Compute features
	Mat features = descriptor(image, validSegments);

	// Instantiate and train SVM
	CvSVM *SVM = new CvSVM();
	CvSVMParams *params = new CvSVMParams();
    params->svm_type    = CvSVM::C_SVC;
    params->kernel_type = CvSVM::LINEAR;
    params->term_crit   = cvTermCriteria(CV_TERMCRIT_ITER, 100, 1e-6);
    SVM->train(features, labelsMat, Mat(), Mat(), *params);

	return SVM;
}

cv::Mat classification::predict(cv::Mat image, Segmentation segmentation, CvSVM *classifier) {
	Mat map(image.rows, image.cols, CV_8UC1);

	for (auto&& s : segmentation.getSegments()) {
		Mat sample = description_vis::GCH(image, s);
		int theClass = (int) classifier->predict(sample);
		map = map | theClass * (s / 255);
	}

	return map;
}