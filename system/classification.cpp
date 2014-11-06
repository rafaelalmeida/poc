#include "classification.h"

using namespace cv;
using namespace std;

using namespace classification;
using namespace segmentation;

cv::Mat classification::predict(cv::Mat image, Segmentation segmentation, CvSVM *classifier) {
	Mat map(image.rows, image.cols, CV_8UC1);

	for (auto&& s : segmentation.getSegments()) {
		Mat sample = description_vis::GCH(image, s);
		int theClass = (int) classifier->predict(sample);
		map = map | theClass * (s / 255);
	}

	return map;
}

Classifier::Classifier(ClassifierEngine engine, cv::Mat vis, Descriptor *descriptor) {
	this->_type = VIS;
	this->_vis = vis;
	this->_descriptor = descriptor;
}

Classifier::Classifier(ClassifierEngine engine, LWIRImage *lwir,  Descriptor *descriptor) {
	this->_type = LWIR;
	this->_lwir = lwir;
	this->_descriptor = descriptor;
}

void Classifier::train(CoverMap training) {
	// Extract regions from training map
	list<Mat> masks = segmentation::getColorBlobs(training.asMat());

	// Recover region labels
	list<float> labels;
	list<Mat> validSegments;
	for (auto mask : masks) {
		float label = training.getRegionClass(mask);
		if (label != 0) { // Disconsider unclassified regions
			labels.push_back(label);
			validSegments.push_back(mask);
		}
	}

	// Create labels training matrix
	Mat labelsMat(labels.size(), 1, CV_32FC1);
	int c = 0;
	for (auto label : labels) {
		labelsMat.at<float>(c) = label;
		c++;
	}

	// Compute features
	Mat features;
	if (_type == VIS) {
		features = _descriptor->describe(_vis, validSegments);
	}
	else {
		features = _descriptor->describe(*_lwir, validSegments);
	}

	// Instantiate and train SVM
	CvSVM *SVM = new CvSVM();
	CvSVMParams *params = new CvSVMParams();
    params->svm_type    = CvSVM::C_SVC;
    params->kernel_type = CvSVM::LINEAR;
    params->term_crit   = cvTermCriteria(CV_TERMCRIT_ITER, 100, 1e-6);

    SVM->train(features, labelsMat, Mat(), Mat(), *params);
}

Mat Classifier::classify() {
	Mat ret;
	return ret;
}
