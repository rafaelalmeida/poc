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

Classifier::Classifier(ClassifierEngine engine, cv::Mat vis, Descriptor *descriptor)
    : _engine(engine),
      _vis(vis),
      _descriptor(descriptor),
      _type(VIS) {
}

Classifier::Classifier(ClassifierEngine engine, LWIRImage *lwir, Descriptor *descriptor)
    : _engine(engine),
      _lwir(lwir),
      _descriptor(descriptor),
      _type(LWIR) {
}

Classifier::~Classifier() {
	delete _descriptor;
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
	else if (_type == LWIR) {
		features = _descriptor->describe(*_lwir, validSegments);
	}
	else {
		assert(false && "Unknown classifier type");
	}

	// Train the correct classifier
	if (_engine == SVM) {
		CvSVMParams params;
	    params.svm_type    = CvSVM::C_SVC;
	    params.kernel_type = CvSVM::LINEAR;
	    params.term_crit   = cvTermCriteria(CV_TERMCRIT_ITER, 100, 1e-6);

	    _svm.train(features, labelsMat, Mat(), Mat(), params);
	}
	else {
		assert(false && "Unknown classifier engine");
	}
}

Mat Classifier::classify(cv::Mat mask) {
	// Compute features
	Mat features;
	if (_type == VIS) {
		features = _descriptor->describe(_vis, mask);
	}
	else if (_type == LWIR) {
		features = _descriptor->describe(*_lwir, mask);
	}
	else {
		assert(false && "Unknown classifier type");
	}

	// Predict the class
	float theClass;
	if (_engine == SVM) {
		theClass = _svm.predict(features);
	}
	else {
		assert(false && "Unknown classifier engine");
	}

	// Fill the matrix
	Mat classification = Mat::zeros(mask.size(), CV_8UC1);
	classification += theClass * (mask / 255);

	return classification;
}
