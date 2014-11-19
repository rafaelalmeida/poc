#include "classification.h"

using namespace cv;
using namespace std;

using namespace classification;
using namespace segmentation;

Classifier::Classifier(ClassifierEngine engine, cv::Mat vis, 
	Descriptor *descriptor)
	: _engine(engine),
	  _vis(vis),
	  _descriptor(descriptor),
	  _type(VIS) {
}

Classifier::Classifier(ClassifierEngine engine, LWIRImage *lwir, 
	Descriptor *descriptor)
	: _engine(engine),
	  _lwir(lwir),
	  _descriptor(descriptor),
	  _type(LWIR) {
}

void Classifier::train(Mat labels, Segmentation trainingSegments) {
	// Get valid segments
	list<SparseMat> validSegments = trainingSegments.getSegments();

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

		_svm.train(features, labels, Mat(), Mat(), params);
	}
	else if (_engine == NBC) {
		_nbc.train(features, labels, Mat(), Mat(), false);
	}
	else if (_engine == KNN) {
		_knn.train(features, labels);
	}
	else if (_engine == DTREE) {
		_dtree.train(features, CV_ROW_SAMPLE, labels);
	}
	else if (_engine == GBT) {
		_gbtrees.train(features, CV_ROW_SAMPLE, labels);
	}
	else if (_engine == RTREES) {
		_rtrees.train(features, CV_ROW_SAMPLE, labels);
	}
	else if (_engine == ERTREES) {
		_ertrees.train(features, CV_ROW_SAMPLE, labels);
	}
	else if (_engine == MLP) {
		Mat layers(3, 1, CV_32S);
		layers.at<int>(0, 0) = features.cols;
		layers.at<int>(1, 0) = 16;
		layers.at<int>(2, 0) = CLASS_COUNT;

		// Since MLPs can't handle categorical data, we encode it as a vector v
		// of CLASS_COUNT elements. If the class is c, then v(c) = 1 and all
		// other elements of v are 0.
		Mat encodedLabels = Mat::zeros(labels.rows, CLASS_COUNT, CV_32F);
		for (int i = 0; i < labels.rows; i++) {
			int c = (int) (labels.at<float>(i, 0) - 1); // Classes are 1-index
			encodedLabels.at<float>(i, c) = 1.0;
		}

		_mlp.create(layers, CvANN_MLP::SIGMOID_SYM);
		_mlp.train(features, encodedLabels, Mat());
	}
	else {
		assert(false && "Unknown classifier engine");
	}
}

Mat Classifier::classify(cv::SparseMat mask) {
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
	else if (_engine == NBC) {
		theClass = _nbc.predict(features);
	}
	else if (_engine == KNN) {
		theClass = _knn.find_nearest(features, KNN_K);
	}
	else if (_engine == DTREE) {
		theClass = _dtree.predict(features)->value;
	}
	else if (_engine == GBT) {
		theClass = _gbtrees.predict(features);
	}
	else if (_engine == RTREES) {
		theClass = _rtrees.predict(features);
	}
	else if (_engine == ERTREES) {
		theClass = _ertrees.predict(features);
	}
	else if (_engine == MLP) {
		Mat classification;
		_mlp.predict(features, classification);

		// Find the class with highest weight
		int maxIndex[2];
		minMaxIdx(classification, NULL, NULL, NULL, maxIndex);

		// The predicted class is the 1-based (since we don't include the 
		// 0 (unclassified) class) index of the output layer element with 
		// maximal value.
		theClass = (float) (maxIndex[1] + 1);
	}
	else {
		assert(false && "Unknown classifier engine");
	}

	// Fill the matrix
	Mat denseMask = densify(mask);
	Mat classification = Mat::zeros(denseMask.size(), CV_8UC1);
	classification += theClass * (denseMask / 255);

	return classification;
}

string Classifier::getID() {
	string id;

	// Place the classifier identification
	if (_engine == SVM) {
		id += "SVM-";
	}
	else if (_engine == NBC) {
		id += "NBC-";
	}
	else if (_engine == KNN) {
		id += "KNN-";
	}
	else if (_engine == DTREE) {
		id += "DTREE-";
	}
	else if (_engine == GBT) {
		id += "GBT-";
	}
	else if (_engine == RTREES) {
		id += "RTREES-";
	}
	else if (_engine == ERTREES) {
		id += "ERTREES-";
	}
	else if (_engine == MLP) {
		id += "MLP-";
	}

	// Place the descriptor identification
	id += _descriptor->getID();

	return id;
}

ClassifierType Classifier::getType() {
	return _type;
}
