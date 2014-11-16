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

Classifier::~Classifier() {
	delete _descriptor;
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

	// Place the descriptor identification
	id += _descriptor->getID();

	return id;
}

ClassifierType Classifier::getType() {
	return _type;
}
