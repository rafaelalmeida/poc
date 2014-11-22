#include "classification.h"

using namespace cv;
using namespace std;

using namespace classification;
using namespace segmentation;

Classifier::Classifier(ClassifierEngine engine, VISImage *vis, 
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

void Classifier::train(Mat labels, Segmentation S) {
	// Get features
	Mat features = S.getDescription(_descriptor->getID());

	// Train the correct classifier
	_swatchTraining.start();
	if (_engine == SVM) {
		CvSVMParams params;
		params.svm_type    = CvSVM::C_SVC;
		params.kernel_type = CvSVM::RBF;
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
		CvDTreeParams params; // Initialize with default parameters
		params.cv_folds = 1; // Disable k-fold cross-validation
		params.min_sample_count = DTREE_MIN_SAMPLE_SIZE;

		// Include the class probabilities, if available
		vector<float> *prob = NULL;
		float probF[CLASS_COUNT];

		if (_trainingMap != NULL) {
			prob = _trainingMap->classProbabilities();

			int idx = 0;
			for (auto p : *prob) {
				probF[idx++] = p;
			}

			params.priors = probF;
		}

		_dtree.train(features, CV_ROW_SAMPLE, labels, Mat(), Mat(), Mat(), 
			Mat(), params);

		// Free class probabilities vector
		if (prob != NULL) {
			delete prob;
		}
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
		layers.at<int>(1, 0) = MLP_HIDDEN_LAYER_SIZE;
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

	_swatchTraining.stop();
}

Mat Classifier::classify(Region region) {
	// Get features
	Mat features = region.getDescription(_descriptor->getID());

	// Predict the class
	_swatchClassification.start();

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

	_swatchClassification.stop();

	// Fill the matrix
	Mat mask = densify(region.getMask());
	Mat classification = Mat::zeros(mask.size(), CV_8UC1);
	classification += theClass * (mask / 255);

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

ImageType Classifier::getType() {
	return _type;
}

double Classifier::getDescriptionTime() {
	return _swatchDescription.read();
}

double Classifier::getTrainingTime() {
	return _swatchTraining.read();
}

double Classifier::getClassificationTime() {
	return _swatchClassification.read();
}

void Classifier::setTrainingMap(ThematicMap *map) {
	_trainingMap = map;
}
