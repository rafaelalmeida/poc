#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& s, CoverMap training) 
    : _consensusType(t),
	  _segmentation(s),
	  _training(training) {

	assert(training.asMat().type() == CV_8UC1);
}

Ensemble::~Ensemble() {
	for (auto c : classifiers) {
		delete c;
	}
}

void Ensemble::addClassifier(Classifier *c) {
	classifiers.push_back(c);
}

void Ensemble::setLogger(Logger *logger) {
	_logger = logger;
}

void Ensemble::train() {
	// Extract regions from training map
	list<Mat> masks = segmentation::getColorBlobs(_training.asMat());

	// Recover region labels
	list<float> labels;
	list<Mat> validSegments;
	for (auto mask : masks) {
		float label = _training.getRegionClass(mask);
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

	int i = 1;
	int n = classifiers.size();

	for (auto c : classifiers) {
		cerr << "training classifier " << i << " of " << n << "     \r" << flush;

		c->train(labelsMat, Segmentation(validSegments));
		i++;
	}

	cerr << "training classifiers... done     " << endl;
}

CoverMap Ensemble::classify() {
	// Some initialization
	Size mapSize = _segmentation.getMapSize();
	list<Mat> classifiedSegments;

	// Saves individual classifications for debugging
	_classifications.resize(classifiers.size());
	for (auto& c : _classifications) {
		c = Mat::zeros(mapSize, CV_8UC1);
	}

	// Runs all classifiers for each segment
	int c = 1;
	int n = _segmentation.getSegments().size();

	for (auto mask : _segmentation.getSegments()) {
		// Show progress
		cerr << "classifying segment " << c << " of " << n << "     \r" << flush;

		// Gets the opinion of all classifiers for this segment
		vector<Mat> opinions;
		opinions.reserve(classifiers.size());

		int i = 0;
		for (auto c : classifiers) {
			Mat classification = c->classify(mask);
			opinions.push_back(classification);

			// Save this segment individually for debugging
			_classifications[i] += classification;

			i++;
		}

		// Initializes the consensus matrix for this segment
		Mat consensus = Mat::zeros(mapSize, CV_8UC1);

		if (_consensusType == MAJORITY_VOTING) {
			consensus += opinions[0];
		}
		else {
			assert(false && "Unknown consensus type");
		}

		classifiedSegments.push_back(consensus);
		c++;
	}

	cerr << "classyfing segments... done   " << endl;

	// Builds the map from the segments
	Mat theMap = Mat::zeros(mapSize, CV_8UC1);
	for (auto seg : classifiedSegments) {
		theMap += seg;
	}

	return CoverMap(theMap);
}

vector<Mat> Ensemble::individualClassifications() {
	return _classifications;
}
