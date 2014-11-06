#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& s, CoverMap training) :
	_consensusType(t),
	_segmentation(s),
	_training(training) {}

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
	for (auto c : classifiers) {
		c->train(_training);
	}
}

CoverMap Ensemble::classify() {
	Size mapSize = _segmentation.getMapSize();
	list<Mat> classifiedSegments;

	// Runs all classifiers for each segment
	int i = 1;
	int n = _segmentation.getSegments().size();
	for (auto mask : _segmentation.getSegments()) {
		// Gets the opinion of all classifiers for this segment
		vector<Mat> classifications;
		classifications.reserve(classifiers.size());
		for (auto c : classifiers) {
			Mat classification = c->classify(mask);
			classifications.push_back(classification);
		}

		// Initializes the consensus matrix for this segment
		Mat consensus = Mat::zeros(mapSize, CV_8UC1);

		if (_consensusType == MAJORITY_VOTING) {
			consensus += classifications[0];
		}
		else {
			assert(false && "Unknown consensus type");
		}

		classifiedSegments.push_back(consensus);
		i++;
	}

	// Builds the map from the segments
	Mat theMap = Mat::zeros(mapSize, CV_8UC1);
	for (auto seg : classifiedSegments) {
		theMap += seg;
	}

	return CoverMap(theMap);
}
