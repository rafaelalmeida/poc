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
	// Some initialization
	Size mapSize = _segmentation.getMapSize();
	list<Mat> classifiedSegments;

	// Saves individual classifications for debugging
	vector<Mat> classifications(classifiers.size());
	for (auto& c : classifications) {
		c = Mat::zeros(mapSize, CV_8UC1);
	}

	// Runs all classifiers for each segment
	for (auto mask : _segmentation.getSegments()) {
		// Gets the opinion of all classifiers for this segment
		vector<Mat> opinions;
		opinions.reserve(classifiers.size());

		int i = 0;
		for (auto c : classifiers) {
			Mat classification = c->classify(mask);
			opinions.push_back(classification);

			// Save this segment individually for debugging
			classifications[i] += classification;

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
	}

	// Saves debugging images, if applies
	int i = 1;
	for (auto& c : classifications) {
		char path[16];
		sprintf(path, "classification_%d", i);

		_logger->saveImage(path, CoverMap(c).coloredMap());

		i++;
	}

	// Builds the map from the segments
	Mat theMap = Mat::zeros(mapSize, CV_8UC1);
	for (auto seg : classifiedSegments) {
		theMap += seg;
	}

	return CoverMap(theMap);
}
