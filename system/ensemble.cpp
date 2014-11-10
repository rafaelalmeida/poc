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
	_classifications.clear();
	_classifications.resize(classifiers.size());

	for (auto& c : _classifications) {
		c = Mat::zeros(mapSize, CV_8UC1);
	}

	// Runs all classifiers for each segment
	int i = 0;
	int n = _segmentation.getSegments().size();

	for (auto& c : classifiers) {
		for (auto mask : _segmentation.getSegments()) {
			Mat classification = c->classify(mask);
			_classifications[i] += classification;
		}

		i++;
	}

	// Get consensus
	Mat consensus = Mat::zeros(mapSize, CV_8UC1);

	if (_consensusType == MAJORITY_VOTING) {
		for (auto mask : _segmentation.getSegments()) {
			Counter<int> counter;

			for (auto op: _classifications) {
				int theClass = (int) segmentation::getSegmentLabel(op, mask);
				counter.inc(theClass);
			}

			consensus += counter.top() * (mask / 255);
		}
	}
	else {
		assert(false && "Unknown consensus type");
	}

	return CoverMap(consensus);
}

vector<Mat> Ensemble::individualClassifications() {
	return _classifications;
}
