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
	list<SparseMat> masks = segmentation::getColorBlobs(_training.asMat());

	// Recover region labels
	list<float> labels;
	list<SparseMat> validSegments;
	for (auto mask : masks) {
		float label = _training.getRegionClass(densify(mask));

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

void Ensemble::doClassify(Classifier* C, Segmentation S, int *cursor, 
		int *classifiedSegments, int totalToClassify) {

	Mat classification = Mat::zeros(S.getMapSize(), CV_8UC1);

	int i = 1, n = S.getSegments().size();
	for (auto mask : S.getSegments()) {
		float progress = 100.0 * (*classifiedSegments) / totalToClassify;
		cerr << "\rclassification: " << 
			setprecision(2) << 
			progress << "%        \r" << flush;

		classification += C->classify(mask);
		i++;
		(*classifiedSegments)++;
	}

	if (!_parallel) {
		_mutex.lock();
	}

	auto p = make_pair(C->getID(), classification);
	_classifications[*cursor] = p;
	(*cursor)++;

	if (!_parallel) {
		_mutex.unlock();
	}
}

CoverMap Ensemble::classify() {
	// Some initialization
	Size mapSize = _segmentation.getMapSize();

	// Initializes individual classifications
	_classifications.clear();
	_classifications.resize(classifiers.size());

	int currentCursor = 0;
	int totalClassified = 0;
	int totalToClassify = _segmentation.getSegments().size() * 
		classifiers.size();

	// Runs all classifiers
	if (_parallel) {
		// Creates the threads
		list<thread> threads;
		for (auto& c : classifiers) {
			threads.push_back(thread(&Ensemble::doClassify, this, c, 
					_segmentation, &currentCursor, &totalClassified,
					totalToClassify));
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}	
	}
	else {
		// Run serially
		for (auto& c : classifiers) {
			this->doClassify(c, _segmentation, &currentCursor, 
				&totalClassified, totalToClassify);
		}
	}

	// Get consensus
	Mat consensus = Mat::zeros(mapSize, CV_8UC1);

	if (_consensusType == MAJORITY_VOTING) {
		for (int row = 0; row < mapSize.height; row++) {
			for (int col = 0; col < mapSize.width; col++) {
				Counter<unsigned char> counter;

				for (auto op : _classifications) {
					counter.inc(op.second.at<unsigned char>(row, col));
				}

				consensus.at<unsigned char>(row, col) = counter.top();
			}
		}
	}
	else {
		assert(false && "Unknown consensus type");
	}

	return CoverMap(consensus);
}

vector<pair<string, Mat> > Ensemble::individualClassifications() {
	return _classifications;
}

void Ensemble::setThreads(int n) {
	numThreads = n;
}

void Ensemble::setParallel(bool p) {
	_parallel = p;
}
