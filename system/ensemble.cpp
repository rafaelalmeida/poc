#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& s, CoverMap trainingVIS,
	CoverMap trainingLWIR)

    : _consensusType(t),
	  _segmentation(s),
	  _trainingVIS(trainingVIS),
	  _trainingLWIR(trainingLWIR) {

	assert(trainingVIS.asMat().type() == CV_8UC1);
	assert(trainingLWIR.asMat().type() == CV_8UC1);
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
	// Extract regions from training map, for VIS image
	list<SparseMat> masksVIS = segmentation::getColorBlobs(
		_trainingVIS.asMat());

	// Recover region labels - VIS
	list<float> labelsVIS;
	list<SparseMat> validSegmentsVIS;
	for (auto mask : masksVIS) {
		float label = _trainingVIS.getRegionClass(densify(mask));

		if (label != 0) { // Disconsider unclassified regions
			labelsVIS.push_back(label);
			validSegmentsVIS.push_back(mask);
		}
	}

	// Recover region labels - LWIR
	list<float> labelsLWIR;
	list<SparseMat> validSegmentsLWIR;
	Mat trainingLWIRMat = _trainingLWIR.asMat();
	Size lwirSize = trainingLWIRMat.size();

	for (int row = 0; row < lwirSize.height; row++) {
		for (int col = 0; col < lwirSize.width; col++) {
			float label = (float) trainingLWIRMat.at<unsigned char>(row, col);

			if (label != 0) {
				labelsLWIR.push_back(label);

				Mat mask = Mat::zeros(lwirSize, CV_8UC1);
				mask.at<unsigned char>(row, col) = 255;

				validSegmentsLWIR.push_back(SparseMat(mask));
			}
		}
	}

	// Create labels training matrix - VIS
	Mat labelsMatVIS(labelsVIS.size(), 1, CV_32FC1);
	int c = 0;
	for (auto label : labelsVIS) {
		labelsMatVIS.at<float>(c) = label;
		c++;
	}

	// Create labels training matrix - LWIR
	Mat labelsMatLWIR(labelsLWIR.size(), 1, CV_32FC1);
	c = 0;
	for (auto label : labelsLWIR) {
		labelsMatLWIR.at<float>(c) = label;
		c++;
	}

	int i = 1;
	int n = classifiers.size();

	for (auto c : classifiers) {
		cerr << "training classifier " << i << " of " << n << "     \r" << flush;

		if (c->getType() == VIS) {
			c->train(labelsMatVIS, Segmentation(validSegmentsVIS));
		}
		else {
			c->train(labelsMatLWIR, Segmentation(validSegmentsLWIR));
		}

		i++;
	}

	cerr << "training classifiers... done     " << endl;
}

void Ensemble::doClassify(Classifier* C, Segmentation S, int *cursor, 
		int *classifiedSegments, int totalToClassify) {

	Mat classification = Mat::zeros(S.getMapSize(), CV_8UC1);

	// Pixelize segmentation if this is a LWIR classifier, since they 
	// are classified pixel by pixel
	if (_pixelizeLWIR && C->getType() == LWIR) {
		S = S.pixelize();
	}

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

void Ensemble::setLWIRPixelize(bool p) {
	_pixelizeLWIR = p;
}
