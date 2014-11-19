#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& segmentationVIS, 
	Segmentation& segmentationLWIR, ThematicMap trainingVIS,
	ThematicMap trainingLWIR)

    : _consensusType(t),
	  _segmentationVIS(segmentationVIS),
	  _segmentationLWIR(segmentationLWIR),
	  _trainingVIS(trainingVIS),
	  _trainingLWIR(trainingLWIR) {

	// Verify training matrixes types
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
	cerr << "enumerating regions from VIS training map..." << endl;
	list<pair<SparseMat, int> > trainingRegionsVIS = 
		_trainingVIS.enumerateRegions();
	
	list<SparseMat> validSegmentsVIS;
	for (auto p : trainingRegionsVIS) {
		validSegmentsVIS.push_back(p.first);
	}

	// Recover region labels - LWIR
	list<float> labelsLWIR;
	list<SparseMat> validSegmentsLWIR;
	Mat trainingLWIRMat = _trainingLWIR.asMat();
	Size lwirSize = trainingLWIRMat.size();

	int currentPixel = 1;
	int totalPixels = lwirSize.height * lwirSize.width;

	for (int row = 0; row < lwirSize.height; row++) {
		for (int col = 0; col < lwirSize.width; col++) {
			// Report progress
			float progress = 100.0 * currentPixel / totalPixels;
			cerr << "recovering region labels from LWIR training map... " << 
				progress << "%      \r";

			float label = (float) trainingLWIRMat.at<unsigned char>(row, col);

			if (label != 0) {
				// Save label
				labelsLWIR.push_back(label);

				// Save segment mask
				int sizes[2] = {lwirSize.height, lwirSize.width};
    			SparseMat mask(2, sizes, CV_8UC1);
    			*(mask.ptr(row, col, true)) = 255;
				validSegmentsLWIR.push_back(mask);
			}

			currentPixel++;
		}
	}

	cerr << endl;

	// Create labels training matrix - VIS
	cerr << "creating VIS training matrix..." << endl;
	Mat labelsMatVIS(trainingRegionsVIS.size(), 1, CV_32FC1);
	int c = 0;
	for (auto region : trainingRegionsVIS) {
		labelsMatVIS.at<float>(c) = (float) region.second;
		c++;
	}

	// Create labels training matrix - LWIR
	cerr << "creating LWIR training matrix..." << endl;
	Mat labelsMatLWIR(labelsLWIR.size(), 1, CV_32FC1);
	c = 0;
	for (auto label : labelsLWIR) {
		labelsMatLWIR.at<float>(c) = label;
		c++;
	}

	int i = 1;
	int n = classifiers.size();

	for (auto c : classifiers) {
		cerr << "training classifier " << i << " of " << n << " (" << 
			c->getID() << ")      \r" << flush;

		if (c->getType() == ClassifierType::VIS) {
			c->train(labelsMatVIS, Segmentation(validSegmentsVIS));
		}
		else {
			c->train(labelsMatLWIR, Segmentation(validSegmentsLWIR));
		}

		i++;
	}

	cerr << "training classifiers... done     " << endl;
}

void Ensemble::doClassify(Classifier* C, Size mapSize, Segmentation S, 
	int *cursor, int *classifiedSegments, int totalToClassify, int idx) {

	Mat classification = Mat::zeros(S.getMapSize(), CV_8UC1);

	// Classify all segments
	int n = S.segmentCount();
	for (auto mask : S.getSegments()) {
		// Report progress from time to time
		if (*classifiedSegments % 10 == 0) {

			_consoleMutex.lock();

			float progress = 100.0 * (*classifiedSegments + 1) / 
				totalToClassify;
			cerr << "\rclassification: " << 
				progress << "%        \r" << flush;

			_consoleMutex.unlock();
		}

		// Execute classification
		Mat classifiedSegment = C->classify(mask);
		classification += classifiedSegment;
		_worklog[idx]++;
		(*classifiedSegments)++;
	}

	// Resize classification map if necessary
	if (classification.size() != mapSize) {
		Mat R;
		resize(classification, R, mapSize, 0, 0, TRAINING_INTERPOLATION_MODE);
		classification = R;
	}

	// Will update shared memory, lock it
	if (!_parallel) _mutex.lock();

	auto p = make_pair(C->getID(), classification);
	_classifications[*cursor] = p;
	(*cursor)++;

	// Finished updating shared memory, unlock it
	if (!_parallel) _mutex.unlock();
}

ThematicMap Ensemble::classify() {
	// We consider the VIS size as the real map size, since the LWIR 
	// classifications will be resized into the VIS size for correspondence
	Size mapSize = _segmentationVIS.getMapSize();

	// Initializes individual classifications
	_classifications.clear();
	_classifications.resize(classifiers.size());

	// Initialize counters (will be shared among threads)
	int currentCursor = 0;
	int totalClassified = 0;

	// Determine total segments to classify
	int segmentsPerVISImage = _segmentationVIS.segmentCount();
	int segmentsPerLWIRImage = _segmentationLWIR.segmentCount();
	int totalToClassify = 0;
	for (auto& c : classifiers) {
		if (_pixelizeLWIR && c->getType() == ClassifierType::LWIR) {
			totalToClassify += segmentsPerLWIRImage;
		}
		else {
			totalToClassify += segmentsPerVISImage;
		}
	}

	// Runs all classifiers
	if (_parallel) {
		// Creates the threads
		list<thread> threads;

		// Initialize work count per thread
		int threadID = 0;
		_worklog.resize(classifiers.size());
		fill(_worklog.begin(), _worklog.end(), 0);

		for (auto& c : classifiers) {
			// Determine which segmentation to send to classifier
			Segmentation S;
			if (_pixelizeLWIR && c->getType() == ClassifierType::LWIR) {
				S = _segmentationLWIR;
			}
			else {
				S = _segmentationVIS;
			}

			// Start the thread
			threads.push_back(thread(&Ensemble::doClassify, this, c, mapSize,
					S, &currentCursor, &totalClassified,
					totalToClassify, threadID));

			threadID++;
		}

		// Wait for all threads to finish
		for (auto& t : threads) {
			t.join();
		}
	}
	else {
		// Run serially
		for (auto& c : classifiers) {
			// Determine which segmentation to send to classifier
			Segmentation S;
			if (_pixelizeLWIR && c->getType() == ClassifierType::LWIR) {
				S = _segmentationLWIR;
			}
			else {
				S = _segmentationVIS;
			}

			// Run classifier
			_worklog.resize(1);
			_worklog[0] = 0;

			this->doClassify(c, mapSize, S, &currentCursor, 
				&totalClassified, totalToClassify, 0);
		}
	}

	// Report classification complete
	cerr << "classification: done           " << endl;

	// Report work done by each thread
	cerr << "thread jobs done: ";
	for (auto l : _worklog) {
		cerr << l << " ";
	}
	cerr << endl;

	// Get consensus
	Mat consensus = Mat::zeros(mapSize, CV_8UC1);

	if (_consensusType == MAJORITY_VOTING) {
		for (int row = 0; row < mapSize.height; row++) {
			for (int col = 0; col < mapSize.width; col++) {
				// Determine winning class
				Counter<unsigned char> counter;
				for (auto op : _classifications) {
					counter.inc(op.second.at<unsigned char>(row, col));
				}

				// Set the winning class on the consensus map
				consensus.at<unsigned char>(row, col) = counter.top();
			}
		}
	}
	else {
		assert(false && "Unknown consensus type");
	}

	return ThematicMap(consensus);
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
