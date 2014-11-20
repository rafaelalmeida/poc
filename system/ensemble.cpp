#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& vis, Segmentation& lwir, 
			Segmentation& visTraining, Segmentation& lwirTraining,
			ThematicMap trainingVIS, ThematicMap trainingLWIR)

    : _consensusType(t),
	  _segmentationVIS(vis),
	  _segmentationLWIR(lwir),
	  _segmentationVISTraining(visTraining),
	  _segmentationLWIRTraining(lwirTraining),
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

	// Trains all classifiers
	int i = 1;
	int n = classifiers.size();

	list<thread> threads;
	for (auto c : classifiers) {
		// Determine which information to send
		Mat lbl;
		Segmentation S;
		if (c->getType() == ImageType::VIS) {
			lbl = labelsMatVIS;
			S = _segmentationVISTraining;
		}
		else {
			lbl = labelsMatLWIR;
			S = _segmentationLWIRTraining;
		}

		// Run or start the threads
		if (_parallel) {
			threads.push_back(thread(&Ensemble::doTrain, this, c, lbl,
				S, &i, n));	
		}
		else {
			this->doTrain(c, lbl, S, &i, n);
		}
	}

	// Syncronization
	if (_parallel) {
		for (auto& t : threads) {
			t.join();
		}
	}

	cerr << "training classifiers... done             " << endl;
}

void Ensemble::doTrain(Classifier *C, Mat labels, Segmentation S, int *trained, 
	int totalToTrain) {

	// Report progress
	if (_parallel) _consoleMutex.lock();

	cerr << "training " << *trained << " of " << totalToTrain << 
		" classifiers (" << C->getID() << ")      \r" << flush;
	
	if (_parallel) _consoleMutex.unlock();

	// Execute training
	C->train(labels, S);
	(*trained)++;

	// Log time taken
	_totalTimeDescription += C->getDescriptionTime();
	_totalTimeTraining += C->getTrainingTime();
}

void Ensemble::doClassify(Classifier* C, Size mapSize, Segmentation S, 
	int *cursor, int *classifiedSegments, int totalToClassify, int idx) {

	Mat classification = Mat::zeros(S.getMapSize(), CV_8UC1);

	// Classify all segments
	int n = S.regionCount();
	for (auto region : S.getRegions()) {
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
		Mat classifiedSegment = C->classify(region);
		classification += classifiedSegment;
		_worklog[idx]++;
		(*classifiedSegments)++;

		// Log time taken
		_totalTimeDescription += C->getDescriptionTime();
		_totalTimeClassification += C->getClassificationTime();
	}

	// Resize classification map if necessary
	if (classification.size() != mapSize) {
		Mat R;
		resize(classification, R, mapSize, 0, 0, TRAINING_INTERPOLATION_MODE);
		classification = R;
	}

	// Will update shared memory, lock it
	if (!_parallel) _mutex.lock();

	// Save individual classification
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
	int segmentsPerVISImage = _segmentationVIS.regionCount();
	int segmentsPerLWIRImage = _segmentationLWIR.regionCount();
	int totalToClassify = 0;
	for (auto& c : classifiers) {
		if (_pixelizeLWIR && c->getType() == ImageType::LWIR) {
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
			if (_pixelizeLWIR && c->getType() == ImageType::LWIR) {
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
			if (_pixelizeLWIR && c->getType() == ImageType::LWIR) {
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

void Ensemble::resetAllTimers() {
	_totalTimeDescription = 0;
	_totalTimeTraining = 0;
	_totalTimeClassification = 0;
}

double Ensemble::getTotalDescriptionTime() {
	return _totalTimeDescription;
}

double Ensemble::getTotalTrainingTime() {
	return _totalTimeTraining;
}

double Ensemble::getTotalClassificationTime() {
	return _totalTimeClassification;
}
