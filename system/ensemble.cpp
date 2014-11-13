#include "ensemble.h"

using namespace cv;

Ensemble::Ensemble(ConsensusType t, Segmentation& s, ThematicMap trainingVIS,
	ThematicMap trainingLWIR)

    : _consensusType(t),
	  _segmentation(s),
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
	cerr << "extracting region masks from VIS training map..." << endl;
	list<SparseMat> masksVIS = segmentation::getColorBlobs(
		_trainingVIS.asMat());

	// Recover region labels - VIS
	cerr << "recovering region labels from VIS training map..." << endl;
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
				labelsLWIR.push_back(label);

				Mat mask = Mat::zeros(lwirSize, CV_8UC1);
				mask.at<unsigned char>(row, col) = 255;

				validSegmentsLWIR.push_back(SparseMat(mask));
			}

			currentPixel++;
		}
	}

	cerr << endl;

	// Create labels training matrix - VIS
	cerr << "creating VIS training matrix..." << endl;
	Mat labelsMatVIS(labelsVIS.size(), 1, CV_32FC1);
	int c = 0;
	for (auto label : labelsVIS) {
		labelsMatVIS.at<float>(c) = label;
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

	// Classify all segments
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

	// Resize classification map if necessary
	if (classification.size() != _trainingVIS.size()) {
		Mat R;
		resize(classification, R, _trainingVIS.size(), 0, 0, 
			TRAINING_INTERPOLATION_MODE);

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
	// Some initialization
	Size mapSize = _segmentation.getMapSize();

	// Initializes individual classifications
	_classifications.clear();
	_classifications.resize(classifiers.size());

	// Initialize counters (will be shared among threads)
	int currentCursor = 0;
	int totalClassified = 0;

	// Create VIS and LWIR segmentations
	Segmentation segmentationVIS = _segmentation;
	Segmentation segmentationLWIR = _pixelizeLWIR ? segmentationVIS.pixelize() 
		: segmentationVIS;

	// Determine total segments to classify
	int segmentsPerVISImage = segmentationVIS.segmentCount();
	int segmentsPerLWIRImage = segmentationLWIR.segmentCount();
	int totalToClassify = 0;
	for (auto& c : classifiers) {
		if (_pixelizeLWIR && c->getType() == LWIR) {
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
		for (auto& c : classifiers) {
			// Determine which segmentation to send to classifier
			Segmentation S;
			if (_pixelizeLWIR && c->getType() == LWIR) {
				S = segmentationLWIR;
			}
			else {
				S = segmentationVIS;
			}

			// Start the thread
			threads.push_back(thread(&Ensemble::doClassify, this, c, 
					S, &currentCursor, &totalClassified,
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
