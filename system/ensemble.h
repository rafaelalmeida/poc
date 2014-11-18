#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <list>
#include <stdio.h>
#include <thread>
#include <mutex>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "classification.h"
#include "common.h"
#include "logging.h"
#include "segmentation.h"

using namespace classification;
using namespace segmentation;

enum ConsensusType {
	MAJORITY_VOTING
};

/**
 * Defines a classifier ensemble. Classifiers are automatically destroyed
 * when Ensemble is destroyed.
 */
class Ensemble {
	// Members
	std::vector<Classifier*> classifiers;
	ConsensusType _consensusType;
	Segmentation& _segmentationVIS;
	Segmentation& _segmentationLWIR;

	// Flag to see if segmentation is ignored for LWIR images, and the image
	// is classified pixel-by-pixel instead
	bool _pixelizeLWIR = true;

	ThematicMap _trainingVIS;
	ThematicMap _trainingLWIR;

	Logger *_logger = NULL;

	std::vector<std::pair<std::string, cv::Mat> > _classifications;

	bool _parallel = false;
	int numThreads;
	mutex _mutex;
	mutex _consoleMutex;
	vector<int> _worklog;

	public:
		// Constructors
		Ensemble(ConsensusType t, Segmentation& segmentationVIS, 
			Segmentation& segmentationLWIR, ThematicMap trainingVIS, 
			ThematicMap trainingLWIR);

		// Destructors
		~Ensemble();

		// Methods
		ThematicMap classify();
		std::vector<std::pair<std::string, cv::Mat> > individualClassifications();
		void addClassifier(Classifier* c);
		void doClassify(Classifier* C, cv::Size mapSize, Segmentation S, 
			int *cursor, int *classifiedSegments, int totalToClassify, 
			int threadID);
		void setLogger(Logger *logger);
		void setThreads(int n);
		void train();
		void setLWIRPixelize(bool p);
		void setParallel(bool p);
};

#endif
