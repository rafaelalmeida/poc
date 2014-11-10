#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <list>
#include <stdio.h>
#include <thread>

#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "classification.h"
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
	Segmentation& _segmentation;

	CoverMap _training;

	Logger *_logger = NULL;

	std::vector<std::pair<std::string, cv::Mat> > _classifications;

	bool _parallel = false;
	int numThreads;
	mutex _mutex;

	public:
		// Constructors
		Ensemble(ConsensusType t, Segmentation& s, CoverMap training);

		// Destructors
		~Ensemble();

		// Methods
		CoverMap classify();
		std::vector<std::pair<std::string, cv::Mat> > individualClassifications();
		void addClassifier(Classifier* c);
		void doClassify(Classifier* C, Segmentation S, int *cursor, 
			int *classifiedSegments, int totalToClassify);
		void setLogger(Logger *logger);
		void setThreads(int n);
		void train();
		void setParallel(bool p);
};

#endif
