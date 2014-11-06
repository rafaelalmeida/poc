#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

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
	std::list<Classifier*> classifiers;
	ConsensusType _consensusType;
	Segmentation& _segmentation;

	CoverMap _training;

	Logger *_logger = NULL;

	std::vector<cv::Mat> _classifications;

	public:
		// Constructors
		Ensemble(ConsensusType t, Segmentation& s, CoverMap training);

		// Destructors
		~Ensemble();

		// Methods
		void addClassifier(Classifier* c);
		void train();
		CoverMap classify();
		void setLogger(Logger *logger);
		std::vector<cv::Mat> individualClassifications();
};

#endif
