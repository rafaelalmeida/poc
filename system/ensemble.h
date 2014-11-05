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
#include "segmentation.h"

using namespace classification;
using namespace segmentation;

enum ConsensusType {
	MAJORITY_VOTING
};

class Ensemble {
	// Members
	std::list<Classifier> classifiers;
	ConsensusType _consensusType;
	Segmentation& _segmentation;

	public:
		// Constructors
		Ensemble(ConsensusType t, Segmentation& s);

		// Methods
		void addClassifier(Classifier c);
};

#endif