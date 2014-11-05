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

using namespace classification;

class Ensemble {
	// Members
	std::list<Classifier> classifiers;

	public:
		// Methods
		void addClassifier(Classifier& c);
};

#endif