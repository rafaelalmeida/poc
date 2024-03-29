#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <opencv2/opencv.hpp>

#include "classification.h"
#include "description.h"
#include "models.h"
#include "segmentation.h"
#include "utils.h"

using namespace segmentation;

namespace classification {
	// Classifier engines
	enum ClassifierEngine {
		SVM,
		NBC,
		KNN,
		DTREE,
		GBT,
		RTREES,
		ERTREES,
		MLP
	};

	/**
	 * Defines a classifier. Descriptors are NOT automatically destroyed when
	 * Classifier is destroyed, you must destroy them manually because you can
	 * reuse them across classifiers.
	 */
	class Classifier {
		// Images this classifier will classify
		VISImage *_vis = NULL;
		LWIRImage *_lwir = NULL;

		// The optional training map. This can be used to compute class
		// probabilities for the DTREE classifier.
		ThematicMap *_trainingMap = NULL;

		cv::Size _mapSize;

		// The descriptor it will use
		Descriptor *_descriptor;

		// Some meta information
		ImageType _type;
		ClassifierEngine _engine;

		// Basic instrumentation
		Stopwatch _swatchDescription;
		Stopwatch _swatchTraining;
		Stopwatch _swatchClassification;

		// Classifier object members
		CvSVM _svm;
		CvNormalBayesClassifier _nbc;
		CvKNearest _knn;
		CvDTree _dtree;
		CvGBTrees _gbtrees;
		CvRTrees _rtrees;
		CvERTrees _ertrees;
		CvANN_MLP _mlp;

		public:
			// Constructors
			Classifier(ClassifierEngine engine, VISImage *vis, 
				Descriptor *descriptor);
			Classifier(ClassifierEngine engine, LWIRImage *lwir, 
				Descriptor *descriptor);

			// Methods
			void train(cv::Mat labels, Segmentation S);
			cv::Mat classify(Region region);
			string getID();
			ImageType getType();
			void setTrainingMap(ThematicMap *map);

			// Instrumentation methods
			double getDescriptionTime();
			double getTrainingTime();
			double getClassificationTime();
	};
}

#endif
