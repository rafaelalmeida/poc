#ifndef CLASSIFICATION_H
#define CLASSIFICATION_H

#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/ml/ml.hpp>

#include "classification.h"
#include "description.h"
#include "description_vis.h"
#include "description_lwir.h"
#include "models.h"
#include "segmentation.h"
#include "utils.h"

using namespace segmentation;

namespace classification {
	cv::Mat predict(cv::Mat image, segmentation::Segmentation segmentation, CvSVM *classifier);

	enum ClassifierType {
		VIS,
		LWIR
	};

	enum ClassifierEngine {
		SVM
	};

	/**
	 * Defines a classifier. Descriptors are automatically destroyed when
	 * Classifier is destroyed.
	 */
	class Classifier {
		// Members
		cv::Mat _vis;
		LWIRImage *_lwir = NULL;

		Descriptor *_descriptor;

		ClassifierType _type;
		ClassifierEngine _engine;

		CvSVM _svm;

		public:
			// Constructors
			Classifier(ClassifierEngine engine, cv::Mat vis, Descriptor *descriptor);
			Classifier(ClassifierEngine engine, LWIRImage *lwir, Descriptor *descriptor);

			// Destructors
			~Classifier();

			// Methods
			void train(CoverMap training);
			cv::Mat classify();
	};
}

#endif
