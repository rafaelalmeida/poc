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
#include "description_vis.h"
#include "description_lwir.h"
#include "models.h"
#include "segmentation.h"

using namespace segmentation;

namespace classification {
	CvSVM *trainSVM(cv::Mat image, cv::Mat trainingMap, cv::Mat (*descriptor)(cv::Mat, const std::list<cv::Mat>));
	cv::Mat predict(cv::Mat image, segmentation::Segmentation segmentation, CvSVM *classifier);

	enum ClassifierType {
		VIS,
		LWIR
	};

	enum ClassifierEngine {
		SVM
	};

	class Classifier {
		// Members
		cv::Mat _vis;
		LWIRImage _lwir;
		Segmentation _segmentation;

		ClassifierType _type;
		ClassifierEngine _engine;

		public:
			// Constructors
			Classifier(cv::Mat vis, Segmentation segments, cv::Mat (*descriptor)(cv::Mat, cv::Mat));
			Classifier(LWIRImage lwir, Segmentation segments, cv::Mat (*descriptor)(LWIRImage, cv::Mat));

			// Methods
			void train();
			cv::Mat classify();
	};
}

#endif