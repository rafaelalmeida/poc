#ifndef SYSTEM_H
#define SYSTEM_H

#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <opencv2/opencv.hpp>

#include "classification.h"
#include "common.h"
#include "config.h"
#include "ensemble.h"
#include "gdal_driver.h"
#include "logging.h"
#include "models.h"
#include "segmentation.h"
#include "statistics.h"
#include "utils.h"

// Class to manage the whole classifier system
class System {
	Configuration _conf;
	Logger _logger;

	// Data members
	VISImage vis;
	LWIRImage lwir;
	Mat training;

	VISImage visTest;
	LWIRImage lwirTest;

	ThematicMap trainingMapVIS;
	ThematicMap trainingMapLWIR;

	Segmentation segmentationVIS;
	Segmentation segmentationLWIR;
	
	Segmentation segmentationVISTest;
	Segmentation segmentationLWIRTest;

	vector<Descriptor*> descriptors;

	Segmentation visTrainingSegmentation;
	Segmentation lwirTrainingSegmentation;

	// Initialize time loggers
	double descriptionTime = 0;
	double trainingTime = 0;
	double classificationTime = 0;

	public:
		// Constructors
		System(Configuration config);

		// Main entry point
		void run();

	private:
		// System phases
		void loadAndSetupData();
		void segment();
		void describe();

		// Other procedures
		void rescale(VISImage& vis, LWIRImage& lwir, float scaleVIS, 
			float scaleLWIR, InterpolationMode interpolationMode);
		void setupClassifiers(Ensemble& e, VISImage& vis, LWIRImage& lwir, 
			vector<Descriptor*> descriptors);
		void printClassHistogram(ThematicMap& M);
		void describeAll(list<Segmentation*>& segmentations, 
			vector<Descriptor*>& descriptors, bool parallel);
		void doDescribe(Segmentation *S, Descriptor *D, int *numDone, int n, 
			mutex *mtx);
		void reportDescriptionProcess(list<Segmentation*> segmentations, 
			int *numDone, int n, mutex *mtx);
		void log(const char *msg);
};

#endif // SYSTEM_H
