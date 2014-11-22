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

	public:
		// Constructors
		System(Configuration config);

		// Main entry point
		void run();

	private:
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
