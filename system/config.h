#ifndef CONFIG_H
#define CONFIG_H

#include "utils.h"
#include "constants.h"

#define MAX_PATH 256

typedef struct {
	// File paths
	char pathVIS[MAX_PATH];
	char pathLWIR[MAX_PATH];
	char pathTraining[MAX_PATH];

	// ROI specification
	int roiX;
	int roiY;
	int roiWidth;
	int roiHeight;

	// Verbosity
	bool verbose = false;

	// Segmentation options
	SegmentationMode segmentationMode;
} Configuration;

namespace config {
	void parse(char **argv, int argc, Configuration &config);
}

#endif
