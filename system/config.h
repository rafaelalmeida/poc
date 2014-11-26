#ifndef CONFIG_H
#define CONFIG_H

#include "common.h"
#include "utils.h"
#include "constants.h"

#define MAX_PATH 256

typedef struct {
	// File paths
	char pathVIS[MAX_PATH];
	char pathLWIR[MAX_PATH];
	char pathTraining[MAX_PATH];

	// Test set file paths
	bool hasTestSet = false;
	char pathVISTest[MAX_PATH];
	char pathLWIRTest[MAX_PATH];

	// ROI specification
	int roiX;
	int roiY;
	int roiWidth;
	int roiHeight;

	// Verbosity
	bool verbose = false;

	// Segmentation options
	SegmentationMode segmentationMode;

	// Segmentation parameters - SLIC
	int slicRegionSize = DEFAULT_SLIC_REGION_SIZE;
	int slicMinRegionSize = DEFAULT_SLIC_MIN_REGION_SIZE;
	float slicRegularization = DEFAULT_SLIC_REGULARIZATION;
	bool slicAutoScaleParameters = false;

	// Segmentation paramaters - GRID
	int gridTileSize = 10;

	// Log options
	bool logEnabled = false;
	char logPath[MAX_PATH];

	// Sampling and resampling
	float scaleVIS = 1.0;
	float scaleLWIR = 1.0;
	InterpolationMode interpolationMode = NEAREST_NEIGHBOR;

	// Parallelization
	bool parallel = false;

} Configuration;

namespace config {
	void parse(char **argv, int argc, Configuration &config);
}

#endif
