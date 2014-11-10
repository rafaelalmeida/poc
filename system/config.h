#ifndef CONFIG_H
#define CONFIG_H

#include "utils.h"
#include "constants.h"

#define MAX_PATH 256

enum SamplingMode {
	UPSAMPLE_LWIR,
	DOWNSAMPLE_VIS
};

enum ResamplingMethod {
	NEAREST_NEIGHBOR,
	LINEAR,
	CUBIC
};

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

	// Segmentation paramaters - GRID
	int gridTileSize = 10;

	// Log options
	bool logEnabled = false;
	char logPath[MAX_PATH];

	// Sampling and resampling
	SamplingMode samplingMode;
	ResamplingMethod resamplingMethod = NEAREST_NEIGHBOR;

	// Parallelization
	bool parallel = false;

} Configuration;

namespace config {
	void parse(char **argv, int argc, Configuration &config);
}

#endif
