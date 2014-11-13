#include "config.h"

void config::parse(char **argv, int argc, Configuration &config) {
	int c = 1;

	while (c < argc) {
		// File paths
		if (streq("--vis", argv[c])) {
			strcpy(config.pathVIS, argv[++c]);
		}
		else if (streq("--lwir", argv[c])) {
			strcpy(config.pathLWIR, argv[++c]);
		}
		else if (streq("--training", argv[c])) {
			strcpy(config.pathTraining, argv[++c]);
		}

		// ROI
		else if (streq("--roiX", argv[c])) {
			config.roiX = atoi(argv[++c]);
		}
		else if (streq("--roiY", argv[c])) {
			config.roiY = atoi(argv[++c]);
		}
		else if (streq("--roiW", argv[c])) {
			config.roiWidth = atoi(argv[++c]);
		}
		else if (streq("--roiH", argv[c])) {
			config.roiHeight = atoi(argv[++c]);
		}
		else if (streq("--disable-roi", argv[c])) {
			// Resets the ROI even if it was filled before,
			// for easy command-line manipulation
			config.roiX = 0;
			config.roiY = 0;
			config.roiWidth = 0;
			config.roiHeight = 0;
		}

		// Verbosity
		else if (streq("--verbose", argv[c])) {
			config.verbose = true;
		}

		// Segmentation mode
		else if (streq("--segmentation-mode", argv[c])) {
			char theMode[512];
			strcpy(theMode, argv[++c]);
			if (streq("GRID", theMode)) {
				config.segmentationMode = GRID;
			}
		}

		// Segmentation parameters - GRID
		else if (streq("--grid-tile-size", argv[c])) {
			config.gridTileSize = atoi(argv[++c]);
		}

		// Log
		else if (streq("--log-path", argv[c])) {
			strcpy(config.logPath, argv[++c]);
			config.logEnabled = true;
		}

		// Scaling
		else if (streq("--scale-vis", argv[c])) {
			config.scaleVIS = atof(argv[++c]);
		}
		else if (streq("--scale-lwir", argv[c])) {
			config.scaleLWIR = atof(argv[++c]);
		}

		// Resampling method
		else if (streq("--resampling-method", argv[c])) {
			char theMethod[512];
			strcpy(theMethod, argv[++c]);
			if (streq("NEAREST_NEIGHBOR", theMethod)) {
				config.resamplingMethod = NEAREST_NEIGHBOR;
			}
			else if (streq("LINEAR" , theMethod)) {
				config.resamplingMethod = LINEAR;
			}
			else if (streq("CUBIC" , theMethod)) {
				config.resamplingMethod = CUBIC;
			}
		}

		// Parallelization
		else if (streq("--parallel", argv[c])) {
			config.parallel = true;
		}

		c++;
	}
}
