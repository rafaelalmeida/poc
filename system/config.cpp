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

		// File paths - test
		else if (streq("--vis-test", argv[c])) {
			config.hasTestSet = true;
			strcpy(config.pathVISTest, argv[++c]);
		}
		else if (streq("--lwir-test", argv[c])) {
			config.hasTestSet = true;
			strcpy(config.pathLWIRTest, argv[++c]);
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
			else if (streq("SLIC", theMode)) {
				config.segmentationMode = SLIC;
			}
			else {
				cerr << "Unrecognized segmentation mode: " << theMode << endl;
				exit(EXIT_FAILURE);
			}
		}

		// Segmentation parameters - GRID
		else if (streq("--grid-tile-size", argv[c])) {
			config.gridTileSize = atoi(argv[++c]);
		}

		// Segmentation parameters - SLIC
		else if (streq("--slic-region-size", argv[c])) {
			config.slicRegionSize = atoi(argv[++c]);
		}
		else if (streq("--slic-min-region-size", argv[c])) {
			config.slicMinRegionSize = atoi(argv[++c]);
		}
		else if (streq("--slic-regularization", argv[c])) {
			config.slicRegularization = atof(argv[++c]);
		}
		else if (streq("--slic-auto-scale-parameters", argv[c])) {
			config.slicAutoScaleParameters = true;
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
			if (streq("NEAREST", theMethod)) {
				config.interpolationMode = NEAREST_NEIGHBOR;
			}
			else if (streq("LINEAR" , theMethod)) {
				config.interpolationMode = LINEAR;
			}
			else if (streq("CUBIC" , theMethod)) {
				config.interpolationMode = CUBIC;
			}
			else {
				cerr << "Unrecognized resampling method: " << theMethod << 
					endl;
				exit(EXIT_FAILURE);
			}
		}

		// Parallelization
		else if (streq("--parallel", argv[c])) {
			config.parallel = true;
		}

		// Default
		else {
			cerr << "Unrecognized option: " << argv[c] << endl;
			exit(EXIT_FAILURE);
		}

		c++;
	}
}
