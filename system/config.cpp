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

		c++;
	}
}
