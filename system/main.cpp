#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <stdio.h>

#include <opencv2/core/core.hpp>
#include <opencv2/gpu/gpu.hpp>

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

#include "system.h"

using namespace std;
using namespace cv;

using namespace segmentation;
using namespace classification;

// Main function
int main(int argc, char **argv) {
	// Parse configuration
	Configuration conf;
	config::parse(argv, argc, conf);

	// Start the system
	System system(conf);
	system.run();

	return 0;
}
