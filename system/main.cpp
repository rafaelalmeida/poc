#include "classification.h"
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
