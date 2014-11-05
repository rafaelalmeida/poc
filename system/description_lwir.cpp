#include "description_lwir.h"

using namespace std;
using namespace cv;

Mat description_lwir::SIG(LWIRImage lwir, Mat mask) {
	return lwir.spectralSignature(mask);
}
