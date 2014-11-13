#include "segmentation.h"

using namespace cv;
using namespace std;

using namespace segmentation;

Segmentation segmentation::segmentVIS_SLIC(Mat M) {
	// Adapted from https://github.com/davidstutz/vlfeat-slic-example/blob/
	// 79bc823000ec1bc7fb204b0613578e13d3dff1ae/vlfeat_slic_cli/main.cpp

	assert(M.channels() == 3);

	// Create a binary mask to show where we have missing data and ignore
	// those pixels later
	Mat gray(M.size(), CV_8UC1);
	cvtColor(M, gray, CV_BGR2GRAY);

	Mat bin(gray.size(), CV_8UC1);
	threshold(gray, bin, 1, 255, CV_THRESH_BINARY);

	// Convert image to one-dimensional array.
	log("converting image format...");
	float* image = new float[M.rows*M.cols*M.channels()];
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.cols; ++j) {
			image[j + M.cols*i + M.cols*M.rows*0] = M.at<cv::Vec3b>(i, j)[0];
			image[j + M.cols*i + M.cols*M.rows*1] = M.at<cv::Vec3b>(i, j)[1];
			image[j + M.cols*i + M.cols*M.rows*2] = M.at<cv::Vec3b>(i, j)[2];
		}
	}

	// The algorithm will store the final segmentation in a one-dimensional 
	// array.
	vl_uint32* segmentation = new vl_uint32[M.rows*M.cols];
	vl_size height = M.rows;
	vl_size width = M.cols;
	vl_size channels = M.channels();

	// The region size defines the number of superpixels obtained.
	// Regularization describes a trade-off between the color term and the
	// spatial term.
	vl_size region = 100;
	float regularization = 1000.;
	vl_size minRegion = 50;

	log("running segmentation...");
	vl_slic_segment(segmentation, image, width, height, channels, region, 
		regularization, minRegion);
            
    // Convert segmentation.
    log("converting segmentation...");
    int** labels = new int*[M.rows];
    for (int i = 0; i < M.rows; ++i) {
        labels[i] = new int[M.cols];
                
        for (int j = 0; j < M.cols; ++j) {
            labels[i][j] = (int) segmentation[j + M.cols*i];
        }
    }

    // Free memory used by float image and segmentation (we already converted
    // it above)
    delete [] image;
    delete [] segmentation;

    // Finds the total number of regions. The elements of the labels
    // matrix are the region labels starting from zero.
    int maxLabel = 0;
    for (int i = 0; i < M.rows; i++) {
    	for (int j = 0; j < M.cols; j++) {
    		int l = labels[i][j];
    		if (l > maxLabel) {
    			maxLabel = l;
    		}
    	}
    }

    int totalRegions = maxLabel + 1;

    // Create the segment masks
    vector<SparseMat> masks;
    for (int i = 0; i < totalRegions; i++) {
    	int sizes[2] = {M.rows, M.cols};
    	SparseMat SM(2, sizes, CV_8UC1);
    	masks.push_back(SM);
    }

    // Fill the segment masks where original data is non missing
    Mat nonMissingPixels = segmentation::makeNonMissingDataMask(M);
    for (int i = 0; i < M.rows; i++) {
    	for (int j = 0; j < M.cols; j++) {
    		if (nonMissingPixels.at<unsigned char>(i, j) == 255) {
    			int label = labels[i][j];
    			*(masks[label].ptr(i, j, true)) = 255;
    		}
    	}
    }

    // Delete the labels array
    for (int i = 0; i < M.rows; i++) {
    	delete [] labels[i];
    }
    delete [] labels;

    // Convert the mask vector into a list, filtering empty segments
    // (segments that were created by SLIC but on which all pixels are missing)
    list<SparseMat> masksList;
    for (auto m : masks) {
    	if (m.nzcount() > 0) {
    		masksList.push_back(m);
    	}
    }

    // Create the segmentation and return it
    return Segmentation(masksList);
}

Segmentation segmentation::segmentVISGrid(cv::Mat M, int tileSize) {
	list<SparseMat> segments;

	int regionsPerLine = (M.cols / tileSize) + (M.cols % tileSize);
	int regionsPerColumn = (M.rows / tileSize) + (M.rows % tileSize);

	Mat nonMissing = makeNonMissingDataMask(M);

	for (int i = 0; i < regionsPerLine; i++) {
		for (int j = 0; j < regionsPerColumn; j++) {
			Mat segment = Mat::zeros(M.size(), CV_8UC1);
			Rect region(i*tileSize, j*tileSize, tileSize, tileSize);

			rectangle(segment, region, Scalar(255), CV_FILLED);

			Mat filteredSegment = segment & nonMissing;
			if (countNonZero(filteredSegment) > 0) {
				segments.push_back(SparseMat(filteredSegment));
			}
		}
	}

	return Segmentation(segments);
}

Segmentation segmentation::segmentLWIRPixelated(LWIRImage& lwir, Mat vis) {
	Mat nonMissing = makeNonMissingDataMask(vis);
	Mat M;
	resize(nonMissing, M, lwir.size(), 0, 0, INTER_NEAREST);

	list<SparseMat> pixels;
	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			if (M.at<unsigned char>(row, col) == 255) {
				int sizes[2] = {M.rows, M.cols};
				SparseMat pixel(2, sizes, CV_8UC1);

				*(pixel.ptr(row, col, true)) = 255;
				pixels.push_back(pixel);
			}
		}
	}

	return Segmentation(pixels);
}

// Create a binary mask to show where we have missing data and ignore
// those pixels later
Mat segmentation::makeNonMissingDataMask(Mat vis) {
	Mat gray(vis.size(), CV_8UC1);
	cvtColor(vis, gray, CV_BGR2GRAY);

	Mat bin(gray.size(), CV_8UC1);
	threshold(gray, bin, 1, 255, CV_THRESH_BINARY);

	// TODO? fill holes in mask to account for pixels inside non missing region
	// that have zero value

	return bin;
}

Segmentation segmentation::segmentVISWatershed(Mat M, int tileSize) {
	// Build some rough markers using the grid tile centers as seed points
	int regionsPerLine = (M.cols / tileSize);
	int regionsPerColumn = (M.rows / tileSize);

	Mat markers(M.size(), CV_32S);
	int currentLabel = 1;

	for (int i = 0; i < regionsPerLine; i++) {
		for (int j = 0; j < regionsPerColumn; j++) {
			Point P(i*tileSize + tileSize/2, j*tileSize + tileSize/2);
			markers.at<int>(i, j) = currentLabel++;
		}
	}

	// Notes the total number of regions
	int totalRegions = currentLabel - 1;

	// Run the watershed
	watershed(M, markers);

	// Creates a list of region masks
	list<SparseMat> regionMasks;
	for (int i = 1; i <= totalRegions; i++) {
		float progress = (float) 100 * i / totalRegions;
		cerr << progress << "%          \r" << flush;

		Mat thisMask(markers.size(), CV_8UC1);
		for (int row = 0; row < markers.rows; row++) {
			for (int col = 0; col < markers.cols; col++) {
				if (markers.at<int>(row, col) == i) {
					thisMask.at<unsigned char>(row, col) = 255;
				}
			}
		}

		regionMasks.push_back(SparseMat(thisMask));
	}

	return Segmentation(regionMasks);
}

float segmentation::getSegmentLabel(Mat classificationMap, Mat mask) {
	assert(classificationMap.type() == CV_8UC1);
	assert(mask.type() == CV_8UC1);
	assert((classificationMap.rows == mask.rows) && (classificationMap.cols = mask.cols));

	return (float) mean(classificationMap, mask)[0];
}

list<SparseMat> segmentation::getColorBlobs(Mat posterized) {
	list<SparseMat> segments;

	Mat clone = posterized.clone();
	Mat mask(posterized.rows+2, posterized.cols+2, CV_8UC1, Scalar::all(0));

	for (int y = 0; y < clone.rows; y++) {
		for (int x = 0; x < clone.cols; x++) {
			// Report progress
			if ((y*clone.rows + x) % 1000 == 0) {
				int progress = 100.0 * ((y*clone.cols + x)) / 
					(clone.rows * clone.cols);

				cerr << "finding blobs... " << progress << 
					"%...     " << "\r" << flush;
			}

			// Find connected components
			if (mask.at<uchar>(y+1, x+1) == 0) {
				Mat previousMask = mask.clone();

				floodFill(clone, mask, Point(x,y), Scalar(0, 0, 0), 0, 
						  Scalar::all(0), Scalar::all(0));

				Mat difference = previousMask ^ mask;
				Mat segmentMask = 255 * difference(Range(1, difference.rows - 1),
												   Range(1, difference.cols - 1));

				segments.push_back(SparseMat(segmentMask));
			}
		}
	}

	cerr << "finding blobs... done" << endl;

	return segments;
}

Segmentation::Segmentation(list<SparseMat> masks) {
	this->_masks = masks;
}

list<SparseMat> Segmentation::getSegments() {
	return _masks;
}

cv::Mat Segmentation::representation() {
	assert(_masks.size() > 0);

	Mat repr = Mat::zeros(getMapSize(), CV_8UC3);
	for (auto mask : _masks) {
		Mat clone = densify(mask);
		vector<vector<Point> > contours;
		findContours(clone, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);

		drawContours(repr, contours, -1, Scalar(0, 0, 255));
	}

	return repr;
}

cv::Size Segmentation::getMapSize() {
	return densify(_masks.front()).size();
}

vector<SparseMat> segmentation::pixelSegmentation(Mat image) {
	vector<SparseMat> segments;
	segments.reserve(image.rows * image.cols);

	for (int row = 0; row < image.rows; row++) {
		for (int col = 0; image.cols; col++) {
			Mat mask(image.size(), CV_8UC1);
			mask.at<unsigned char>(row, col) = 255;

			SparseMat sparseMask(mask);

			segments.push_back(sparseMask);
		}
	}

	return segments;
}

Segmentation Segmentation::pixelize() {
	Mat fullMask = Mat::zeros(this->getMapSize(), CV_8UC1);
	for (auto m : _masks) {
		fullMask += densify(m);
	}

	list<SparseMat> pixelatedSegments;

	for (int row = 0; row < fullMask.rows; row++) {
		for (int col = 0; col < fullMask.cols; col++) {
			if (fullMask.at<unsigned char>(row, col) != 0) {
				Mat mask = Mat::zeros(fullMask.size(), CV_8UC1);
				mask.at<unsigned char>(row, col) = 255;

				pixelatedSegments.push_back(SparseMat(mask));
			}
		}
	}

	return Segmentation(pixelatedSegments);
}

int Segmentation::segmentCount() {
	return _masks.size();
}