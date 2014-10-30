#include <stdio.h>
#include <iostream>
#include <cstdint>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

vector<Mat> readLWIR(const char*);
Mat readVIS(const char*);
Mat readTrainingData(const char*);
list<Mat> segment(Mat M);
vector<float> statSummary(vector<float> samples);
Mat floatImageTo8UC3Image(Mat floatImage);
void showImage(const char *winname, Mat img, int delay=0);

int main() {
	Mat training = readTrainingData("data/subset/TrainingMap_ENVI_RAW_format.raw");
	Mat vis = readVIS("data/subset/TelopsDatasetCityVisible_20cm_Subset.img");
	vector<Mat> matx = readLWIR("data/subset/TelopsDatasetCityLWIR_Subset.img");
	return 0;

	Mat M = matx[0];
	Mat avg = Mat::zeros(M.rows, M.cols, M.type());

	for (int c = 0; c < matx.size(); c++) {
		avg += matx[c];
	}
	avg = avg / matx.size();

	list<Mat> segments = segment(avg);

	Mat repr = floatImageTo8UC3Image(avg);

	imwrite("lwir.png", repr);
	return 0;

	// Plot contours
	for (list<Mat>::iterator i = segments.begin(); i != segments.end(); ++i)
	{
		vector<vector<Point> > contours;
  		vector<Vec4i> hierarchy;

  		findContours(*i, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
  		drawContours(repr, contours, -1, Scalar(0, 0, 255), 1);
	}

	showImage("win", repr);

	// Calculate and show stats
	vector<float> sizes(segments.size());
	int c = 0;
	for (list<Mat>::iterator i = segments.begin(); i != segments.end(); ++i)
	{
		sizes[c] = countNonZero(*i);
		c++;
	}

	vector<float> stats = statSummary(sizes);
	for (std::vector<float>::iterator i = stats.begin(); i != stats.end(); ++i)
	{
		cout << *i << " ";	
	}
	cout << "\n";

	return 0;
}

Mat floatImageTo8UC3Image(Mat floatImage) {
	Mat M = floatImage;

	double max, min;
	minMaxIdx(M, &min, &max);
	Mat ret(M.rows, M.cols, CV_8UC3);
	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			uchar val = (uchar) (255 * M.at<float>(row, col) / max);
			ret.at<Vec3b>(row, col) = Vec3b(val, val, val);
		}
	}

	return ret;
}

void showImage(const char *winname, Mat img, int delay) {
	namedWindow(winname);
	imshow(winname, img);
	waitKey(delay);
}

list<Mat> segment(Mat M) {
	// Convert to 8UC3 image
	Mat img8uc3 = floatImageTo8UC3Image(M);

	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(img8uc3, res, 9, 9, 0);

	// Separates the segments
	list<Mat> segments;
	RNG rng = theRNG();
	Mat mask(res.rows+2, res.cols+2, CV_8UC1, Scalar::all(0));
	for (int y = 0; y < res.rows; y++)
	{
		for (int x = 0; x < res.cols; x++)
		{
			if (mask.at<uchar>(y+1, x+1) == 0)
			{
				Mat previousMask = mask.clone();

				Scalar newVal(rng(256), rng(256), rng(256));
				floodFill(res, mask, Point(x,y), newVal, 0, Scalar::all(1), Scalar::all(1));

				Mat difference = previousMask ^ mask;
				segments.push_back(difference);
			}
		}
	}

	return segments;
}

vector<float> statSummary(vector<float> samples) {
	// Sort array to calculate some metrics
	sort(samples.begin(), samples.end());

	// Calculate median
	float median;
	if (samples.size() > 0) {
		if (samples.size() % 2 == 1) {
			median = samples[samples.size() / 2];
		}
		else {
			median = (samples[samples.size() / 2 - 1] + 
			          samples[samples.size() / 2]) / 2;
		}
	}
	else {
		median = 0;
	}

	// Calculate quartiles, min and max
	float firstQuartile = samples[samples.size() / 4];
	float thirdQuartile = samples[3 * samples.size() / 4];
	float min = samples[0];
	float max = samples[samples.size() - 1];

	// Calculate mean
	float mean = 0;
	for (vector<float>::iterator i = samples.begin(); i != samples.end(); ++i) {
		mean += *i;
	}
	mean = mean / samples.size();

	// Calculate population standard deviation
	float variance = 0;
	for (vector<float>::iterator i = samples.begin(); i != samples.end(); ++i) {
		float d = (*i - mean);
		variance += d*d;
	}
	variance = variance / samples.size();
	float stddev = sqrt(variance);

	// Organize stats
	vector<float> stats;
	stats.push_back(samples.size());
	stats.push_back(mean);
	stats.push_back(stddev);
	stats.push_back(median);
	stats.push_back(firstQuartile);
	stats.push_back(thirdQuartile);
	stats.push_back(min);
	stats.push_back(max);

	return stats;
}

vector<Mat> readLWIR(const char* path) {
	GDALAllRegister();
	GDALDataset *data = (GDALDataset*) GDALOpen(path, GA_ReadOnly);

	int bands = data->GetRasterCount();
	int width = data->GetRasterXSize();
	int height = data->GetRasterYSize();

	vector<Mat> matx(bands);

	for (int channel = 0; channel < bands; channel++) {
		Mat M(height, width, CV_32FC1);

		GDALRasterBand *band = data->GetRasterBand(channel + 1);
		float *scanline = new float[width];

		for (int row = 0; row < height; row++) {
			band->RasterIO(GF_Read, 0, row, width, 1, scanline, width, 1, GDT_Float32, 0, 0);

			for (int col = 0; col < width; col++) {
				M.at<float>(row, col) = scanline[col];
			}
		}

		matx[channel] = M;
		delete scanline;
	}

	return matx;
}

Mat readVIS(const char* path) {
	GDALAllRegister();
	GDALDataset *data = (GDALDataset*) GDALOpen(path, GA_ReadOnly);

	int bands = data->GetRasterCount();
	int width = data->GetRasterXSize();
	int height = data->GetRasterYSize();

	assert(bands == 3);

	Mat M(height, width, CV_8UC3);
	GDALRasterBand *bandRed = data->GetRasterBand(1);
	GDALRasterBand *bandGreen = data->GetRasterBand(2);
	GDALRasterBand *bandBlue = data->GetRasterBand(3);

	uint16_t *scanlineRed = new uint16_t[width];
	uint16_t *scanlineGreen = new uint16_t[width];
	uint16_t *scanlineBlue = new uint16_t[width];

	for (int row = 0; row < height; row++) {
		bandRed->RasterIO(GF_Read, 0, row, width, 1, scanlineRed, width, 1, GDT_UInt16, 0, 0);
		bandGreen->RasterIO(GF_Read, 0, row, width, 1, scanlineGreen, width, 1, GDT_UInt16, 0, 0);
		bandBlue->RasterIO(GF_Read, 0, row, width, 1, scanlineBlue, width, 1, GDT_UInt16, 0, 0);

		for (int col = 0; col < width; col++) {
			M.at<Vec3b>(row, col) = Vec3b(scanlineBlue[col], scanlineGreen[col], scanlineRed[col]);
		}
	}

	delete scanlineRed;
	delete scanlineGreen;
	delete scanlineBlue;

	return M;
}

Mat readTrainingData(const char* path) {
	GDALAllRegister();
	GDALDataset *data = (GDALDataset*) GDALOpen(path, GA_ReadOnly);

	int bands = data->GetRasterCount();
	int width = data->GetRasterXSize();
	int height = data->GetRasterYSize();

	assert(bands == 1);
	GDALRasterBand *band = data->GetRasterBand(1);

	Mat M(height, width, CV_8UC1);
	unsigned char *scanline = new unsigned char[width];
	for (int row = 0; row < height; row++) {
		band->RasterIO(GF_Read, 0, row, width, 1, scanline, width, 1, GDT_Byte, 0, 0);

		for (int col = 0; col < width; col++) {
			M.at<unsigned char>(row, col) = scanline[col];
		}
	}

	return M;
}
