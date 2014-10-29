#include <stdio.h>
#include <iostream>
#include <list>
#include <gdal_priv.h>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

vector<Mat> readLWIR(const char*);
list<Mat> segment(Mat M);
vector<float> statSummary(vector<float> samples);

int main() {
	vector<Mat> matx = readLWIR("data/subset/TelopsDatasetCityLWIR_Subset.img");

	Mat M = matx[0];
	Mat avg = Mat::zeros(M.rows, M.cols, M.type());

	for (int c = 0; c < matx.size(); c++) {
		avg += matx[c];
	}
	avg = avg / matx.size();

	list<Mat> segments = segment(avg);
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

list<Mat> segment(Mat M) {
	// Convert to 8UC3 image
	double max, min;
	minMaxIdx(M, &min, &max);
	Mat tri(M.rows, M.cols, CV_8UC3);
	for (int row = 0; row < M.rows; row++) {
		for (int col = 0; col < M.cols; col++) {
			uchar val = (uchar) (255 * M.at<float>(row, col) / max);
			tri.at<Vec3b>(row, col) = Vec3b(val, val, val);
		}
	}

	// Executes the filtering
	Mat res;
	pyrMeanShiftFiltering(tri, res, 9, 9, 0);

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
	sort(samples.begin(), samples.end());

	float median = samples[samples.size() / 2];
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
	}

	return matx;
}
