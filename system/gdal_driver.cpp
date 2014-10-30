#include "gdal_driver.h"

using namespace cv;

std::vector<cv::Mat> gdal_driver::loadLWIR(const char* path) {
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

Mat gdal_driver::loadVIS(const char* path) {
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

Mat gdal_driver::loadTrainingData(const char* path) {
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