#ifndef COMMON_H
#define COMMON_H

enum ResamplingMethod {
	NEAREST_NEIGHBOR,
	LINEAR,
	CUBIC
};

enum SegmentationMode {
	GRID,
	SLIC
};

// The training matrixes need to be resampled by nearest neighbor 
// interpolation, because it uses gray levels as class labels, so other
// interpolation types will produce values with wrong semantics.
#define TRAINING_INTERPOLATION_MODE cv::INTER_NEAREST

// Hardcoded (sorry) number of bands to keep on dimensionality reduction.
#define LWIR_BANDS_TO_KEEP_ON_PCA 5

// Empirically determined and hardcoded (sorry) SLIC parameters
#define SLIC_REGION_SIZE 30
#define SLIC_REGULARIZATION 1000.
#define SLIC_MIN_REGION_SIZE 10

#endif // COMMON_H
