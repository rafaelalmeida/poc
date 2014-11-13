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

// Some (kind of) reasonable defaults for the SLIC segmentation
#define DEFAULT_SLIC_REGION_SIZE 100
#define DEFAULT_SLIC_MIN_REGION_SIZE 30
#define DEFAULT_SLIC_REGULARIZATION 500.

#endif // COMMON_H
