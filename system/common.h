#ifndef COMMON_H
#define COMMON_H

enum ResamplingMethod {
	NEAREST_NEIGHBOR,
	LINEAR,
	CUBIC
};

enum SegmentationMode {
	GRID
};

// The training matrixes need to be resampled by nearest neighbor 
// interpolation, because it uses gray levels as class labels, so other
// interpolation types will produce values with wrong semantics.
#define TRAINING_INTERPOLATION_MODE cv::INTER_NEAREST

#endif // COMMON_H
