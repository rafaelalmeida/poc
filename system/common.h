#ifndef COMMON_H
#define COMMON_H

// Interpolation modes
enum InterpolationMode {
	NEAREST_NEIGHBOR,
	LINEAR,
	CUBIC
};

// Segmentation modes
enum SegmentationMode {
	GRID,
	SLIC
};

// Image types
enum ImageType {
	VIS,
	LWIR
};

// Number of classes in the thematic map
#define CLASS_COUNT 7

// The training matrixes need to be resampled by nearest neighbor 
// interpolation, because it uses gray levels as class labels, so other
// interpolation types will produce values with wrong semantics.
#define TRAINING_INTERPOLATION_MODE cv::INTER_NEAREST

// Hardcoded (sorry) number of bands to keep on dimensionality reduction.
#define LWIR_BANDS_TO_KEEP_ON_PCA 5

// Some (kind of) reasonable defaults for the SLIC segmentation
#define DEFAULT_SLIC_REGION_SIZE 100
#define DEFAULT_SLIC_MIN_REGION_SIZE 30
#define DEFAULT_SLIC_REGULARIZATION 500.0

// Random seed
#define RANDOM_SEED 13

// Number of k-fold splits
#define K_FOLDS 5

// Number of neighbors to consider in KNN classifier (hardcoded for now)
#define KNN_K 5

// MLP hidden layer size
#define MLP_HIDDEN_LAYER_SIZE 55

// Minimum sample size for node splitting in decision tree classifier
#define DTREE_MIN_SAMPLE_SIZE 5

// ANSI escape sequence to clear the console screen
#define ANSI_CLEAR_SCREEN "\x1b[2J\x1b[1;1H"

// Macro to raise an error
#define FATAL_ERROR(msg) cerr << "Fatal error: " << msg << " (" << \
	__FILE__ << ":"<< __LINE__ << ")" << endl; exit(1);

#endif // COMMON_H
