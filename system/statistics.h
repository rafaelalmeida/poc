#ifndef STAT_H
#define STAT_H

#include <vector>
#include <set>
#include <cmath>
#include <iostream>

#include <opencv2/opencv.hpp>

#include "utils.h"

#define MAX_CATEGORIES 256;

// A K-fold index fold. That is, a pair (T, V) where each member is a 
// collection of indexes that make up the training and validation sets.
typedef std::pair<std::vector<int>, std::vector<int> > KFoldIndexes;

// Class to help in k-fold splits
class KFolder {
	// Number of folds
	int k;

	// The seed used by the random number generator. By saving this explicitly,
	// we can use the same KFolder instance to split many different objects
	// (with the same size) and guarantee that the splits correspond.
	int seed;

	public:
		KFolder(int k, int seed);

		std::list<KFoldIndexes> makeFolds(int collectionSize);
};

namespace statistics {
	double agreement(cv::Mat A, cv::Mat B);
	double kappa(cv::Mat A, cv::Mat B);
	
	std::vector<float> moments(std::vector<float> samples, int maxOrder);
	
	std::vector<float> summary(std::vector<float> samples);
}

#endif