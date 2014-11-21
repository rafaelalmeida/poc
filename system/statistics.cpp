#include "statistics.h"

using namespace std;

std::vector<float> statistics::summary(std::vector<float> samples) {
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

double statistics::agreement(cv::Mat A, cv::Mat B) {
	assert(A.rows == B.rows && A.cols == B.cols && A.type() == B.type());

	int agreementCount = 0;
	int pixelsConsidered = 0;

	// Loop through the image to count agreements
	for (int row = 0; row < A.rows; row++) {
		for (int col = 0; col < A.cols; col++) {
			unsigned char labelA = A.at<unsigned char>(col, row);
			unsigned char labelB = B.at<unsigned char>(col, row);

			// Disconsider pixels marked as "not classified" (category id 0)
			if (labelA != 0 && labelB != 0) {
				// Update agreement count
				if (labelA == labelB) {
					agreementCount++;
				}

				// Update pixels considered
				pixelsConsidered++;
			}
		}
	}

	// Calculates Pr(a), the relative agreement observed
	double pr_a = 1.0 * agreementCount / pixelsConsidered;

	return pr_a;
}

double statistics::kappa(cv::Mat A, cv::Mat B) {
	assert(A.rows == B.rows && A.cols == B.cols && A.type() == B.type());

	Counter<unsigned char> labelsA;
	Counter<unsigned char> labelsB;

	int pixelsConsidered = 0;

	// Loop through the image to count the occurence of each category
	for (int row = 0; row < A.rows; row++) {
		for (int col = 0; col < A.cols; col++) {
			unsigned char labelA = A.at<unsigned char>(col, row);
			unsigned char labelB = B.at<unsigned char>(col, row);

			// Disconsider pixels marked as "not classified" (category id 0)
			if (labelA != 0 && labelB != 0) {
				labelsA.inc(labelA);
				labelsB.inc(labelB);

				// Update pixels considered
				pixelsConsidered++;
			}
		}
	}

	// Calculates Pr(e), the hypothetical probability of chance agreement
	double pr_e = 0.0;
	for (auto l : labelsA.getCounts()) {
		double probA = 1.0 * l.second / pixelsConsidered;
		double probB;

		int countB = labelsB.getCount(l.first);
		probB = (float) countB / pixelsConsidered;

		pr_e += probA * probB;
	}

	// Calculates Pr(a), the relative agreement observed
	double pr_a = statistics::agreement(A, B);

	// Calculates the Cohen's Kappa
	double kappa = (pr_a - pr_e) / (1 - pr_e);

	return kappa;
}

std::vector<float> statistics::moments(std::vector<float> samples, int maxOrder) {
	vector<float> moments(maxOrder);

	int n = samples.size();
	for (int i = 1; i <= maxOrder; i++) {
		float m = 0;
		for (auto x : samples) {
			m += pow(x, (float) i) / n;
		}

		moments[i-1] = m;
	}

	return moments;
}

KFolder::KFolder(int k, int seed) : k(k), seed(seed) {}

list<KFoldIndexes> KFolder::makeFolds(int collectionSize) {
	vector<int> folds(collectionSize);

	// Assign the folds
	int foldToGo = 0;
	for (int i = 0; i < collectionSize; i++) {
		folds[i] = foldToGo++;

		// Circle back to first fold
		if (foldToGo == k) {
			foldToGo = 0;
		}
	}
	
	// Shuffle the assigned folds
	auto engine = default_random_engine(seed);
	shuffle(begin(folds), end(folds), engine);

	// Create the fold (TRAINING, VALIDATION) pairs
	list<KFoldIndexes> pairs;
	for (int i = 0; i < k; i++) {
		// Create the vectors that will hold the indexes
		vector<int> trainingIdx;
		vector<int> validationIdx;
		trainingIdx.reserve(collectionSize);
		validationIdx.reserve(collectionSize);

		// Populate the vectors. Each fold is a pair, where the validation 
		// set is composed of the elements in the i-th folder, and the training
		// set by elements in the other folds.
		for (int j = 0; j < collectionSize; j++) {
			if (folds[j] == i) {
				validationIdx.push_back(j);
			}
			else {
				trainingIdx.push_back(j);
			}
		}

		// Add the pair to the list
		pairs.push_back(make_pair(trainingIdx, validationIdx));
	}

	// Return the pairs
	return pairs;
}
