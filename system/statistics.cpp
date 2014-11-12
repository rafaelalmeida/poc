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

double statistics::kappa(cv::Mat A, cv::Mat B) {
	assert(A.rows == B.rows && A.cols == B.cols && A.type() == B.type());

	typedef map<unsigned char, int> LabelCountMap;
	LabelCountMap labelsA;
	LabelCountMap labelsB;

	int agreementCount = 0;
	int pixelsConsidered = 0;

	// Loop through the image to count the occurence of each category
	// and the agreements
	for (int row = 0; row < A.rows; row++) {
		for (int col = 0; col < A.cols; col++) {
			unsigned char labelA = A.at<unsigned char>(col, row);
			unsigned char labelB = B.at<unsigned char>(col, row);

			// Disconsider pixels marked as "not classified" (category id 0)
			if (labelA != (unsigned char)0 && labelB != (unsigned char)0) {
				// Update A labels count
				if (labelsA.count(labelA) == 0) {
					labelsA[labelA] = 1;
				}
				else {
					labelsA[labelA]++;
				}

				// Update B labels count
				if (labelsB.count(labelB) == 0) {
					labelsB[labelB] = 1;
				}
				else {
					labelsB[labelB]++;
				}

				// Update agreement count
				if (labelA == labelB) {
					agreementCount++;
				}

				// Update pixels considered
				pixelsConsidered++;
			}
		}
	}

	// Calculates Pr(e), the hypothetical probability of chance agreement
	double pr_e = 0.0;
	for (LabelCountMap::iterator it = labelsA.begin(); it != labelsA.end(); ++it) {
		double probA = 1.0 * it->second / pixelsConsidered;
		double probB;

		LabelCountMap::iterator valB = labelsB.find(it->first);
		if (valB == labelsB.end()) {
			probB = 0;
		}
		else {
			probB = 1.0 * valB->second / pixelsConsidered;
		}

		pr_e += probA * probB;
	}

	// Calculates Pr(a), the relative agreement observed
	double pr_a = 1.0 * agreementCount / pixelsConsidered;

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
