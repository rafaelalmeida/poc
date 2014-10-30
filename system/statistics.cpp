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