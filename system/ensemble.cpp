#include "ensemble.h"

void Ensemble::addClassifier(Classifier c) {
	classifiers.push_back(c);
}

Ensemble::Ensemble(ConsensusType t, Segmentation& s, CoverMap training) :
	_consensusType(t),
	_segmentation(s),
	_training(training) {}

void Ensemble::train() {
	for (auto& c : classifiers) {
		c.train(_training);
	}
}

void Ensemble::classify() {
	for (auto& c : classifiers) {
		c.classify();
	}
}