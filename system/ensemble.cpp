#include "ensemble.h"

void Ensemble::addClassifier(Classifier c) {
	classifiers.push_back(c);
}

Ensemble::Ensemble(ConsensusType t, Segmentation& s) :
	_consensusType(t),
	_segmentation(s) {

}
