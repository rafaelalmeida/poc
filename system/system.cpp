#include "system.h"

using namespace std;
using namespace cv;

using namespace segmentation;
using namespace classification;

System::System(Configuration config) 
	: _conf(config),
	  _logger(config.logPath) {}

// =============================================================================
// MAIN ENTRY POINT IMPLEMENTATION
// =============================================================================

void System::run() {
	// Start timer
	Stopwatch swatchMain;
	swatchMain.start();

	// Load and setup data
	loadAndSetupData();

	// Reduce dimensionality of LWIR image
	log("reducing LWIR dimensionality...");
	lwir.reduceDimensionality(LWIR_BANDS_TO_KEEP_ON_PCA);

	// Run segmentation phase
	segment();

	// Create descriptors
	describe();

	// Create k-fold splits
	log("creating k-fold splits...");
	KFolder folder(K_FOLDS, RANDOM_SEED);
	vector<pair<ThematicMap, ThematicMap> > splits = trainingMapVIS.split(
		folder);
	vector<pair<Segmentation, Segmentation> > visSegmentationSplits = 
		visTrainingSegmentation.split(folder);
	vector<pair<Segmentation, Segmentation> > lwirSegmentationSplits = 
		lwirTrainingSegmentation.split(folder);

	// Open result file
	ofstream *results = _logger.makeFile("results.txt");
	*results << "MAP / AGREEMENT / KAPPA" << endl;

	// Run k-fold cross validation
	float bestKappa = -numeric_limits<float>::max();
	int bestFold = -1; // Sentinel
	vector<pair<string, Mat> > bestClassifications;
	vector<ThematicMap> foldValidationMaps;
	ThematicMap bestConsensus;

	cerr << "running k-fold cross-validation (k = " << K_FOLDS << ")" << endl;
	for (int fold = 0; fold < K_FOLDS; fold++) {
		// Report progress
		cerr << "running fold " << (fold+1) << " of " << K_FOLDS << endl;

		// Get training and validation maps
		ThematicMap T = splits[fold].first;
		ThematicMap V = splits[fold].second;

		// Resize this fold's training map to the correct sizes
		ThematicMap T_VIS = T.clone();
		ThematicMap T_LWIR = T.clone();
		T_VIS.resize(vis.size());
		T_LWIR.resize(lwir.size());

		// Get VIS and LWIR training segmentations
		Segmentation TS_VIS = visSegmentationSplits[fold].first;
		Segmentation TS_LWIR = segmentationLWIR.cloneWithMask(
			T_LWIR.getFullMask());

		// Show classes counts in this fold
		cerr << "Amount of pixels per class in fold: " << endl;
		printClassHistogram(T_VIS);

		// Build the ensemble
		Ensemble E(MAJORITY_VOTING, segmentationVIS, segmentationLWIR,
			TS_VIS, TS_LWIR, T_VIS, T_LWIR);
		E.setParallel(_conf.parallel);
		setupClassifiers(E, vis, lwir, descriptors);

		// Train the ensemble
		E.train();

		// Run the classification
		ThematicMap C = E.classify();

		// Use the current fold for validation
		foldValidationMaps.push_back(V.asMat());

		// Calculate performance metrics
		Mat X = C.asMat();
		float a = statistics::agreement(V.asMat(), X);
		float k = statistics::kappa(V.asMat(), X);

		// Show statistics
		cerr << "agreement = " << a << ", kappa = " << k << endl;

		// Register time taken
		descriptionTime += E.getTotalDescriptionTime();
		trainingTime += E.getTotalTrainingTime();
		classificationTime += E.getTotalClassificationTime();

		// See if there is improvement
		if (k > bestKappa) {
			bestKappa = k;
			bestFold = fold;
			bestClassifications = E.individualClassifications();
			bestConsensus = C;
		}

		// Separate visually from other folds
		cerr << endl;
	}

	cerr << "best kappa = " << bestKappa << " on fold " << (bestFold) << 
		endl;

	// Process best individual classifications of best performing ensemble
	log("calculating individual statistics...");
	Mat G = foldValidationMaps[bestFold].asMat();
	for (auto c : bestClassifications) {
		// Save results
		Mat X = c.second;

		float agreement = statistics::agreement(G, X);
		float kappa = statistics::kappa(G, X);

		*results << c.first << " " << agreement << " " << kappa 
			<< endl;
	}

	// Calculate consensus kappa
	log("calculating consensus statistics...");

	Mat C = bestConsensus.asMat();
	float agreement = statistics::agreement(G, C);
	float kappa = statistics::kappa(G, C);

	*results << "MAJORITY " << agreement << " " << kappa << endl;

	// Stop timer and log total execution time
	swatchMain.stop();
	double totalTime = swatchMain.read();
	*results << endl << "Total wall time: " << totalTime << endl;

	// Log detailed execution times
	*results << endl << 
		"Description time: " << descriptionTime << endl <<
		"Training time: " << trainingTime << endl <<
		"Classification time: " << classificationTime << endl;

	// Close result file
	results->close();
	delete results;

	// Destroy descriptors
	for (auto d : descriptors) {
		delete d;
	}

	// Print kappa and accuracy on stdout
	cerr << endl;
	cout << agreement << " " << kappa << endl;
}

// =============================================================================
// SYSTEM PHASES IMPLEMENTATION
// =============================================================================

void System::loadAndSetupData() {
	// Load images
	log("loading VIS image...");
	vis = gdal_driver::loadVIS(_conf.pathVIS);
	log("loading LWIR image...");
	lwir = gdal_driver::loadLWIR(_conf.pathLWIR);
	log("loading training data...");
	training = gdal_driver::loadTrainingData(_conf.pathTraining);

	// Rescale images if necessary
	if (_conf.scaleVIS != 1.0 || _conf.scaleLWIR != 1.0) {
		log("rescaling images...");
		rescale(vis, lwir, _conf.scaleVIS, _conf.scaleLWIR, 
			_conf.interpolationMode);
	}

	log("creating training thematic maps...");

	// Create training map objects
	trainingMapVIS = ThematicMap(training);
	trainingMapLWIR = ThematicMap(training);

	// Scale the training map
	trainingMapVIS.resize(vis.size());
	trainingMapLWIR.resize(lwir.size());

	// Set ROI if there is one
	Rect roi;
	if (_conf.roiWidth > 0 && _conf.roiHeight > 0) {
		roi = Rect(_conf.roiX, _conf.roiY, _conf.roiWidth, _conf.roiHeight);
		vis.setRoi(roi);
		training = training(roi);

		log("applying ROIs to LWIR...");
		lwir.setRoi(roi);
	}
}

void System::segment() {
	// Create segmentation - VIS
	log("creating segmentation...");
	if (_conf.segmentationMode == GRID) {
		segmentationVIS = segmentation::segmentVISGrid(vis.asMat(), 
			_conf.gridTileSize);
	}
	else if (_conf.segmentationMode == SLIC) {
		// Try to guess rescaled parameters for SLIC, if desired
		if (_conf.slicAutoScaleParameters) {
			float s = _conf.scaleVIS;

			_conf.slicRegionSize *= s;
			_conf.slicMinRegionSize *= s;
		}

		segmentationVIS = segmentation::segmentVIS_SLIC(vis.asMat(), 
			_conf.slicRegionSize, _conf.slicMinRegionSize, 
			_conf.slicRegularization);
	}
	else {
		FATAL_ERROR("Unsupported segmentation mode");
	}

	// Create segmentation - LWIR
	segmentationLWIR = segmentation::segmentLWIRPixelated(
		lwir, vis);

	// Save segmentation representation
	_logger.saveImage("segmentation", segmentationVIS.representation());

	// Assign images to segmentations
	segmentationVIS.setImage(&vis);
	segmentationLWIR.setImage(&lwir);
}

void System::describe() {
	descriptors.push_back(new GCHDescriptor("GCH"));
	descriptors.push_back(new BICDescriptor("BIC"));
	descriptors.push_back(new LCHDescriptor("LCH"));
	descriptors.push_back(new SIGDescriptor("SIG"));
	descriptors.push_back(new REDUCEDSIGDescriptor("RSIG"));
	descriptors.push_back(new MOMENTSDescriptor("MMT"));

	// Describe the images
	cerr << "describing images..." << endl;
	visTrainingSegmentation = Segmentation(trainingMapVIS);
	visTrainingSegmentation.setImage(&vis);

	list<Segmentation*> segmentations = {&segmentationVIS, &segmentationLWIR,
		&visTrainingSegmentation};

	describeAll(segmentations, descriptors, _conf.parallel);

	// Describe LWIR training set - reuse the main LWIR segmentation, which
	// is already described, and take a subset of it.
	lwirTrainingSegmentation = segmentationLWIR.cloneWithMask(
		trainingMapLWIR.getFullMask());
}

// =============================================================================
// SYSTEM PHASES IMPLEMENTATION
// =============================================================================

void System::rescale(VISImage& vis, LWIRImage& lwir, float scaleVIS, 
	float scaleLWIR, InterpolationMode interpolationMode) {

	// Scale LWIR
	lwir.rescale(scaleLWIR, interpolationMode);

	// Scale VIS
	vis.rescale(scaleVIS, interpolationMode);
}

void System::setupClassifiers(Ensemble& ensemble, VISImage& vis, 
	LWIRImage& lwir, vector<Descriptor*> descriptors) {

	// List classifier engines
	vector<ClassifierEngine> engines = {
		ClassifierEngine::SVM,
		//ClassifierEngine::KNN,
		//ClassifierEngine::DTREE,
		//ClassifierEngine::MLP
	};

	// Create the classifier <-> descriptor pairs and add them to the ensemble
	for (auto engine : engines) {
		for (auto descriptor : descriptors) {
			Classifier *classifier;
			if (descriptor->getType() == ImageType::VIS) {
				classifier = new Classifier(engine, &vis, descriptor);
			}
			else {
				classifier = new Classifier(engine, &lwir, descriptor);
			}

			ensemble.addClassifier(classifier);
		}
	}
}

void System::printClassHistogram(ThematicMap& M) {
	auto counts = M.getClassesCounts();
	for (int i = 0; i < counts.size(); i++) {
		cerr << counts[(unsigned char) i] << " ";
	}
	cerr << endl;
}

void System::describeAll(list<Segmentation*>& segmentations, 
	vector<Descriptor*>& descriptors, bool parallel) {

	// Count how many descriptions will be performed
	int total = 0;
	for (auto& s : segmentations) {
		for (auto& d : descriptors) {
			if (s->getImage()->getType() == d->getType()) {
				total++;
			}
		}
	}

	// Start the description threads
	list<thread> threads;
	int numDone = 0;
	mutex mtx;

	for (auto& s : segmentations) {
		for (auto& d : descriptors) {
			if (s->getImage()->getType() == d->getType()) {
				if (parallel) {
					threads.push_back(thread(&System::doDescribe, this, s, d, 
						&numDone, total, &mtx));
				}
				else {
					doDescribe(s, d, &numDone, total, NULL);
					s->describe(d);
				}
			}
		}
	}

	// If parallel, start the progress reporting procedure
	thread *progressThread = NULL;
	if (parallel) {
		progressThread = new thread(&System::reportDescriptionProcess, this, 
			segmentations, &numDone, total, &mtx);
	}

	// Synchronize
	if (parallel) {
		for (auto& t : threads) {
			t.join();
		}

		if (progressThread != NULL) {
			progressThread->join();
			delete progressThread;
		}
	}

	// Finish progress reporting
	cerr << ANSI_CLEAR_SCREEN << std::flush;
	cerr << "describing: done           " << endl;
}

void System::doDescribe(Segmentation *S, Descriptor *D, int *numDone, int n, 
	mutex *mtx) {

	// If the description is not parallel, we report at the beginning of
	// procedure. If it is, we just update the statistics, and another
	// function will take care of reporting.
	if (mtx == NULL) {
		cerr << "describing: " << *numDone << " of " << n << 
		"       \r" << flush;
	}

	S->describe(D);

	if (mtx != NULL) mtx->lock();
	(*numDone)++;
	if (mtx != NULL) mtx->unlock();
}

void System::reportDescriptionProcess(list<Segmentation*> segmentations, 
	int *numDone, int n, mutex *mtx) {

	while (*numDone < n) {
		mtx->lock();
		/*cerr << "describing: " << *numDone << " of " << n << 
			"       \r" << flush;*/

		// Clear the screen with an ANSI escape sequence
		cerr << ANSI_CLEAR_SCREEN << std::flush;

		int i = 0, n = segmentations.size();
		for (auto S : segmentations) {
			cerr << "Segmentation " << i << " of " << n << endl;
			cerr << "===================" << endl;
			for (auto dp : S->getDescriptionCounts().getCounts()) {
				float progress = 100.0 * dp.second / S->regionCount();
				cerr << dp.first << ": " << progress << "%" << endl;
			}
			cerr << endl;

			i++;
		}

		mtx->unlock();

		this_thread::sleep_for(std::chrono::milliseconds(100));
	}
}

void System::log(const char *msg) {
	if (_conf.verbose) {
		cerr << msg << endl;
	}
}