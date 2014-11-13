#/usr/bin/bash

# run.sh - runs a sample instance for debugging purposes

./bin/classify \
	--verbose \
	--vis "data/subset/TelopsDatasetCityVisible_20cm_Subset.img" \
	--lwir "data/subset/TelopsDatasetCityLWIR_Subset.img" \
	--training "data/subset/TrainingMap_ENVI_RAW_format.raw" \
	--scale-vis 1.0 \
	--scale-lwir 1.0 \
	--segmentation-mode SLIC \
	--slic-region-size 100 \
	--slic-min-region-size 75 \
	--slic-regularization 250 \
	--resampling-method NEAREST \
	--log-path "scratch/logs" \
	--parallel \

# TEST ROI - HIGH RES
# --roiX 940 \
# --roiY 2619 \
# --roiW 549 \
# --roiH 675 \

# TEST ROI - LOW RES
# --roiX 285 \
# --roiY 384 \
# --roiW 111 \
# --roiH 122 \