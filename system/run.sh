#/usr/bin/bash

# run.sh - runs a sample instance for debugging purposes

./bin/classify \
	--verbose \
	--vis "data/subset/TelopsDatasetCityVisible_20cm_Subset.img" \
	--lwir "data/subset/TelopsDatasetCityLWIR_Subset.img" \
	--training "data/subset/TrainingMap_ENVI_RAW_format.raw" \
	--scale-vis 0.1 \
	--scale-lwir 0.1 \
	--segmentation-mode SLIC \
	--slic-region-size 200 \
	--slic-min-region-size 100 \
	--slic-regularization 1000 \
	--slic-auto-scale-parameters \
	--resampling-method NEAREST \
	--log-path "scratch/logs" \
	#--parallel \

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