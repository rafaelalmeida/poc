#/usr/bin/bash

# run.sh - runs a sample instance for debugging purposes

./bin/classify \
	--verbose \
	--vis "data/subset/TelopsDatasetCityVisible_20cm_Subset.img" \
	--lwir "data/subset/TelopsDatasetCityLWIR_Subset.img" \
	--training "data/subset/TrainingMap_ENVI_RAW_format.raw" \
	--roiX 285 \
	--roiY 384 \
	--roiW 111 \
	--roiH 122 \
	--segmentation-mode GRID \
	--sampling-mode DOWNSAMPLE_VIS \
	--log-path "scratch/logs" \

# --roiX 940 \
# --roiY 2619 \
# --roiW 549 \
# --roiH 675 \