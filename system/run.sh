#/usr/bin/bash

# run.sh - runs a sample instance for debugging purposes

export LD_LIBRARY_PATH=bin/

if [ -z $1 ]; then
	SCALE_VIS=0.1
else
	SCALE_VIS=$1
fi

if [ -z $2 ]; then
	SCALE_LWIR=0.1
else
	SCALE_LWIR=$2
fi

lldb -- ./bin/classify \
	--verbose \
	--vis "data/subset/TelopsDatasetCityVisible_20cm_Subset.img" \
	--lwir "data/subset/TelopsDatasetCityLWIR_Subset.img" \
	--training "data/subset/TrainingMap_ENVI_RAW_format.raw" \
	--vis-test "data/full/TelopsDatasetCityVisible.img" \
	--lwir-test "data/full/TelopsDatasetCityLWIR.img" \
	--scale-vis $SCALE_VIS \
	--scale-lwir $SCALE_LWIR \
	--segmentation-mode SLIC \
	--slic-region-size 200 \
	--slic-min-region-size 100 \
	--slic-regularization 1000 \
	--slic-auto-scale-parameters \
	--resampling-method NEAREST \
	--log-path "scratch/logs" \
	#--parallel