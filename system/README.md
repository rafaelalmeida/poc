## A Multiclassifier Approach for Different Sensors Image Fusion

This system implements the approach described in my honors project, titled
above.

### Dependencies
* OpenCV
* GDAL
* VLFeat

### Setup

After installing the dependencies and cloning the repository, one needs to
add the IEEE GRSS 2014 Data Fusion Contest subset data into the 
`/system/data/subset` folder. One way to do this is:

```
mkdir data/subset
cd data/subset
wget https://www.dropbox.com/s/kkxogu2vjhbzf36/IEEE_GRSS_IADFTC_Contest2014_subset.zip
unzip IEEE_GRSS_IADFTC_Contest2014_subset.zip
```

To add the full set:
```
mkdir data/fullset
cd data/fullset
https://www.dropbox.com/s/41xykmsmj4qbxde/IEEE_GRSS_IADFTC_Contest2014_full_data_set.zip
unzip IEEE_GRSS_IADFTC_Contest2014_full_data_set.zip
```

### Building and running

To build, simply `make` from the `/system` directory. To run a sample,
simply `make run`.