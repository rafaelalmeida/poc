# -*- coding: utf8 -*-

import cv2
from osgeo import gdal
import numpy as np

def read_dataset(filename):
	return gdal.Open(filename).ReadAsArray()

def read_vis_data(filename):
	return np.rollaxis(read_dataset(filename), 0, 3)
