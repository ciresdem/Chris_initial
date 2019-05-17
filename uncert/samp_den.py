import sys
import numpy as np
import gdal
#from gdalconst import *

#open raster
raster = sys.argv[1]

#loop thru bands of raster and append each band of data to 'layers'
#note that 'ReadAsArray()' returns a numpy array
ds = gdal.Open(raster)
ds_array = ds.ReadAsArray()
constrained = np.sum(ds_array)
print constrained, " cells contrained by measurements"
total_cells = ds_array.size
print total_cells, " total cells in raster"
sampling_d = (constrained/total_cells)*100
print sampling_d, "% sampling density"
print sampling_d
