import numpy as np
from osgeo import gdal
import sys
import warnings

#Version Change
#Changed so that if the interpolation distance (95th percentile) is less than 2, make it 2. To make sure I get interpolation errors as a function of distance.
# FUTURE WORK: 
# Make Option so that IF all cells have measurements (i.e., no interpolation), it skips interpolation process.

raster=sys.argv[1]
ds = gdal.Open(raster)
myarray = np.array(ds.GetRasterBand(1).ReadAsArray())
x_dim=myarray.shape[0]
myarray[myarray==0]=np.nan
print myarray
myarray_flat=myarray.flatten()
#print myarray_pos
#p = np.nanpercentile(myarray_flat, 95)

p = np.nanpercentile(myarray_flat, 95)

if p==p:
	print "valid percentile"
	percentile=p
	if percentile < 2:
		print "percentile is less than 2, changing it to 2"
		percentile=2
else:
	print "percentile is nan because every cell has a measurement, changing it to 1"
	percentile=1

print "95th percentile of distance to nearest measurement is", percentile
print percentile









