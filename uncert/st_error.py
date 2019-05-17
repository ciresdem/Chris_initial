
#VERSION CHANGES
#Changed back to 1 sigma
#Created weight raster
#Fixed standard error of the mean to divide by N, not N-1.

import numpy as np
#from numpy import genfromtxt
from numpy import copy
import sys
import numpy as np
import gdal
from gdalconst import *
import osr
import math
from collections import Counter

def GetGeoInfo(FileName):
    SourceDS = gdal.Open(FileName, GA_ReadOnly)
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()
    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())
    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    return xsize, ysize, GeoT, Projection, DataType

def CreateGeoTiff(Name, Array, driver,
                  xsize, ysize, GeoT, Projection, DataType):
    if DataType == 'Float32':
        DataType = gdal.GDT_Float32
    NewFileName = Name + '.tif'
    DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projection.ExportToWkt())
    DataSet.GetRasterBand(1).WriteArray(Array)
    return NewFileName

def comp_VAR(elevs, SDs, weights):
    '''Calculate the pooled sample variance, following https://arxiv.org/ftp/arxiv/papers/1007/1007.1012.pdf
       Inputs are:
         elevs, the elevation value
         SDs, the measurement uncertainty at 1 standard deviation (from metadata), i.e., population standard deviation
         weights, the weighting value given to the datalist that contains the elevation
       Result is the pooled sample variance
    '''
    # calculate total number of samples, N, and grand mean, GM
    N = len(elevs)  # total number of samples
    GM = 0.0
    for i in range(N):
        GM += elevs[i] * weights[i]
    #print "elev is", elevs[i]
    #print "weight is", weights[i]
    #print "elev*weight is", G
    #print "sum of weights is", sum(weights)
    GM /= sum(weights)  # grand mean
    #print "Elevation Mean is", GM

    #Calculating weighted mean of the measurement variances
    mean_var = 0.0
    for i in range(N):
    	mean_var += (SDs[i]**2) * float(weights[i])/float(sum(weights))

    #Calculating weighted variance of the means (where means = elevation values)
    var_mean = 0.0
    for i in range(N):
    	var_mean += ((elevs[i]-GM)**2) * float(weights[i])/float(sum(weights))

    #Calculating Pooled Sample Variance
    var_samp = ((mean_var+var_mean)*(float(N)/(N-1)))

    return var_samp, GM, sum(weights)/float(N)

#Get data (.csv) from source_uncert.sh
data=sys.argv[1]
tile_name=sys.argv[2]
raster=sys.argv[3]

print "starting data import into array"
my_data = np.loadtxt(data, delimiter=',')
print "created array from data"

#See if there are more than one data points in my_array
try:
    print(my_data.shape[0])
    #If this works, then more than one data point
    data_unique_ids = np.unique(my_data[:,0])
    #print "unique ids are", data_unique_ids
    #Create empty dictionaries
    #d_elev = {}
    #d_st_dev = {}
    d_st_error = {}
    d_weights= {}

    for data_unique_id in data_unique_ids:
        cell_vals = my_data[np.ix_(my_data[:,0] == data_unique_id, (1,2,3))]
        elevs=cell_vals[:,0].tolist()
        weights=cell_vals[:,1].tolist()
        SDs=cell_vals[:,2].tolist()
        if len(elevs) == 1:
            #print "ONLY 1 DATA POINT IN CELL"
            st_error = SDs[0]
            d_st_error[data_unique_id] = round(st_error, 3)
            avg_weight = weights[0]
            d_weights[data_unique_id] = avg_weight

            #mean_elev = elevs[0]
            #st_dev = SDs[0]
            #d_elev[data_unique_id] = round(mean_elev, 3)
            #d_st_dev[data_unique_id] = round(st_dev, 3)
        else:
            #print "MULTIPLE DATA POINTS IN CELL"
            var_samp, mean_elev, avg_weight = comp_VAR(elevs, SDs, weights)
            st_dev_samp = math.sqrt(var_samp)
            st_error_samp = float(st_dev_samp)/(math.sqrt((len(elevs))))
            d_st_error[data_unique_id] = round(st_error_samp, 3)
            d_weights[data_unique_id] = avg_weight

            #d_elev[data_unique_id] = round(mean_elev, 3)
            #d_st_dev[data_unique_id] = round(st_dev_samp, 3)
except IndexError:
    #print("only 1 data point in ARRAY")
    d_st_error = {}
    d_weights = {}
    data_unique_id=my_data[0]
    st_error = my_data[3]
    d_st_error[data_unique_id] = round(st_error, 3)
    avg_weight = my_data[2]
    d_weights[data_unique_id] = avg_weight



#Get raster info from uncert_source.sh (extents,cell size, etc.)
#Make raster of unique ids
xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
driver = gdal.GetDriverByName('GTiff')
total_cells = xsize*ysize
unique_ids = np.arange(total_cells)
unique_ids_neg = np.negative(unique_ids)
unique_array = unique_ids_neg.reshape(ysize,xsize).astype(float)
#print unique_array
print "created unique id array"

d={}
#make copy of unique_id array for each metric of interest (mean elev, st dev, st error, etc)
st_error_array = np.copy(unique_array)
for k, v in d.iteritems(): st_error_array[unique_array==k] = v
print "mapped unique ids to st_error array"
for k, v in d_st_error.iteritems(): st_error_array[unique_array==k] = v
print "mapped standard errors"

CreateGeoTiff('st_error',st_error_array, driver, xsize, ysize, GeoT, Projection, DataType)
print "Created Standard Error geotiff"

weights_array = np.copy(st_error_array)
#print "mapped unique ids to weights array"
for k, v in d_weights.iteritems(): weights_array[unique_array==k] = v
print "mapped weights"

CreateGeoTiff('weights',weights_array, driver, xsize, ysize, GeoT, Projection, DataType)
print "Created Weights geotiff"


#st_dev_array = np.copy(st_error_array)
#print "mapped unique ids to st_array array"
#for k, v in d_st_dev.iteritems(): st_dev_array[unique_array==k] = v
#print "mapped standard deviations"

#mean_elev_array = np.copy(st_error_array)
#print "mapped unique ids to mean_elev array"
#for k, v in d_elev.iteritems(): mean_elev_array[unique_array==k] = v
#print "mapped mean elevs"

#CreateGeoTiff('mean_elev',mean_elev_array, driver, xsize, ysize, GeoT, Projection, DataType)
#print "Created Mean Elevation geotiff"

#CreateGeoTiff('st_dev',st_dev_array, driver, xsize, ysize, GeoT, Projection, DataType)
#print "Created Standard Deviation geotiff"

#Exit python script
sys.exit()
