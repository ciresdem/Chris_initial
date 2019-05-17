
#VERSION CHANGES
#changed unique id to negative values
#Fixed standard error of the mean to divide by N, not N-1.


import numpy as np
#from numpy import genfromtxt
from numpy import copy
import sys
import numpy as np
import gdal
from gdalconst import *
import osr
import glob
import math
from collections import Counter
import time
import random


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



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
    #print "elev*weight is", GM
    #print "sum of weights is", sum(weights)
    GM /= sum(weights)  # grand mean
    #print "Elevation Mean is", GM

    #Calculating weighted mean of the measurement variances
    mean_var = 0.0
    for i in range(N):
    	mean_var += ((SDs[i])**2) * float(weights[i])/float(sum(weights))

    #Calculating weighted variance of the means (where means = elevation values)
    var_mean = 0.0
    for i in range(N):
    	var_mean += ((elevs[i]-GM)**2) * float(weights[i])/float(sum(weights))

    #Calculating Pooled Sample Variance
    var_samp = ((mean_var+var_mean)*(float(N)/(N-1)))

    return var_samp, GM

#Get data (.csv) from source_uncert.sh
data=sys.argv[1]
tile_name=sys.argv[2]
raster=sys.argv[3]
num_measurements=sys.argv[4]

print "starting data import into array"
my_data = np.loadtxt(data, delimiter=',')
print "created array from data"

data_unique_ids = np.unique(my_data[:,0])
#print "unique ids are", data_unique_ids

#Create empty dictionaries
#d_elev = {}
#d_st_dev = {}
d_st_error = {}

for data_unique_id in data_unique_ids:
    cell_vals = my_data[np.ix_(my_data[:,0] == data_unique_id, (1,2,3))]
    elevs=cell_vals[:,0].tolist()
    weights=cell_vals[:,1].tolist()
    SDs=cell_vals[:,2].tolist()
    if len(elevs) == 1:
        #mean_elev = elevs[0]
        #st_dev = SDs[0]
        st_error = SDs[0]
        #d_elev[data_unique_id] = round(mean_elev, 3)
        #d_st_dev[data_unique_id] = round(st_dev, 3)
        d_st_error[data_unique_id] = round(st_error, 3)
        #print "ONLY 1 DATA POINT"
        #print "id is", data_unique_id
        #print "st_error is st dev ", st_error
    else:
    	#Calculate Sample St Error
        #print "data_unique_id is", data_unique_id
        var_samp, mean_elev = comp_VAR(elevs, SDs, weights)

        st_dev_samp = math.sqrt(var_samp)
        st_error_samp = float(st_dev_samp)/(math.sqrt((len(elevs))))

        #d_elev[data_unique_id] = round(mean_elev, 3)
        #d_st_dev[data_unique_id] = round(st_dev_samp, 3)
        d_st_error[data_unique_id] = round(st_error_samp, 3)

        #Plot the elevations in cell
        #x_val = [1]*len(elevs)

        if len(elevs) >= float(num_measurements):
            x = random.sample(xrange(2*len(elevs)), len(elevs))
            y = elevs
            plt_data=plt.scatter(x,y, zorder=4, label='Measurements')
            plt.show
            error_bars = SDs
            plt_data_uncert=plt.errorbar(x,y,yerr=error_bars, c="red", linestyle="None", zorder=3, label='Measurement Uncertainty')
            plt_st_dev=plt.axhspan(mean_elev-st_dev_samp, mean_elev+st_dev_samp, facecolor='0.5', alpha=0.5, zorder=1, label='St. Dev')
            plt_mean=plt.axhline(y=mean_elev,xmin=0,xmax=100,c="black",linewidth=2,zorder=5, linestyle='--', label='Mean Elevation')
            plt_st_error=plt.axhspan(mean_elev-st_error_samp, mean_elev+st_error_samp, facecolor='0.5', alpha=1, zorder=2, label='St. Error')

            #plt.axes.get_xaxis().set_ticks([])
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # labels along the bottom edge are off
            #handles,labels = ax.get_legend_handles_labels()
            #handles = [handles[0], handles[2], handles[1]]
            #labels = [labels[0], labels[2], labels[1]]

            #plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], loc='upper right')
            plt.legend([plt_data, plt_data_uncert, plt_mean, plt_st_dev, plt_st_error], ['Measurement', 'Measurement Uncert.', 'Mean Elevation', 'St. Deviation', 'St. Error'], loc='upper center', bbox_to_anchor=(0.5, 1.14), fancybox=True, shadow=True, ncol=3, fontsize=12)
            
            #plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
            #plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], loc=8)
            #plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            #plt.legend(loc='upper right')
            plt.ylabel('Elevation (m)')
            #fig,ax = plt.subplots()
            # Plot the data
            #data_line = ax.plot(x,y, label='Data', marker='o')
            # Plot the average line
            #mean_line = ax.plot(x,y_mean, label='Mean', linestyle='--'
            #fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
            #ax.plot(x_val, elevs)
            png_name = str(tile_name)+"_"+str(data_unique_id)+".png"
            plt.savefig(png_name)   # save the figure to file
            plt.close()
            #print "created plot"


        #Calculate Population St Error
        #var_pop = var_samp*(float(len(elevs)-1)/(len(elevs)))
        #st_dev_pop = math.sqrt(var_pop)
        #st_error_pop = float(st_dev_pop)/math.sqrt((len(elevs)))
        #d_st_dev[data_unique_id] = round(st_dev_pop, 3)
        #d_st_error[data_unique_id] = round(st_error_pop, 3)
        #print "MULTIPLE DATA POINTS"
        #print "id is", data_unique_id
        #print "st_error was calculated as ", round(st_error_pop, 3)

#Get raster info from uncert_source.sh (extents,cell size, etc.)

#Make raster of nans
xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
driver = gdal.GetDriverByName('GTiff')
total_cells = xsize*ysize
unique_ids = np.arange(total_cells)
unique_ids_neg = np.negative(unique_ids)
unique_array = unique_ids_neg.reshape(ysize,xsize).astype(float)
print "created unique id array"

d={}
#make copy of unique_id array for each metric of interest (mean elev, st dev, st error, etc)
st_error_array = np.copy(unique_array)
for k, v in d.iteritems(): st_error_array[unique_array==k] = v
print "mapped unique ids to st_error array"
for k, v in d_st_error.iteritems(): st_error_array[unique_array==k] = v
print "mapped standard errors"

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

CreateGeoTiff('st_error',st_error_array, driver, xsize, ysize, GeoT, Projection, DataType)
print "Created Standard Error geotiff"


#Exit python script
sys.exit()
