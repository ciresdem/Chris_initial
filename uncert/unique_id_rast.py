import numpy as np
import gdal
from gdalconst import *
import osr
import glob
import sys

#Verision Change: Made unique IDs negative so that I can eliminate all unique ids where I don't have data by being less than 0.

def GetGeoInfo(FileName):
    SourceDS = gdal.Open(FileName, GA_ReadOnly)
    xsize = SourceDS.RasterXSize
    ysize = SourceDS.RasterYSize
    GeoT = SourceDS.GetGeoTransform()

    Projection = osr.SpatialReference()
    Projection.ImportFromWkt(SourceDS.GetProjectionRef())

    # Projection = SourceDS.GetProjection()

    DataType = SourceDS.GetRasterBand(1).DataType
    DataType = gdal.GetDataTypeName(DataType)
    return xsize, ysize, GeoT, Projection, DataType


def CreateGeoTiff(Name, Array, driver,
                  xsize, ysize, GeoT, Projection, DataType):
    if DataType == 'Float32':
        DataType = gdal.GDT_Float32
    NewFileName = Name + '.tif'
    # Set up the dataset
    DataSet = driver.Create(NewFileName, xsize, ysize, 1, DataType)
    # the '1' is for band 1.

    DataSet.SetGeoTransform(GeoT)
    DataSet.SetProjection(Projection.ExportToWkt())

    # DataSet.SetProjection(Projection)

    # Write the array
    DataSet.GetRasterBand(1).WriteArray(Array)
    return NewFileName


raster=sys.argv[1]

xsize, ysize, GeoT, Projection, DataType = GetGeoInfo(raster)
# Set up the GTiff driver
driver = gdal.GetDriverByName('GTiff')

#print "number of x columns is", xsize
#print "number of y columns is", ysize

total_cells = xsize*ysize
unique_ids = np.arange(total_cells)
unique_ids_neg = np.negative(unique_ids)
#print "number of cells is", len(unique_ids_neg)
#print unique_ids_neg
unique_array = unique_ids_neg.reshape(ysize,xsize).astype(int)
#print "created array"
#print unique_array

CreateGeoTiff('unique_id',unique_array, driver, xsize, ysize, GeoT, Projection, GDT_Int32)

print "Created geotiff"




