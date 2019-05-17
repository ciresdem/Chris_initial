# Chris_initial
Initial code for Matt to work with
3 subdirs:

-processing
create_dem.sh = script to create DEM
create_bs.sh = script to create bathy surface
qc_archive.sh = script to qc an area of a DEM by copying all data in that area to one location. Could also be useful for copying all data in a DEM tile to one location for archiving purposes.

-spatial_meta
spatial_meta.sh = script to create polygon shp of source data. Merges in two ways; all datasets in a single DEM tile and single dataset across multiple DEM tiles.

-uncert
a whole bunch of scripts used to estimate DEM uncertainty at invidiual cell-level. These scripts were used to generate results in my Journal of Coastal Research (JCR) article (https://www.jcronline.org/doi/pdf/10.2112/JCOASTRES-D-17-00211.1). The main two scripts are:
master_uncert.sh
tile_analysis.sh

