#!/bin/sh -e
function help () {
echo "create_dem- Script that creates a DEM for multiple DEM tiles and smooths the bathymetry based on user defined parameter to reduce artifacts. IMPORTANT, make sure cell size is 0.00003086420 for 1/9th DEMs in name_cell_extents, otherwise program will just take bathy surface and smooth"
	echo "Usage: $0 name_cell_extents datalist bs_path smooth_factor "
	echo "* name_cell_extents: <csv file with name,target spatial resolution in decimal degrees,tile_exents in W,E,S,N>"
	echo "* datalist: <master datalist file that points to individual datasets datalists>"
	echo "* bs_path: <path to bathy surface tifs for 1/3rd DEMs in format /path/to/bs/tifs NOTE: bathy surface must have same name from name_cell_extents $basename_"_bs.tif">"
	echo "* smooth_factor: <Smoothing factor to reduce bathymetric artifacts, recommended value between 5 and 10 >"
}

#see if 4 parameters were provided
#show help if not
if [ ${#@} == 4 ]; 
then
mkdir -p cmd
mkdir -p save_mb1
mkdir -p save_datalists
mkdir -p tifs
mkdir -p tifs/smoothed

name_cell_extents=$1
datalist_orig=$2
bs_path=$3
smooth_factor=$4

# Get Tile Name, Cellsize, and Extents from name_cell_extents.csv
IFS=,
sed -n '/^ *[^#]/p' $name_cell_extents |
while read -r line
do
basename_=$(echo $line | awk '{print $1}')
target_res=$(echo $line | awk '{print $2}')
west_quarter=$(echo $line | awk '{print $3}')
east_quarter=$(echo $line | awk '{print $4}')
south_quarter=$(echo $line | awk '{print $5}')
north_quarter=$(echo $line | awk '{print $6}')

echo
echo "Tile Name is" $basename_
echo "Cellsize in degrees is" $target_res
echo "West is" $west_quarter
echo "East is" $east_quarter
echo "South is" $south_quarter
echo "North is" $north_quarter
echo
echo "Original Datalist is" $datalist_orig
echo "Bathy Surface Path is" $bs_path
echo "Bathy Smooth Factor is" $smooth_factor
echo

# If the DEM is bathy 1_3 arc-second, just take bathy surface and smooth
if [ $target_res = 0.00003086420 ]; then
	echo "DEM is 1/9th arc-sec"
	echo "Starting Gridding Process"

	#if mb1 file exists for tile, use that.
	if [ -f $"save_mb1/"$basename_"_DEM.mb-1" ]; then
		echo "Mb1 file exists, using as datalist"
		cp "save_mb1/"$basename_"_DEM.mb-1" $basename_"_DEM.datalist"
		datalist=$(echo $basename_"_DEM.datalist")
	else
		echo "MB1 file doesn't exist, using orig datalist"
		cp $datalist_orig $basename_"_DEM.datalist"
	fi
	datalist=$(echo $basename_"_DEM.datalist")

	#############################################################################
	#############################################################################
	#############################################################################
	######################      DERIVED VARIABLES     ###########################
	#############################################################################
	#############################################################################
	#############################################################################
	#Expand DEM extents by 6 cells to provide overlap between tiles
	six_cells_target=$(echo "$target_res * 6" | bc -l)
	#echo six_cells_target is $six_cells_target
	west=$(echo "$west_quarter - $six_cells_target" | bc -l)
	north=$(echo "$north_quarter + $six_cells_target" | bc -l)
	east=$(echo "$east_quarter + $six_cells_target" | bc -l)
	south=$(echo "$south_quarter - $six_cells_target " | bc -l)

	#Take in a half-cell on all sides so that grid-registered raster edge aligns exactly on desired extent
	half_cell=$(echo "$target_res / 2" | bc -l)
	echo half_cell is $half_cell
	west_reduced=$(echo "$west + $half_cell" | bc -l)
	north_reduced=$(echo "$north - $half_cell" | bc -l)
	east_reduced=$(echo "$east - $half_cell" | bc -l)
	south_reduced=$(echo "$south + $half_cell" | bc -l)

	echo "West_reduced is" $west_reduced
	echo "East_reduced is" $east_reduced
	echo "South_reduced is" $south_reduced
	echo "North_reduced is" $north_reduced

	#Determine number of rows and columns with the desired cell size, rounding up to nearest integer.
	#i.e., 1_9 arc-second
	x_diff=$(echo "$east - $west" | bc -l)
	y_diff=$(echo "$north - $south" | bc -l)
	x_dim=$(echo "$x_diff / $target_res" | bc -l)
	y_dim=$(echo "$y_diff / $target_res" | bc -l)
	x_dim_int=$(echo "($x_dim+0.5)/1" | bc)
	y_dim_int=$(echo "($y_dim+0.5)/1" | bc)

	echo -- Creating interpolated DEM for tile $basename_
	dem_name=$basename_"_DEM"
	grid_dem=$dem_name".grd"
	mb_range="-R$west_reduced/$east_reduced/$south_reduced/$north_reduced"
	echo mb_range is $mb_range

	echo --Running mbgrid...
	mbgrid -I$datalist -O$dem_name \
	$mb_range \
	-A2 -D$x_dim_int/$y_dim_int -G3 -N \
	-C810000000/3 -S0 -F1 -T0.25 -X0.1

	# echo --Running mbgrid with no interpolation to test...
	# mbgrid -I$datalist -O$dem_name \
	# $mb_range \
	# -A2 -D$x_dim_int/$y_dim_int -G3 -N \
	# -C0 -S0 -F1

	echo -- Converting to tif
	gdal_translate $grid_dem -a_srs EPSG:4269 -a_nodata -99999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" "tifs/"$basename_"_DEM.tif"

	rm $grid_dem

	echo -- Smoothing Bathy in DEM
	smooth_dem_bathy.py "tifs/"$basename_"_DEM.tif" -s $smooth_factor
	mv "tifs/"$basename_"_DEM_smooth_"$smooth_factor".tif" "tifs/smoothed/"$basename_"_DEM_smooth_"$smooth_factor".tif"

	mv $basename_"_DEM.mb-1" save_mb1/$basename_"_DEM.mb-1"
	mv $datalist save_datalists/$datalist

	if [ -f $basename_"_DEM.grd.cmd" ]; then
		echo "cmd file exists, move to subdir"
		mv $basename_"_DEM.grd.cmd" cmd/$basename_"_DEM.grd.cmd"
	else
		echo "cmd file didn't exist"
	fi

else
	echo "DEM is 1/3 arc-sec Bathy Only"
	echo "Copying Bathy surface to DEM"
	#bs_path=/media/sf_external_hd/ga_sc_nc/data/bathy/bathy_surf/tifs
	cp $bs_path"/"$basename_"_bs.tif" "tifs/"$basename_"_DEM_tmp.tif"

	#Resampling to 1/9th arc-sec to run in arcpy mosaic later
	gdalwarp "tifs/"$basename_"_DEM_tmp.tif" -r cubicspline -tr 0.00003086420 0.00003086420 -t_srs EPSG:4269 "tifs/"$basename_"_DEM.tif" -overwrite
	rm "tifs/"$basename_"_DEM_tmp.tif"
	
	echo "Smoothing DEM"
	smooth_dem_bathy.py "tifs/"$basename_"_DEM.tif" -s $smooth_factor
	mv "tifs/"$basename_"_DEM_smooth_"$smooth_factor".tif" "tifs/smoothed/"$basename_"_DEM_smooth_"$smooth_factor".tif"
	echo
fi
done


else
	help
fi



