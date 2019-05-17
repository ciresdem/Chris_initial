#!/bin/sh
tile_extents=name_cell_extents_qc.csv

cwd=$(pwd)

# Get Tile Name, Cellsize, and Extents from name_cell_extents.csv
IFS=,
sed -n '/^ *[^#]/p' $tile_extents |
while read -r line
do
basename_=$(echo $line | awk '{print $1}')
descrip=$(echo $line | awk '{print $2}')
cellsize_degrees=$(echo $line | awk '{print $3}')
west_quarter=$(echo $line | awk '{print $4}')
east_quarter=$(echo $line | awk '{print $5}')
south_quarter=$(echo $line | awk '{print $6}')
north_quarter=$(echo $line | awk '{print $7}')
datalist="ga_sc_nc_qc.datalist"
#datalist="dummy.datalist"

full_name=$basename_"_"$descrip

# #if mb1 file exists for tile, use that.
# if [ -f $"save_mb1/"$basename_"_DEM.mb-1" ]; then
# 	echo "Mb1 file exists, using as datalist"
# 	cp "save_mb1/"$basename_"_DEM.mb-1" $basename_".datalist"
# 	datalist=$(echo $basename_".datalist")
# else
# 	echo "MB1 file doesn't exist, using orig datalist"
# 	datalist="e_fl.datalist"
# fi

#cp $basename_".datalist" $basename_"_save.datalist"
#change ned weight values from 0.000000 to 0.000001 
#sed -i -e 's/ 0.000000/ 0.000001/g' $basename_".datalist" 

echo "Datalist is" $datalist
echo

echo
echo "Tile Name is" $basename_
echo "Description is" $descrip
echo "Full name is" $full_name
echo "Cellsize in degrees is" $cellsize_degrees
echo "West is" $west_quarter
echo "East is" $east_quarter
echo "South is" $south_quarter
echo "North is" $north_quarter
echo "Datalist is" $datalist


#############################################################################
#############################################################################
#############################################################################
######################      DERIVED VARIABLES     ###########################
#############################################################################
#############################################################################
#############################################################################
# six_cells_target=$(echo "$cellsize_degrees * 6" | bc -l)
# #half_cell_target=$(echo "$target_res / 2" | bc -l)
# #echo six_cells_target is $six_cells_target

# west=$(echo "$west_quarter - $six_cells_target" | bc -l)
# north=$(echo "$north_quarter + $six_cells_target" | bc -l)
# east=$(echo "$east_quarter + $six_cells_target" | bc -l)
# south=$(echo "$south_quarter - $six_cells_target " | bc -l)

# #Take in a half-cell on all sides so that grid-registered raster edge aligns exactly on desired extent
# half_cell=$(echo "$cellsize_degrees / 2" | bc -l)
# echo half_cell is $half_cell
# west_reduced=$(echo "$west + $half_cell" | bc -l)
# north_reduced=$(echo "$north - $half_cell" | bc -l)
# east_reduced=$(echo "$east - $half_cell" | bc -l)
# south_reduced=$(echo "$south + $half_cell" | bc -l)

# echo "West_reduced is" $west_reduced
# echo "East_reduced is" $east_reduced
# echo "South_reduced is" $south_reduced
# echo "North_reduced is" $north_reduced

# #Determine number of rows and columns with the desired cell size, rounding up to nearest integer.
# #i.e., 1_9 arc-second
# x_diff=$(echo "$east - $west" | bc -l)
# y_diff=$(echo "$north - $south" | bc -l)
# x_dim=$(echo "$x_diff / $cellsize_degrees" | bc -l)
# y_dim=$(echo "$y_diff / $cellsize_degrees" | bc -l)
# x_dim_int=$(echo "($x_dim+0.5)/1" | bc)
# y_dim_int=$(echo "($y_dim+0.5)/1" | bc)


#echo -- Creating interpolated DEM for tile $basename_
grid_dem=$full_name".grd"
#mb_range="-R$west_reduced/$east_reduced/$south_reduced/$north_reduced"
mb_range="-R$west_quarter/$east_quarter/$south_quarter/$north_quarter"
echo mb_range is $mb_range

# Run mbgrid
#echo --Running mbgrid...
# mbgrid -I$datalist -O$full_name \
# $mb_range \
# -A2 -D$x_dim_int/$y_dim_int -G3 -N \
# -C810000000/3 -S0 -F1 -T0.25 -X0.1

echo --Running mbgrid with no interpolation...
mbgrid -I$datalist -O$full_name $mb_range -E$cellsize_degrees"/"$cellsize_degrees"/degrees!" -A2 -G3 -N -C0 -S0 -F1

#echo -- Converting to tif
gdal_translate $grid_dem -a_srs EPSG:4269 -a_nodata -99999 $full_name".tif"
rm $grid_dem
rm $grid_dem".cmd"


# create subdir and copy data used in gridding to subdir
mb1=$full_name".mb-1"
#make dir
new_dir=$basename_"/"$descrip
mkdir -p $new_dir

#mkdir -p "media/sf_external_hd/ga_sc_nc/data/qc/"$basename_
#mkdir -p "media/sf_external_hd/ga_sc_nc/data/qc/"$basename_"/"$descrip
mv $full_name".tif" $new_dir"/"$full_name".tif" 
mv $mb1 $new_dir"/"$mb1

cd $new_dir

mb1=$full_name".mb-1"
# Get Tile Name, Cellsize, and Extents from name_cell_extents.csv
IFS= 
sed -n '/^ *[^#]/p' $mb1 |
while read -r line
do
name_tmp=$(echo $line | awk '{print $1}')
name="${name_tmp:2}"
base=$(basename $name)
echo copying $base
dirname=$(dirname "$name_tmp") 
create_folder_tmp=$(echo $dirname | awk -Fdata '{print $NF}')
create_folder="${create_folder_tmp:1}"
mkdir -p $create_folder
cp $name $create_folder"/"$base
name_xyc="${name::-4}.xyc"
echo "mv $name $name_xyc" >> xyz2xyc.sh
done

cd ..
cd ..
IFS=,

done






