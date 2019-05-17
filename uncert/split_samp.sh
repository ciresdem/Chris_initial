#!/bin/sh -e

#Verison Change:
#Fixed issue so that it gets errors with distances less than or equal to max interpolation distance (before it was just less than)

basename_=$1
short_name=$2
datalist=$3
cellsize_degrees=$4
num_loops=$5
max_int_dist=$6
ss_samp_den=$7
buff_cells=$8

#From csv
tile_num=$9
west_tile=${10}
east_tile=${11}
south_tile=${12}
north_tile=${13}
x_dim_int_tile=${14}
y_dim_int_tile=${15}
x_win_tile=${16}
y_win_tile=${17}

tile_lat=${18} 
tile_lon=${19} 
tile_zone=${20} 
tile_samp_den=${21} 
tile_zMin=${22} 
tile_zMax=${23} 

#Derived
mb_range_tile="-R$west_tile/$east_tile/$south_tile/$north_tile"
#dem_name_true="tile_"$basename_"_"$tile_num"_dem_true"
dem_name_source="tile_"$basename_"_"$tile_num"_dem_source"

#reduce number of decimal places for ss_samp_den so name in csv isn't really long by rounding and chopping.
ss_samp_den_round_tmp=$(echo "scale=4;  $ss_samp_den / 1 " | bc -l)
ss_samp_den_round=$(printf "%.3f\n" $ss_samp_den_round_tmp)


echo Making outer polygon
gdaltindex $dem_name_source"_interp_outer.shp" $dem_name_source".tif"

echo Making inner polygon
gdal_translate -srcwin $buff_cells $buff_cells $x_win_tile $y_win_tile $dem_name_source".tif" $dem_name_source"_inner.tif"
gdaltindex $dem_name_source"_interp_inner.shp" $dem_name_source"_inner.tif"

echo Converting inner raster to xyz
gdal_translate -of XYZ $dem_name_source"_inner.tif" $dem_name_source"_true_elevs_nan.xyz"
echo -- Removing NaN values
grep -v "\nan" $dem_name_source"_true_elevs_nan.xyz" > $dem_name_source"_true_elevs.xyz"

echo creating ring polygon to extract buffer measurements
python ./ring_poly.py $dem_name_source"_interp_inner.shp" $dem_name_source"_interp_outer.shp" $dem_name_source"_ring.shp"

#create projection file for ring
cp $dem_name_source"_interp_inner.prj" $dem_name_source"_ring.prj"

#get raster values in ring
gdalwarp -cutline $dem_name_source"_ring.shp" -crop_to_cutline $dem_name_source".tif" $dem_name_source"_ring.tif" -overwrite
gdal_translate -of XYZ $dem_name_source"_ring.tif" $dem_name_source"_ring_nan.xyz"
echo -- Removing NaN values
grep -v "\nan" $dem_name_source"_ring_nan.xyz" > $dem_name_source"_ring.xyz"


#Do split-sample anaylsis in loop
for ((n=0;n<$num_loops;n++)); 
do
echo "Starting Loop $n"

#get number of lines
total_pts=$(wc -l < $dem_name_source"_true_elevs.xyz")
total_pts_nan=$(wc -l < $dem_name_source"_true_elevs_nan.xyz")

#initial sample density
intial_samp_den=$(echo "($total_pts/$total_pts_nan)" | bc -l)
intial_samp_pct=$(echo "($intial_samp_den * 100)" | bc -l)
echo intial sampling density is $intial_samp_pct
echo desired sampling density is $ss_samp_den

#desired sample desnity
ss_samp_dec=$(echo "(($ss_samp_den)/100)" | bc -l)

#determine number of points needed for provided sampling density
sampling_desired=$(echo "(($total_pts_nan * $ss_samp_dec))" | bc )
sampling_desired_tmp=$( printf "%.0f" $sampling_desired)
sampling_desired_int=$(echo "($sampling_desired_tmp + 1)" | bc)

echo number of initial elevations is $total_pts
echo number of elevations in sampling desired is $sampling_desired_int
#echo total number of cells in grid is $total_pts_nan

#if initial sampling density is <= desired sampling density, use all the points except 1?
#if not, remove more points so that sampling density equals desired amount
if (( $(echo "$total_pts <= $sampling_desired_int" | bc) )); then
    echo initial sampling density is less than or equal to desired amount 
    echo removing only 1 data point to do splitsample
    sampling=$(echo "($total_pts - 1)" | bc)
else 
    echo initial sampling density is greater than desired amount
    echo reducing sampling further
    sampling=$sampling_desired_int
fi

#
echo orig number of true elevations was $total_pts
echo keeping $sampling of these elevations

# Randomly sample 
shuf -n $sampling $dem_name_source"_true_elevs.xyz" > $dem_name_source"_samples_tmp.xyz"

num_elevs_interp=$(wc -l < $dem_name_source"_samples_tmp.xyz")
echo randomly sampled $num_elevs_interp elevations to use for interpolation

echo adding buffer points to guide interpolation
cat $dem_name_source"_samples_tmp.xyz" $dem_name_source"_ring.xyz" > $dem_name_source"_samples.xyz"


#exit 1


#Create datalist with sample points
dem_int_name="tile_"$basename_"_"$tile_num"_dem_int"
grid_dem_int=$dem_int_name".grd"

echo $dem_name_source"_samples.xyz 168" > $dem_int_name"_interp.datalist"



#make executable
chmod +x $dem_int_name"_interp.datalist"
#reduce grid dimension to avoid edge effect, but use -X to use those data points

# Run mbgrid do on larger area in order for points to match up directly
echo --Running mbgrid...
mbgrid -I$dem_int_name"_interp.datalist" -O$dem_int_name \
        $mb_range_tile \
        -A2 -D$x_dim_int_tile/$y_dim_int_tile -G3 -N \
        -C3/100000000000000 -S0 -F1 -X0.05

grdinfo $dem_int_name".grd"


#clip to inner boundary
gdal_translate $grid_dem_int -a_srs EPSG:4326 $dem_int_name".tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win_tile $y_win_tile $dem_int_name".tif" $dem_int_name"_inner.tif"




gdal_translate -of XYZ $dem_int_name"_inner.tif" $dem_int_name"_interp_nan.xyz"
echo -- Removing NaN values
grep -v "\nan" $dem_int_name"_interp_nan.xyz" > $dem_int_name"_interp.xyz"

#delete files 
#rm $basename_"_scratch"/$dem_int_name".grd"
#rm $basename_"_scratch"/$dem_int_name".tif"
#rm $basename_"_scratch"/$dem_int_name"_inner.tif"


#calulate errors at locations of withheld

#remove lines in samples.xyz in both interpolated and true elevations
awk -F' ' '
NR==FNR {
  a[$1,$2]=$2
  next
}
{
  line=$0
  gsub(/\"/,"")
  gsub(/ *; */,";")
  if (a[$1,$2]!=$2) {
    print line
    line=""
  }
}' $dem_name_source"_samples.xyz" $dem_int_name"_interp.xyz" > $dem_int_name"_interp_no_samples.xyz"


awk -F' ' '
NR==FNR {
  a[$1,$2]=$2
  next
}
{
  line=$0
  gsub(/\"/,"")
  gsub(/ *; */,";")
  if (a[$1,$2]!=$2) {
    print line
    line=""
  }
}' $dem_name_source"_samples.xyz" $dem_name_source"_true_elevs.xyz" > $dem_name_source"_true_no_samples.xyz"






#keep only lines w/ true elevations (i.e., ignore no data values in true elevations that were interpolated)
awk -F' ' '
NR==FNR {
  a[$1,$2]=$2
  next
}
{
  line=$0
  gsub(/\"/,"")
  gsub(/ *; */,";")
  if (a[$1,$2]==$2) {
    print line
    line=""
  }
}' $dem_name_source"_true_no_samples.xyz" $dem_int_name"_interp_no_samples.xyz" > $dem_int_name"_interp_no_samples_no_nans.xyz"



#awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' interp_no_samples.xyz true_elevs.xyz > interp_no_samples_no_nans.xyz
#awk 'FNR==NR{a[$1,$2]=$0;next}{if(b=a[$1,$2]){print b}}' "tb_test_scratch"/interp_no_samples.xyz "tb_test_scratch"/$dem_name_source"_true_elevs.xyz" > "tb_test_scratch"/interp_no_samples_no_nans.xyz

#if number of lines isn't the same, then quit and print error message


num_int_lines=$(wc -l < "$dem_int_name"_interp_no_samples_no_nans.xyz"")
echo number of lines in interpolated file is $num_int_lines

num_true_lines=$(wc -l < "$dem_name_source"_true_no_samples.xyz"")
echo number of lines in true file is $num_true_lines

echo finding matching lat/lons and calulating errors
echo negative errors means interpoled value is lower than true value

awk 'NR==FNR {h[$1 $2] = $3; next} {print $1,$2,$3,h[$1 $2]}' $dem_name_source"_true_no_samples.xyz" "$dem_int_name"_interp_no_samples_no_nans.xyz | awk '{print $1, $2, $3-$4}' > $dem_int_name"_int_errors.xyz"

#TEST, this should work. See what's happening after this.
#paste true_elevs.xyz interp.xyz | awk '{print $1, $2, $6-$3}' > int_errors.xyz  

#calculate distance grid
# Run mbgrid w no interpolation to get location of sample measurements (use larger extent)
echo --Running mbgrid...
samples_name=$basename_"_"$tile_num"_samples"
grid_samples_name=$samples_name".grd"
mbgrid -I$dem_int_name"_interp.datalist" -O$samples_name \
        $mb_range_tile \
        -A2 -D$x_dim_int_tile/$y_dim_int_tile -G3 -N \
        -C0 -S0 -F1


gdal_translate $grid_samples_name -a_srs EPSG:4326 $samples_name".tif"
#gdal_translate tb_test_scratch/tile_tb_test_1_dem_true.grd -a_srs EPSG:4326 tb_test_scratch/tile_tb_test_1_dem_true.tif

gdal_calc.py -A $samples_name".tif" --outfile=$samples_name"_1_0.tif"  --calc="1*(A>-999999999999999999999)" --overwrite

gdal_proximity.py $samples_name"_1_0.tif" $samples_name"_dist.tif"
#gdal_proximity.py tb_test_scratch/tile_tb_test_1_dem_true.tif tb_test_scratch/tile_tb_test_1_dem_true_distance.tif 

echo max_int_dist is $max_int_dist
#extract distance and keep sign of error
#remove errors at distances farther that farthest interpolation distance in grid (i.e., 140 cells in tile C3)
gdal_query.py -delimiter " " -s_format "0,1,2" -d_format "zg" -d_nodata "-99999" $samples_name"_dist.tif" $dem_int_name"_int_errors.xyz" | awk -F'[ ]' -v dist=$max_int_dist '{if ($2 <= dist) {printf "%.2f,%.2f\n", $1,$2}}' > $dem_int_name"_int_errors_dist.xyz"
#add to previous loops
cat $dem_int_name"_int_errors_dist.xyz" >> int_errors_dist_total.csv

echo "Finshed Loop $n";

done

echo Finished All Loops

#Add tile info to csv
echo $tile_num","$tile_lat","$tile_lon","$tile_zone","$tile_samp_den","$tile_zMin","$tile_zMax >> training_tile_info.csv



#rm $dem_name_source"_interp_inner.shp"
#rm $dem_name_source"_interp_inner.dbf"
#rm $dem_name_source"_interp_inner.prj"
#rm $dem_name_source"_interp_inner.shx"

rm $dem_name_source"_interp_outer.shp"
rm $dem_name_source"_interp_outer.dbf"
rm $dem_name_source"_interp_outer.prj"
rm $dem_name_source"_interp_outer.shx"

rm $dem_name_source"_ring.shp"
rm $dem_name_source"_ring.dbf"
rm $dem_name_source"_ring.prj"
rm $dem_name_source"_ring.shx"

rm $dem_name_source"_ring.xyz"
rm $dem_name_source"_ring_nan.xyz"
rm $dem_name_source"_samples.xyz"
rm $dem_name_source"_samples_tmp.xyz"
rm $dem_name_source"_true_elevs.xyz"
rm $dem_name_source"_true_elevs_nan.xyz"
rm $dem_name_source"_true_no_samples.xyz"

rm $dem_name_source"_ring.tif"
rm $dem_name_source"_inner.tif"

#rm $dem_int_name"_errors.xyz"
rm $dem_int_name"_interp.xyz"
rm $dem_int_name"_interp_nan.xyz"
rm $dem_int_name"_interp_no_samples.xyz"
rm $dem_int_name"_interp_no_samples_no_nans.xyz"

rm $dem_int_name"_int_errors.xyz"
rm $dem_int_name"_int_errors_dist.xyz"

rm $grid_samples_name

rm $dem_int_name".grd"
rm $dem_int_name".tif"


rm $samples_name".tif"
rm $samples_name".grd.cmd"
rm $samples_name".mb-1"
rm $samples_name".tif.aux.xml"
rm $samples_name"_1_0.tif"

rm $samples_name"_dist.tif"

rm $dem_name_source".tif"
