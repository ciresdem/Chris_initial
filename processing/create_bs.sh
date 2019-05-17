#!/bin/sh -e
function help () {
echo "create_bs- Script that creates a bathy surface at 1 arc-sec for multiple DEM tiles, resamples it back to the target resolution, and converts it xyz. All inputs must be in same, local directory."
	echo "Usage: $0 name_cell_extents datalist coastline bs_res "
	echo "* name_cell_extents: <csv file with name,target spatial resolution in decimal degrees,tile_exents in W,E,S,N>"
	echo "* datalist: <master datalist file that points to individual datasets datalists>"
	echo "* coastline: <coastline shapefile for clipping. Don't include .shp extension >"
}

#see if 3 parameters were provided
#show help if not
if [ ${#@} == 3 ]; 
then
mkdir -p xyz
mkdir -p topo_guide
mkdir -p cmd
mkdir -p save_mb1
mkdir -p save_datalists
mkdir -p tifs
mkdir -p coast_shp
mkdir -p nan_grds

name_cell_extents=$1
datalist_orig=$2
coastline_full=$3

#other pre-set parameters
bs_res=0.00027777777
bs_extend=2

#############################################################################
#############################################################################
#############################################################################
######################      DERIVED VARIABLES     ###########################
#############################################################################
#############################################################################
#############################################################################

#This is used to eventually resample back to target resolution of 1/9 or 1/3rd
ninth_clip_factor=$(echo "((8100+(9*$bs_extend*2))-8112) / 2" | bc )
third_clip_factor=$(echo "((2700+(3*$bs_extend*2))-2712) / 2" | bc )
bs_align_cells=$(echo "$bs_res * $bs_extend" | bc -l)

# Get Tile Name, Cellsize, and Extents from name_cell_extents.csv
IFS=,
sed -n '/^ *[^#]/p' $name_cell_extents |
while read -r line
do
basename_=$(echo $line | awk '{print $1}')
target_res=$(echo $line | awk '{print $2}')
west_tmp=$(echo $line | awk '{print $3}')
east_tmp=$(echo $line | awk '{print $4}')
south_tmp=$(echo $line | awk '{print $5}')
north_tmp=$(echo $line | awk '{print $6}')

#Add on additional cells at bathy_surf resolution to ensure complete coverage of each tile and alignment with 1/9th and 1/3rd res
west=$(echo "$west_tmp - $bs_align_cells" | bc -l)
north=$(echo "$north_tmp + $bs_align_cells" | bc -l)
east=$(echo "$east_tmp + $bs_align_cells" | bc -l)
south=$(echo "$south_tmp - $bs_align_cells" | bc -l)

#if mb1 file already exists for tile, use that. This speeds up processing time if input data files didn't change.
if [ -f $"save_mb1/"$basename_".mb-1" ]; then
	echo "Mb1 file exists, using as datalist"
	cp "save_mb1/"$basename_".mb-1" $basename_".datalist"
	datalist=$(echo $basename_".datalist")
else
	echo "MB1 file doesn't exist, using orig datalist"
	cp $datalist_orig $basename_".datalist"
fi

datalist=$(echo $basename_".datalist")

echo
echo "Tile Name is" $basename_
echo "Cellsize in degrees is" $cellsize_degrees
echo "West is" $west_quarter
echo "East is" $east_quarter
echo "South is" $south_quarter
echo "North is" $north_quarter
echo "Datalist is" $datalist

#Take in a half-cell on all sides so that grid-registered raster edge aligns exactly on desired extent
half_cell=$(echo "$bs_res / 2" | bc -l)
west_reduced=$(echo "$west + $half_cell" | bc -l)
north_reduced=$(echo "$north - $half_cell" | bc -l)
east_reduced=$(echo "$east - $half_cell" | bc -l)
south_reduced=$(echo "$south + $half_cell" | bc -l)

#Determine number of rows and columns with the desired cell size, rounding up to nearest integer.
#i.e., 1_9 arc-second
x_diff=$(echo "$east - $west" | bc -l)
y_diff=$(echo "$north - $south" | bc -l)
x_dim=$(echo "$x_diff / $bs_res" | bc -l)
y_dim=$(echo "$y_diff / $bs_res" | bc -l)
x_dim_int=$(echo "($x_dim+0.5)/1" | bc)
y_dim_int=$(echo "($y_dim+0.5)/1" | bc)

#Target Resolution
#Take orig extents, add on 6 cells at target resolution. All DEM tiles have 6 cell overlap.
six_cells_target=$(echo "$target_res * 6" | bc -l)

west_grdsamp=$(echo "$west_tmp - $six_cells_target" | bc -l)
north_grdsamp=$(echo "$north_tmp + $six_cells_target" | bc -l)
east_grdsamp=$(echo "$east_tmp + $six_cells_target" | bc -l)
south_grdsamp=$(echo "$south_tmp - $six_cells_target " | bc -l)

mb_range_grdsamp="-R$west_grdsamp/$east_grdsamp/$south_grdsamp/$north_grdsamp"

#echo grdsamp range is $mb_range_grdsamp
x_diff_grdsamp=$(echo "$east_grdsamp - $west_grdsamp" | bc -l)
y_diff_grdsamp=$(echo "$north_grdsamp - $south_grdsamp" | bc -l)

x_dim_grdsamp=$(echo "$x_diff_grdsamp / $target_res" | bc -l)
y_dim_grdsamp=$(echo "$y_diff_grdsamp / $target_res" | bc -l)
x_dim_int_grdsamp=$(echo "($x_dim_grdsamp+0.5)/1" | bc)
y_dim_int_grdsamp=$(echo "($y_dim_grdsamp+0.5)/1" | bc)
#echo x_dim_int_grdsamp is $x_dim_int_grdsamp
#echo y_dim_int_grdsamp is $y_dim_int_grdsamp

#bathy surf to target res factor
bs_target_factor=$(echo "$bs_res / $target_res" | bc -l)
bs_target_factor_int=$(echo "($bs_target_factor+0.5)/1" | bc)

#echo bs_target_factor_int is $bs_target_factor_int


#############################################################################
#############################################################################
#############################################################################
######################      Topo Guide     		  ###########################
#############################################################################
#############################################################################
#############################################################################

echo -- Creating Topo Guide... 
#This adds in values of zero to constain interpolation in inland areas without data.
dem_name=$basename_"_topo_guide"
grid_dem=$dem_name".grd"
mb_range="-R$west_reduced/$east_reduced/$south_reduced/$north_reduced"
#create empty datalist
touch dummy.datalist
#echo mb_range is $mb_range
# Run mbgrid
#echo --Running mbgrid...
mbgrid -Idummy.datalist -O$dem_name \
$mb_range \
-A2 -D$x_dim_int/$y_dim_int -G3 -N \
-C0/0 -S0 -F1 -T0.25 -X0.1

# Get Grid Area
echo -- Reclasifying any NaNs to 0s
gdal_calc.py -A $grid_dem --outfile=$grid_dem"_fix_nan.tif" --calc="nan_to_num(A)" --overwrite

echo -- Clipping coastline shp to grid extents
# get the extent of the raster
x_min_tmp=`gmt grdinfo $grid_dem | grep -e "x_min" | awk '{print $3}'`
x_max_tmp=`gmt grdinfo $grid_dem | grep -e "x_max" | awk '{print $5}'`
y_min_tmp=`gmt grdinfo $grid_dem | grep -e "y_min" | awk '{print $3}'`
y_max_tmp=`gmt grdinfo $grid_dem | grep -e "y_max" | awk '{print $5}'`

#Add on 6 more cells just to make sure there is no edge effects when burnining in shp.
x_min=$(echo "$x_min_tmp - $six_cells_target" | bc -l)
x_max=$(echo "$x_max_tmp + $six_cells_target" | bc -l)
y_min=$(echo "$y_min_tmp - $six_cells_target" | bc -l)
y_max=$(echo "$y_max_tmp + $six_cells_target" | bc -l)

ogr2ogr $basename_"_coast.shp" $coastline_full".shp" -clipsrc $x_min $y_min $x_max $y_max

echo -- Setting Topo to -0.1 Prior to Gridding
gdal_rasterize -burn -0.1 -l $basename_"_coast" $basename_"_coast.shp" $grid_dem"_fix_nan.tif"

echo -- Converting to xyz
gdal_translate -of XYZ $grid_dem"_fix_nan.tif" $grid_dem"_rc_tmp.xyz"

echo -- Keeping only negative values
awk '{if ($3 < 0.0) {printf "%.7f %.7f %.2f\n", $1,$2,$3}}' $grid_dem"_rc_tmp.xyz" > $grid_dem"_rc_tmp2.xyz"

echo -- Randomly Sample 100 xyz pnts
shuf -n 100 $grid_dem"_rc_tmp2.xyz" > $grid_dem"_rc.xyz"

rm $grid_dem"_rc_tmp2.xyz"

[ -e topo_guide/$grid_dem"_rc.xyz" ] && rm topo_guide/$grid_dem"_rc.xyz"

mv $grid_dem"_rc.xyz" topo_guide/$grid_dem"_rc.xyz"

rm $grid_dem"_fix_nan.tif"
rm $grid_dem"_rc_tmp.xyz"
rm $grid_dem
rm $dem_name".mb-1"
rm $grid_dem".cmd"

echo -- Creating Datalist for Topo Guide
cd topo_guide/
ls *.xyz > temp
awk '{print $1, 168}' temp > topo_guide.datalist
rm temp
mbdatalist -F-1 -Itopo_guide.datalist -O -V
#
echo
echo "All done"
cd .. 


#############################################################################
#############################################################################
#############################################################################
######################      BATHY SURFACE    		  #######################
#############################################################################
#############################################################################
#############################################################################


#echo -- Creating interpolated DEM for tile $basename_
dem_name=$basename_
grid_dem=$dem_name".grd"
mb_range="-R$west_reduced/$east_reduced/$south_reduced/$north_reduced"
echo mb_range is $mb_range
# Run mbgrid
echo --Running mbgrid...
mbgrid -I$datalist -O$dem_name \
$mb_range \
-A2 -D$x_dim_int/$y_dim_int -G3 -N \
-C810000000/3 -S0 -F1 -T0.25 -X0.1

#Check to see if any valid data
z_min=`gmt grdinfo $grid_dem | grep -e "z_min" | awk '{print $3}'`
z_max=`gmt grdinfo $grid_dem | grep -e "z_max" | awk '{print $3}'`
echo "z_min is" $z_min
echo "z_max is" $z_max




if [ "$z_min" = "NaN" ] && [ "$z_max" = "NaN" ]
then
      	echo "Tile has no data, moving to nan_grds folder..."
      	echo -- Reclasifying any NaNs to 0s
		gdal_calc.py -A $grid_dem --outfile=$grid_dem"_zeroes.tif" --calc="nan_to_num(A)" --overwrite
      	rm $grid_dem
      	mv $grid_dem"_zeroes.tif" nan_grds/$grid_dem"_zeroes.tif"
      	rm $basename_".mb-1"
		rm $basename_"_coast.shp" 
		rm $basename_"_coast.dbf" 
		rm $basename_"_coast.prj" 
		rm $basename_"_coast.shx"
		rm $datalist
else
      	echo "tile has data, keeping..."

      	echo -- Reclasifying any NaNs to 0s
      	gdal_calc.py -A $grid_dem --outfile=$grid_dem"_fix_nan.tif" --calc="nan_to_num(A)" --overwrite
      	
      	echo -- Changing any 0 or above to -0.1m
		gdal_calc.py -A $grid_dem"_fix_nan.tif" --outfile=$grid_dem"_rc.tif" --calc="-0.1*(A >= 0.0)+A*(A < 0.0)" --format=GTiff --overwrite
		rm $grid_dem
		rm $grid_dem"_fix_nan.tif"

		echo -- Resampling to target resolution of 1/9th arc-second or 1/3rd arc-second. 
		# If 1/9th arc-second, I need to take smaller area.
		# If 1/3rd arc-second, I can use entire area.
		echo -- Target Res is $target_res
		echo -- bs_target_factor_int is $bs_target_factor_int
		if [ "$bs_target_factor_int" = 9 ]
		then 
			echo "Resampling to 1/9th Resolution"
			[ -e $grid_dem"_rc_tr_tmp.tif" ] && rm $grid_dem"_rc_tr_tmp.tif"
			gdalwarp $grid_dem"_rc.tif" -r cubicspline -tr $target_res $target_res -t_srs EPSG:4269 $grid_dem"_rc_tr_tmp.tif" -overwrite
			
			#subset image
			[ -e $grid_dem"_rc_tr_tmp2.tif" ] && rm $grid_dem"_rc_tr_tmp2.tif"
			gdal_translate $grid_dem"_rc_tr_tmp.tif" -tr $target_res $target_res -srcwin $ninth_clip_factor $ninth_clip_factor $x_dim_int_grdsamp $y_dim_int_grdsamp -a_srs EPSG:4269 -a_nodata -9999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" $grid_dem"_rc_tr.tif"
			
		else
			echo "Resamping to 1/3rd Resolution"
			#don't need to take smaller area, take full extent
			[ -e $grid_dem"_rc_tr_tmp.tif" ] && rm $grid_dem"_rc_tr_tmp.tif"
			gdalwarp $grid_dem"_rc.tif" -r cubicspline -tr $target_res $target_res -t_srs EPSG:4269 $grid_dem"_rc_tr_tmp.tif" -overwrite
			
			#subset image
			[ -e $grid_dem"_rc_tr_tmp2.tif" ] && rm $grid_dem"_rc_tr_tmp2.tif"
			gdal_translate $grid_dem"_rc_tr_tmp.tif" -tr $target_res $target_res -srcwin $third_clip_factor $third_clip_factor $x_dim_int_grdsamp $y_dim_int_grdsamp -a_srs EPSG:4269 -a_nodata -9999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" $grid_dem"_rc_tr.tif"
			
		fi

		rm $grid_dem"_rc.tif"
		rm $grid_dem"_rc_tr_tmp.tif"

		echo -- Masking out Topo
		gdal_rasterize -burn 1 -l $basename_"_coast" $basename_"_coast.shp" $grid_dem"_rc_tr.tif"

		echo -- Compressing tif
		gdal_translate $grid_dem"_rc_tr.tif" -a_srs EPSG:4269 -a_nodata 999999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" $basename_"_target_res.tif"
		
		rm $grid_dem"_rc_tr.tif"

		#If bathy mask file exists, mask out additional topo
		if [ -f $basename_"_DEM_bathy_1_0.tif" ]; then
			echo "Bathy Mask exists, masking out additional topo areas"
			echo -- Masking out additional topo from DEM
			gdal_calc.py -A $basename_"_target_res.tif" -B $basename_"_DEM_bathy_1_0.tif" --outfile=$basename_"_target_res_final_tmp.tif" --calc="A*B"
			gdal_translate $basename_"_target_res_final_tmp.tif" -a_srs EPSG:4269 -a_nodata -99999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" $basename_"_bs.tif"
			rm $basename_"_target_res_final_tmp.tif"
		else
			echo "Bathy Mask doesn't exist, skipping..."
			mv $basename_"_target_res.tif" $basename_"_bs.tif"
		fi

		echo -- Tiling Bathy Surface and Converting to XYZ
		#Input parameters
		input_file=${basename_}"_bs.tif"

		cp $input_file tifs/$input_file
		mv $basename_".mb-1" save_mb1/$basename_".mb-1"
		mv $basename_"_coast.shp" coast_shp/$basename_"_coast.shp"
		mv $basename_"_coast.dbf" coast_shp/$basename_"_coast.dbf"
		mv $basename_"_coast.prj" coast_shp/$basename_"_coast.prj"
		mv $basename_"_coast.shx" coast_shp/$basename_"_coast.shx"
		mv $datalist save_datalists/$datalist


		#8112 and 2712 are both divisible by 8, create 64 tiles
		#or by 12, create 144 tiles
		tile_x_div=12
		tile_y_div=12

		###################################################################################
		###################################################################################
		#tile_x_pct_tmp=$(echo "100 / $tile_x_pct" | bc -l)
		#tile_y_pct_tmp=$(echo "100 / $tile_y_pct" | bc -l)

		#get input grid dimensions
		x_dim=`gdalinfo $input_file | grep -e "Size is" | awk '{print $3}' | sed 's/.$//'`
		y_dim=`gdalinfo $input_file | grep -e "Size is" | awk '{print $4}'`

		#calculate tile grid dimensions
		tile_dim_x_tmp=$(echo "$x_dim / $tile_x_div" | bc -l)
		tile_dim_x_int=$(echo "($tile_dim_x_tmp+0.5)/1" | bc)

		tile_dim_y_tmp=$(echo "$y_dim / $tile_y_div" | bc -l)
		tile_dim_y_int=$(echo "($tile_dim_y_tmp+0.5)/1" | bc)

		#OR, you can just manually input the number of rows/cols you want
		#tile_dim_x_int=ENTER NUMBER OF COLUMNS HERE
		#tile_dim_y_int=ENTER NUMBER OF ROWS HERE
		#Example; to create 100 tiles when original dimensions are 475100, 61450:
		
		#tile_dim_x_int=1014
		#tile_dim_y_int=1014

		echo input grid is $input_file	

		echo input grid x_dim is $x_dim
		echo input grid y_dim is $y_dim

		echo tile x_dim is $tile_dim_x_int
		echo tile y_dim is $tile_dim_y_int

		echo
		echo -- Starting Tile Analysis
		echo

		#initiate tile names with numbers, starting with 1
		tile_name="1"
		#remove file extension to get basename from input file
		input_name=${input_file::-4}
		#starting point for tiling
		xoff=0
		yoff=0

		while [ "$(bc <<< "$xoff < $x_dim")" == "1"  ]; do
		    yoff=0
		    while [ "$(bc <<< "$yoff < $y_dim")" == "1"  ]; do
		    tile_name_full=$input_name"_tile_"$tile_name".tif"
		    echo creating tile $tile_name_full
		    echo xoff is $xoff
		    echo yoff is $yoff
		    echo tile_dim_x_int is $tile_dim_x_int
		    echo tile_dim_y_int $tile_dim_y_int
		    gdal_translate -of GTiff -a_nodata 999999 -srcwin $xoff $yoff $tile_dim_x_int $tile_dim_y_int $input_file $tile_name_full -stats
		    z_min=`gmt grdinfo $tile_name_full | grep -e "z_min" | awk '{print $3}'`
		    #valid_check=
		    #check to see if tile has data, and detele if no valid values.
		    #valid_check=`gdalinfo $tile_name_full | grep -e "STATISTICS_MAXIMUM"`
		    echo "z_min is" $z_min
		    if (( $(echo "$z_min > 0" | bc -l) ));
			then
			      echo "tile has no data, deleting..."
			      rm $tile_name_full
			else
			      echo "tile has data, keeping..."
			      echo -- Converting to xyz, only keeping negative values
			      gdal_translate -of XYZ $tile_name_full $tile_name_full"_tmp.xyz"
			      awk '{if ($3 < -0.00) {printf "%.7f %.7f %.2f\n", $1,$2,$3}}' $tile_name_full"_tmp.xyz" > $tile_name_full".xyz"
			      echo -- Converted to xyz
			      rm $tile_name_full
			      rm $tile_name_full"_tmp.xyz" 
			      mv $tile_name_full".xyz" xyz/$tile_name_full".xyz"
			fi
			yoff=$(echo "$yoff+$tile_dim_y_int" | bc)
		    tile_name=$((tile_name+1))
		    done
		  xoff=$(echo "$xoff+$tile_dim_x_int" | bc)
		done
		rm $input_file
		echo -- Creating Datalist for Bathy Surface XYZ
		cd xyz/
		ls *.xyz > temp
		awk '{print $1, 168}' temp > bathy_surf.datalist
		rm temp
		mbdatalist -F-1 -Ibathy_surf.datalist -O -V
		#
		echo
		echo "All done"
		cd ..
fi
#else
#	echo file $grid_dem does not exist. Skipping...
#fi

if [ -f $basename_".grd.cmd" ]; then
	echo "cmd file exists, move to subdir"
	mv $basename_".grd.cmd" cmd/$basename_".grd.cmd"
else
	echo "cmd file didn't exist"
fi

done



else
	help
fi

