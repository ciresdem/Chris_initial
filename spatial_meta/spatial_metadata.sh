#!/bin/sh -e
function help () {
echo "xyz2poly- Script that creates polygon shps for each xyz dataset and for each DEM tile"
	echo "Usage: $0 name_cell_extents datalist "
	echo "* name_cell_extents: <csv file with name,target spatial resolution in decimal degrees,tile_exents in W,E,S,N>"
	echo "* datalist: <master datalist file that points to individual datasets datalists>"
}

#################################################################
#see if 2 parameters were provided
#show help if not
if [ ${#@} == 2 ]; 
then
name_cell_extents=$1
datalist=$2

#Split Factor
#This splits the tile into smaller pieces to improve processing speed.
#rows/cols divisible by number below, both 1/9th (8112) and 1/3rd (2712) are divisible by 12
#also by 8 or 6 or 4
split_factor=8

#create folder for entire study area
mkdir -p all_tiles

#create folders for each tile
IFS=,
sed -n '/^ *[^#]/p' $name_cell_extents |
while read -r entire_line_tiles
	do
	tile_name=$(echo $entire_line_tiles | awk '{print $1}')
	mkdir -p $tile_name
	done

#creat datalist for each data set
sed -n '/^ *[^#]/p' $datalist |
while read -r entire_line
	do
	line1=$(echo $entire_line | awk '{print $1}')
	line2=$(basename $line1)
	#dataset_long=${line2::-9}
	#dataset=${dataset_long::10}
	dataset=${line2::-9}
	echo $entire_line > $dataset.datalist
	shp_count=0

	#make directory for each dataset
	mkdir -p $dataset

	mv $dataset.datalist $dataset/$dataset.datalist
	cp $name_cell_extents $dataset/$name_cell_extents
	

	#Go into dataset and make shp for each tile
	cd $dataset/

	IFS=,
	sed -n '/^ *[^#]/p' $name_cell_extents |
	while read -r entire_line_tiles
		do
		#Start Loop to Grid a Data set for Each Tile
		tile_name=$(echo $entire_line_tiles | awk '{print $1}')
		cellsize_degrees=$(echo $entire_line_tiles | awk '{print $2}')
		west_quarter=$(echo $entire_line_tiles | awk '{print $3}')
		east_quarter=$(echo $entire_line_tiles | awk '{print $4}')
		south_quarter=$(echo $entire_line_tiles | awk '{print $5}')
		north_quarter=$(echo $entire_line_tiles | awk '{print $6}')

		#Expand DEM extents by 6 cells to provide overlap between tiles
		six_cells_target=$(echo "$cellsize_degrees * 6" | bc -l)
		#echo six_cells_target is $six_cells_target
		west=$(echo "$west_quarter - $six_cells_target" | bc -l)
		north=$(echo "$north_quarter + $six_cells_target" | bc -l)
		east=$(echo "$east_quarter + $six_cells_target" | bc -l)
		south=$(echo "$south_quarter - $six_cells_target " | bc -l)

		#make subdirectory for each dataset for each tile to split up
		dir_tiles=split_tiles_$tile_name
		mkdir -p $dir_tiles

		#Take in a half-cell on all sides so that grid-registered raster edge aligns exactly on desired extent
		half_cell=$(echo "$cellsize_degrees / 2" | bc -l)
		#echo half_cell is $half_cell
		west_reduced=$(echo "$west + $half_cell" | bc -l)
		north_reduced=$(echo "$north - $half_cell" | bc -l)
		east_reduced=$(echo "$east - $half_cell" | bc -l)
		south_reduced=$(echo "$south + $half_cell" | bc -l)
		#echo "West_reduced is" $west_reduced
		#echo "East_reduced is" $east_reduced
		#echo "South_reduced is" $south_reduced
		#echo "North_reduced is" $north_reduced

		#Determine number of rows and columns with the desired cell size, rounding up to nearest integer.
		#i.e., 1_9 arc-second
		x_diff=$(echo "$east - $west" | bc -l)
		y_diff=$(echo "$north - $south" | bc -l)
		x_dim=$(echo "$x_diff / $cellsize_degrees" | bc -l)
		y_dim=$(echo "$y_diff / $cellsize_degrees" | bc -l)
		x_dim_int=$(echo "($x_dim+0.5)/1" | bc)
		y_dim_int=$(echo "($y_dim+0.5)/1" | bc)

		echo "Creating a Shapefile for Data set: " $dataset 
		echo "For Tile " $tile_name

		dem_name=$dataset"_"$tile_name
		grid_dem=$dem_name".grd"
		mb_range="-R$west_reduced/$east_reduced/$south_reduced/$north_reduced"

		#echo mb_range is $mb_range
		# Run mbgrid
		#echo --Running mbgrid...
		mbgrid -I$dataset.datalist -O$dem_name \
		$mb_range \
		-A2 -D$x_dim_int/$y_dim_int -G3 -N \
		-C0 -S0 -F1

		z_min=`gmt grdinfo $grid_dem | grep -e "z_min" | awk '{print $3}'`
		z_max=`gmt grdinfo $grid_dem | grep -e "z_max" | awk '{print $5}'`

		echo "z_min is" $z_min
		echo "z_max is" $z_max
		#exit 1

		if [ "$z_min" = "NaN" ] && [ "$z_max" = "NaN" ]
		then
		      	echo "tile has no data, deleting..."
		      	rm $grid_dem
		else
		      	echo "tile has data, converting to tif..."
		      	gdal_translate $grid_dem -a_srs EPSG:4269 -a_nodata -9999 -co "COMPRESS=DEFLATE" -co "PREDICTOR=3" -co "TILED=YES" $dem_name".tif"
		      	clip_dir=$(pwd)
		      	gdal_calc.py -A $dem_name".tif" --outfile=$dem_name"_1s_nan.tif"  --calc="1*(A>-999999999999999999999)" --NoDataValue=0 --overwrite
				echo Created source data grid to 1s and NANs
				rm $dem_name".tif"
				rm $grid_dem
				
				#get input grid dimensions
				x_dim=`gdalinfo $dem_name"_1s_nan.tif" | grep -e "Size is" | awk '{print $3}' | sed 's/.$//'`
				y_dim=`gdalinfo $dem_name"_1s_nan.tif" | grep -e "Size is" | awk '{print $4}'`

				#calculate tile grid dimensions
				tile_dim_x_int=$(echo "$x_dim / $split_factor" | bc -l)
				tile_dim_y_int=$(echo "$y_dim / $split_factor" | bc -l)

				echo input grid is $dem_name"_1s_nan.tif"

				echo input grid x_dim is $x_dim
				echo input grid y_dim is $y_dim

				echo tile x_dim is $tile_dim_x_int
				echo tile y_dim is $tile_dim_y_int

				echo
				echo -- Splitting up grid into smaller pieces
				echo

				#initiate tile names with numbers, starting with 1
				split_tile_name="1"
				#remove file extension to get basename from input file
				full_name=${dem_name}"_1s_nan.tif"
				input_name=${full_name::-4}
				#starting point for tiling
				xoff=0
				yoff=0

				while [ "$(bc <<< "$xoff < $x_dim")" == "1"  ]; do
					#Start Splitting
				    yoff=0
				    while [ "$(bc <<< "$yoff < $y_dim")" == "1"  ]; do
				    split_tile_name_full=$input_name"_tile_"$split_tile_name".tif"
				    #echo creating tile $split_tile_name_full
				    #echo xoff is $xoff
				    #echo yoff is $yoff
				    #echo tile_dim_x_int is $tile_dim_x_int
				    #echo tile_dim_y_int $tile_dim_y_int
				    gdal_translate -of GTiff -srcwin $xoff $yoff $tile_dim_x_int $tile_dim_y_int $dem_name"_1s_nan.tif" $split_tile_name_full -stats

				    z_max=`gmt grdinfo $split_tile_name_full | grep -e "z_min" | awk '{print $3}'`
				    if [ "$z_max" = "0" ]
				    then
				        echo "tile has no data, deleting..."
				        rm $split_tile_name_full
				    else
				        echo "tile has data, keeping..."
				        echo "creating shapefile"
				        gdal_polygonize.py $split_tile_name_full -f "ESRI Shapefile" $dir_tiles"/"$split_tile_name_full"_tmp.shp" 
				    	echo Created polygon from raster
				    	rm $split_tile_name_full
				    fi
					yoff=$(echo "$yoff+$tile_dim_y_int" | bc)
				    split_tile_name=$((split_tile_name+1))
				done
				  xoff=$(echo "$xoff+$tile_dim_x_int" | bc)
				#Finishing Splitting
				done
				rm $dem_name"_1s_nan.tif"

			cd $dir_tiles/
			#merge all split tile shps together
			#if there are shps, merge them
			shp_count=$(ls *shp -l | wc -l)
			echo "number of shps is" $shp_count
			if [ $shp_count -gt 0 ];
			then
				echo "Merging all tile shps for " $dataset "and tile" $tile_name
				mkdir -p $dataset"_split_merge/"
				ogrmerge.py -single -o $dataset"_split_merge/"$dem_name".shp" *.shp

				cd $dataset"_split_merge/"
				echo "Dissolving shp"
				ogr2ogr $dem_name"_diss.shp" $dem_name".shp" -dialect sqlite -sql "SELECT DN,ST_Union(geometry) as geometry FROM \"$dem_name\" GROUP BY DN"
				
				echo "Removing DN field"
				ogrinfo $dem_name"_diss.shp" -sql "ALTER TABLE $dem_name"_diss" DROP COLUMN DN"

				cd ..

			    #copy shapefile to directory with specific tile and to the main data set folder
		    	cd ..
		    	cd ..
				#echo path is
				#echo $PWD
				echo "Copying shp to specific tile"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shp" $tile_name/$dem_name".shp"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.dbf" $tile_name/$dem_name".dbf"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shx" $tile_name/$dem_name".shx"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.prj" $tile_name/$dem_name".prj"

		    	echo "Copying shp to specific data set"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shp" $dataset/$dem_name".shp"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.dbf" $dataset/$dem_name".dbf"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shx" $dataset/$dem_name".shx"
		    	cp $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.prj" $dataset/$dem_name".prj"

		    	rm $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shp"
				rm $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.dbf"
				rm $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.shx"
				rm $dataset/$dir_tiles/$dataset"_split_merge/"$dem_name"_diss.prj" 

			else
				echo "No shps, skipping merge..."
				#Move back up to dataset level
				cd ..
			fi
	    	cd $dataset/
	    	rm -rf $dir_tiles
	    	
		fi
	done
	#Finish gridding all tiles for data set, merging shp, then moving on to next dataset.

	#After going through all tiles in dataset, merge all tile shps into 1.
	#if there are shps, merge them
	shp_count=$(ls *shp -l | wc -l)
	echo "number of shps is" $shp_count
	if [ $shp_count -gt 0 ];
	then
		echo "Merging all tile shps for " $dataset
		mkdir -p $dataset"_merge/"
		ogrmerge.py -single -o $dataset"_merge/"$dataset"_merge.shp" *.shp -src_layer_field_name tile
		#Copy shapefile to single directory
		cd .. 
		cp $dataset/$dataset"_merge/"$dataset"_merge.shp" "all_tiles/"$dataset"_merge.shp"
		cp $dataset/$dataset"_merge/"$dataset"_merge.dbf" "all_tiles/"$dataset"_merge.dbf"
		cp $dataset/$dataset"_merge/"$dataset"_merge.shx" "all_tiles/"$dataset"_merge.shx"
		cp $dataset/$dataset"_merge/"$dataset"_merge.prj" "all_tiles/"$dataset"_merge.prj"
		cd $dataset/
	else
		echo "No shps, skipping merge..."
	fi
	echo
	echo "Going on to next dataset"
	cd ..

	rm -rf $dataset
done
	#Finished gridding all datasets


#Post-Processing Merge for Each Tile
#merge all shps for each tile
IFS=,
sed -n '/^ *[^#]/p' $name_cell_extents |
while read -r entire_line_tiles
	do
	tile_name=$(echo $entire_line_tiles | awk '{print $1}')
	echo "Merging shps for tile" $tile_name
	cd $tile_name"/"
	shp_count2=$(ls *shp -l | wc -l)
	echo "number of shps is" $shp_count2
	if [ $shp_count2 -gt 0 ];
	then
		echo "Merging all tile shps for " $tile_name
		rm -rf $tile_name"_merge/"
		mkdir -p $tile_name"_merge/"
		ogrmerge.py -single -o $tile_name"_merge/"$tile_name"_merge.shp" *.shp -src_layer_field_name dataset
		#Copy shapefile to single directory
	else
		echo "No shps, skipping merge..."
	fi
	cd ..
	done


else
	help
fi
