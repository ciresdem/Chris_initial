#!/bin/sh

#Version Changes:
# Fixed Issue with 5th percentile of sampling density when there is only 1 tile. Now will default to that 1 tile sampling density. 

# Get Tile Name, Cellsize, and Extents from tile_extents_gridding.txt
name_cell_extents=name_cell_extents_JCR.txt
sed -n '/^ *[^#]/p' $name_cell_extents |
while read -r line
do
basename_=$(echo $line | awk '{print $1}')
cellsize_degrees=$(echo $line | awk '{print $2}')
west=$(echo $line | awk '{print $3}')
east=$(echo $line | awk '{print $4}')
south=$(echo $line | awk '{print $5}')
north=$(echo $line | awk '{print $6}')

echo
echo
echo
echo
echo
echo
echo
echo
echo


####################################################################
####################################################################
########################                  ##########################
########################    ENTER USER    ##########################
########################    PARAMETERS    ##########################
########################       HERE       ##########################
########################                  ##########################
####################################################################
####################################################################

#Change other parameters below

#########################################################################################
#########     NUMBER OF SPLIT-SAMPLE LOOPS TO QUANTIFY INTERPOLATION ERROR    ###########
#########################################################################################
#num_loops=500
#num_loops=250
#num_loops=100
num_loops=50
#num_loops=3
####################################################################
####################     TILE size                       ###########
####################################################################
#Correction factor for tile size
#i.e. 1=95th percentile of distance to nearest measurement
#2= 2 * 95th perctile of distance to nearest measurement = 2 times the tile size dimensions, 4 times larger
tile_factor=4
#tile_factor=1
####################################################################
###############     Training                   #####################
####################################################################
#Final number of training tiles per zone (i.e., Bathy, BathyTopo, Topo), so max will be num_training * 3.
num_training=25
#num_training=3
#minimum percentile of sampling density to be considered a potential training
#50th percentile, or median sampling density is a good starting point
#If there is little data in zone, such as bathy, increase to higher percentile to only use tiles with more data for training
#may need to tweek, see what percentile results in a certain area of data (i.e. strip of multibeam) being used as training
min_samp_pct_bathy=50
min_samp_pct_bathytopo=50
min_samp_pct_topo=50
####################################################################
########################     DATALIST    ###########################
####################################################################
#datalist="sw_fl_C3_1_sigma_clip.datalist"
#datalist="sw_fl_95_conf.datalist"
#datalist="tb_2015.datalist"
#datalist="nos_test.datalist"
#datalist="sw_fl_equal_w.datalist"
#datalist="sw_fl_equal_w_no_uncert.datalist"
#datalist="no_uncert.datalist"
#datalist="test.datalist"
#datalist="topo_test.datalist"
datalist="sw_fl_1_sigma_clip_JCR_diff_w.datalist"
####################################################################
##############             DEM GENERATION           ################
####################################################################
#make interpolated DEM in addition to uncertainty analysis; 0=no, 1=yes
#DEM_generation=0
DEM_generation=1
####################################################################
###################    VDATUM UNCERTAINTY    #######################
####################################################################
#vdatum uncertainty at 1 sigma 
vdatum_uncert=0.12
####################################################################
###################    DATA BUFFER     #############################
####################################################################
# number of cells to expand grid to more accurately calculate interpolation uncertainty (distance to measurement)
buff_cells=1
####################################################################
##############             PLOTS                    ################
####################################################################
#make plots of measurement uncertainty; 0=no, 1=yes
plots=0
#only make plots where cells have more than this number of measurements
num_measurements=1000
####################################################################
####################################################################
########################                  ##########################
########################    END OF USER   ##########################
########################    PARAMETERS    ##########################
########################                  ##########################
####################################################################
####################################################################

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

####################################################################
######################                        ######################
######################   Derived Variables    ######################
######################                        ######################
####################################################################

###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
echo
if (( $(echo "$DEM_generation == 0" | bc -l) )); then
    echo No interpolated DEM will be generated
    echo Change -- DEM_generation -- variable equal to 1 to do so.
elif (( $(echo "$DEM_generation == 1" | bc -l) )); then
    echo Interpolated DEM will be generated.
    echo Change -- DEM_generation -- variable equal to 0 to only do uncertainty analysis.
else
    echo Invalid variable input for -- DEM_generation -- variable. Exiting
    exit 1
fi

echo

if (( $(echo "$plots == 0" | bc -l) )); then
    echo No plots of measurement uncertainty will be made.
    echo Change -- plots -- variable equal to 1 to do so. 
    echo Using --st_error.py--
    source_uncert_py="st_error.py"
elif (( $(echo "$plots == 1" | bc -l) )); then
    echo Plots will be made for all raster cells with atleast $num_measurements measurements.
    echo Change -- plots -- variable equal to 0 to not generate plots. 
    echo Change num_measurements variable for different number than $num_measurements. 
    echo Using --st_error_plots.py--
    source_uncert_py="st_error_plots.py"
else
    echo Invalid variable input for -- Plots -- variable. Exiting
    exit 1
fi

echo

west_orig=$west
east_orig=$east
south_orig=$south
north_orig=$north

#initiate tile names with numbers, starting with 1
tile_name="1"
#get data buffer in geographical space, based on cell size
data_buff=$(echo "$cellsize_degrees * $buff_cells" | bc -l)
#echo data_buff is $data_buff

#Take in a half-cell on all sides so that grid-registered raster edge aligns exactly on desired extent
half_cell=$(echo "$cellsize_degrees / 2" | bc -l)
#echo half_cell is $half_cell

west_reduced=$(echo "$west + $half_cell" | bc -l)
north_reduced=$(echo "$north - $half_cell" | bc -l)
east_reduced=$(echo "$east - $half_cell" | bc -l)
south_reduced=$(echo "$south + $half_cell" | bc -l)

#keep copy of final mb-system bounds
west_buff_final=$(echo "$west_reduced - $data_buff" | bc -l)
north_buff_final=$(echo "$north_reduced + $data_buff" | bc -l)
east_buff_final=$(echo "$east_reduced + $data_buff" | bc -l)
south_buff_final=$(echo "$south_reduced - $data_buff" | bc -l)

#Determine number of rows and columns with the desired cell size, rounding up to nearest integer.
#i.e., 1_9 arc-second
x_diff=$(echo "$east - $west + $data_buff + $data_buff " | bc -l)
y_diff=$(echo "$north - $south + $data_buff + $data_buff " | bc -l)
x_dim=$(echo "$x_diff / $cellsize_degrees" | bc -l)
y_dim=$(echo "$y_diff / $cellsize_degrees" | bc -l)
x_dim_int=$(echo "($x_dim+0.5)/1" | bc)
y_dim_int=$(echo "($y_dim+0.5)/1" | bc)

#break point to exit code for testing

#TEMP FOR HORILLO
# JUST SET NUMBER OF X ROWS AND Y ROWS

#x_dim_int=13500
#y_dim_int=16200
#cellsize_degrees=0.000092592595872
#west_buff_final=-82.750000000
#east_buff_final=-81.50000000
#south_buff_final=25.75000000
#north_buff_final=27.25000000

#0.000046296297936
#exit 1

####################################################################
######################          END           ######################
######################   Derived Variables    ######################
######################                        ######################
####################################################################


####################################################################
######################                        ######################
######################   FOLDER MANAGEMENT    ######################
######################                        ######################
####################################################################

echo 
echo "Starting Folder Management"

#delete scratch directory if exists, then make a new one
rm -rf $basename_"_scratch"
mkdir $basename_"_scratch"


#create datalist that ignores Resampled bathy surface
sed -n '/^ *[^#]/p' $datalist |
while read -r line
do
    line_name="$line"
    resamp_flag=$(echo $line_name | awk '{print $7}')
    if  [ "$resamp_flag" == "RESAMP" ]  
    then
      echo $line_name >> $basename_"_scratch"/resamp.datalist
    elif [ "$resamp_flag" == "IGNORE" ]
    then
        echo $line_name >> $basename_"_scratch"/ignore_uncert.datalist  
    else
      echo $line_name >> $basename_"_scratch"/no_resamp.datalist
    fi
done



#Copy Necessary scripts to scratch
cp tile_analysis.sh $basename_"_scratch"/tile_analysis.sh
cp st_error.py $basename_"_scratch"/st_error.py
cp $datalist $basename_"_scratch"/$datalist
cp st_error_plots.py $basename_"_scratch"/st_error_plots.py
cp unique_id_rast.py $basename_"_scratch"/unique_id_rast.py
cp split_samp.sh $basename_"_scratch"/split_samp.sh
cp rast_percentile.py $basename_"_scratch"/rast_percentile.py
cp samp_den.py $basename_"_scratch"/samp_den.py
cp distance_pts.py $basename_"_scratch"/distance_pts.py
cp ring_poly.py $basename_"_scratch"/ring_poly.py
cp error_distance_plots.py $basename_"_scratch"/error_distance_plots.py
cp error_distance_plots_small_distances.py $basename_"_scratch"/error_distance_plots_small_distances.py


#navigate to scratch file so that all results are moved there
main_dir=$PWD
#echo $main_dir
scratch_dir=$basename_"_scratch"
scratch_dir_full=$PWD/$scratch_dir
cd "${scratch_dir_full}"

resamp_datalist="resamp.datalist"
no_resamp_datalist="no_resamp.datalist"




echo
echo
echo
echo
echo
echo
echo
echo "Tile Name is" $basename_
echo "Cellsize in degrees is" $cellsize_degrees
echo "West is" $west
echo "East is" $east
echo "South is" $south
echo "North is" $north
echo
####################################################################
######################                         #####################
######################  END FOLDER MANAGEMENT  #####################
######################                         #####################
####################################################################


####################################################################
######################                         #####################
######################         ANALYSIS        #####################
######################                         #####################
####################################################################

#Grid data with no interpolation to get cells constrained by data (for later use in data mask)
#
echo
echo --Creating Source DEM mask -- no resampled datasets
echo
mb_source_name=$basename_"_mb_source"
grid_mb_source=$mb_source_name".grd"
mb_range="-R$west_buff_final/$east_buff_final/$south_buff_final/$north_buff_final"
mb_cellsize="-E$cellsize_degrees/$cellsize_degrees/degrees!"

echo --Running mbgrid...
mbgrid -I$no_resamp_datalist -O$mb_source_name \
        $mb_range \
        -A2 -D$x_dim_int/$y_dim_int -G3 -N -M \
        -C0 -S0 -F1

gdal_translate $grid_mb_source -a_srs EPSG:4326 $basename_"_mb_source.tif"
gdal_translate $mb_source_name"_num.grd" -a_srs EPSG:4326 $mb_source_name"_num.tif"
gdal_translate $mb_source_name"_sd.grd" -a_srs EPSG:4326 $mb_source_name"_sd.tif"

echo -- Creating source data grid to 1s and 0s
gdal_calc.py -A $basename_"_mb_source.tif" --outfile=$basename_"_all_source_1_0.tif"  --calc="1*(A>-999999999999999999999)" --overwrite

echo -- Calculating distance for entire grid
gdal_proximity.py $basename_"_all_source_1_0.tif" $basename_"_dist.tif"

echo -- Getting 95th percentile of distance to nearest
#python ./rast_percentile.py $basename_"_dist.tif"
max_int_dist=`python ./rast_percentile.py $basename_"_dist.tif" | tail -n 1`

echo -- 95th percentile of distance to nearest measurement is $max_int_dist

#exit 1


echo -- Creating shapefile of extents with buffer
gdaltindex $basename_"_DEM_extent_buff.shp" $basename_"_mb_source.tif"


#determine tile dimensions from max int dist, applying tile_factor to make larger or smaller
max_int_dist_int=$(echo "($max_int_dist+0.5)/1" | bc)

if (( $(echo "$max_int_dist_int > $x_dim_int" | bc -l) )); then
    echo max_int_dist is greater than entire grid size, changing max_int_dist_int to grid dimensions
    max_int_dist_int=$x_dim_int
elif (( $(echo "$max_int_dist_int > $y_dim_int" | bc -l) )); then
    echo max_int_dist is greater than entire grid size, changing max_int_dist_int to grid dimensions
    max_int_dist_int=$y_dim_int
else
    echo max_int_dist is less than entire grid size, continuing... 
fi


tile_dim_x_int=$(echo "scale=0; ((($tile_factor*$max_int_dist_int)+0.5) / 1) " | bc)
tile_dim_y_int=$(echo "scale=0; ((($tile_factor*$max_int_dist_int)+0.5) / 1) " | bc)

# Uncomment below to get dimensions other than 95 percentile of distance to data. 

# if (( $(echo "$tile_dim_x_int < 200" | bc -l) )); then
#   echo tiles are less than 200 rows/cols
#   echo increasing tile dimensions to 200 to improve speed
#   tile_dim_x_int=200
#     tile_dim_y_int=200
#     if (( $(echo "$tile_dim_x_int >= $x_dim_int" | bc -l) )); then
#     echo tiles are larger than grid, making tile dimensions the same dimensions as grid w buffer
#     tile_dim_x_int=$(echo "scale=0; ($x_dim_int) " | bc)
#     tile_dim_y_int=$(echo "scale=0; ($y_dim_int) " | bc)
#   else
#       echo tile dimensions are smaller than entire grid size, continuing with dimensions of 200... 
#   fi
# elif (( $(echo "$tile_dim_x_int >= $x_dim_int" | bc -l) )); then
#     echo tiles are larger than grid, making tile dimensions the same dimensions as grid w buffer
#     tile_dim_x_int=$(echo "scale=0; ($x_dim_int) " | bc)
#     tile_dim_y_int=$(echo "scale=0; ($y_dim_int) " | bc)
# else
#     echo tile dimensions are smaller than entire grid size and are larger than 200, continuing... 
# fi


# # Just Make tiles 200 by 200 cells for now. #THINK ABOUT MORE AND CHANGE TO SOMETHING THAT MAKES SENSE W INTERPOLATION DISTANCES
# # Maybe something like 100 to 300 cells?
# tile_dim_x_int=200
# tile_dim_y_int=200




tile_dim_x_cell=$(echo "scale=6;  $tile_dim_x_int * $cellsize_degrees / 1 " | bc -l)
tile_dim_y_cell=$(echo "scale=6;  $tile_dim_y_int * $cellsize_degrees / 1 " | bc -l)

tile_dim_x_cell_round=$(printf "%.5f\n" $tile_dim_x_cell)
tile_dim_y_cell_round=$(printf "%.5f\n" $tile_dim_y_cell)

echo x_dim_int is $x_dim_int
echo tile_dim_x_int is $tile_dim_x_int


x_dim_ratio=$(echo "($x_dim_int / ($tile_dim_x_int - $buff_cells - $buff_cells))" | bc)
#echo $x_dim_ratio
#always round up
x_dim_ratio_int=$(echo "($x_dim_ratio+1)" | bc)
#echo $x_dim_ratio_int

y_dim_ratio=$(echo "($y_dim_int / ($tile_dim_y_int - $buff_cells - $buff_cells))" | bc)
#always round up
y_dim_ratio_int=$(echo "($y_dim_ratio+1)" | bc)
#echo $y_dim_ratio_int

total_num_tiles=$(echo "($x_dim_ratio_int * $y_dim_ratio_int)" | bc)

echo total_num_tiles is $total_num_tiles

echo

x_win=$(echo "$x_dim_int - $buff_cells - $buff_cells " | bc)
echo Final width for grid analysis is $x_win
y_win=$(echo "$y_dim_int - $buff_cells - $buff_cells " | bc)
echo Final height for grid analysis is $y_win

#exit 1


echo
echo
echo
echo

echo
echo
echo -- Starting Tile Analysis
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo

#starting point for tiling
xoff=0
yoff=0
xsize=$(echo "$tile_dim_x_int+$buff_cells" | bc)
ysize=$(echo "$tile_dim_y_int+$buff_cells" | bc)
grid_xsize=$x_dim_int
grid_ysize=$y_dim_int

echo xsize w buffer is $xsize
echo ysize w buffer is $ysize
echo grid_xsize is $grid_xsize
echo grid_ysize is $grid_ysize

#exit 1

while [ "$(bc <<< "$xoff < $grid_xsize ")" == "1"  ]; do
    yoff=0
    while [ "$(bc <<< "$yoff < $grid_ysize ")" == "1"  ]; do
    echo --Extracting tile $tile_name out of approximately $total_num_tiles
    tile_name_full="tile_"$basename_"_"$tile_name
    dem_name_source=$tile_name_full"_dem_source"
    gdal_translate -srcwin $xoff $yoff $xsize $ysize $basename_"_all_source_1_0.tif" $dem_name_source"_ext.tif"
    xsize_clip=$(echo "$xsize-$buff_cells-$buff_cells" | bc)
    ysize_clip=$(echo "$ysize-$buff_cells-$buff_cells" | bc)
    gdal_translate -srcwin $buff_cells $buff_cells $xsize_clip $ysize_clip $dem_name_source"_ext.tif" $dem_name_source"_1_0_final.tif"
    echo Clipping distance grid to tile extents w buffer
    gdal_translate -srcwin $xoff $yoff $xsize $ysize $basename_"_dist.tif" "tile_"$tile_name"_dist.tif"
  
    west_tmp=`gdalinfo $dem_name_source"_ext.tif" | grep -e "Upper Left"  | awk '{print $4}' | sed 's/.$//'`
    north_tmp=`gdalinfo $dem_name_source"_ext.tif" | grep -e "Upper Left"  | awk '{print $5}' | sed 's/.$//'`
    east_tmp=`gdalinfo $dem_name_source"_ext.tif" | grep -e "Lower Right"  | awk '{print $4}' | sed 's/.$//'`
    south_tmp=`gdalinfo $dem_name_source"_ext.tif" | grep -e "Lower Right"  | awk '{print $5}' | sed 's/.$//'`
    west_tile=$(echo "$west_tmp + $half_cell" | bc -l)
    north_tile=$(echo "$north_tmp - $half_cell" | bc -l)
    east_tile=$(echo "$east_tmp - $half_cell" | bc -l)
    south_tile=$(echo "$south_tmp + $half_cell" | bc -l)
    rm $dem_name_source"_ext.tif"

    mb_range_tile="-R$west_tile/$east_tile/$south_tile/$north_tile"
    echo $mb_range_tile
    dem_name_source=$tile_name_full"_dem_source"
    grid_dem_source=$dem_name_source".grd"
    mbgrid -I$no_resamp_datalist -O$dem_name_source $mb_range_tile -A2 -D$xsize/$ysize -G3 -N -C0 -S0 -F1
    gdal_translate $grid_dem_source -a_srs EPSG:4326 $dem_name_source".tif"
    rm $grid_dem_source

    ./tile_analysis.sh $basename_ $tile_name $west_tile $east_tile $south_tile $north_tile $no_resamp_datalist $vdatum_uncert $source_uncert_py $num_measurements $buff_cells
    yoff=$(echo "$yoff+$ysize-$buff_cells-$buff_cells" | bc)
    tile_name=$((tile_name+1))
    done
  xoff=$(echo "$xoff+$xsize-$buff_cells-$buff_cells" | bc)
done

echo -- FINISHED ALL TILES

#########################################################################
#############################                    ########################
############################# SOURCE UNCERTAINTY ########################
#############################                    ########################
#########################################################################
echo -- Starting Source Data Uncertainty Analysis
# Source Uncertainty

echo
echo
echo
echo
echo
echo

#delete results directory if exists
rm -rf $basename_"_results"
mkdir $basename_"_results"
mkdir $basename_"_results"/$basename_"_training"
mkdir $basename_"_results"/$basename_"_ss_plots"
mkdir $basename_"_results"/$basename_"_ss_data"

echo
echo --Merging Source St Error Rasters
echo
gdal_merge.py -o $basename_"_st_error_merge.tif" st_error_tile_*.tif -a_nodata -9999

cp st_error_tile_1.tif $basename_"_results"/st_error_tile_1.tif
cp $basename_"_dist.tif" $basename_"_results"/$basename_"_distance_buff.tif"


#remove st_error tiles
rm st_error_tile_*.tif
#clip source measurement uncertainty to extents with buffer to eliminate data not needed
gdalwarp -cutline $basename_"_DEM_extent_buff.shp" -crop_to_cutline $basename_"_st_error_merge.tif" $basename_"_st_error_merge_clip.tif" -overwrite

#Create uncertainty surface input for resampled datasets
echo RESAMPLED DATASETS ANALYSIS
resamp_file=$resamp_datalist
if [ -f $resamp_file ] ; then
    #create data masks to mask out resampled values where there are measurements
    gdal_calc.py -A $basename_"_all_source_1_0.tif" --outfile=$basename_"_all_source_0_1.tif" --calc="0*(A==1)" --calc="1*(A==0)" 
    gdal_calc.py -A $basename_"_all_source_0_1.tif" --outfile=$basename_"_all_source_0_1_nan.tif" --calc="A*(A>0)" --NoDataValue=0
    #go through each resampled dataset
    sed -n '/^ *[^#]/p' $resamp_datalist |
    while read -r line
    do
        line_name="$line"
        #echo "Name read from file - $line_name"
        resamp_flag=$(echo $line_name | awk '{print $7}')
        #echo resamp_flag is $resamp_flag
        echo $line_name > ind_data.datalist
        line=$(echo $line_name | awk '{print $1}')
        line2=$(basename $line)
        line_final=${line2::-9}
        #echo "dataset is $line_final"
        ind_name=$basename_"_ind_source_$line_final"
        ind_data="-Iind_data.datalist"
        grid=$ind_name".grd"
        data_file=$ind_name".mb-1"
        echo "gridding $ind_name"

        # Run mbgrid
        #echo --Running mbgrid...
        mbgrid $ind_data -O$ind_name \
            $mb_range \
            -A2 -D$x_dim_int/$y_dim_int -G3 -N \
            -C0 -S0 -F1 

        #Check to see if there is data
        if [ -s $data_file ]
        then 
            echo "--Data in region"
            gdal_translate $grid -a_srs EPSG:4326 $ind_name".tif"
            rm $grid

            echo --Converting elevations to source data uncertainty
            weight_value=$(echo $line_name | awk '{print $3}')
            echo "weight value for $ind_name is $weight_value"
            uncert_value=$(echo $line_name | awk '{print $4}')
            echo "uncertainty value for $ind_name is $uncert_value"
            depth_func_tmp=$(echo $line_name | awk '{print $5}')
            echo "uncertainty value for $ind_name is $depth_func_tmp"
            vdatum_flag=$(echo $line_name | awk '{print $6}')
            echo "vdatum flag for $ind_name is $vdatum_flag"
            resamp_flag=$(echo $line_name | awk '{print $7}')
            echo "resamp flag for $ind_name is $resamp_flag"

            column_val='A'
            depth_func=${depth_func_tmp//D/$column_val}
            echo --updated depth function is $depth_func

            #creating absolute value of data so that depth uncertainty calculation is correct
            gdal_calc.py -A $ind_name".tif" --outfile=$ind_name"_abs.tif"  --calc="abs(A)" --overwrite

            if [ "$uncert_value" != "DEPTH" ] && [ "$vdatum_flag" == "VDATUM" ]
            then
                echo "name is" $ind_name
                echo --Dataset used VDATUM and has a constant uncertainty value 

                #echo "creating uncertainty raster with constant value of $uncert_value"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="uncert_value_"$ind_name"_abs.tif"  --calc="$uncert_value*(A>-999999999999999999999)" --overwrite

                #echo "creating vdatum uncertainty raster with constant value of $vdatum_uncert"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="vdatum_uncert_"$ind_name"_abs.tif"  --calc="$vdatum_uncert*(A>-999999999999999999999)" --overwrite

                #echo --Calculating square root sum of squares of source and vdatum uncertainty
                gdal_calc.py -A "vdatum_uncert_"$ind_name"_abs.tif" -B "uncert_value_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_"$ind_name"_abs.tif"  --calc="((A*A)+(B*B))**(1/2.0)" --overwrite

                #echo Masking out where there is original source measurements
                gdal_calc.py -A $basename_"_all_source_0_1_nan.tif" -B "resamp_uncert_value_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_mask_"$ind_name"_abs.tif"  --calc="A*B"
                #convert to xyz
                gdal_translate -of XYZ "resamp_uncert_value_mask_"$ind_name"_abs.tif" "resamp_uncert_value_mask_xyz_"$ind_name".xyz"

                #remove nan values
                grep -v "nan" "resamp_uncert_value_mask_xyz_"$ind_name".xyz" > "resamp_uncert_value_mask_xyz_"$ind_name"_final.xyz"

                rm "resamp_uncert_value_mask_xyz_"$ind_name".xyz"
                rm $ind_name".tif"
                rm $ind_name"_abs.tif"
                rm "resamp_uncert_value_"$ind_name"_abs.tif"
                rm "resamp_uncert_value_mask_"$ind_name"_abs.tif"


            elif [ "$uncert_value" == "DEPTH" ] && [ "$vdatum_flag" == "VDATUM" ]
            then
                echo "name is" $ind_name
                echo --Dataset used VDATUM and uncertainty is a function of depth

                #echo "creating uncertainty raster with uncertainty as a function of depth according to $depth_func"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="uncert_value_fdepth_"$ind_name"_abs.tif"  --calc="abs($depth_func)" --overwrite
                
                #echo "creating vdatum uncertainty raster with constant value of $vdatum_uncert"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="vdatum_uncert_"$ind_name"_abs.tif"  --calc="$vdatum_uncert*(A>-999999999999999999999)" --overwrite

                echo --Calculating square root sum of squares of source and vdatum uncertainty
                gdal_calc.py -A "vdatum_uncert_"$ind_name"_abs.tif" -B "uncert_value_fdepth_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_"$ind_name"_abs.tif"  --calc="((A*A)+(B*B))**(1/2.0)" --overwrite

                #echo Masking out where there is original source measurements
                gdal_calc.py -A $basename_"_all_source_0_1_nan.tif" -B "resamp_uncert_value_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_mask_"$ind_name"_abs.tif"  --calc="A*B"
                #convert to xyz
                gdal_translate -of XYZ "resamp_uncert_value_mask_"$ind_name"_abs.tif" "resamp_uncert_value_mask_xyz_"$ind_name".xyz"

                #remove nan values
                grep -v "nan" "resamp_uncert_value_mask_xyz_"$ind_name".xyz" > "resamp_uncert_value_mask_xyz_"$ind_name"_final.xyz"
                
                rm "resamp_uncert_value_mask_xyz_"$ind_name".xyz"
                rm $ind_name".tif"
                rm $ind_name"_abs.tif"
                rm "resamp_uncert_value_"$ind_name"_abs.tif"
                rm "resamp_uncert_value_mask_"$ind_name"_abs.tif"

            elif [ "$uncert_value" == "DEPTH" ] && [ "$vdatum_flag" != "VDATUM" ]
            then
                echo "name is" $ind_name
                echo --Dataset DID NOT use VDATUM and uncertainty is a function of depth

                #echo "creating uncertainty raster with uncertainty as a function of depth according to $depth_func"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="resamp_uncert_value_"$ind_name"_abs.tif"  --calc="abs($depth_func)" --overwrite

                #echo Masking out where there is original source measurements
                gdal_calc.py -A $basename_"_all_source_0_1_nan.tif" -B "resamp_uncert_value_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_mask_"$ind_name"_abs.tif"  --calc="A*B"
                #convert to xyz
                gdal_translate -of XYZ "resamp_uncert_value_mask_"$ind_name"_abs.tif" "resamp_uncert_value_mask_xyz_"$ind_name".xyz"

                #remove nan values
                grep -v "nan" "resamp_uncert_value_mask_xyz_"$ind_name".xyz" > "resamp_uncert_value_mask_xyz_"$ind_name"_final.xyz"
                
                rm "resamp_uncert_value_mask_xyz_"$ind_name".xyz"
                rm $ind_name".tif"
                rm $ind_name"_abs.tif"
                rm "resamp_uncert_value_"$ind_name"_abs.tif"
                rm "resamp_uncert_value_mask_"$ind_name"_abs.tif"

            else
                echo "name is" $ind_name
                echo --Dataset DID NOT use VDATUM and has a constant uncertainty value

                #echo "creating uncertainty raster with constant value of $uncert_value"
                gdal_calc.py -A $ind_name"_abs.tif" --outfile="resamp_uncert_value_"$ind_name"_abs.tif"  --calc="$uncert_value*(A>-999999999999999999999)" --overwrite

                #echo Masking out where there is original source measurements
                gdal_calc.py -A $basename_"_all_source_0_1_nan.tif" -B "resamp_uncert_value_"$ind_name"_abs.tif" --outfile="resamp_uncert_value_mask_"$ind_name"_abs.tif"  --calc="A*B"
                #convert to xyz
                gdal_translate -of XYZ "resamp_uncert_value_mask_"$ind_name"_abs.tif" "resamp_uncert_value_mask_xyz_"$ind_name".xyz"

                #remove nan values
                grep -v "nan" "resamp_uncert_value_mask_xyz_"$ind_name".xyz" > "resamp_uncert_value_mask_xyz_"$ind_name"_final.xyz"
                
                rm "resamp_uncert_value_mask_xyz_"$ind_name".xyz"
                rm $ind_name".tif"
                rm $ind_name"_abs.tif"
                rm "resamp_uncert_value_"$ind_name"_abs.tif"
                rm "resamp_uncert_value_mask_"$ind_name"_abs.tif"

            fi
        else
            echo "file is empty"
            [ -e $grid ] && rm $grid

        fi
    done
else
    echo No resampled datasets, SKIPPING ANALYSIS
fi



# Check to see if there were actually resampled data.
num_resamp_ds=$(ls resamp_uncert_value_mask_xyz* | wc -l)

echo "number of resampled datasets is" $num_resamp_ds

if [ $num_resamp_ds -gt 0 ]
then
    echo
    echo Merging Resampled Uncertainty
    echo
    #combine all resampled xyzs and remove null values
    cat resamp_uncert_value_mask_xyz_* > $basename_"_resamp_st_error_merge_mb_source_tmp.xyz"
    #cat C3_bathytopo_test_scratch/resamp_uncert_value_mask_xyz_* > C3_bathytopo_test_scratch/TEST_resamp_st_error_merge_mb_source_tmp.xyz
    #remove no data values
    grep -v "1.17549435082228751e-38" $basename_"_resamp_st_error_merge_mb_source_tmp.xyz" > $basename_"_resamp_st_error_merge_mb_source_tmp2.xyz"
    #blockmedian in case there are multiple resampled datasets
    blockmedian $basename_"_resamp_st_error_merge_mb_source_tmp2.xyz" -R$grid_mb_source | awk -F" " '{print $1,$2,$3}' > $basename_"_resamp_st_error.xyz"
else
    echo No resampled datasets in ROI
fi




if (( $(echo "$DEM_generation == 0" | bc -l) )); then
    echo No interpolated DEM will be generated
    echo Change -- DEM_generation -- variable equal to 1 to do so.
elif (( $(echo "$DEM_generation == 1" | bc -l) )); then
    echo Interpolated DEM will be generated.
    echo Change -- DEM_generation -- variable equal to 0 to only do uncertainty analysis.
    #echo -- Creating interpolated DEM
    dem_name=$basename_"_DEM"
    grid_dem=$dem_name".grd"
    mb_range="-R$west_buff_final/$east_buff_final/$south_buff_final/$north_buff_final"
    mb_cellsize="-E$cellsize_degrees/$cellsize_degrees/degrees!"
    # Run mbgrid
    #echo --Running mbgrid...
    mbgrid -I$datalist -O$dem_name \
            $mb_range \
            -A2 -D$x_dim_int/$y_dim_int -G3 -N \
            -C810000000/3 -S0 -F1 -T0.25 -X0.01

    gdal_translate $grid_dem -a_srs EPSG:4326 $basename_"_DEM.tif"
    gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_DEM.tif" $basename_"_results"/$basename_"_DEM.tif"
    #echo Moving files to results folder
    mv $basename_"_DEM.mb-1" $basename_"_results"/$basename_"_DEM.mb-1"
    mv $basename_"_DEM.grd.cmd" $basename_"_results"/$basename_"_DEM.grd.cmd"
    rm $basename_"_results"/$basename_"_DEM.tif.aux.xml"
else
    echo Invalid variable input for --DEM_generation-- variable. Exiting.
    exit 1
fi

##################################################################################
#############################                        #############################
#############################     SOURCE + VDATUM    #############################
#############################         SURFACE        #############################
#############################                        #############################
##################################################################################

#echo -- Converting weighted average source tif to xyz
gdal_translate -of XYZ $basename_"_st_error_merge_clip.tif" $basename_"_source_st_error_merge_nan.xyz" 

#echo -- Removing NaN values
grep -v "\-9999" $basename_"_source_st_error_merge_nan.xyz" > $basename_"_source_st_error_merge.xyz" 

resamp_file2=$resamp_datalist
if [ -f $resamp_file2 ] ; then
    echo merging resampled uncertainty with source data uncertainty
    cat $basename_"_source_st_error_merge.xyz" $basename_"_resamp_st_error.xyz"  > $basename_"_st_error_all.xyz"
else
    echo no resampled dataset
    mv $basename_"_source_st_error_merge.xyz" $basename_"_st_error_all.xyz" 
fi


echo
echo -- Creating interpolated surface of source and vdatum uncertainty
echo
surface $basename_"_st_error_all.xyz" -G$basename_"_s_v.grd" -R$grid_mb_source -T0.25 -Ll0.0
gdal_translate $basename_"_s_v.grd" -a_srs EPSG:4326 $basename_"_s_v_uncert.tif" 



#########################################################################
#############################               #############################
############################# INTERPOLATION #############################
#############################               #############################
#########################################################################

echo -- Starting Interpolation Uncertainty Analysis
echo
echo
echo
echo
echo
echo


#Calculate 5th percentile of Sampling Density for all tiles.
#This will be sampling density used in split-sample process

#if there is only one tile, use its sampling density.


ss_samp_den=`awk -F, '{print $5}' tile_info.csv | sort -nr | awk 'BEGIN{c=0} length($0){a[c]=$0;c++}END{p5=(c/100*5); p5=p5%1?int(p5)+1:p5; print a[c-p5-1]}'`
echo The 5th percentile of the sampling density for all tiles is $ss_samp_den

if (( $(echo "$ss_samp_den < 0.001" | bc -l) )); then
    echo ss_samp_den is less than 0.001, changing to 0.001
    ss_samp_den=0.001
elif [ -z "$ss_samp_den" ]; then
    echo ss_samp_den was not computed. Most likely due to there being only 1 tile. Setting sampling density to that of only tile. 
    ss_samp_den=`awk -F, '{print $5}' tile_info.csv`
else
    :
fi

echo The Sampling Density for split-sample process will be $ss_samp_den


#Do in Loop with different zones (Bathy, BathyTopo, Topo)
for zone in Bathy BathyTopo Topo; do
  echo Starting Analyis on $zone
  echo Subsetting tile_info.csv
  awk -F, -v var="$zone" '$4==var' tile_info.csv > "tile_info_"$zone"_all.csv"

  if [ $zone == Bathy ]; then
        echo zone is Bathy
        train_percentile=$min_samp_pct_bathy
        echo train_percentile is $train_percentile
        min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var="$train_percentile" 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        if [ -z "$min_samp_var" ]
        then
              echo min_samp_var can not be determined from given percentile...
              echo tiles may all have the same density.
              echo running again with 50th percentile
              min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var=50 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        else
              :
        fi
        echo min_samp_var is $min_samp_var
  elif [ $zone == BathyTopo ]; then
        echo zone is BathyTopo
        train_percentile=$min_samp_pct_bathytopo
        echo train_percentile is $train_percentile
        min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var="$train_percentile" 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        if [ -z "$min_samp_var" ]
        then
              echo min_samp_var can not be determined from given percentile...
              echo tiles may all have the same density.
              echo running again with 50th percentile
              min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var=50 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        else
              :
        fi
        echo min_samp_var is $min_samp_var
  else
        echo zone is Topo
        train_percentile=$min_samp_pct_topo
        echo train_percentile is $train_percentile
        min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var="$train_percentile" 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        if [ -z "$min_samp_var" ]
        then
              echo min_samp_var can not be determined from given percentile...
              echo tiles may all have the same density.
              echo running again with 50th percentile
              min_samp_var=`awk -F, '{print $5}' "tile_info_"$zone"_all.csv" | sort -nr | awk -v var=50 'BEGIN{c=0} length($0){a[c]=$0;c++}END{pvar=(c/100*var); pvar=pvar%1?int(pvar)+1:pvar; print a[c-pvar-1]}'`
        else
              :
        fi
        echo min_samp_var is $min_samp_var
  fi

  #FIND PERSPECTIVE TRAINING SITES.
  #Any tiles that have atleast the median sampling density for that zone can be potential training
  echo Finding Perspective Training tiles for zone $zone
  awk -F, -v var="$min_samp_var" '$5>=var' "tile_info_"$zone"_all.csv" > "tile_info_"$zone"_samples.csv"

  #randomly sample and calculate cumalative distance between every pair of points
  #Find 5 training zones the farthest apart from the 10 possible training sites
  #Randomly choose 5 training sites 100 times. Export results to csv with all info, and max cumaltive distance
  #Choose 5 with largest cumulative distance. 
  echo Selecting Best training site locations from perspective tiles for zone $zone
  #only do if there are perspective tiles for that zone
  if [[ -s "tile_info_"$zone"_samples.csv" ]] ; then
    echo "$zone has data."
    echo Proceeding with analysis

    max_cum=0
    i=0

    #if there is only one training tile available
    #testVar=( tile_info_"$zone"_samples.csv" )
    #num_lines_train=${#testVar[@]}
    #
    #
    #NEW, removed this below and based it only on num_training. 
    #Num_lines_train states the number of possible training zones, not the number of desired training zones. 
    num_lines_train=$(wc -l < "tile_info_"$zone"_samples.csv")
    echo number of possible training zones is $num_lines_train
    #

    if [ "$num_training" -lt "2" ] || [ "$num_lines_train" -lt "2" ] ;
    then
        echo Only 1 training, so using this one.
        #chmod 777 tile_info_"$zone"_samples.csv"
        cat "tile_info_"$zone"_samples.csv" > "tile_info_"$zone"_rand.csv"
        cat "tile_info_"$zone"_rand.csv" > "tile_info_"$zone"_training_tmp.csv"
        #chmod 777 tile_info_"$zone"_rand.csv"
    else 
        # do 500 simulations to get random combinations of training and try to maximize distance between them.
        echo More than 1 training zone, finding optimal training zones to maximize distance
        # NEW 
        # start off with some training zones in case it messes up?
        shuf -n $num_training "tile_info_"$zone"_samples.csv" > "tile_info_"$zone"_rand.csv"
        cat "tile_info_"$zone"_rand.csv" > "tile_info_"$zone"_training_tmp.csv"
        echo num_training is $num_training
        # Now try to get training zones with more distance between them. 
        while [ $i -lt 1000 ]
        do
            shuf -n $num_training "tile_info_"$zone"_samples.csv" > "tile_info_"$zone"_rand.csv"
            chmod 777 "tile_info_"$zone"_rand.csv"
            #head tile_info_"$zone"_rand.csv"
            #store variable
            max_cum_new=`python ./distance_pts.py "tile_info_"$zone"_rand.csv" | tail -n 1`
            echo old max_cum distance is $max_cum
            echo new max_cum distance is $max_cum_new
            if (( $(echo "$max_cum_new > $max_cum" | bc -l) )); then
                echo Found larger cumalative distance of $max_cum_new
                cat "tile_info_"$zone"_rand.csv" > "tile_info_"$zone"_training_tmp.csv"
                max_cum=$max_cum_new
            else
                :
            fi
            i=$((i + 1))
        done
    fi


    #Sort so that training tiles go in order 1...N
    sort -t, -nk1 "tile_info_"$zone"_training_tmp.csv" > "tile_info_"$zone"_training.csv"
    rm "tile_info_"$zone"_training_tmp.csv"

    cat "tile_info_"$zone"_training.csv" | while read line
    do
       #tile_num=$( awk -F, '{print $1}' $line )
       tile_num=`echo $line | awk -F, '{print $1}'`
       west_tile=`echo $line | awk -F, '{print $8}'`
       east_tile=`echo $line | awk -F, '{print $9}'`
       south_tile=`echo $line | awk -F, '{print $10}'`
       north_tile=`echo $line | awk -F, '{print $11}'`
       x_dim_int_tile=`echo $line | awk -F, '{print $12}'`
       y_dim_int_tile=`echo $line | awk -F, '{print $13}'`
       x_win_tile=`echo $line | awk -F, '{print $14}'`
       y_win_tile=`echo $line | awk -F, '{print $15}'`
       #
       #Additional variables that I decided I wanted to keep for tile_training
       tile_lat=`echo $line | awk -F, '{print $2}'`
       tile_lon=`echo $line | awk -F, '{print $3}'`
       tile_zone=`echo $line | awk -F, '{print $4}'`
       tile_samp_den=`echo $line | awk -F, '{print $5}'`
       tile_zMin=`echo $line | awk -F, '{print $6}'`
       tile_zMax=`echo $line | awk -F, '{print $7}'`

       #create empty file for results
       #cat > int_errors_dist_total.csv
       echo Running Split-Sample Process
       #echo $tile_num $west_tile $east_tile $south_tile $north_tile $x_dim_int_tile $y_dim_int_tile $x_win_tile $y_win_tile
       ./split_samp.sh $basename_ $tile_zone $no_resamp_datalist $cellsize_degrees $num_loops $max_int_dist_int $ss_samp_den $buff_cells $tile_num $west_tile $east_tile $south_tile $north_tile $x_dim_int_tile $y_dim_int_tile $x_win_tile $y_win_tile $tile_lat $tile_lon $tile_zone $tile_samp_den $tile_zMin $tile_zMax
       # SPLIT SAMPLE CREATES training_tile_info.csv

       # combine errors/distance csv
       cat int_errors_dist_total.csv >> int_errors_dist_total_all.csv
       rm -rf int_errors_dist_total.csv
       #
       #rm int_errors_dist_total.csv
       #exit 1
    done
  else
    echo "$zone is empty."
    #create empty text file that will be merged later
    touch "tile_info_"$zone"_training.csv"
  fi ;

done
echo "Finished zone analysis"



#exit 1

#Derive Equation for all the errors/dist.csv
#create backup of data with 50% of the lines and final backup with 1% of the lines


error_dist_num_lines=$(wc -l < int_errors_dist_total_all.csv)

if (( $(echo "$error_dist_num_lines > 50000000" | bc -l) )); then
    echo Quantified too many interpolation errors, file could cause memory issues in python.
    echo
    echo Creating file with random percentages of interpolation errors to avoid memory issue.
    echo
    echo
    echo IMPORTANT MESSAGE: Reduce the num_loops variable and/or num_training variable to save time in future.
    echo
    echo
    shuf -n 50000000 int_errors_dist_total_all.csv > int_errors_dist_total_py.csv
else
    echo Shouldnt encounter memory issue. Reduce the num_loops variable and/or num_training variable if you get KILLED message.
    echo KILLED message means memory issue in python. 
    #rename so file is the same for both orig and if taking random 50%
    mv int_errors_dist_total_all.csv int_errors_dist_total_py.csv
fi


echo Max interpolation distance of interest is $max_int_dist_int

#if max_interpolation_distance is less than 5, don't bin data. Just derive equation by calculating st dev for each unique distance to nearest measurement.
if (( $(echo "$max_int_dist_int < 5" | bc -l) )); then
    echo Max interpolation distance of interest is less than 5
    echo Not binning data.
    echo
    #co_eff=`python ./error_distance_plots.py int_errors_dist_total_py.csv "int_uncert_BIN_TEST" $ss_samp_den $max_int_dist_int | tail -n 1`
    rm error_distance_plots.py
    echo Didnt use, so deleted error_distance_plots.py
    echo Running error_distance_plots_small_distances.py...
    co_eff=`python ./error_distance_plots_small_distances.py int_errors_dist_total_py.csv "int_uncert" $ss_samp_den $max_int_dist_int | tail -n 1`

else
    echo Max interpolation distance is 5 cells or greater
    echo Binning data.
    rm error_distance_plots_small_distances.py
    echo Didnt use, so deleted error_distance_plots_small_distances.py
    echo Running error_distance_plots.py
    co_eff=`python ./error_distance_plots.py int_errors_dist_total_py.csv "int_uncert" $ss_samp_den $max_int_dist_int | tail -n 1`
fi



echo $co_eff
co_eff0=$(cut -d "(" -f 2 <<< $co_eff | awk -F, '{print $1}') 
co_eff1=$(cut -d "(" -f 2 <<< $co_eff | awk -F, '{print $2}' | sed 's/^.//')
co_eff2=$(cut -d "(" -f 2 <<< $co_eff | awk -F, '{print $3}' | sed 's/.$//' | sed 's/^.//')

echo COEFFICIENTS:
echo $co_eff0
echo $co_eff1
echo $co_eff2

equation="($co_eff0+($co_eff1*(A**$co_eff2)))"

echo -- Calculating Interpolation Uncertainty
echo "Applying equation: " $equation

gdal_calc.py -A $basename_"_dist.tif" --outfile=$basename_"_i_uncert.tif" --calc="$equation" --overwrite


echo --Calculating Total Uncertainty
gdal_calc.py -A $basename_"_s_v_uncert.tif" -B $basename_"_i_uncert.tif" --outfile=$basename_"_s_v_i_uncert.tif"  --calc="((A*A)+(B*B))**(1/2.0)" --overwrite

echo --Calculating Total Uncertainty at 95% confidence
gdal_calc.py -A $basename_"_s_v_i_uncert.tif" --outfile=$basename_"_s_v_i_uncert_95conf.tif"  --calc="A*1.96" --overwrite

echo -- Starting Training Tile Info Management

echo -- Reduce tile_info.csv to important columns
#reduce tile_info.csv to same columns
awk -F, '{print $1,$2,$3,$4,$5,$6,$7}' OFS=, tile_info.csv > tile_info_reduced.csv


#Note:training_tile_info.csv is generated from split_sample.sh and thus will not be generated if there are no interpolation errors

if [ ! -f training_tile_info.csv ]
then
  echo FATAL ERROR.
  echo No Tiles Suitable for Training. 
  echo Program DID NOT run sucessfully. Try again in area of more data. 
  exit 1
fi

echo -- Determining training and non-training tiles
awk -F, -vOFS=, '(NR==FNR){a[$1]=$0; next}
                   {
                      if(a[$1]){print $0,"training"}
                      else{print $0,"non_training"}
                   }' training_tile_info.csv  tile_info_reduced.csv > training_id.csv


#Join based on ID and put coefficients for derived equations in tile_info.csv
awk -F, -v OFS=, 'NR==FNR{a[$1]=$8","$9","$10;next}{print $0,a[$1]}' training_tile_info.csv training_id.csv > tile_info_join.csv


#Get training shps to save
cat tile_info_Bathy_training.csv tile_info_BathyTopo_training.csv tile_info_Topo_training.csv > tile_info_training_all.csv

cat tile_info_training_all.csv | while read line
do
   training_tile_id=`echo $line | awk -F, '{print $1}'`
   dem_name_source_train="tile_"$basename_"_"$training_tile_id"_dem_source"
   zone_train=`echo $line | awk -F, '{print $4}'`
   mv $dem_name_source_train"_interp_inner.shp" $basename_"_results"/$basename_"_training"/$zone_train"_train_tile_"$training_tile_id".shp"
   mv $dem_name_source_train"_interp_inner.dbf" $basename_"_results"/$basename_"_training"/$zone_train"_train_tile_"$training_tile_id".dbf"
   mv $dem_name_source_train"_interp_inner.prj" $basename_"_results"/$basename_"_training"/$zone_train"_train_tile_"$training_tile_id".prj"
   mv $dem_name_source_train"_interp_inner.shx" $basename_"_results"/$basename_"_training"/$zone_train"_train_tile_"$training_tile_id".shx"
done

#remove tile source data used to make polygons for training
rm tile_*dist.tif
rm weights_tile*.tif
rm tile_$basename_*.tif

echo
echo "Clipping to Final Extents"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_dist.tif" $basename_"_results"/$basename_"_distance.tif"

gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $mb_source_name"_num.tif" $basename_"_results"/$mb_source_name"_num.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $mb_source_name"_sd.tif" $basename_"_results"/$mb_source_name"_sd.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_all_source_1_0.tif" $basename_"_results"/$basename_"_all_source_1_0.tif"

#gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $mb_source_name"_st_error.tif" -a_srs EPSG:4326 $basename_"_results"/$mb_source_name"_st_error.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_st_error_merge_clip.tif" -a_srs EPSG:4326 $basename_"_results"/$basename_"_st_error.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_i_uncert.tif" $basename_"_results"/$basename_"_i_uncert.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_s_v_i_uncert.tif" $basename_"_results"/$basename_"_s_v_i_uncert.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_s_v_i_uncert_95conf.tif" $basename_"_results"/$basename_"_s_v_i_uncert_95conf.tif"
gdal_translate -srcwin $buff_cells $buff_cells $x_win $y_win $basename_"_s_v_uncert.tif" $basename_"_results"/$basename_"_s_v_uncert.tif"




#TEST TO SEE WHATS GOING ON W SOURCE UNCERTAINTY SURFACE
#mv $basename_"_resamp_st_error.xyz" $basename_"_results"/$basename_"_resamp_st_error.xyz"
#mv $basename_"_source_st_error_merge.xyz" $basename_"_results"/$basename_"_source_st_error_merge.xyz"
mv $basename_"_st_error_all.xyz" $basename_"_results"/$basename_"_st_error_all.xyz"





echo
#Create shapefile of final DEM extent
gdaltindex $basename_"_results"/$basename_"_DEM_extent.shp" $basename_"_results"/$basename_"_distance.tif"

mv training_tile_info.csv $basename_"_results"/training_tile_info.csv
mv int_uncert_best_fit.png $basename_"_results"/$basename_"_ss_plots"/$basename_"_int_uncert_best_fit.png"
mv int_uncert_scatter.png $basename_"_results"/$basename_"_ss_plots"/$basename_"_int_uncert_scatter.png"
mv int_errors_dist_total_py.csv $basename_"_results"/$basename_"_ss_data"/$basename_"_int_errors_dist.csv"

#rm *"int_uncert_tile_"*
#remove extra xml files
rm -rf $basename_"_results"/$basename_"_s_v_uncert.tif.aux.xml"
rm -rf $basename_"_results"/$basename_"_mb_source_st_error.tif.aux.xml"
rm -rf $basename_"_results"/$basename_"_mb_source_sd.tif.aux.xml"
rm -rf $basename_"_results"/$basename_"_mb_source_num.tif.aux.xml"

#
# rm $basename_"_mb_source.grd"
# rm $basename_"_st_error_merge.tif"
# rm $basename_"_s_v.grd"
# rm $basename_"_s_v.tif"
# rm $basename_"_mb_source.tif"
# rm $basename_"_DEM.tif"
# rm $basename_"_s_v_i_uncert.tif"
# rm $basename_"_s_v_i_uncert_95conf.tif"
# rm $basename_"_i_uncert.tif"
# rm $basename_"_all_source_1_0.tif"
# rm $basename_"_all_source_0_1.tif"
# rm $basename_"_st_error_merge_clip.tif"
# rm $basename_"_all_source_0_1_nan.tif"
# rm $basename_"_dist.tif"

# rm id_elev_w_u.csv
# rm id_elev_w_u_tmp.csv
# rm tile_*int_errors.xyz
# [ -e $basename_".grd" ] && rm $basename_".grd"


rm *.tif
rm *.xyz
rm *.grd
rm *.csv
rm *.tif.aux.xml

cd ".."
rm -rf $basename_"_results"
cd $basename_"_scratch"

# move results up a directory
mv $basename_"_results"/ ..




#rm unique_id.tif
echo 
echo Finished Analysis!
echo

#
# REMOVE SCRATCH FOLDER IF NOT NEEDED
#rm -rf mydir 
#

# Go up a directory to start next tile
cd ".."



done