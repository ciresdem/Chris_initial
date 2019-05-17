#!/bin/sh -e
#To Do:

# NEGATIVE ERROR MEANS INTERPOLATED SURFACE IS DEEPER THAN TRUE SURFACE

#Version Changes:
# Created variable, samples_threshold, so that if a raster cell has over threshold, uncertainty is set to 0.001. 
# This speeds up computation in cases where the are thousands of measurements in a coarse cell.
# When creating DEM from stack, these areas of dense data will be higher resolution anyways, so doesn't matter if they are set to arbitrarily low uncertainty.

# Fixed issue with samples threshold producing empty file. 
# AFTER CHANGING SCRATCH

####################################################################
####################################################################
########################                  ##########################
########################  GET PARAMETERS  ##########################
########################       FROM       ##########################
######################## MASTER_UNCERT.SH ##########################
########################              	  ##########################
########################                  ##########################
####################################################################
####################################################################

basename_=$1
tile_name=$2
west_tile_buff=$3
east_tile_buff=$4
south_tile_buff=$5
north_tile_buff=$6
datalist=$7
vdatum_uncert=$8
source_uncert_py=$9
num_measurements=${10}
buff_cells=${11}
# New Variable (add to master_uncert.sh later)
samples_threshold=100

####################################################################
####################################################################
########################                  ##########################
########################  END PARAMETERS  ##########################
########################       FROM       ##########################
######################## MASTER_UNCERT.SH ##########################
########################                  ##########################
########################                  ##########################
####################################################################
####################################################################

####################################################################
####################################################################
########################                  ##########################
########################    Additional    ##########################
########################    Parameters    ##########################
########################                  ##########################
####################################################################
####################################################################

#make data buffer to have data points on outside (DO THIS IN MASTER)
#buff_cells=5
#Field Sep
IFS=''
#current_dir
current_dir="$PWD"

####################################################################
######################                        ######################
######################   FOLDER MANAGEMENT    ######################
######################                        ######################
####################################################################

#Delete csvs if exists
if [ -f id_elev_w_u.csv ] ; then
    rm id_elev_w_u.csv
fi

if [ -f id_elev_w_u_final.csv ] ; then
    rm id_elev_w_u_final.csv
fi

if [ -f id_elev_w_u_tmp.csv ] ; then
    rm id_elev_w_u_tmp.csv
fi

if [ -f id_elev_w_u_tmp2.csv ] ; then
    rm id_elev_w_u_tmp2.csv
fi

if [ -f files_used_list.txt ] ; then
    rm files_used_list.txt
fi

##################################################################
#############################        #############################
#############################   DEM  #############################
#############################        #############################
##################################################################

# echo "Creating Source DEM w No Interpolation to get data files in grid bounds"

# dem_name_source=$name"_dem_source"
# grid_dem_source=$dem_name_source".grd"

# echo --Running mbgrid...
# mbgrid -I$datalist -O$dem_name_source $mb_range_tile -A2 -D$x_dim_int_tile/$y_dim_int_tile -G3 -N -C0 -S0 -F1
# gdal_translate $grid_dem_source -a_srs EPSG:4326 $dem_name_source".tif"


dem_name_source="tile_"$basename_"_"$tile_name"_dem_source"


# echo SOURCE DATA UNCERTAINTY ANALYIS
# #Create Unique ID raster based on DEM w buffer
python ./unique_id_rast.py $dem_name_source".tif"

#echo "clipping grid to remove edge effect"
echo Getting number of rows and columns
x_dim_tile_buff="$(gdalinfo $dem_name_source".tif" | grep -oP "is\s+\K\w+" | head -1)"
echo Number of columns--width is $x_dim_tile_buff
y_dim_tile_buff="$(gdalinfo $dem_name_source".tif" | grep -oP ",\s+\K\w+" | head -1)"
echo Number of rows--height is $y_dim_tile_buff


x_dim_tile=$(echo "$x_dim_tile_buff - $buff_cells - $buff_cells" | bc)
echo final extent width is $x_dim_tile
y_dim_tile=$(echo "$y_dim_tile_buff - $buff_cells - $buff_cells" | bc)
echo final extent height is $y_dim_tile


# gdal_translate -srcwin $buff_cells $buff_cells $x_dim_tile $y_dim_tile $dem_name_source".tif" $dem_name_source"_final.tif"

# #echo "Creating data mask"
gdal_calc.py -A $dem_name_source".tif"  --outfile=$dem_name_source"_mask.tif"  --calc="1*(A>-999999999999999999999)" --overwrite




#################################################################################
#############################                               #############################
#############################   Source Uncertainty  #############################
#############################                               #############################
#################################################################################
dem_name_source="tile_"$basename_"_"$tile_name"_dem_source"

echo -- Getting datafiles that were used to create DEM from $dem_name_source".tif"
files_used_list=$dem_name_source".mb-1"

echo $files_used_list
#stat -s $files_used_list
file_size=$( stat -c %s $files_used_list)
echo size of file is $file_size

#check to see if there is any data:
if [ $file_size = 0 ]
then
    echo No data, skipping...
else
    echo TILE HAS DATA, STARTING ANALYSIS
    #Remove R: (First 2 characters) of each line
    cat $files_used_list | sed 's/^..//' > $tile_name"_dem_tmp2.mb-1"
    files_used_list=$tile_name"_dem_tmp2.mb-1"


    #remove lines with # sign and go through each line of original datalist
    sed -n '/^ *[^#]/p' $datalist |
    while read line_tmp || [ -n "$line_tmp" ];
      do
      echo line_tmp is $line_tmp
      #Get name of each dataset
      DIR=$(dirname "${line_tmp}")
      echo dir is $DIR
      #echo directory is $DIR
      line=$(basename $line_tmp)
      echo line is $line

      #Get list of files in mb-1 that start with DIR (get all the files in the dataset)
      grep "$DIR" $files_used_list | awk '{print $1}' > files_used_list.txt
      cat files_used_list.txt
      chmod +x files_used_list.txt

    #Read info from each dataset on datalist
    cat files_used_list.txt | while read input_name
      do          
      #add the weight value to the new datalist from the original datalist
      weight_value=$(echo $line | awk '{print $3}')
      echo "weight value for $input_name is $weight_value"
      #use the original uncertainty value (1 SD) to use in the python script
      st_dev=$(echo $line | awk '{print $4}')
      echo "standard deviation for $input_name is $st_dev"
      depth_func_tmp=$(echo $line | awk '{print $5}')
      #Replace D in function to column 4, which will be absolute value of depth
      column_val='$4'
      depth_func=${depth_func_tmp//D/$column_val}
      echo "uncertainty as a function of depth for $input_name is $depth_func"
      vdatum_flag=$(echo $line | awk '{print $6}')
      echo "vdatum flag for $input_name is $vdatum_flag"
      resamp_flag=$(echo $line | awk '{print $7}')
      echo "resamp flag for $input_name is $resamp_flag"

      delim_tmp=$(head -1 $input_name)
      delim="Not Recognized"
      echo "$delim_tmp" | grep -q ' ' && delim=" " 
      echo "$delim_tmp" | grep -q ',' && delim=","
      echo "$delim_tmp" | grep -q '|' && delim="|"

      if  [ "$st_dev" != "DEPTH" ] && [ "$vdatum_flag" == "VDATUM" ]
      then
        echo "name is" $input_name
        echo --Dataset used VDATUM and has a constant uncertainty value

        st_dev=$(echo $line | awk '{print $4}')
        echo "standard deviation for $input_name is $st_dev"

        #echo vdatum uncert is $vdatum_uncert

        #combine vdatum and constant uncertainty value using root sum of squares
        combined_st_dev=$(echo "sqrt (($st_dev*$st_dev)+($vdatum_uncert*$vdatum_uncert))" | bc -l)
        #echo --combined standard deviation is, $combined_st_dev

        #Extract unique_ids
        if [ "$delim" == " " ]
          then
          echo Delimiter is Space
          gdal_query.py -delimiter " " -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS=" " OFS=, w=$weight_value u=$combined_st_dev > id_elev_w_u_tmp.csv
        elif [ "$delim" == "," ]
          then
          echo Delimiter is Comma
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS="," OFS=, w=$weight_value u=$combined_st_dev > id_elev_w_u_tmp.csv
        elif [ "$delim" == "|" ]
          then
          echo Delimiter is Pipe
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS="|" OFS=, w=$weight_value u=$combined_st_dev > id_elev_w_u_tmp.csv
        else
          echo IMPORTANT: ERROR with File $i
          echo Delimiter is not recognized as space, comma, or pipe.
          echo Fix delimiter for file $i to space, commma, or pipe and re-run program!
          echo Exiting...
          echo
          exit 1
        fi

        #remove lines with -99999
        sed '/-99999.0/d' id_elev_w_u_tmp.csv | awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,$4}' >> id_elev_w_u.csv

      elif [ "$st_dev" == "DEPTH" ] && [ "$vdatum_flag" == "VDATUM" ]
      then
        echo "name is" $input_name
        echo --Dataset used VDATUM and uncertainty is a function of depth according to $depth_func

        #Extract unique_ids
        
        if [ "$delim" == " " ]
          then
          echo Delimiter is Space
          gdal_query.py -delimiter " " -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS=" " OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        elif [ "$delim" == "," ]
          then
          echo Delimiter is Comma
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS="," OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        elif [ "$delim" == "|" ]
          then
          echo Delimiter is Pipe
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS="|" OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        else
          echo IMPORTANT: ERROR with File $i
          echo Delimiter is not recognized as space, comma, or pipe.
          echo Fix delimiter for file $i to space, commma, or pipe and re-run program!
          echo Exiting...
          echo
          exit 1
        fi

        #remove lines with -99999 and FIXED TO TAKE ABSOLUTE VALUE OF ELEVATION TO USE IN CALCULATING UNCERTAINTY DEPTH
        sed '/-99999.0/d' id_elev_w_u_tmp.csv | awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,$2<0?$2*-1:$2}' > id_elev_w_u_tmp_1.csv

        #Calcute uncertainty as a function of depth
        awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,'"$depth_func"'}' id_elev_w_u_tmp_1.csv > id_elev_w_u_tmp2.csv

        #Incorporate vdatum uncertainty
        awk -v vdatum_uncert="$vdatum_uncert" -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,sqrt(($4*$4)+(vdatum_uncert*vdatum_uncert))}' id_elev_w_u_tmp2.csv  >> id_elev_w_u.csv


      elif [ "$st_dev" == "DEPTH" ] && [ "$vdatum_flag" != "VDATUM" ]
      then
        echo "name is" $input_name
        echo --Dataset DID NOT USE VDATUM and uncertainty is a function of depth according to $depth_func
        #Extract unique_ids
        
        if [ "$delim" == " " ]
          then
          echo Delimiter is Space
          gdal_query.py -delimiter " " -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS=" " OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        elif [ "$delim" == "," ]
          then
          echo Delimiter is Comma
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS="," OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        elif [ "$delim" == "|" ]
          then
          echo Delimiter is Pipe
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w}' FS="|" OFS=, w=$weight_value > id_elev_w_u_tmp.csv
        else
          echo IMPORTANT: ERROR with File $i
          echo Delimiter is not recognized as space, comma, or pipe.
          echo Fix delimiter for file $i to space, commma, or pipe and re-run program!
          echo Exiting...
          echo
          exit 1
        fi

        #remove lines with -99999 and FIXED TO TAKE ABSOLUTE VALUE OF ELEVATION TO USE IN CALCULATING UNCERTAINTY DEPTH
        sed '/-99999.0/d' id_elev_w_u_tmp.csv | awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,$2<0?$2*-1:$2}' > id_elev_w_u_tmp_1.csv

        #Calcute uncertainty as a function of depth
        awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,'"$depth_func"'}' id_elev_w_u_tmp_1.csv > id_elev_w_u.csv

      else
        echo "name is" $input_name
        echo --Dataset DID NOT USE VDATUM and has a constant uncertainty value
        
        #Extract unique_ids
        if [ "$delim" == " " ]
          then
          echo Delimiter is Space
          gdal_query.py -delimiter " " -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS=" " OFS=, w=$weight_value u=$st_dev > id_elev_w_u_tmp.csv
        elif [ "$delim" == "," ]
          then
          echo Delimiter is Comma
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS="," OFS=, w=$weight_value u=$st_dev > id_elev_w_u_tmp.csv
        elif [ "$delim" == "|" ]
          then
          echo Delimiter is Pipe
          gdal_query.py -delimiter "$delim" -s_format "0,1,2" -d_format "gz" -d_nodata "-99999" unique_id.tif $input_name | awk '{print $1,$2,w,u}' FS="|" OFS=, w=$weight_value u=$st_dev > id_elev_w_u_tmp.csv
        else
          echo IMPORTANT: ERROR with File $i
          echo Delimiter is not recognized as space, comma, or pipe.
          echo Fix delimiter for file $i to space, commma, or pipe and re-run program!
          echo Exiting...
          echo
          exit 1
        fi        

        #remove lines with -99999
        sed '/-99999.0/d' id_elev_w_u_tmp.csv | awk -F "," '{printf "%.0f,%.3f,%.10f,%.3f\n", $1,$2,$3,$4}' >> id_elev_w_u.csv 
      fi
      done
  done
  
  #check again if data, run script
  file_size2=$( stat -c %s id_elev_w_u.csv )
  echo size of file is $file_size2
  #num_id_elev_w_u=$(ls id_elev_w_u.csv | wc -l)
  #echo number of lines in num_id_elev_w_u is $num_id

  #check to see if there is any data:
  if [ $file_size2 = 0 ]
  then
      echo No data, skipping...
  else
  	echo Data Found, running source_uncert_py...

  #id_elev_w_u.csv
  # Determine number of measurements in file
  #note: num_measurements is a user defined variable in master_uncert.sh that set the number of measurements needed to create a figure of cell-level uncertainty in st_error_plot.py
  # num_resamp_ds=$(ls resamp_uncert_value_mask_xyz* | wc -l)

  # change delimiter
  #tr ',' ' ' < id_elev_w_u.csv > id_elev_w_u2.csv

  echo Determining is there are cells with samples over threshold to reduce computational cost.

  echo Counting samples for each unique id
  # Count the number of samples for each unique id
  awk -F ',' '{print $1}' id_elev_w_u.csv | sort | uniq -c > id_elev_w_u_num_samples.csv

  # Make file of unique ids with samples over a set limit (i.e., 100)
  echo Samples Threshold is $samples_threshold
  awk -v samples_threshold_var="$samples_threshold" '$1 > samples_threshold_var' id_elev_w_u_num_samples.csv | awk '{print $2}'> id_elev_w_u_over_thresh.csv
  #awk '$1 > 100' id_elev_w_u_num_samples.csv | awk '{print $2}'> id_elev_w_u_over_thresh.csv

  # See if there are samples over threshold
  # If so, then remove and replace with dummy 0.001 uncertainty value
  # If not, then proceed to python script.
  over_thresh_size=$( stat -c %s id_elev_w_u_over_thresh.csv )

  #check to see if there is any data over threshold:
  if [ $over_thresh_size = 0 ]
  then
      echo No data over threshold, proceeding to source_uncert_py...
      cat id_elev_w_u.csv > id_elev_w_u_final.csv
  else
    echo Data over threshold. 
    echo Removing lines with ids over threshold
    awk -F ',' 'NR==FNR{a[$1]=$1;next} !($1 in a){print $0}' id_elev_w_u_over_thresh.csv id_elev_w_u.csv > id_elev_w_u_under_thresh.csv

    echo Setting uncertainty to zero for ids over threshold
    #create file with over_thresh ids with zeros for elevation, weight (0.001), uncertainty (0.001). Uncertainty should be the only one that matters
    awk 'BEGIN{FS=",";OFS=","}{print $1,0,0.001,0.001}' id_elev_w_u_over_thresh.csv > id_elev_w_u_over_thresh_zero_uncert.csv
  
    echo Merging back ids over threshold with ids under threshold
    # Add back in line to fine with uncertainty of zero and the average weight. For now, set average weight to zero. It appears I don't do anything with weights raster in master_uncert.sh anyways.
    # With only one line, source_uncert_py will just set uncertainty for that ID to zero.
    cat id_elev_w_u_under_thresh.csv id_elev_w_u_over_thresh_zero_uncert.csv > id_elev_w_u_final.csv
  fi 
  
  echo Starting Standard Error Analysis with source_uncert_py
  python ./$source_uncert_py id_elev_w_u_final.csv $tile_name $dem_name_source".tif" $num_measurements

  #if png exists:
  #mv *$basename.png plots

  #exit 1

  #apply data masks to remove values with no data (and thus have the value of the unique id of that cell still)
  echo -- Removing unique_id values where there is no data
  #echo --Applying Data Masks
  #old method using data mask
  #gdal_calc.py -A $dem_name_source"_mask.tif" -B "st_error.tif" --outfile="tmp_st_error_"$tile_name".tif"  --calc="A*B" --overwrite
  
  #Changed method: This will ensure no unique id values, whereas previous method, the data mask would occasionally not work on a few cells
  #remove any negative values (areas w/o data that will be unique id value)
  gdal_calc.py -A "st_error.tif" --outfile="tmp_st_error_"$tile_name".tif" --calc="A*(A>0)" --NoDataValue=0
  gdal_calc.py -A "weights.tif" --outfile="tmp_weights_"$tile_name".tif" --calc="A*(A>0)" --NoDataValue=0


  #Clip to remove edge effect
  echo --Clipping grid to remove edge effect
  gdal_translate -srcwin $buff_cells $buff_cells $x_dim_tile $y_dim_tile "tmp_st_error_"$tile_name".tif" "st_error_tile_"$tile_name".tif"
  gdal_translate -srcwin $buff_cells $buff_cells $x_dim_tile $y_dim_tile "tmp_weights_"$tile_name".tif" "weights_tile_"$tile_name".tif"
  

  #remove files not needed
  rm "tmp_st_error_"$tile_name".tif"
  rm unique_id.tif
  rm st_error.tif
  rm "tmp_weights_"$tile_name".tif"
  rm weights.tif
  
  fi

fi


#exit 1
echo 
echo
echo FINISHED SOURCE UNCERTAINTY FOR TILE $tile_name
echo
echo



echo 

echo STARTING INTERPOLATION UNCERTAINTY ANALYSIS FOR TILE $tile_name


echo --Determine center lat and lon 
tile_lon=`gdalinfo $dem_name_source".tif" | grep -e "Center"  | awk '{print $3}' | sed 's/.$//'`
tile_lat=`gdalinfo $dem_name_source".tif" | grep -e "Center"  | awk '{print $4}' | sed 's/.$//'`

echo lat is $tile_lat
echo lon is $tile_lon

echo --Determining zone using min and max elevation
#all negative = bathy
#negative and positive = bathytopo
#all positive = topo

zMin=`gdalinfo -mm $dem_name_source".tif" | sed -ne 's/.*Computed Min\/Max=//p'| tr -d ' ' | cut -d "," -f 1`
zMax=`gdalinfo -mm $dem_name_source".tif" | sed -ne 's/.*Computed Min\/Max=//p'| tr -d ' ' | cut -d "," -f 2` 


if (( $(echo "$zMin < 0.0000" | bc -l) )) && (( $(echo "$zMax < 0.0000" | bc -l) )); then
    zone="Bathy"
elif (( $(echo "$zMin > 0.0000" | bc -l) )) && (( $(echo "$zMax > 0.0000" | bc -l) )); then
    zone="Topo"
else
    zone="BathyTopo"
fi


if [ -z "$zMin" ]
then
      echo variable is empty, tile has no data
      zMin=99999
      zMax=99999
      zone="NoData"
else
      :
fi

echo zMin is $zMin
echo zMax is $zMax
echo zone is $zone




#echo --Coverting source raster to 1s and 0s.
#gdal_calc.py -A $basename_"_scratch"/$dem_name_source".tif" --outfile=$basename_"_scratch"/$dem_name_source"_1_0_buff.tif"  --calc="1*(A>-999999999999999999999)" --overwrite

#echo --Removing Buffer before calculating sampling density
#gdal_translate -srcwin $buff_cells $buff_cells $x_dim_tile $y_dim_tile $basename_"_scratch"/$dem_name_source"_1_0_buff.tif" $basename_"_scratch"/$dem_name_source"_1_0.tif"

#make shp for every tile to clip distance and calculate tile interpolation uncer
#echo --Creating Polygon for tile w/o buffer
#gdaltindex "tile_"$tile_name"_extent.shp" $dem_name_source"_1_0_final.tif"

echo --Calculating Sampling Density...
sampling_d_sh=`python samp_den.py $dem_name_source"_mask.tif" | tail -n 1`

echo -- Sampling density from shell script is $sampling_d_sh


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
echo
echo
echo
echo
echo


if (( $(echo "$zMin < 0.0000" | bc -l) )) && (( $(echo "$zMax < 0.0000" | bc -l) )); then
    zone="Bathy"
elif (( $(echo "$zMin == 99999" | bc -l) )); then
  zone="NoData"    
elif (( $(echo "$zMin > 0.0000" | bc -l) )) && (( $(echo "$zMax > 0.0000" | bc -l) )); then
    zone="Topo"
else
    zone="BathyTopo"
fi

echo tile number is $tile_name
echo zMin is $zMin
echo zMax is $zMax
echo zone is $zone


echo Adding info to csv
echo $tile_name","$tile_lat","$tile_lon","$zone","$sampling_d_sh","$zMin","$zMax","$west_tile_buff","$east_tile_buff","$south_tile_buff","$north_tile_buff","$x_dim_tile_buff","$y_dim_tile_buff","$x_dim_tile","$y_dim_tile >> "tile_info.csv"


rm $dem_name_source"_1_0_final.tif"
rm $dem_name_source"_mask.tif"
