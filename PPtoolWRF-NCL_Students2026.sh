#!/bin/bash
#SBATCH -p serial
#SBATCH -J PPTOOL
##SBATCH -C node_type=64core256G
#SBATCH --mem=128G

ulimit -s unlimited
source /etc/profile
source /home/hpc/46/mazzeo/.aliases.bash

module add anaconda3/4.3.1
source activate /storage/hpc/46/mazzeo/envs/nco_students

###########################################################################
#  Author : AMZ, 2024 
#
###########################################################################
# Select Name of experiment out"
ExpName_In='HW2023'
# Select wrf domain
dom=02 #d02 d03 d04
# Select Name of experiment out"
ExpName_Out=${ExpName_In}_ct_v3.2
DomOut=5K  # 9K 3K 1K
# Select path where are input data"
DataDir_in=/mmfs1/scratch/hpc/46/mazzeo/data/WRF-Chem/${ExpName_In}/emis_test/ct_v3.2/

# Select path where write output data"
DataDir_out=/mmfs1/scratch/hpc/46/mazzeo/data/WRF-Chem/${ExpName_In}/wrfout_pptool_${DomOut}_weath
mkdir -p ${DataDir_out}

# Select file to NCL processes
NCLFile=./1File_PPtoolWRFChem_updated.ncl
###########################################################################
# Select time data
FrecData=PP
# First day of input data
FirstDay=20
# First month of input data
FirstMonth=05
# Last day of input data
LastDay=21
# Last month of input data
LastMonth=05 #${FirstMonth}
# Year of input data
Year=2023
# Select wich variable is processed (true or false)
###########################################################################

OutFile2D=${FrecData}_WRF_${ExpName_Out}-2D_${DomOut}_${FirstDay}-${FirstMonth}-${LastDay}-${LastMonth}.nc

echo '--------------------------------------------- '
echo '               :>> NCL Options <<:'
echo '   o     File selected for NCL process:>> '${NCLFile}
echo ' '
echo '               :>> Input Options <<:'
echo '   o     Expriment selected:>> '${ExpName_In}
echo '   o     Domain selected:>> d'${dom}
echo '   o     File selected:>> '$InFile
echo '   o     Path selected:>> '${DataDir_in}
echo ' '
echo '               :>> Output Options <<:'
echo '   o     Expriment selected:>> '${ExpName_Out}
echo '   o     Domain selected:>> '${DomOut}
echo '   o     Path selected:>> '${DataDir_out}
echo '--------------------------------------------- '
###########################################################################

                 ######################################
                 #   Processing WRFchem Variables     #
                 ######################################
rm ./wrfout_*
echo ' ' 
echo 'Linking original wrfout files in outdir/tmp directory'		
for mon in "05" ; do #"07" "08" "09" ; do
   if [ "${mon}" = "05" ]; then
     days=$(seq -w 20 21)
   elif [ "${mon}" = "06" ]; then
     days=$(seq -w 01 30)
   elif [ "${mon}" = "07" ] || [ "${mon}" = "08" ] ; then
     days=$(seq -w 01 31)
   elif [ "${mon}" = "09" ]; then
     days=$(seq -w 01 01)
   fi
   for day in ${days}; do
      echo '--------------------------------------------- '
      echo ':>>   Processing WRF Variables with NCL   <<:'
      echo '      Day: '${Year}'-'${mon}'-'${day}
      echo '--------------------------------------------- '
      toto_ncl=toto-ncl_${mon}${day}.nc
      ln -s ${DataDir_in}/wrfout_d${dom}_${Year}-${mon}-${day}_00:00:00 wrfout_d${dom}_${Year}-${mon}-${day}_0000 
      InFile=wrfout_d${dom}_${Year}-${mon}-${day}_0000
      ncl 'file_in="'${InFile}'"' 'file_out="'${toto_ncl}'"' 'Year="'${Year}'"' 'Month="'${mon}'"' 'Day="'${day}'"' ${NCLFile}
      mv toto-ncl_${mon}${day}.nc ${DataDir_out}/.
      echo ' '
      echo ':Processing WRF Variables with NCL:>> DONE!'
      echo '--------------------------------------------- '

      echo '--------------------------------------------- '
      echo ':>>   Processing WRF Variables with NCO   <<:'
      echo '--------------------------------------------- '
      echo '   :>>  Transferring 2D WRF Variables ...'
      echo '  '
      echo ' - Processing Time variable'
      ncks -O -v Time ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing grid scale precipitation variable'
      ncks -A -v precip_g ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Latitude variable'
      ncks -A -v lat ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Longitude variable'
      ncks -A -v lon ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Terrain Height variable'
      ncks -A -v Z_sfc ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Temperature 2m variable'
      ncks -A -v T_2m ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Surface Pressure variable'
      ncks -A -v p_sfc ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Relative Humidity 2m variable'
      ncks -A -v rh_2m ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing U-Wind 10m variable'
      ncks -A -v u_10m_tr ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing V-Wind 10m variable'
      ncks -A -v v_10m_tr ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Wind Speed 10m variable'
      ncks -A -v ws_10m ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Wind Direction variable'
      ncks -A -v wd_10m ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing PBLH variable'
      ncks -A -v pblh ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Short Wave Radiaton downward variable'
      ncks -A -v SW_d ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Large Wave Radiation downward variable'
      ncks -A -v LW_d ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Leaf Area Index variable'
      ncks -A -v LAI ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Soil Moisture variable'
      ncks -A -v SoilMoist ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' - Processing Soil Temperature variable'
      ncks -A -v SoilTemp ${DataDir_out}/toto-ncl_${mon}${day}.nc ${DataDir_out}/toto_${mon}${day}.nc
      echo ' - DONE!'
      echo ' --------------------------------------------- '   
   done
   echo ' --------------------------------------------- '
   echo '   :>> Merging daily 2D Chem Variables ...'
   echo '   o   Output File selected:>> '$OutFile2D
   echo ' '
   echo 'Merging days of month :'${mon}	
   ncrcat -O ${DataDir_out}/toto_${mon}??.nc ${DataDir_out}/toto_${mon}.nc
done    
echo 'Merging all months :'
ncrcat -O ${DataDir_out}/toto_??.nc ${DataDir_out}/$OutFile2D
echo ' --------------------------------------------- '

#Change names of variables and dimenssions
echo '   o   Rename variables'
ncrename -v .ChemPM10,PM10 -v .ChemPM25,PM25 -v .ChemOzon,O3 -v .ChemNO2,NO2 -v .ChemCO,CO -v .ChemNO,NO -v .ChemSO2,SO2 ${DataDir_out}/$OutFile2D
ncrename -v .Z_sfc,HG -v .T_2m,TEMP -v .p_sfc,P -v .rh_2m,RH -v .u_10m_tr,WSPEEDU -v .v_10m_tr,WSPEEDV -v .ws_10m,WSPEED ${DataDir_out}/$OutFile2D
ncrename -v .wd_10m,VD -v .SW_d,RADSW -v .LW_d,RADLW -d .south_north,lat -d .west_east,lon ${DataDir_out}/$OutFile2D
ncrename -v .pblh,PBLH -v .precip_g,PRAIN ${DataDir_out}/$OutFile2D
ncrename -v .SoilMoist,SMOIS -v .SoilTemp,TSLB ${DataDir_out}/$OutFile2D
rm -rf ${DataDir_out}/toto*
rm -rf wrfout_d*
echo ':Processing WRF Variables with NCO:>> DONE!'
echo '--------------------------------------------- '
