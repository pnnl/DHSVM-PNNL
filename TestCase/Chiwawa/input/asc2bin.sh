#!/bin/bash

#--------------------------------------------------------------------------------
# This is a bash shell script to convert ascii file to binary file
# please run the script in ./input directory
# use "dos2unix" commant first to remove the PC '/r/n' newline character
# written by Hongxiang Yan at PNNL on Mar 6, 2017
#--------------------------------------------------------------------------------

find . -type f -exec dos2unix {} \;  #change all files from PC to Linux format

# 1. remove the header file (6 lines)
tail -n +7 dem.asc       > new_dem.asc
tail -n +7 mask.asc      > new_mask.asc
tail -n +7 soilclass.asc > new_soilclass.asc
tail -n +7 vegclass.asc  > new_vegclass.asc
# tail -n +7 flowdir.asc   > new_flowdir.asc
tail -n +7 soildepth.asc > new_soildepth.asc

#2. save the header file
head -6 mask.asc > headerfile.asc

#3. do the conversion 
../programs/myconvert asc float new_dem.asc       dem.bin        5 100    #row, column
../programs/myconvert asc char  new_mask.asc      mask.bin       5 100
../programs/myconvert asc char  new_soilclass.asc soilclass.bin  5 100
../programs/myconvert asc int  new_vegclass.asc  vegclass.bin   5 100
../programs/myconvert asc float new_soildepth.asc soildepth.bin  5 100

#4. remove the temporal files
rm new_dem.asc new_mask.asc new_vegclass.asc new_soilclass.asc new_soildepth.asc 
