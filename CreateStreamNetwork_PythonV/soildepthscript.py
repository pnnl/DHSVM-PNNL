# SUMMARY:      soildepth.py
# Called by: 	createstreamnetwork.py
# USAGE:        Part of python version createstreamnetwork. Computes
#		the aspect and length of stream arc within each DEM 
#		cell.
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts soildepth.aml as part of DHSVM

# This aml creates a soildepth file for DHSVM based on local slope (determined from DEM),
# upstream source area, and elevation. There are a number of variables that need to be set:
#
# mindepth - the minimum depth of the soil (this is a floor)
# maxdepth - the maximum depth of the soil (this will never be exceeded)
# wtslope - the relative weighting for the slope
# wtsource - the relative weighting for the source area
# wtelev - the relative weighting for the elevation
# maxslope - anything greater than this will create the slope function = 1
# maxsource - anything greater than this will create the source function = 1
# maxelev - anything greater than this will create the elevation function = 1
# powslope - raise the slope fraction by this power
# powsource - raise the source area fraction by this power
# powelev - raise the elevation fraction by this power
#
#
# Last Change: 2017-08-10 by zduan

import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import os

def soildepthfun(flowacc, elev, mindepth, maxdepth, soildepth):

# the below variables can/should be modified by the user
    wtslope = 0.7
    print('The relative weighting for slope is '+ str(wtslope))
    wtsource = 0.0
    print('The relative weighting for source area is '+str(wtsource))
    wtelev = 0.3
    print('The relative weighting for elevation is '+ str(wtelev))
    maxslope = 30.
    maxsource = 100000.
    maxelev = 1500
    powslope = .25
    powsource = 1.
    powelev = .75
# end of variables can/should be modified by the user

    print('All value read in')
    totalwt = wtslope + wtsource + wtelev
    print(str(totalwt))
     
    if totalwt <> 1.0:
        print('the weights must add up to 1.')
    else:
        print('start creating file')
     
        #env.mask = elev
        #env.cellSize = elev
        
        slopegrid=Slope(elev,"DEGREE",1)
        slopegrid.save("slopegrid")
        
        tempsrc = Con(Raster(flowacc)>maxsource,maxsource,flowacc)
        tempelev = Con(Raster(elev)>maxelev,maxelev,elev)
        tempslope = Con(Raster("slopegrid")>maxslope,maxslope,slopegrid)

        print('calculate soil  map')  
        soild = mindepth + (maxdepth - mindepth)*\
                (wtslope * (1.0 - Power((tempslope/ maxslope), powslope))+ \
                 (wtsource * Power( (tempsrc/ maxsource), powsource))+\
                 wtelev *(1.0 - Power( (tempelev / maxelev), powelev)))
   
        soild.save(soildepth)
        
        print('delete temp files')  
        arcpy.Delete_management(tempsrc)
        arcpy.Delete_management(tempelev)       
        arcpy.Delete_management(tempslope)
        arcpy.Delete_management(slopegrid)

        #return soildepth



              
