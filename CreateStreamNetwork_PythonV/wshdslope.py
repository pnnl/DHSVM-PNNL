# SUMMARY:      wshdslop.py
# USAGE:        Part of python version createstreamnetwork. Compute 
#		watershed slope and aspect using four neighbours.  
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts createstreamnetwork.aml as part of DHSVM.
# 		Reason for 4 neighbour: make the computation of surface 
# 		and subsurface flow pathways consistent with the digital 
# 		elevation model networks (DEMON) described by Costa-Cabral
#		and Burges (1994)
# Last Change: 2017-08-10 by zduan


import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import os
import math
import numpy as np


def wshdslopefun(dem, slope, aspect, flowacc):
    # to make life easier, create an arry from dem to do calculation
    rad2deg= 180.0 / math.pi

    # Convert DEM to numpy array with no data value set to -9999
    demarray = arcpy.RasterToNumPyArray(dem,nodata_to_value=-9999)
    demraster=arcpy.sa.Raster(dem)
  #  arcpy.env.outputCoordinateSystem = demraster

    # Get raster properties
    elevcellx = arcpy.GetRasterProperties_management(dem, "CELLSIZEX")
    deltax=float(elevcellx.getOutput(0))
    elevcelly = arcpy.GetRasterProperties_management(dem, "CELLSIZEY")
    deltay=float(elevcelly.getOutput(0))
    
   # arcpy.env.cellSize = demraster
    lowerLeft = arcpy.Point(demraster.extent.XMin,demraster.extent.YMin)
    

    # Skip wshdoutlet.aml, do calculation directly
    # use max flow accumulation point instead of lowest elevation point as outlet
    # Convert flow accumulation to numpy array with no data value set to -9999
    accarr = arcpy.RasterToNumPyArray(flowacc,nodata_to_value=-9999)
    minloc=np.unravel_index(accarr.argmax(), accarr.shape)
    accarr[minloc[0]][minloc[1]]=0
    outlet=arcpy.NumPyArrayToRaster(accarr,lowerLeft,deltax,deltay)
    # raster layer outlet was not used?
    
    
    size=demarray.shape  # size=[Number of Rows, Number of Columns]

    # Extend boundary of dem input by add empty row/column before/after first/last row/column
    demarr=np.zeros((size[0]+2,size[1]+2))
    slopearr=np.zeros((size[0],size[1]))
    flag=np.zeros((size[0],size[1]))
    tmp=np.zeros((size[0],size[1]))
    # Change array value to nodata value in DEM array
    demarr[:][:]=-9999
    for ii in range(1,size[0]+1):
        demarr[ii][1:(size[1]+1)]=demarray[ii-1][:]

    # Calculate slope for extended DEM array
    for i in range(1,size[0]+1):  # loop through i-row
        for j in range(1,size[1]+1): # loop through j-col

            #Force boundary cells have aspect towards inner basin
            #(-1,-1)--(0,-1)--(1,-1)
            #   |       |        |
            #   |       |        | 
            #(-1, 0)--(0, 0)--(1, 0)
            #   |       |        |
            #   |       |        |
            #(-1, 1)--(0, 1)--(1, 1)
        

            # calculate dx
            if (demarr[i][j+1]==-9999) and (demarr[i][j-1]==-9999):
                dzdx=0.0
            elif demarr[i][j-1]==-9999: #aspect(0-180) or(0,pi),dzdz>=0
                dzdx=max((demarr[i][j]-demarr[i][j+1])/(deltax),0.0)
            elif demarr[i][j+1]==-9999: #aspect(180-360) or(-pi,0),dzdz<=0
                dzdx=min((demarr[i][j-1]-demarr[i][j])/(deltax),0.0)
            else:
                dzdx=(demarr[i][j+1]-demarr[i][j-1])/(2*deltax)

            # calculate dy
            if (demarr[i+1][j]==-9999) and (demarr[i-1][j]==-9999):
                dzdy=0.0
            elif demarr[i-1][j]==-9999: #aspec (pi/2,pi)||((-pi,-pi/2)), dzdy<=0
                dzdy=min((demarr[i+1][j]-demarr[i][j])/(deltay),0.0)
            elif demarr[i+1][j]==-9999: #aspec (0, pi/2)||((-pi/2,0)), dzdy>=0
                dzdy=max((demarr[i][j]-demarr[i-1][j])/(deltay),0.0)
            else:
                dzdy=(demarr[i+1][j]-demarr[i-1][j])/(2*deltay)
                
            slopearr[i-1][j-1]=math.sqrt(math.pow(dzdx,2)+math.pow(dzdy,2))
            if ((dzdx==0) and (dzdy==0)):
                #since atam2(0,0) is undefined
                flag[i-1][j-1]=1        
            else:
                flag[i-1][j-1]=0
                tmp[i-1][j-1]=int(math.atan2(dzdx,dzdy) * rad2deg +360)

    flag
    tmp
    FlagRaster = arcpy.NumPyArrayToRaster(flag,lowerLeft,deltax,deltay)
    TmpRaster = arcpy.NumPyArrayToRaster(tmp,lowerLeft,deltax,deltay)

    tmp2=Aspect(dem) 
    junk=Con(FlagRaster==1,tmp2,TmpRaster)
    tmpaspect=Con(junk>360,junk%360,junk)

    tmpslope = arcpy.NumPyArrayToRaster(slopearr,lowerLeft,deltax,deltay)
    sloperas=SetNull(dem==-9999,tmpslope)
    aspectras=SetNull(dem==-9999,tmpaspect)
    
    sloperas.save(slope)
    aspectras.save(aspect) 
    
