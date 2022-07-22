# SUMMARY:      rowcolmap.py
# USAGE:        Part of python version createstreamnetwork. 
#               Original Description: It is a simple script used to produce a
#               polygon representation of a grid. Each polygon in outcover
#               represents a single cell, and is assigned row and column 
#               numbers.
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts roadaspect.aml as part of DHSVM.
#
# Last Change: 2017-08-10 by zduan


import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import math
import numpy as np

def rowcolmapfun(elev,aspect,slope,MaskExt):
    # Get raster properties
    #aspectcellx = arcpy.GetRasterProperties_management(aspect, "CELLSIZEX")
    #deltax=float(aspectcellx.getOutput(0))
    #aspectcelly = arcpy.GetRasterProperties_management(aspect, "CELLSIZEY")
    #deltay=float(aspectcelly.getOutput(0))
    env.mask = MaskExt

    env.cellSize=elev
    arcpy.env.extent=elev
    #.outputCoordinateSystem = elev
    
    demraster=arcpy.sa.Raster(elev)
    RightRaster=CreateConstantRaster (1, "INTEGER", env.cellSize, demraster.extent) 
    #UpRaster=CreateConstantRaster (64, "INTEGER", env.cellSize, demraster.extent)
    UpRaster=CreateConstantRaster (4, "INTEGER", env.cellSize, demraster.extent)
    

    col=FlowAccumulation(RightRaster)
    col.save("colraster")
    
    rownum=FlowAccumulation(UpRaster)
    rownum.save("rowraster")

    # Create Grid-Points and Grid-Polygon
	
    if arcpy.Exists("cb"):
        arcpy.Delete_management("cb")
        
    print('joining row map and col map')	
    rowcolcomb=Combine(["colraster","rowraster"])
    rowcolcomb.save("cb")

    if arcpy.Exists("rowcolpoly.shp"):
        arcpy.Delete_management("rowcolpoly.shp")

    arcpy.RasterToPolygon_conversion(rowcolcomb, "rowcolpoly.shp", "NO_SIMPLIFY", "VALUE")

    
    print('Joinning table, this may take a while...')

    print('Copying fields...')
    joinfields = ['VALUE', 'COLRASTER', 'ROWRASTER']
    joindict = {}
    with arcpy.da.SearchCursor("cb", joinfields) as rows:
        for row in rows:
            joinval = row[0]
            val1 = row[1]
            val2 = row[2]
            joindict[joinval]=[val1, val2]
    del row, rows

    print('Pasting fields...')

    arcpy.AddField_management("rowcolpoly.shp", "ROW", "LONG")
    arcpy.AddField_management("rowcolpoly.shp", "COL", "LONG")

    ## Step 3
    targfields = ['GRIDCODE', 'COL', 'ROW']
    with arcpy.da.UpdateCursor("rowcolpoly.shp", targfields) as recs:
        for rec in recs:
            keyval = rec[0]

            ## Step 4
            if joindict.has_key(keyval):
                rec[1] = joindict[keyval][0]
                rec[2] = joindict[keyval][1]
            else:
                rec[1] = 0
                rec[2] = 0
            recs.updateRow(rec)
    del rec, recs
    
    print('Join table finished')
    
    
    arcpy.Delete_management(col)
    arcpy.Delete_management(rownum)
    arcpy.Delete_management(rowcolcomb)
    #arcpy.Delete_management(tmp)
    arcpy.Delete_management(RightRaster)
    arcpy.Delete_management(UpRaster)
    #arcpy.Delete_management("tmpraster")
    

    
    

