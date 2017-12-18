#Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting
import sys
import os
import math
import numpy as np
import csv 


#-------------------------------------------------------------------#
#--------------------------- WorkSpace  ----------------------------#    
#-------------------------------------------------------------------#
env.workspace = "C:\\Users\\username\\Documents\\foldername"   
path = "C:/Users/username/Documents/foldername/"

#-------------------------------------------------------------------#
###########           Setup Input          
#-------------------------------------------------------------------#
elev = "dem"                          # name of DEM GRID file
wshed = "mask"                        # name of MASK file
soildepth = "soild"                   # name of soil depth file
streamfile = "streamfile"             # name of stream arc file
key = 'MASK'                          # Enter 'MASK' or 'MOUTH'
source = 4860000                      # Min source area to initiate stream
mindepth = 0.76                       # Minimum Soil Depth                        
maxdepth = 2.01                       # Maximum Soil Depth

#-------------------------------------------------------------------#
#------------------------   End of Edits ---------------------------#
#---------------  Spatial Analysis License Required   --------------# 
#-------------------------------------------------------------------#


#Check if inputs are valid
if not arcpy.Exists(elev):
    sys.exit("elevation input not valid, exit program")

if not arcpy.Exists(wshed):
    sys.exit("mask/pour point input not valid, exit program")

# Set the cell size environment using a raster dataset.
env.outputCoordinateSystem = elev
env.extent = elev

env.cellSize = elev
arcpy.env.overwriteOutput = True

arcpy.CheckOutExtension("Spatial")

# Create a mask of full extention for rowcolmap
if arcpy.Exists("maskras"):
    arcpy.Delete_management("maskras")
MaskExt=CreateConstantRaster (1, "INTEGER", env.cellSize, env.extent)
MaskExt.save("maskras")

# Execute FlowDirection
if arcpy.Exists("flow_dir"):
    arcpy.Delete_management("flow_dir")

print('Flow direction')
flowdir = FlowDirection(elev, "", "")

# Execute Watershed
if key=='MASK':
    print('Mask file provided')
    wshd = wshed
    env.mask = wshd
    flowdir = FlowDirection(elev, "", "")
elif key=='MOUTH':
    print('Mask file not provided, computing watershed')
    flowdir = FlowDirection(elev, "", "")
    wshd = Watershed(flowdir , wshed, "")
    env.mask = wshd
else:
    print('Input key not valid')

flowdir.save("flow_dir")

# Execute flow accumulation

if arcpy.Exists("flow_acc"):
    arcpy.Delete_management("flow_acc")
env.mask = wshd    
print('Flow accumulation')    
flowacc = FlowAccumulation(flowdir, "", "INTEGER")
flowacc.save("flow_acc")

detlaxResult = arcpy.GetRasterProperties_management(elev,"CELLSIZEX")
deltax = detlaxResult.getOutput(0)
sourcepix = source / (float(deltax) * float(deltax))

temp = Con(flowacc > sourcepix, 1)  

rivg= temp * wshd

if arcpy.Exists("/output.gdb"):
    print('geodatabase already exist')
else:
    print('geodatabase does NOT exist, create geodatabse')
    streamgdb = "output.gdb"
    arcpy.CreateFileGDB_management(path, streamgdb)

streamlink=StreamLink (rivg, flowdir)
if arcpy.Exists(path+"/output.gdb/"+streamfile):
    arcpy.Delete_management(path+"/output.gdb/"+streamfile)
    print('deleted previous stream file and create a new one')

streamnet=path+"/output.gdb/"+streamfile
print('creating stream shapefile')
StreamToFeature(streamlink, flowdir, streamnet, "NO_SIMPLIFY")                                                        

# Find contributing area to each cell along stream network 

if arcpy.Exists(rivg):
    arcpy.Delete_management(rivg)


arcpy.PolylineToRaster_conversion (streamnet, "arcid", rivg, "MAXIMUM_LENGTH","NONE", env.cellSize)

local=Watershed(flowdir, rivg, "VALUE")

if arcpy.Exists("local"):
    arcpy.Delete_management("local")

local.save("local")
arcpy.AddField_management (streamnet, "local", "LONG")
arcpy.JoinField_management(streamnet,"arcid",local,"Value","Count")

fields=['COUNT','local']

with arcpy.da.UpdateCursor(streamnet, fields) as cursor:
    for row in cursor:
        if row[0] is None:
            row[1] = 0
        else:
            row[1]=row[0]
        cursor.updateRow(row)
del row

arcpy.DeleteField_management(streamnet, "COUNT")

print('local file created')
arcpy.Delete_management(rivg)

#-------------------------------------------------------------------#
#### create node point coverage and find elevations  
#-------------------------------------------------------------------#

if arcpy.Exists(path+"output.gdb/nodestart"):
    arcpy.Delete_management(path+"/output.gdb/nodestart")
    print('start node file already exists, delete and create new')

nodestart=path+"output.gdb/nodestart"
arcpy.FeatureVerticesToPoints_management(streamnet, nodestart, "START")

if arcpy.Exists(path+"output.gdb/nodeend"):
    arcpy.Delete_management(path+"/output.gdb/nodeend")
    print('end node file already exists, delete and create new')

nodeend=path+"output.gdb/nodeend"
arcpy.FeatureVerticesToPoints_management(streamnet, nodeend, "END")

# Find Contributing area at end of each arc

env.mask="maskras"
tmpacc=Con(IsNull("flow_acc")==1,int(flowacc.maximum),"flow_acc")

ExtractMultiValuesToPoints(nodestart, [[tmpacc, "MAXGRID"]])

arcpy.AddField_management (streamnet, "downarc", "LONG")

ExtractMultiValuesToPoints(nodestart, [[elev, "SELEV"]], "NONE")
elevras=Raster(elev)
tmpelev=Con(IsNull(elevras)==1,int(elevras.minimum),elev)
ExtractMultiValuesToPoints(nodeend, [[tmpelev, "EELEV"]], "NONE")

env.mask = wshd  

# easier to add elevation infomation with 3D extention 
#arcpy.CheckOutExtension ("3D")
#arcpy.AddSurfaceInformation_3d (nodes, elev, "Z", "CONFLATE_NEAREST" , "", "", "", "")

arcpy.JoinField_management(streamnet,"arcid",nodestart,"arcid","SELEV")
arcpy.JoinField_management (streamnet, "arcid", nodeend, "arcid", "EELEV")
arcpy.JoinField_management (streamnet, "arcid", nodestart, "arcid", "MAXGRID")

arcpy.AddField_management (streamnet, "uparc", "LONG")
arcpy.AddField_management (streamnet, "dz", "FLOAT", 12, 3)
arcpy.AddField_management (streamnet, "slope", "FLOAT", 12, 5)
arcpy.AddField_management (streamnet, "meanmsq",  "FLOAT")
arcpy.AddField_management (streamnet, "segorder", "LONG")
arcpy.AddField_management (streamnet, "chanclass", "SHORT", 8, "")
arcpy.AddField_management (streamnet, "hyddepth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnet, "hydwidth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnet, "effwidth", "FLOAT", 8, 2)
arcpy.AddField_management (streamnet, "effdepth", "FLOAT", 8, 2)

print('Calculating Slope of channel segment')
arcpy.CalculateField_management (streamnet, "dz", "abs(!SELEV! - !EELEV!)", "PYTHON_9.3", "")

expression = "clacSlope(float(!dz!),float(!Shape_Length!))"
codeblock = """def clacSlope(dz,length):
    if (dz/length)>0.00001:
        return dz/length
    else:
        return 0.00001"""

arcpy.CalculateField_management (streamnet, "slope", expression, "PYTHON_9.3", codeblock)

arcpy.CalculateField_management (streamnet, "segorder", "0", "PYTHON_9.3")
arcpy.CalculateField_management (streamnet, "uparc", "0", "PYTHON_9.3")
arcpy.CalculateField_management (streamnet, "downarc", "-1", "PYTHON_9.3")
arcpy.CalculateField_management (streamnet, "meanmsq", "0.0", "PYTHON_9.3")

#-------------------------------------------------------------------#
# Calculate Segment Order       
#-------------------------------------------------------------------#
arr=arcpy.da.TableToNumPyArray(streamnet,('from_node','to_node','segorder','local','MAXGRID','meanmsq','uparc','downarc','arcid'))

print('Looking for downstream arc')
# Calculate downstream arc 
for jj, ii in enumerate(arr['to_node']):
    arr2=[]
    for i, j in enumerate(arr['from_node']):
        if j == ii:
            arr['downarc'][jj]=arr['arcid'][i]


print('Looking for upstream arc')
# Calculate upstream arc based on max conributing area
for jj, ii in enumerate(arr['from_node']):
    arr2=[]
    for i, j in enumerate(arr['to_node']):
        if j == ii:
            arr2=np.append(arr2,i)
    #print arr2
    if not len(arr2):
        arr['uparc'][jj]=-1
        arr['segorder'][jj]=1
        arr['meanmsq'][jj]=arr['local'][jj]/ 2 * float(env.cellSize) * float(env.cellSize)
    else:
        arr3=arr2.astype(int)
        loc=np.argmax(arr['MAXGRID'][arr3]+arr['local'][arr3])
        arr['uparc'][jj]=arr['arcid'][arr3[loc]]
        arr['meanmsq'][jj]=(arr['MAXGRID'][jj]+arr['local'][jj]/ 2) * float(env.cellSize) * float(env.cellSize)
    
# Calculate segorder
order=1
a=99
print('Calculating segment order')
while a > 0:
    a=0
    for jj, ii in enumerate(arr['segorder']):
        if ii==order:
            a+=1
            for i, j in enumerate(arr['arcid']):
                if j == arr['downarc'][jj]:
                    arr['segorder'][i]=max(order+1,arr['segorder'][i])            
    order+=1
   
               
# Append the array to an existing table
arcpy.da.ExtendTable(streamnet,  "arcid", arr, "arcid", append_only=False)

#-------------------------------------------------------------------#
###########          Chanel Hydraulic Classes          
#-------------------------------------------------------------------#

from channelclass import channelclassfun
print('Assign channl class')
channelclassfun(streamnet)

#-------------------------------------------------------------------#
###########           SOIL DEPTH          
#-------------------------------------------------------------------#

if arcpy.Exists(soildepth):
    print('Soil depth file already exists, delete and create new')
    arcpy.Delete_management(soildepth)
else:
    print('soil depth file not provided, creating map')

print('making the soil depth with mindepth of '+ str(mindepth)+' and maxdepth of '+ str(maxdepth))
from soildepthscript import soildepthfun

soildepthfun("flow_acc", elev, mindepth, maxdepth, soildepth)

#-------------------------------------------------------------------#
#######   Shallowest Soil Detpth for channel segment 
#-------------------------------------------------------------------#
stmlineras="stmlineras"

#arcpy.CheckOutExtension ("3D")
#arcpy.AddSurfaceInformation_3d (streamnet, soildepth, "Z_MIN", "CONFLATE_NEAREST" , "", "", "", "") 

if arcpy.Exists(stmlineras):
    print('Table already exists, delete and create new')
    arcpy.Delete_management(stmlineras)

arcpy.PolylineToRaster_conversion (streamnet, "arcid", stmlineras, "", "", env.cellSize)
soild_table="soild_zonal"

if arcpy.Exists(soild_table):
    arcpy.Delete_management(soild_table)
    print('zonal statistics table already exists, delete and create new')

ZonalStatisticsAsTable (stmlineras, "VALUE", soildepth, soild_table, "DATA", "MINIMUM")

arcpy.JoinField_management(streamnet,"arcid",soild_table ,"VALUE","MIN")

########################################
# When channel lenght is smaller than grid cell size, there're chances that the segment will not be
# represented in the stream to raster conversion, therefore no seg-depth value generated for the arc.
# To avoid such problem, find min soil depth at all feature points of an arc to replace null values

all_nodes=path+"output.gdb/all_nodes"
arcpy.FeatureVerticesToPoints_management(streamnet, all_nodes, "ALL")
ExtractMultiValuesToPoints(all_nodes, [[soildepth, "soildepth"]], "BILINEAR")

arcpy.Statistics_analysis(all_nodes,"all_nodes_statistics","soildepth MIN","arcid")
arcpy.JoinField_management(streamnet,"arcid","all_nodes_statistics","arcid","MIN_soildepth")

arcpy.AddField_management (streamnet, "segdepth", "FLOAT", 8, 2)
fields=['MIN_soildepth','MIN','segdepth']

with arcpy.da.UpdateCursor(streamnet, fields) as cursor:
    for row in cursor:
        if row[1] is None:
            row[2]=row[0]
        elif (row[0] < row[1]):
            row[2]=row[0]
        elif (row[0] > row[1]):
            row[2]=row[1]
        else:
            row[2]=row[0]
        cursor.updateRow(row)

arcpy.DeleteField_management(streamnet, "MIN")
arcpy.DeleteField_management(streamnet, "MIN_soild")
arcpy.Delete_management(all_nodes)
arcpy.Delete_management("all_nodes_statistics")

arcpy.CalculateField_management (streamnet, "effdepth", "0.95*float(!segdepth!)", "PYTHON_9.3","")

#-------------------------------------------------------------------#
#######     Greate stream network file for DHSVM 
#-------------------------------------------------------------------#
    
if os.path.exists("stream.network.dat"):
    os.remove("stream.network.dat")
    print('stream.network.dat sucessfully deleted')

print('creating new stream.network.dat file')


# Select certain fields from table to a new numpy array

arr_export=arcpy.da.TableToNumPyArray(streamnet,('arcid','segorder','slope','Shape_Length','chanclass','downarc'))

np.savetxt(path+'stream.network.dat', arr_export,fmt="%5d %3d %11.5f %17.5f %3d %7d")

#-------------------------------------------------------------------#
#######     Create stream map file for DHSVM
#-------------------------------------------------------------------#

print('running wshdslope')

from wshdslope import wshdslopefun
    
slope="gridslope_4d"
aspect="gridaspect"

wshdslopefun(elev, slope, aspect, flowacc)

from rowcolmap import rowcolmapfun

print('running rowcolmap')

rowcolmapfun(elev,aspect,slope,"maskras")

print('creating streammap')

if arcpy.Exists('/output.gdb/outcover.shp'):
    arcpy.Delete_management('/output.gdb/outcover.shp')

outcover="./output.gdb/outcover"
arcpy.Intersect_analysis([streamnet,"rowcolpoly.shp"], outcover, "", "", "")

print('running roadaspect')

from roadaspect import roadaspectfun 

roadaspectfun(outcover)

print('generating stream map file')

from streammapfile import streammapfilefun

streammapfilefun(outcover,path)

#-------------------------------------------------------------------#
#######     Clean Up Tmp files 
#-------------------------------------------------------------------#

print('Starting Cleanning Up')

arcpy.Delete_management(flowdir)
arcpy.Delete_management(flowacc)
arcpy.Delete_management(slope)
arcpy.Delete_management(aspect)
arcpy.Delete_management(tmpelev)


arcpy.Delete_management(nodestart)
arcpy.Delete_management(nodeend)


arcpy.Delete_management(tmpacc)
arcpy.Delete_management(stmlineras)







