# SUMMARY:      channelclass.py
# USAGE:        Part of python version createstreamnetwork. Classify
#		channel based on slope and average contributing area.
#		Requires mannual modification of channel classification   
#		criteria.
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts createstreamnetwork.aml as part of DHSVM
#
# Last Change: 2017-08-10 by zduan

import arcpy
from arcpy import env
from arcpy.sa import *
import arcgisscripting

def channelclassfun(streamnet):
    
    fields=['slope','meanmsq','chanclass','hyddepth','hydwidth','effwidth']
    
    with arcpy.da.UpdateCursor(streamnet, fields) as cursor:

        # row[0] - 'slope'
        # row[1] - 'meanmsq' - mean contributing area to channel segment
        # row[2] - 'chanclass' - channel class reference number
        # row[3] - 'hyddepth'
        # row[4] - 'hydwidth'
        # row[5] - 'effwitdth'
        
        for row in cursor:
            if (row[0] <= 0.002 and row[1] <= 1000000):
                row[2] = 1
                row[3] = 0.5
                row[4] = 0.03
                row[5] = 0.06
            elif (row[0] <= 0.002 and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 2
                row[3] = 1.0
                row[4] = 0.03
                row[5] = 0.09
            elif (row[0] <= 0.002 and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 3
                row[3] = 2.0
                row[4] = 0.03
                row[5] = 0.12
            elif (row[0] <= 0.002 and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 4
                row[3] = 3.0
                row[4] = 0.03
                row[5] = 0.15
            elif (row[0] <= 0.002 and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 5
                row[3] = 4.0
                row[4] = 0.03
                row[5] = 0.18
            elif (row[0] <= 0.002 and row[1] > 40000000) :
                row[2] = 6
                row[3] = 4.5
                row[4] = 0.03
                row[5] = 0.21
            elif ((row[0] > 0.002 and row[0] <= 0.1) and row[1] <= 1000000):
                row[2] = 7
                row[3] = 0.5
                row[4] = 0.05
                row[5] = 0.1
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 8
                row[3] = 1.0
                row[4] = 0.05
                row[5] = 0.15
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 9
                row[3] = 2.0
                row[4] = 0.05
                row[5] = 0.2
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 10
                row[3] = 3.0
                row[4] = 0.05
                row[5] = 0.25
            elif ((row[0] > 0.002 and row[0] <= 0.1) and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 11
                row[3] = 4.0
                row[4] = 0.05
                row[5] = 0.3
            elif ((row[0] > 0.002 and row[0] <= 0.1) and row[1] > 40000000):
                row[2] = 12
                row[3] = 4.5
                row[4] = 0.05
                row[5] = 0.35
            elif (row[0] > 0.01 and row[1] <= 1000000):
                row[2] = 13
                row[3] = 0.5
                row[4] = 0.1
                row[5] = 0.2
            elif (row[0] > 0.01 and (row[1] > 1000000 and row[1] <= 10000000)):
                row[2] = 14
                row[3] = 1.0
                row[4] = 0.1
                row[5] = 0.3
            elif (row[0] > 0.01 and (row[1] > 10000000 and row[1] <= 20000000)):
                row[2] = 15
                row[3] = 2.0
                row[4] = 0.1
                row[5] = 0.4
            elif (row[0] > 0.01 and (row[1] > 20000000 and row[1] <= 30000000)):
                row[2] = 16
                row[3] = 3.0
                row[4] = 0.1
                row[5] = 0.5
            elif (row[0] > 0.01 and (row[1] > 30000000 and row[1] <= 40000000)):
                row[2] = 17
                row[3] = 4.0
                row[4] = 0.1
                row[5] = 0.6
            elif (row[0] > 0.01 and row[1] > 40000000) :
                row[2] = 18
                row[3] = 4.5
                row[4] = 0.1
                row[5] = 0.7

            # Update the cursor with the updated list
            cursor.updateRow(row)
