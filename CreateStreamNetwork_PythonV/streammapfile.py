# SUMMARY:      streammapfile.py
# USAGE:        Part of python version createstreamnetwork. Assign
#		row and col information to road/channel coverage   
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts as part of DHSVM
#
# Last Change: 2017-08-10 by zduan


import arcpy
import numpy as np
import os
import datetime

def streammapfilefun(outcover,path):

    if os.path.exists("stream.map.dat"):
        os.remove("stream.map.dat")
        print('stream.map.dat sucessfully deleted')

    print('creating new stream.map.dat file')
        
    info = '###### This file has been automatically generated #####'
    info += '\n######             EDIT WITH CARE!!!              #####'
    info += '\n# Generated: '+str(datetime.datetime.now())
    info += '\n# Created by streammapfile.py'
    info += '\n# Workspace '+str(path)
    info += '\n#                   Segment  Cut/Bank     Cut     Segment'
    info += '\n#  Col  Row  ID      Length   Height     Width     Aspect   SINK?'
    info += '\n#                     (m)      (m)        (m)       (d)    (optional)'
    info += '\n# \n'

    arr_export=arcpy.da.TableToNumPyArray(outcover,('COL','ROW','arcid','Shape_Length','effdepth','effwidth','rd_aspect'))

##    with open(path+'stream.map.dat', 'w') as f:
##          f.write(info)
##          np.savetxt(f, arr_export,fmt="%5d %4d %4d %11.4f %10.4f %9.4f %10.4f")

    with open(path+'stream.map.dat', 'w') as f:
        f.write(info)
        np.savetxt(f, arr_export,fmt="%5d %5d %5d %11.4f %10.4f %9.4f %10.4f")
