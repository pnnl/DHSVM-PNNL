# SUMMARY:      roadaspect.py
# USAGE:        Part of python version createstreamnetwork. Computes
#		the aspect and length of stream arc within each DEM 
#		cell.
# ORG:          Pacific Northwest National Laboratory
# E-MAIL:       zhuoran.duan@pnnl.gov
# ORIG-DATE:    Apr-2017
# DESCRIPTION:  Python version of original createstreamnetwork aml
# DESCRIP-END.
# COMMENTS:     This python script is created based on original 
#		AML scripts roadaspect.aml as part of DHSVM
#
# Last Change: 2017-08-10 by zduan


import arcpy
import math

def roadaspectfun(outcover):
    arcpy.AddField_management (outcover, "rd_aspect", "FLOAT", 12, 5)
    arcpy.AddField_management (outcover, "rd_st_len", "FLOAT", 12, 5)
    arcpy.AddField_management (outcover, "rd_efflen", "FLOAT", 12, 5)
    
    expression1 = "GetGeographicalDegrees( !SHAPE! )"
    
    codeblock1 = """def GetGeographicalDegrees(shape):				
                      radian = math.atan2(shape.firstpoint.x - shape.lastpoint.x ,
                                          shape.firstpoint.y - shape.lastpoint.y)
		      rad2deg= 180.0 / math.pi
                      degrees = round((radian * rad2deg + 360)%360)
                      return degrees """

    arcpy.CalculateField_management (outcover, "rd_aspect", expression1, "PYTHON_9.3", codeblock1)

    expression2 = "CalcStLen( !SHAPE!)"
    codeblock2 = """def CalcStLen(shape):
                      len = math.sqrt(math.pow(shape.lastpoint.y - shape.firstpoint.y,2)+
                          math.pow(shape.lastpoint.x - shape.firstpoint.x,2))
                      return len """

    arcpy.CalculateField_management (outcover, "rd_st_len", expression2, "PYTHON_9.3", codeblock2)

    expression3 = "CalcEffectLen( !SHAPE!,!rd_st_len! ,!rd_aspect! )"
    codeblock3 = """def CalcEffectLen(shape,len,rd_aspect):
                      eff_len = math.fabs( len * math.sin(rd_aspect* math.pi /180 )) #convert from degree to radius
                      return eff_len """

    arcpy.CalculateField_management (outcover, "rd_efflen", expression3, "PYTHON_9.3", codeblock3)
    
