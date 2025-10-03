import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine
import os
import sys
sys.path.append(os.path.abspath(".."))

from Utils.Generalization import check_lines_connectivity, check_proper_division

arcpy.env.overwriteOutput = True

if __name__ == '__main__':

    # Get input parameters from user
    sorted_segments = arcpy.GetParameterAsText(0)
    points_classification = arcpy.GetParameterAsText(1)
    output_table = arcpy.GetParameterAsText(2) 

    check_proper_division(sorted_segments, points_classification)
    check_lines_connectivity(sorted_segments, output_table)
