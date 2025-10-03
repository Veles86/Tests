
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, CreateTable_management as CreateTable, ListFields
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField,EnableEditorTracking, SelectLayerByLocation as SelectByLocation
from arcpy.analysis import Statistics, Frequency

from collections import Counter
from collections import defaultdict


from arcpy.edit import FlipLine
import os
import sys

from Utils.Generalization import check_self_loops_new
from Utils.Helpers import count_features

arcpy.env.overwriteOutput = True




if __name__ == '__main__':

    # Get input parameters from user
    #sorted_segments = arcpy.GetParameterAsText(0) 
    #point_coordinates = arcpy.GetParameterAsText(1) 

    sorted_segments = r"E:\git\LevellingThesis\Project\Thesis.gdb\common_line_segments_recent_luba_sorted_for_generalization"
    points_coordinates = r"E:\git\LevellingThesis\Project\Thesis.gdb\common_line_segments_recent_luba_points_data"

    check_self_loops_new(sorted_segments, points_coordinates)

   
