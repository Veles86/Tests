from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace, count_features
from Utils.Lines import turn_segments_by_min_point_name_old_format

def turn_lines_by_smallest_point(lines_fc):

    turn_segments_by_min_point_name_old_format(lines_fc)





if __name__ == "__main__":

    lines_fc = GetParameterAsText(0)
    turn_lines_by_smallest_point(lines_fc)