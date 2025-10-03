import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from Utils.Helpers import extract_coordinates_from_external_source,start_editing, stop_editing, reverse_line_segments, get_first_segment_with_point, add_segment, delete_segment,get_last_segment_end_point_code
from Utils.Helpers import get_segments_end_point_names,get_segments_end_point_codes, check_self_loops, create_points_dictionary, count_features, mmyy_int_to_yyyymm_str, fill_midpoints_and_dates, fill_dH_and_dist_differences


arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True




def main_function(p1,p2,p3) -> None:
    # Main script execution




if __name__ == '__main__':

    # Get input parameters from user
    p1 = arcpy.GetParameterAsText(0) 
    p2 = arcpy.GetParameterAsText(1) 
    p3 = arcpy.GetParameterAsText(2) 

  




    main_function(p1,p2,p3)


