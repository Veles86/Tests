
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
script_dir = os.path.dirname(os.path.abspath(__file__))
utils_path = os.path.join(script_dir, '..', 'utils')
sys.path.append(os.path.abspath(".."))

from Utils.Points import create_points_dictionary, create_points_dictionary, get_distance_between_points,extract_coordinates_from_external_source, get_point_names

from Utils.Lines import reverse_segments_new, reverse_line_segments_without_copy,reverse_line_segments_from_tuple, add_segment_from_tuple, get_single_segment_fc, fill_dH_and_dist_differences, reverse_line_segments, get_first_segment_with_point, add_segment, delete_segment,get_last_segment_end_point_code, get_last_segment
from Utils.Helpers import count_features, mmyy_int_to_yyyymm_str,start_editing, stop_editing, get_feature_layer_workspace
from Utils.Generalization import check_self_loops_new


arcpy.env.overwriteOutput = True

def break_self_loop(sorted_segments, point_coordinates, line_num,start_point) -> int:
    """
    Recieving a set of segments that creating a self-loop and breaking it at one point so the would be loop turns to 2 lines
    The breakage point will be set according to 2 criterias: the participating point's ranks and their distance from the start point:
    If one (and only one) of the participating points (without the starting point) has rank 1 or 2, the breakage point will be that point
    otherwise, the breakage point will be the farthest point from the starting point


    Parameters
    ----------
    sorted_segments : FeatureLayer
        Feature layer of sorted segments
    point_coordinates : FeatureLayer
        Feature layer of control points
    line_num : int
        Line number of the generalized line
    start_point : int
        Starting point of the generalized line

    Returns
    -------
    breakage_point_code : int
        Code of the breakage point

    
    """
    breakage_point_code = None

    lines_segments = SelectByAttribute(sorted_segments, where_clause=f"line_num = {line_num}", selection_type="NEW_SELECTION")
    num_of_segments = count_features(lines_segments)
    participating_points = extract_coordinates_from_external_source(lines_segments, point_coordinates)
    points_dictionary = create_points_dictionary(participating_points)

    high_ranked_points = []
    farthest_point = None
    max_distance = 0
    start_point_data = points_dictionary[start_point]
    for point_code in points_dictionary:
        if point_code != start_point:
            point_data = points_dictionary[point_code]
            if point_data["HRank"] in [1, 2, 11, 12]:
                high_ranked_points.append(point_code)
            temp_dist = get_distance_between_points(start_point_data["X"], start_point_data["Y"], point_data["X"], point_data["Y"])
            if temp_dist > max_distance:
                max_distance = temp_dist
                farthest_point = point_code



    if len(high_ranked_points) ==1:
        breakage_point_code = high_ranked_points[0]
    else:
        breakage_point_code = farthest_point


   
    break_objectID = 0
    break_num = 0

    with SearchCursor(lines_segments, ["OBJECTID", "UniqueRowNum", "StartPointCode","EndPointCode"]) as cursor:
        for row in cursor:
            if row[3] == breakage_point_code:
                break_objectID = row[0]
                break_num = row[1]
                break

    
    increase_line_num = False
    SelectByAttribute(in_layer_or_view=sorted_segments,selection_type="CLEAR_SELECTION")
    with UpdateCursor(sorted_segments, ["OBJECTID", "UniqueRowNum","line_num"]) as cursor:
        for row in cursor:

            if increase_line_num:
                row[2] = row[2]+1
                cursor.updateRow(row)
            elif row[0] == break_objectID and row[1] == break_num:
                increase_line_num = True

    return breakage_point_code


def check_self_loops_new(sorted_segments, point_coordinates) -> None:
    """
    Checking if the sorted segments data will have self-loops and breaking them at one point so the would be loop turns to 2 lines

    Parameters
    ----------
    sorted_segments : FeatureLayer
        Feature layer of sorted segments
    point_coordinates : FeatureLayer
        Feature layer of control points

    
    """
    num_of_lines = 0
    num_of_self_loops = 0
    points_dictionary = create_points_dictionary(point_coordinates)
    self_looping_points = []
    breakage_points = []
    self_looping_lines = []




    temp_statistics = Statistics(in_table=sorted_segments,out_table=r"memory\temp_statistics",statistics_fields="StartPointCode FIRST;EndPointCode LAST",case_field="line_num",concatenation_separator="")
    num_of_lines = count_features(temp_statistics)
    self_looping_lines = SelectByAttribute(in_layer_or_view=temp_statistics,selection_type="NEW_SELECTION",where_clause="FIRST_StartPointCode = LAST_EndPointCode",invert_where_clause=None)
    num_of_self_loops = count_features(self_looping_lines)
    AddMessage(f"The initial data has {num_of_lines} lines and {num_of_self_loops} self-looping lines")
    if num_of_self_loops > 0:
        AddMessage(f"Breaking self-loops...")
        num_of_fixed_lines = 0
        with SearchCursor(self_looping_lines, ["line_num", "FIRST_StartPointCode"]) as cursor:
            for row in cursor:
                self_looping_lines.append(row[0]+num_of_fixed_lines)
                start_point = row[1]
                self_looping_points.append(start_point)
                AddMessage(f"line {row[0]+num_of_fixed_lines} is self-looping")
                temp_point = break_self_loop(sorted_segments, point_coordinates, row[0]+num_of_fixed_lines ,start_point)
                breakage_points.append(temp_point)
                num_of_fixed_lines += 1

        
        for i in range(num_of_self_loops):
            current_line = self_looping_lines[i] + num_of_fixed_lines
            self_loop_point_code = self_looping_points[i]
            self_loop_point_data = points_dictionary[self_loop_point_code]
            breake_point_code = breakage_points[i]
            breake_point_data = points_dictionary[breake_point_code]
            AddMessage(f"   The self loop around point {self_loop_point_data['Name']} ({self_loop_point_code}) was broken at point {breake_point_data['Name']} ({breake_point_code})")

    SelectByAttribute(in_layer_or_view=sorted_segments,selection_type="CLEAR_SELECTION")

    del points_dictionary   


if __name__ == '__main__':

    # Get input parameters from user
    #sorted_segments = arcpy.GetParameterAsText(0) 
    #point_coordinates = arcpy.GetParameterAsText(1) 

    sorted_segments = r"E:\git\LevellingThesis\Project\Thesis.gdb\common_line_segments_recent_luba_sorted_for_generalization"
    points_coordinates = r"E:\git\LevellingThesis\Project\Thesis.gdb\common_line_segments_recent_luba_points_data"

    check_self_loops_new(sorted_segments, points_coordinates)

    DeleteField(in_table=sorted_segments,drop_field="new_order;reverse_flag;ORIG_FID",method="DELETE_FIELDS")

