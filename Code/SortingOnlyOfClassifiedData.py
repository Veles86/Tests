import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, GetParameterAsText
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField


from Utils.Generalization import sort_segments_matlab
from Utils.Generalization import classify_points, sort_segments, generalize_sorted_segments, fill_missing_data_for_generalized_lines
from Utils.Lines import fill_dH_and_dist_differences
from Utils.Helpers import count_features

arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True

def points_classification_to_tuple(points_classification_input) -> tuple:
    with SearchCursor(points_classification_input, ["Code","Type"]) as cursor:
        junction_points_codes = []
        regular_points_codes = []
        dead_end_points_codes = []
        for row in cursor:
            if row[1] == 'junction':
                junction_points_codes.append(row[0])
            elif row[1] == 'regular':
                regular_points_codes.append(row[0])
            elif row[1] == 'dead-end':
                dead_end_points_codes.append(row[0])

    return junction_points_codes,regular_points_codes,dead_end_points_codes



if __name__ == '__main__':





    METHOD = "ARCGIS PRO"
    #METHOD = "DEBUGGING"

    if METHOD == "DEBUGGING":

        gdb_path = fr"E:\git\LevellingThesis\Project\Thesis.gdb"
        levelling_segments_input = gdb_path + fr"\common_line_segments_recent_luba"
        source_points_input =  gdb_path + fr"\common_line_segments_recent_luba_points_data"
        points_classification_input = gdb_path + fr"\common_line_segments_recent_luba_points_classification"
        sorted_segments_output = gdb_path + fr"\segments_sorted_for_generalization"

        AddWarning(f"The script is running in DEBUGGING mode. Please ensure that the input data is correct and the output paths are valid.")

    elif METHOD == "ARCGIS PRO":
    # Get input parameters from user
        levelling_segments_input = GetParameterAsText(0)  # Input feature class of selected segments
        source_points_input = GetParameterAsText(1)  # Input feature class of control points
        points_classification_input = GetParameterAsText(2)  # Output feature class for classified points
        sorted_segments_output = GetParameterAsText(3)  # Output feature class for sorted segments
    else:
        AddError("Invalid method specified. Please choose either 'DEBUGGING' or 'ARCGIS PRO'.")
        raise ValueError("Invalid method specified. Please choose either 'DEBUGGING' or 'ARCGIS PRO'.")
    
    
    junction_points_codes, regular_points_codes, dead_end_points_codes = points_classification_to_tuple(points_classification_input)
    sort_segments_matlab(levelling_segments_input, source_points_input, junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output, print_messages = True)

   


