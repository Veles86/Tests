import arcpy
from arcpy import GetParameterAsText,ListFields, Describe
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, GetParameterAsText
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField


import os
import sys
sys.path.append(os.path.abspath(".."))

from Utils.Configs import CONTROL_POINTS_FULL_DATA


def extract_used_points_data(lines_segments_input, source_points_input, points_data_output) -> None:

    """Extracts unique control points used in the provided lines or segments feature class and saves them to a new feature class.
    
    Args:
        lines_segments_input (str): The input feature class of selected lines or segments.
        source_points_input (str): The input feature class of control points.
        points_data_output (str): The output feature class for the extracted points.
    
    Returns:
        None
    
    """

    if not source_points_input:
        source_points_input = CONTROL_POINTS_FULL_DATA

    fields = ListFields(lines_segments_input)
    field_names = [f.name for f in fields]

    if "MidPointCodes" in field_names:
        generalized_lines = True
    else:
        generalized_lines = False

    # Create a set to store unique codes
    unique_codes = set()


    if generalized_lines:
        with SearchCursor(lines_segments_input, ["StartPointCode", "EndPointCode", "MidPointCodes"]) as cursor:
            for row in cursor:
                unique_codes.add(row[0])
                unique_codes.add(row[1])

                if row[2]:
                    string_of_codes = row[2]
                    codes_list = [code.strip() for code in string_of_codes.split(",")]
                    for code in codes_list:
                        unique_codes.add(int(code))

    else: # If not generalized lines hence simple segments, use only StartPointCode and EndPointCode

        with SearchCursor(lines_segments_input, ["StartPointCode", "EndPointCode"]) as cursor:
            for row in cursor:
                unique_codes.add(row[0])
                unique_codes.add(row[1])

    # Convert set to a list and format for SQL IN clause
    code_list = ", ".join(map(str, unique_codes))

    # Create a query to filter ControlPoints
    query = f"Code IN ({code_list})"



    filtered_layer = SelectByAttribute(source_points_input, "NEW_SELECTION", query)

    CopyFeatures(filtered_layer, points_data_output)
    SelectByAttribute(in_layer_or_view=filtered_layer,selection_type="CLEAR_SELECTION")





if __name__ == '__main__':

    lines_segments_input = GetParameterAsText(0)  # Input feature class of selected segments
    source_points_input = GetParameterAsText(1)  # Input feature class of control points
    points_data_output = GetParameterAsText(2)  # Output feature class for generalized polylines

    extract_used_points_data(lines_segments_input, source_points_input, points_data_output)