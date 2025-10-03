import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, GetParameterAsText
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField


from Utils.Generalization import sort_segments_matlab
from Utils.Generalization import classify_points, sort_segments, generalize_sorted_segments, fill_missing_data_for_generalized_lines
from Utils.Lines import fill_dH_and_dist_differences, reverse_generalized_lines


arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True


def generalize_segments(levelling_segments_input, source_points_input, points_classification_output, sorted_segments_output, generalized_lines_output) -> None:
    # Main script execution

    AddMessage("Starting the generalization process...")

    # Step 1: Classify the edge points (junction, regular, dead-end) using codeA and codeB, and create output points feature class
    junction_points_codes, regular_points_codes, dead_end_points_codes = classify_points(levelling_segments_input, source_points_input, points_classification_output, print_messages = True)
    
    # dict = create_points_dictionary(point_coordinates)
    # Step 2: Sort and order the levelling segments between junction points with addition of a field with the final generalized line number
    sorted_segments = sort_segments_matlab(levelling_segments_input, source_points_input, junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output, print_messages = True)

    
    # Step 3: Generalize the sorted segments between junction points and save the output feature class
    generalized_lines = generalize_sorted_segments(sorted_segments, generalized_lines_output, print_messages=True)


    #lines_to_reverse = SelectByAttribute(in_layer_or_view=generalized_lines, where_clause="LineDirection = -1", selection_type="NEW_SELECTION")
    #reverse_generalized_lines(lines_to_reverse)
    #SelectByAttribute(in_layer_or_view=generalized_lines, selection_type="CLEAR_SELECTION")
    #AddMessage("Turned all lines in the lines' code direction")


    # Step 4: Check for unique "UniqueRowNum" values in the generalized_lines table
    unique_row_nums = set()
    duplicates = set()
    
    with SearchCursor(generalized_lines, ["UniqueRowNum"]) as cursor:
        for row in cursor:
            if row[0] in unique_row_nums:
                duplicates.add(row[0])
            else:
                unique_row_nums.add(row[0])
    
    if duplicates:
        AddWarning(f"Duplicate UniqueRowNum values found: {duplicates}")


    #TODO
    #fill_missing_data_for_generalized_lines(generalized_lines, source_points_input)

    AddMessage(f"Generalization complete.")
    
    

if __name__ == '__main__':


    #METHOD = "DEBUGGING"

    METHOD = "ARCGIS PRO"

    if METHOD == "DEBUGGING":

        gdb_path = fr"E:\git\LevellingThesis\Data\_Testing\GeneralizationQA.gdb"
        levelling_segments_input = gdb_path + fr"\Luba_segments"
        source_points_input = gdb_path + fr"\Luba_points_data"
        points_classification_output = gdb_path + fr"\Luba_points_classification"
        sorted_segments_output = gdb_path + fr"\Luba_sorted_for_generalization"
        generalized_lines_output = gdb_path + fr"\Luba_generalized"

        AddWarning(f"The script is running in DEBUGGING mode. Please ensure that the input data is correct and the output paths are valid.")

    elif METHOD == "ARCGIS PRO":
        # Get input parameters from user
        levelling_segments_input = GetParameterAsText(0)  # Input feature class of selected segments
        source_points_input = GetParameterAsText(1)  # Input feature class of control points
        points_classification_output = GetParameterAsText(2)  # Output feature class for classified points
        sorted_segments_output = GetParameterAsText(3)  # Output feature class for sorted segments
        generalized_lines_output = GetParameterAsText(4)  # Output feature class for generalized polylines
    else:
        AddError("Invalid method specified. Please choose either 'DEBUGGING' or 'ARCGIS PRO'.")
        raise ValueError("Invalid method specified. Please choose either 'DEBUGGING' or 'ARCGIS PRO'.")


    generalize_segments(levelling_segments_input, source_points_input, points_classification_output, sorted_segments_output, generalized_lines_output)


