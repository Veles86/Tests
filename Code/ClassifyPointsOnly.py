import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, GetParameterAsText
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField


from Utils.Generalization import sort_segments_matlab
from Utils.Generalization import classify_points, sort_segments, generalize_sorted_segments, fill_missing_data_for_generalized_lines
from Utils.Lines import fill_dH_and_dist_differences


arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True




def update_regular_points_to_junction_by_input(input_classification_results, points_classification_output) -> None:

    # Create an empty list to store the values from the "Code" field
    classification_count_dict = {}

    # Use a SearchCursor to iterate through the records in input_classification_results
    with SearchCursor(input_classification_results, ["Code", "Count"]) as cursor:
        for row in cursor:
            # Add the value from the "Code" field as the key and "Count" field as the value to the dictionary
            classification_count_dict[row[0]] = row[1]


    count = 0
    with UpdateCursor(points_classification_output, ["Code","Count","Type"]) as cursor:
        for row in cursor:
            if row[0] in classification_count_dict:
                new_count = classification_count_dict[row[0]] + row[1]
                if row[1] == 2 and new_count > 2:
                    count += 1
                    row[1] = new_count
                    row[2] = "junction"
                    cursor.updateRow(row)
                
                

        AddMessage(f"Updated the original classification data and turned {count} regular points into junctions.")

if __name__ == '__main__':



    # Get input parameters from userEnforced additional junctions o
    levelling_segments_input = GetParameterAsText(0)  # Input feature class of selected segments
    source_points_input = GetParameterAsText(1)  # Input feature class of control points
    input_classification_results = GetParameterAsText(2)  # Output feature class for classified points
    points_classification_output = GetParameterAsText(3)  # Output feature class for sorted segments


    classify_points(levelling_segments_input, source_points_input, points_classification_output, print_messages = True)
    if input_classification_results:
        update_regular_points_to_junction_by_input(input_classification_results, points_classification_output)

