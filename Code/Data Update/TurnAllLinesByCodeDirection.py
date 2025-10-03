from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace, count_features
from Utils.Lines import reverse_segments,reverse_segments_new, reverse_generalized_lines

def turn_lines_by_code_direction(lines_fc):

    # Check if SegmentDirection field exists
    fields = [f.name for f in ListFields(lines_fc)]
    if "SegmentDirection" in fields:
        AddMessage("Reversing segments with SegmentDirection = -1")
        lines_to_reverse = SelectByAttribute(lines_fc, selection_type="SUBSET_SELECTION", where_clause="SegmentDirection = -1")
        num_before = count_features(lines_to_reverse)
        reverse_segments_new(lines_to_reverse)
        SelectByAttribute(lines_to_reverse, selection_type="CLEAR_SELECTION")

        AddMessage(f"Reversed {num_before} segments.")
    elif "LineDirection" in fields:
        AddMessage("Reversing lines with LineDirection = -1")
        lines_to_reverse = SelectByAttribute(lines_fc, selection_type="SUBSET_SELECTION", where_clause="LineDirection = -1")
        num_before = count_features(lines_to_reverse)

        reverse_generalized_lines(lines_to_reverse)
        SelectByAttribute(lines_to_reverse, selection_type="CLEAR_SELECTION")

        AddMessage(f"Reversed {num_before} lines.")

    else:
        AddError("No field found to indicate the direction of the lines or segments.")
        return
    





if __name__ == "__main__":

    lines_fc = GetParameterAsText(0)
    turn_lines_by_code_direction(lines_fc)