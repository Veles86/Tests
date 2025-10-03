from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField
import os
import sys
sys.path.append(os.path.abspath(".."))
#from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace, count_features
#from Utils.Gravimetry import get_theoretic_normal_gravity
from Utils.Gravimetry import calculate_theoretic_normal_gravity

"""
def update_theoretical_grav_value(points_fc):

    fields = [field.name for field in ListFields(points_fc)]
    if "Latitude" not in fields or "NormG" not in fields:
        AddError("Required fields 'Latitude' or 'NormG' do not exist in the feature class.")
        return
    

    # Select features where "Phi" is not null
    where_clause = "Latitude IS NOT NULL"
    selected_features = SelectByAttribute(points_fc, "SUBSET_SELECTION", where_clause)
    
    count = count_features(selected_features)
    if count == 0:
        AddWarning("No features with non-null 'Phi' values found.")
        return
    else:
        AddMessage(f"{count} features will be updated with theoretical G values.")
    
    wc = get_feature_layer_workspace(points_fc)
    editor = start_editing(wc)
    counter = 0
    with UpdateCursor(points_fc, ["NormG", "Latitude"]) as cursor:
        for row in cursor:
            row[0] = get_theoretic_normal_gravity(row[1])
            cursor.updateRow(row)
            counter += 1
    stop_editing(editor)
    AddMessage(f"Updated {counter} features with theoretical G values.")
"""
if __name__ == "__main__":
    # --- 4. Parameters when running as a script or script tool -------------
    # For a script tool, you might set these as parameters in ArcGIS Pro:
    points_fc = GetParameterAsText(0)

    calculate_theoretic_normal_gravity(points_fc)
    #update_theoretical_grav_value(points_fc)