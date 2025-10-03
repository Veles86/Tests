from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField
import os
import sys
sys.path.append(os.path.abspath(".."))
#from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace, count_features
#from Utils.Gravimetry import get_theoretic_normal_gravity
from Utils.Lines import enumerate_lines_by_code


if __name__ == "__main__":
    # --- 4. Parameters when running as a script or script tool -------------
    # For a script tool, you might set these as parameters in ArcGIS Pro:
    lines_fc = GetParameterAsText(0)

    enumerate_lines_by_code(lines_fc)
