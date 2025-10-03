from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields, Array, Polyline
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor

import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace

import arcpy




def reconstruct_segments_geometry(points_fc, segments_fc):
    # Create a dictionary to map point codes to their geometry
    point_geom_dict = {}
    with arcpy.da.SearchCursor(points_fc, ["Code", "SHAPE@XY"]) as cursor:
        for code, xy in cursor:
            point_geom_dict[code] = xy

    # Use an update cursor to edit the line geometries
    with arcpy.da.UpdateCursor(segments_fc, ["StartPointCode", "EndPointCode", "SHAPE@"]) as cursor:
        for start_code, end_code, _ in cursor:
            if start_code in point_geom_dict and end_code in point_geom_dict:
                start_xy = point_geom_dict[start_code]
                end_xy = point_geom_dict[end_code]

                # Create a new line geometry
                array = arcpy.Array([
                    arcpy.Point(*start_xy),
                    arcpy.Point(*end_xy)
                ])
                new_line = arcpy.Polyline(array)
                cursor.updateRow([start_code, end_code, new_line])
            else:
                arcpy.AddWarning(f"Missing point(s) for Start: {start_code}, End: {end_code}")


if __name__ == "__main__":

    points_fc = GetParameterAsText(0)
    segments_fc = GetParameterAsText(1)


    reconstruct_segments_geometry(points_fc,segments_fc)
