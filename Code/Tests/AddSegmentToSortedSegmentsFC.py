import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine



arcpy.env.overwriteOutput = True



def add_reversed_segment(segment, target_feature_class):
    with arcpy.da.InsertCursor(target_feature_class, ["SHAPE@"]) as cursor:
        cursor.insertRow([segment])

def print_polyline_coords(polyline, label="Polyline"):
    AddMessage(f"\n{label} - Parts: {polyline.partCount}")
    for i, part in enumerate(polyline):
        AddMessage(f"  Part {i + 1}:")
        for pt in part:
            if pt:
                AddMessage(f"    ({pt.X}, {pt.Y})")

def add_segments_to_sorted_fc(segments_to_add, target_feature_class) -> None:

    with arcpy.da.UpdateCursor(feature_class, ["OID@", "SHAPE@"]) as cursor:
        for oid, geom in cursor:
            AddMessage(f"\nFeature OID: {oid}")
            print_polyline_coords(geom, "Original")

            #reversed_geom = reverse_polyline_geometry(geom)
            print_polyline_coords((geom), "Reversed")
            #print_polyline_coords(reversed_geom, "Reversed")

            # Optionally update
            #cursor.updateRow([oid, reversed_geom])




if __name__ == '__main__':

    # Get input parameters from user
    segments_to_add = arcpy.GetParameterAsText(0) 
    target_feature_class = arcpy.GetParameterAsText(1) 


    add_segments_to_sorted_fc(segments_to_add, target_feature_class)

