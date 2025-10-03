import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine



arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True


def reverse_polyline_geometry(polyline):
    reversed_parts = []
    # Reverse the order of parts AND reverse each part's vertices
    for part in reversed(list(polyline)):
        reversed_part = arcpy.Array([pt for pt in reversed(part) if pt])
        reversed_parts.append(reversed_part)
    return arcpy.Polyline(arcpy.Array(reversed_parts), polyline.spatialReference)



def reverse_polyline_geometry_old(polyline):
    reversed_parts = []
    for part in polyline:
        reversed_part = arcpy.Array([pt for pt in reversed(part) if pt])
        reversed_parts.append(reversed_part)
    return arcpy.Polyline(arcpy.Array(reversed_parts), polyline.spatialReference)

def print_polyline_coords(polyline, label="Polyline"):
    AddMessage(f"\n{label} - Parts: {polyline.partCount}")
    for i, part in enumerate(polyline):
        AddMessage(f"  Part {i + 1}:")
        for pt in part:
            if pt:
                AddMessage(f"    ({pt.X}, {pt.Y})")

def test_polylines(feature_class) -> None:
    FlipLine(feature_class)
    with arcpy.da.UpdateCursor(feature_class, ["OID@", "SHAPE@"]) as cursor:
        for oid, geom in cursor:
            AddMessage(f"\nFeature OID: {oid}")
            print_polyline_coords(geom, "Original")

            reversed_geom = reverse_polyline_geometry(geom)
            #print_polyline_coords((geom), "Reversed")
            print_polyline_coords(reversed_geom, "Reversed")

            # Optionally update
            #cursor.updateRow([oid, reversed_geom])




if __name__ == '__main__':

    # Get input parameters from user
    feature_class = arcpy.GetParameterAsText(0) 


    test_polylines(feature_class)


