import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, ListFields
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine
import os
import sys
sys.path.append(os.path.abspath(".."))

from Utils.Generalization import check_self_loops
from Utils.Points import get_point_height

arcpy.env.overwriteOutput = True

def calculate_approved_height_diff(levelling_lines, points_data):

    # Check if the field 'ApprovedHeightDiff' exists, if not, add it
    fields = [field.name for field in ListFields(levelling_lines)]
    if "ApprovedHeightDiff" not in fields:
        arcpy.AddField_management(levelling_lines, "ApprovedHeightDiff", "DOUBLE")
    
    with arcpy.da.UpdateCursor(levelling_lines, ["StartPointCode", "EndPointCode", "ApprovedHeightDiff"]) as cursor:
        for row in cursor:
            height_A = get_point_height(row[0],points_data)
            height_B = get_point_height(row[1],points_data)

            if height_A and height_B:
                row[2] = height_B - height_A
                cursor.updateRow(row)

            

if __name__ == '__main__':

    # Get input parameters from user
    levelling_lines = arcpy.GetParameterAsText(0) 
    points_data = arcpy.GetParameterAsText(1) 

    calculate_approved_height_diff(levelling_lines, points_data)



