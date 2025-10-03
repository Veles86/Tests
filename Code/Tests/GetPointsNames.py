import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine
import os
import sys
sys.path.append(os.path.abspath(".."))

from Utils.Points import get_point_names


arcpy.env.overwriteOutput = True



if __name__ == '__main__':

    s = "563760, 563759, 569064, 409447, 569229"

    # Convert to list of integers
    codes_list = [int(x.strip()) for x in s.split(",")]
    #codes_list = ['563760', '563759', '569064', '409447', '569229', '569228', '569227', '569226', '413155', '569232', '569231', '569236', '545972']
    #codes_list = [563760,563759,569064,409447,569229,569228,569227,569226,413155,569232,569231,569236,545972]


    names_list = get_point_names(codes_list)
    AddMessage(f"SUCCESS: Point names for codes {codes_list}: {names_list}")


