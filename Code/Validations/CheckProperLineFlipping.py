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

from Utils.Generalization import check_lines_connectivity, check_proper_division
#from Utils.Lines import check_proper_line_flipping
arcpy.env.overwriteOutput = True


def check_proper_line_flipping(levelling_lines, points_data) -> None:
    # Create a search cursor to iterate through levelling_lines
    with SearchCursor(levelling_lines, ['SHAPE@', 'StartPointCode', 'EndPointCode', 'OID@']) as line_cursor:
        # Create a search cursor to iterate through points_data
        with SearchCursor(points_data, ['SHAPE@', 'Code']) as point_cursor:
            # Convert points_data into a dictionary for quick lookup
            points_dict = {row[1]: row[0] for row in point_cursor}

        # Iterate through each line in levelling_lines
        for line_row in line_cursor:
            line_geometry = line_row[0]
            start_point_code = line_row[1]
            end_point_code = line_row[2]

            # Determine the first and last points of the line geometry
            if line_geometry.isMultipart:
                first_point = line_geometry.getPart(0)[0]
                last_part = line_geometry.getPart(line_geometry.partCount - 1)
                last_point = last_part[len(last_part) - 1]
            else:
                first_point = line_geometry.firstPoint
                last_point = line_geometry.lastPoint

            # Check if the start point matches the corresponding point in points_data
            if start_point_code in points_dict:
                if not first_point.equals(points_dict[start_point_code]):
                    AddWarning(f"Start point of line {line_row[3]} does not match StartPointCode: {start_point_code}")
            else:
                AddWarning(f"StartPointCode {start_point_code} not found in points_data")

            # Check if the end point matches the corresponding point in points_data
            if end_point_code in points_dict:
                if not last_point.equals(points_dict[end_point_code]):
                    AddWarning(f"End point of line {line_row[3]} does not match EndPointCode: {end_point_code}")
            else:
                AddWarning(f"EndPointCode {end_point_code} not found in points_data")





if __name__ == '__main__':

    # Get input parameters from user
    levelling_lines = arcpy.GetParameterAsText(0)
    points_data = arcpy.GetParameterAsText(1)
    #output_table = arcpy.GetParameterAsText(2) 

    check_proper_line_flipping(levelling_lines, points_data)