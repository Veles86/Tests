import os
import arcpy
import networkx as nx
import numpy as np
import math
import random
import string
from typing import Literal
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, ListFields
from arcpy.da import SearchCursor, Editor, UpdateCursor, InsertCursor
from arcpy import MakeFeatureLayer_management as MakeFeatureLayer
from arcpy.edit import FlipLine
from arcpy.analysis import Statistics
from arcpy.management import Append, CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount,Dissolve, Sort, DeleteField,FindIdentical
from arcpy import CopyFeatures_management as CopyFeatures
from arcpy import Delete_management as Delete, Describe
from Utils.Helpers import get_feature_layer_workspace,start_editing, stop_editing, count_features
from Utils.Points import get_distance_between_points
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace



def update_is_duplicated_field(levelling_lines) -> None:

    SelectByAttribute(in_layer_or_view=levelling_lines,selection_type="CLEAR_SELECTION")

    fields = [f.name for f in ListFields(levelling_lines)]
    if "LineCode" in fields:
        field_with_the_name = "LineCode"
    elif "SegmentCode" in fields:
        field_with_the_name = "SegmentCode"
    else:
        AddError("No field found to indicate the line or segment code.")
        return

    duplicated_lines = FindIdentical(in_dataset=levelling_lines,out_dataset=r"memory\duplicated_lines",fields=field_with_the_name,
                                    xy_tolerance=None,z_tolerance=0,output_record_option="ONLY_DUPLICATES")
    
    duplicated_oids = []
    with SearchCursor(duplicated_lines, ["IN_FID"]) as cursor:
        for row in cursor:
            duplicated_oids.append(row[0])

    with UpdateCursor(levelling_lines, ["OID@", "IsDuplicated"]) as cursor:
        for row in cursor:
            if row[0] in duplicated_oids:
                row[1] = True
            else:
                row[1] = False
            cursor.updateRow(row)

    Delete(duplicated_lines)








def enumerate_lines_by_code(generalized_lines) -> None:

# TODO: insert it into the function that creates the generalized lines feature class instead of explicit code
# TODO: make it also work on segments
    """
    Numerate lines by code in the generalized lines feature class.
    
    Parameters:
    - generalized_lines (FeatureLayer): Input feature class of generalized lines
    
    Returns:
    - None: The function updates the input feature class in place.
    """
    # Create a dictionary to store line codes and their corresponding unique numbers
    line_code_to_unique_num = {}
    unique_num = 1

    with UpdateCursor(generalized_lines, ["LineCode", "UniqueLineNum"]) as cursor:
        for row in cursor:
            line_code = row[0]
            if line_code not in line_code_to_unique_num:
                line_code_to_unique_num[line_code] = unique_num
                unique_num += 1
            row[1] = line_code_to_unique_num[line_code]
            cursor.updateRow(row)


def reverse_line_segments_without_copy(line_segments_fc) -> FeatureLayer:
    """"
    Reverse the direction of the given line segments

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of line segments to be reversed

    Returns
    -------
    line_segments_fc : FeatureLayer
        Feature layer of reversed line segments
    """
        
    #reversed_line_segment = CopyFeatures(line_segments_fc, fr"memory\reversed_line_segment")
        
    with UpdateCursor(line_segments_fc, ["StartPointCode", "EndPointCode", "StartPointName","EndPointName", "SegmentName", "MeasHeightDiff","SegmentDirection","GravCorrection","CorrHeightDiff"]) as cursor:
        for row in cursor:
            row[4] = row[3] + " to " + row[2]
            row[0], row[1] = row[1], row[0]
            row[2], row[3] = row[3], row[2]
            row[5], row[6] = -row[5], -row[6]
            if row[7] != None:
                row[7] = -row[7]
            if row[8] != None:
                row[8] = -row[8]
            cursor.updateRow(row)

    FlipLine(line_segments_fc)

    return line_segments_fc


def reverse_line_segments_from_tuple(row, field_names):
   # improved version of the function suitable for matlab version of sort_segments
    """
    Reverse the values of a SearchCursor row (tuple) according to geodetic logic.

    Parameters
    ----------
    row : tuple
        Full row including SHAPE@ and fields (excluding OID@).
    field_names : list
        The names of the fields in the row (aligned with the row's indexes).

    Returns
    -------
    tuple
        A reversed version of the row (with geometry and key attributes flipped).
    """
    row = list(row)
    geom = row[0]
    # Reverse geometry using point order
    reversed_geom = arcpy.Polyline(arcpy.Array(list(geom.getPart(0))[::-1]), geom.spatialReference)
    row[0] = reversed_geom

    # Define fields to reverse
    def swap_fields(a, b):
        idx_a = field_names.index(a)
        idx_b = field_names.index(b)
        row[idx_a], row[idx_b] = row[idx_b], row[idx_a]

    def negate_field(f):
        idx = field_names.index(f)
        if row[idx] is not None:
            row[idx] = -row[idx]

    # Swap point codes and names
    swap_fields("StartPointCode", "EndPointCode")
    swap_fields("StartPointName", "EndPointName")

    # Update SegmentName
    try:
        idx_start = field_names.index("StartPointName")
        idx_end = field_names.index("EndPointName")
        idx_segname = field_names.index("SegmentName")
        row[idx_segname] = f"{row[idx_start]} to {row[idx_end]}"
    except ValueError:
        pass  # Skip if missing fields

    # Negate fields
    for f in ["MeasHeightDiff", "SegmentDirection", "GravCorrection", "CorrHeightDiff"]:
        if f in field_names:
            negate_field(f)

    return tuple(row)



def add_segment_from_tuple(line_segments_fc, row_data, line_number=0):
    # improved version of the function suitable for matlab version of sort_segments
    """
    Inserts a single reversed or original row into output FC using InsertCursor.

    Parameters
    ----------
    line_segments_fc : FeatureLayer
    row_data : tuple (SHAPE@, ... fields)
    line_number : int
    """
    # fields = [f.name for f in arcpy.ListFields(line_segments_fc) if f.type not in ("OID", "Geometry","GlobalID")]
    # if "line_num" not in fields:
    #     AddField(line_segments_fc, "line_num", "SHORT")
    #     fields.append("line_num")

    # insert_fields = ["SHAPE@"] + fields

    # # Build insert row
    # row_data = list(row_data)
    # idx_line = fields.index("line_num")
    # if len(row_data) <= idx_line:
    #     row_data.append(line_number)
    # else:
    #     row_data[idx_line] = line_number

    # with arcpy.da.InsertCursor(line_segments_fc, insert_fields) as cursor:
    #     insert_row = [row_data[0]]  # SHAPE@
    #     insert_row += row_data[1:len(fields)+1]  # match number of fields exactly
    #     cursor.insertRow(insert_row)


    fields = [
        f.name for f in arcpy.ListFields(line_segments_fc)
        if f.type not in ("OID", "Geometry", "GlobalID")
        and f.name.lower() not in ("shape", "shape_length")
    ]


    if "line_num" not in fields:
        AddField(line_segments_fc, "line_num", "SHORT")
        fields.append("line_num")

    insert_fields = ["SHAPE@"] + fields

    row_data = list(row_data)
    idx_line = fields.index("line_num")
    if len(row_data) <= idx_line:
        row_data.append(line_number)
    else:
        row_data[idx_line] = line_number

    insert_row = [row_data[0]] + row_data[1:len(fields)+1]

    with arcpy.da.InsertCursor(line_segments_fc, insert_fields) as cursor:
        cursor.insertRow(insert_row)


def fill_dH_and_dist_differences(levelling_lines) -> None:

    """
    Calculates the columns dH_from_DB, dH_from_file, dH_diff, dist_from_DB, dist_from_file, dist_diff in the levelling_lines feature class
    according the coordinates data and the measurement data in the levelling_lines
    
    """

    with UpdateCursor(levelling_lines, ["dH","dist", "Y1", "X1", "H1", "Y2", "X2", "H2","dH_from_DB", "dH_from_file", "dH_diff", "dist_from_DB", "dist_from_file", "dist_diff"]) as cursor:
        for row in cursor:
            dH = row[0]
            dist = row[1]
            Y1 = row[2]
            X1 = row[3]
            H1 = row[4]
            Y2 = row[5]
            X2 = row[6]
            H2 = row[7]

            if (H1 and H2):
                dH_from_DB = round(H2 - H1, 3)
                dH_from_file = round(dH, 5)
                dH_diff = round(abs(dH_from_file - dH_from_DB),5)
            else:
                dH_from_DB = None
                dH_from_file = round(dH, 5)
                dH_diff = None
                
            dist_from_DB = round(get_distance_between_points(Y1, X1, Y2, X2), 0)
            dist_from_file = round(dist, 0)
            dist_diff = round(abs(dist_from_DB - dist_from_file),0)

            row[8] = dH_from_DB
            row[9] = dH_from_file
            row[10] = dH_diff
            row[11] = dist_from_DB
            row[12] = dist_from_file
            row[13] = dist_diff

            cursor.updateRow(row)


def turn_segments_by_min_point_name_old_format(line_segments_fc) -> None:
    """
    Reverse the direction of the given line segments

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of line segments to be reversed

    Returns
    -------
    None
    """

    #ws = get_feature_layer_workspace(line_segments_fc)
    #editor = start_editing(ws)
    with UpdateCursor(line_segments_fc, ["A", "B", "dH","codeA","codeB","Shape@","OriginalDirection","OBJECTID"]) as cursor:
        for row in cursor:
            nameA = row[0]
            nameB = row[1]

            if nameA > nameB:            
                row[0], row[1] = row[1], row[0]
                row[3], row[4] = row[4], row[3]
                row[6] = -row[6]
                try:
                    row[5] = row[5].reverse()
                except:
                    pass
                try:
                    row[2] = -row[2]
                except:
                    AddWarning(f"Could not reverse the segment with OBJECTID: {row[7]}")

            cursor.updateRow(row)

    #stop_editing(editor)
'''
def reverse_polyline_geometry(polyline):
    reversed_parts = []
    for part in polyline:
        reversed_part = arcpy.Array([pt for pt in reversed(part) if pt])
        reversed_parts.append(reversed_part)
    return arcpy.Polyline(arcpy.Array(reversed_parts), polyline.spatialReference)
'''

def reverse_polyline_geometry(polyline):
    reversed_parts = []
    # Reverse the order of parts AND reverse each part's vertices
    for part in reversed(list(polyline)):
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

def reverse_generalized_lines(lines_fc) -> None:
    """
    Reverse the direction of the given lines

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of segments to be reversed

    Returns
    -------
    None
    """


    #ws = get_feature_layer_workspace(lines_fc)
    #editor = start_editing(ws)
    with UpdateCursor(lines_fc, ["StartPointCode", "EndPointCode", "StartPointName","EndPointName", "LineName",  "MeasHeightDiff", "LineDirection","GravCorrection","CorrHeightDiff","SHAPE@"]) as cursor:
        for row in cursor:
            row[4] = row[3] + " to " + row[2]
            row[0], row[1] = row[1], row[0]
            row[2], row[3] = row[3], row[2]
            row[5], row[6] = -row[5], -row[6]
            if row[7] != None:
                row[7] = -row[7]
            if row[8] != None:
                row[8] = -row[8]
            
            row[9] = reverse_polyline_geometry(row[9]) 


            cursor.updateRow(row)

  
    #stop_editing(editor)

def reverse_segments_new(line_segments_fc) -> None:
    """
    Reverse the direction of the given segments

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of segments to be reversed

    Returns
    -------
    None
    """


    #ws = get_feature_layer_workspace(line_segments_fc)
    #editor = start_editing(ws)
    with UpdateCursor(line_segments_fc, ["StartPointCode", "EndPointCode", "StartPointName","EndPointName", "SegmentName",  "MeasHeightDiff", "SegmentDirection","GravCorrection","CorrHeightDiff"]) as cursor:
        for row in cursor:
            row[4] = row[3] + " to " + row[2]
            row[0], row[1] = row[1], row[0]
            row[2], row[3] = row[3], row[2]
            row[5], row[6] = -row[5], -row[6]
            if row[7] != None:
                row[7] = -row[7]
            if row[8] != None:
                row[8] = -row[8]
            cursor.updateRow(row)

    FlipLine(line_segments_fc)
    #stop_editing(editor)


def reverse_segments(line_segments_fc, line_type: Literal["segments", "lines"] = "segments") -> None:
    """
    Reverse the direction of the given line segments

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of line segments to be reversed

    Returns
    -------
    None
    """

    if line_type == "segments":
        name_type = "SegmentName"
        direction_type = "SegmentDirection"
    elif line_type == "lines":
        name_type = "LineName"
        direction_type = "LineDirection"
    else:
        raise ValueError("Invalid line type. Use 'segments' or 'lines'.")
    ws = get_feature_layer_workspace(line_segments_fc)
    #AddMessage(f"Reversing {line_type} in {ws}")
    editor = start_editing(ws)
    #AddMessage(f"Entered edit in in {ws}")
    with UpdateCursor(line_segments_fc, ["StartPointCode", "EndPointCode", "StartPointName","EndPointName", name_type,  "MeasHeightDiff", direction_type,"GravCorrection","CorrHeightDiff"]) as cursor:
        #AddMessage("Edit session started.")
        for row in cursor:
            #AddMessage("entered for loop")
            row[4] = row[3] + " to " + row[2]
            row[0], row[1] = row[1], row[0]
            row[2], row[3] = row[3], row[2]
            row[5], row[6] = -row[5], -row[6]
            if row[7] != None:
                row[7] = -row[7]
            if row[8] != None:
                row[8] = -row[8]
            #AddMessage("before updateRow")
            cursor.updateRow(row)

    FlipLine(line_segments_fc)
    stop_editing(editor)
    #AddMessage(f"Finished edit in {ws}")


def reverse_line_segments(line_segments_fc) -> FeatureLayer:
    """"
    Reverse the direction of the given line segments

    Parameters  
    ----------
    line_segments_fc : FeatureLayer
        Feature layer of line segments to be reversed

    Returns
    -------
    line_segments_fc : FeatureLayer
        Feature layer of reversed line segments
    """
        
    reversed_line_segment = CopyFeatures(line_segments_fc, fr"memory\reversed_line_segment")
        
    with UpdateCursor(reversed_line_segment, ["StartPointCode", "EndPointCode", "StartPointName","EndPointName", "SegmentName", "MeasHeightDiff","SegmentDirection","GravCorrection","CorrHeightDiff"]) as cursor:
        for row in cursor:
            row[4] = row[3] + " to " + row[2]
            row[0], row[1] = row[1], row[0]
            row[2], row[3] = row[3], row[2]
            row[5], row[6] = -row[5], -row[6]
            if row[7] != None:
                row[7] = -row[7]
            if row[8] != None:
                row[8] = -row[8]
            cursor.updateRow(row)

    FlipLine(reversed_line_segment)

    return reversed_line_segment

def get_first_segment_with_point(line_segments, point_code) -> FeatureLayer:
    """
    Get the first segment with the given point code

    Parameters
    ----------
    line_segments : FeatureLayer
        Feature layer of line segments
    point_code : int
        Code of the point    

    Returns
    -------
    first_row : FeatureLayer
        Feature layer of the first segment with the given point code
    """
    segments_with_point = SelectByAttribute(line_segments, where_clause=f"StartPointCode = {point_code} OR EndPointCode = {point_code} ", selection_type='NEW_SELECTION')


    # Use a SearchCursor to get the first row from the selected features
    with SearchCursor(segments_with_point, ["UniqueRowNum"]) as cursor:
        first_row = next(cursor, None)  # Get the first selected row

    # If a row was found, select only that row by OBJECTID (or another unique field)
    if first_row:
        query = f"UniqueRowNum = {first_row[0]}"
        # Clear the previous selection and select only the first row
        selected_first_row = SelectByAttribute(line_segments,where_clause = query, selection_type= "NEW_SELECTION")
        first_row = CopyFeatures(selected_first_row, fr"memory\selected_first_row")
        SelectByAttribute(in_layer_or_view=line_segments,selection_type="CLEAR_SELECTION")
        return first_row
    else:
        return None
    
def get_segments_end_point_names(segment) -> list:
    """
    Get the end point names of the given segment

    Parameters
    ----------
    segment : FeatureLayer
        Feature layer of line segment

    Returns
    -------
    list
        List of end point names
    """
    point_names = []

    with SearchCursor(segment, ["StartPointName", "EndPointName"]) as cursor:
        for row in cursor:
            if row[0] != None:
                point_names.append(row[0])
            if row[1] != None:
                point_names.append(row[1])

    return point_names

def get_segments_end_point_codes(segment) -> list:
    """
    Get the end point codes of the given segment

    Parameters
    ----------
    segment : FeatureLayer
        Feature layer of line segment

    Returns
    -------
    list
        List of end point codes
    """
    point_codes = []

    with SearchCursor(segment, ["StartPointCode", "EndPointCode"]) as cursor:
        for row in cursor:
            if row[0] != None:
                point_codes.append(row[0])
            if row[1] != None:
                point_codes.append(row[1])

    return point_codes





def get_last_segment(line_segments) -> FeatureLayer:
    """
    Get the last segment in the given feature class

    Parameters
    ----------
    line_segments : FeatureLayer
        Feature layer of line segments    

    Returns
    -------
    segments_with_point : FeatureLayer
        Feature layer of the last line segment
    """

    # Get the spatial reference from the input segments feature class
    sr = Describe(line_segments).spatialReference

    last_row = CreateFeatureclass(out_path="memory", out_name="last_row",geometry_type="POLYLINE",template=line_segments,spatial_reference=sr)

    # List fields from both the input and the created feature class
    #input_fields = [field.name for field in ListFields(line_segments)]
    #output_fields = [field.name for field in ListFields(last_row)]
    #AddMessage(f"input_fields count: {len(input_fields)}")
    #AddMessage(f"input_fields: {input_fields}")
    #AddMessage(f"output_fields count: {len(output_fields)}")
    #AddMessage(f"input_fields: {output_fields}")
    # Ensure only matching fields are used in the cursors
    #common_fields = [field for field in input_fields if field in output_fields]

    input_fields = [field for field in ListFields(line_segments) if field.type not in ["OID", "GlobalID"]]
    output_fields = [field for field in ListFields(last_row) if field.type not in ["OID", "GlobalID"]]

    common_fields = [field.name for field in input_fields if field.name in [f.name for f in output_fields]]



    num_of_rows = int(GetCount(line_segments)[0])

    row_num = 0

    # Use SearchCursor to iterate through all rows
    with SearchCursor(line_segments,common_fields) as search_cursor:
        for row in search_cursor:
            row_num += 1
            if row_num == num_of_rows:
                with InsertCursor(last_row,common_fields) as insert_cursor:
                    insert_cursor.insertRow(row)
                    

    return last_row

  

def get_last_segment_end_point_code(line_segments) -> int:

    """
    Get the last segment code frpm the sorted feature class

    Parameters
    ----------
    line_segments : FeatureLayer
        Feature class of line segments    

    Returns
    -------
    Last segment end point code : int
        
    """

    if int(GetCount(line_segments)[0]) == 1:
        last_segment = line_segments
    else:
        last_segment = get_last_segment(line_segments)

    if last_segment:
        # Use SearchCursor to get the "EndPointCode" field for the row
        with SearchCursor(last_segment, ["EndPointCode"]) as cursor:
            # Use 'next()' to retrieve the first row or handle if the feature class is empty
            last_row = next(cursor, None)
            
            if last_row:
                return last_row[0]  # Since "EndPointCode" is the first field in the cursor
            else:
                return None  # Return None if the cursor is empty
    else:
        AddError("Last segment not found")

def get_single_segment_fc(input_fc, oid, temp_name="temp_segment"):
    """
    Extracts a 1-row in-memory feature class by OID
    """
    where = f"OBJECTID = {oid}"
    selection = arcpy.management.SelectLayerByAttribute(input_fc, "NEW_SELECTION", where)
    return arcpy.management.CopyFeatures(selection, fr"in_memory\{temp_name}_{oid}")

def add_segment(line_segments, segment_for_addition, line_number=0) -> None:
    """"
    Appends the given segment to the given line segments

    Parameters
    ----------
    line_segments : FeatureLayer
        Feature layer of line segments
    segment_for_addition : Feature
        Feature of segment to be added
    line_number : int (optional)
        Line number of the generalized line
    """
    #arcpy.env.preserveGlobalIds = False

    if line_number == 0:
        Append(inputs=segment_for_addition,target=line_segments,schema_type="NO_TEST",field_mapping=None,subtype="",match_fields=None,update_geometry="NOT_UPDATE_GEOMETRY")
    else:
        AddField(segment_for_addition, "line_num", "SHORT")
        with UpdateCursor(segment_for_addition, ["line_num"]) as cursor:
            for row in cursor:
                row[0] = line_number
                cursor.updateRow(row)
        #arcpy.env.preserveGlobalIds = False
        Append(inputs=segment_for_addition,target=line_segments,schema_type="NO_TEST",field_mapping=None,subtype="",match_fields=None,update_geometry="NOT_UPDATE_GEOMETRY")



def delete_segment(line_segments, segment_for_deletion) -> None:
    """"
    Deletes the given segment from the given line segments

    Parameters
    ----------
    line_segments : FeatureLayer
        Feature layer of line segments
    segment_for_deletion : Feature
        Feature of segment to be deleted
    """

    with SearchCursor(segment_for_deletion, ["UniqueRowNum"]) as cursor1:
        for row1 in cursor1:
            query = fr"UniqueRowNum = {row1[0]}"
            with UpdateCursor(line_segments, ["UniqueRowNum"], query) as cursor2:
                for row2 in cursor2:
                    cursor2.deleteRow()

