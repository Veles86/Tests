from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields, Array, Polyline
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor

import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace

import arcpy

def reconstruct_lines_geometry(
    lines_fc, 
    points_fc, 
    code_field_points="Code", 
    codeA_field="codeA", 
    codeB_field="codeB", 
    mid_field="mid_codes"
):
    """
    Rebuilds each line in lines_fc so that its geometry goes from the point with codeA,
    through any intermediate points in mid_codes, and ends at the point with codeB.
    
    Optimized to only load point geometries that are actually needed by lines_fc.

    :param lines_fc: Feature class (lines). Must have codeA_field, codeB_field, mid_field.
    :param points_fc: Feature class (points). Must have code_field_points.
    :param code_field_points: The name of the code field in the points_fc (default 'Code').
    :param codeA_field: The name of the codeA field in the lines_fc.
    :param codeB_field: The name of the codeB field in the lines_fc.
    :param mid_field: The name of the field in the lines_fc that holds intermediate codes 
                      (comma-separated or some other format).
    """

    # ------------------------------------------------------
    # 1. Build a set of all needed codes by scanning the lines FC
    # ------------------------------------------------------
    needed_codes = set()

    with SearchCursor(lines_fc, [codeA_field, codeB_field, mid_field]) as line_cursor:
        for row in line_cursor:
            codeA_val = row[0]
            codeB_val = row[1]
            mid_codes_val = row[2]

            # Add codeA and codeB if they exist
            if codeA_val:
                needed_codes.add(codeA_val)
            if codeB_val:
                needed_codes.add(codeB_val)

            # Parse mid_codes if it exists
            if mid_codes_val:
                # Example: comma-separated list "C100,C101"
                mid_split = [c.strip() for c in mid_codes_val.split(",") if c.strip()]
                needed_codes.update(mid_split)

    # ------------------------------------------------------
    # 2. Read point geometries for only those needed codes
    # ------------------------------------------------------
    point_dict = {}

    # We'll read all points, then filter by needed_codes in Python.
    # If your needed_codes set is extremely large, an "IN" clause might also get huge.
    # For smaller sets, you could do something like:
    #
    #   where_clause = f"{code_field_points} IN ('" + "','".join(needed_codes) + "')"
    #
    # But for extremely large sets, it's often easier to just filter in Python.

    with SearchCursor(points_fc, [code_field_points, "SHAPE@"]) as point_cursor:
        for code_val, pt_geom in point_cursor:
            if code_val in needed_codes:
                # Store the first occurrence, or handle duplicates as needed
                if code_val not in point_dict:
                    point_dict[code_val] = pt_geom

    # ------------------------------------------------------
    # 3. Update each line geometry
    # ------------------------------------------------------
    update_fields = [codeA_field, codeB_field, mid_field, "SHAPE@"]
    with UpdateCursor(lines_fc, update_fields) as u_cursor:
        for row in u_cursor:
            codeA_val = row[0]
            codeB_val = row[1]
            mid_codes_val = row[2]  # e.g. "C100,C101"
            #old_shape = row[3]      # not used in this example, we overwrite

            # Parse mid_codes (if any)
            mid_list = []
            if mid_codes_val:
                mid_list = [c.strip() for c in mid_codes_val.split(",") if c.strip()]

            # Build final code sequence: [codeA, mid_code_1, ..., codeB]
            code_sequence = []
            if codeA_val:
                code_sequence.append(codeA_val)
            code_sequence.extend(mid_list)
            if codeB_val:
                code_sequence.append(codeB_val)

            # Gather point geometries in that order
            point_geometries = []
            for c in code_sequence:
                if c in point_dict:
                    point_geometries.append(point_dict[c])
                else:
                    # If the code is missing in point_dict, skip or handle differently
                    pass

            # Build a new polyline if we have at least two points
            new_polyline = None
            if len(point_geometries) >= 2:
                # We'll assume the points share the same spatial reference as lines_fc
                # or that you can project them if needed
                new_polyline = Polyline(
                    Array(point_geometries),
                    point_geometries[0].spatialReference
                )

            # Overwrite the line shape with the newly constructed geometry
            if new_polyline and not new_polyline.isEmpty:
                row[3] = new_polyline
                u_cursor.updateRow(row)
            else:
                # Optionally keep the original geometry if we didn't build a valid new one
                # row[3] = old_shape
                # u_cursor.updateRow(row)
                pass

    print("Line geometry reconstruction is complete!")





if __name__ == "__main__":

    lines_fc_path = GetParameterAsText(0)
    points_fc_path = GetParameterAsText(1)

    reconstruct_lines_geometry(
        lines_fc=lines_fc_path, 
        points_fc=points_fc_path,
        code_field_points="Code",    # name of code field in point FC
        codeA_field="codeA",         # name of codeA field in line FC
        codeB_field="codeB",         # name of codeB field in line FC
        mid_field="mid_codes"        # name of mid_codes field in line FC
    )
