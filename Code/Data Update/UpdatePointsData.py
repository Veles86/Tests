from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env, ListFields
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor

import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Helpers import start_editing, stop_editing, get_feature_layer_workspace

def update_feature_class_data(main_fc, update_fc, key_field):
    """
    Updates both attributes and geometry from update_fc to main_fc where the key_field matches.
    Excludes GlobalID and OBJECTID from the update.
    """

    # 1. Identify all valid (editable) fields in main_fc and update_fc
    main_fields_all = ListFields(main_fc)
    update_fields_all = ListFields(update_fc)

    # Convert field objects to names, and filter out read-only or system fields
    excluded_fields = {"GlobalID", "OBJECTID"}

    # We'll create lists of field names for each FC that are:
    #   - Not in the excluded list
    #   - Editable by ArcPy
    main_fields = [
        f.name for f in main_fields_all
        if f.editable and f.name not in excluded_fields
    ]
    update_fields = [
        f.name for f in update_fields_all
        if f.editable and f.name not in excluded_fields
    ]

    # 2. Find the common fields (excluding the key_field from the intersection, if desired)
    common_fields = set(main_fields).intersection(set(update_fields))
    # We will not exclude the key field from the dictionary creation entirely (we need to look it up),
    # but we don't treat it as an attribute to "update"
    # so we won't remove it *yet* from the intersection—just handle it carefully later if needed.

    # 3. Build a dictionary from the update_fc
    #    We'll include the key field, geometry (SHAPE@), and the common fields
    fields_for_search = [key_field, "SHAPE@"] + list(common_fields)
    counter = 0
    update_dict = {}
    with SearchCursor(update_fc, fields_for_search) as s_cursor:
        for row in s_cursor:
            code_val = row[0]
            shape_val = row[1]
            attr_vals = row[2:]  # matches the common fields in the list

            attr_dict = dict(zip(common_fields, attr_vals))
            attr_dict["SHAPE@"] = shape_val
            update_dict[code_val] = attr_dict
            counter += 1
	
    AddMessage(f"Read {counter} records from the update feature class.")

    # 4. Use an update cursor on the main_fc to update geometry and attributes
    #    The order is: key_field, SHAPE@, then the same common fields
    fields_for_update = [key_field, "SHAPE@"] + list(common_fields)
    counter = 0

    #ws = get_feature_layer_workspace(main_fc)
    #editor = start_editing(ws)
    with UpdateCursor(main_fc, fields_for_update) as u_cursor:
        for row in u_cursor:
            code_val = row[0]  # key_field

            if code_val in update_dict:
                # Update geometry
                new_shape = update_dict[code_val]["SHAPE@"]
                row[1] = new_shape

                # Update the common attribute fields
                for i, field_name in enumerate(common_fields, start=2):
                    row[i] = update_dict[code_val][field_name]

                u_cursor.updateRow(row)
                counter += 1

    #stop_editing(editor)
    AddMessage(f"Updated {counter} records in the main feature class.")


if __name__ == "__main__":
    # --- 4. Parameters when running as a script or script tool -------------
    # For a script tool, you might set these as parameters in ArcGIS Pro:
    main_fc = GetParameterAsText(0)
    update_fc = GetParameterAsText(1)
    key_field = GetParameterAsText(2)

    update_feature_class_data(main_fc, update_fc, key_field)


