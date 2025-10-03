import arcpy

def compare_layers(layer_a, layer_b, compare_fields, output_table):
    """
    Compare two layers and check if each feature in layer_a has a matching feature in layer_b.
    
    :param layer_a: Path or layer name for Feature Class A
    :param layer_b: Path or layer name for Feature Class B
    :param compare_fields: List of fields to compare. If "SHAPE" is in the list, geometry is compared.
    :param output_table: Path to an output table to store mismatches (Mandatory)
    :return: None (Results are written to the output table)
    """
    mismatches = []

    # --- Ensure the output table is always created ---
    if not output_table:
        output_table = "in_memory\\Mismatches"
        arcpy.AddMessage("No output table specified. Results will be saved to: in_memory\\Mismatches")

    # --- Prepare the field lists ---
    if "SHAPE" in compare_fields:
        fields_b = ["OBJECTID"] + [f for f in compare_fields if f != "SHAPE"] + ["SHAPE@"]
        fields_a = ["OBJECTID"] + [f for f in compare_fields if f != "SHAPE"] + ["SHAPE@"]
    else:
        fields_b = ["OBJECTID"] + compare_fields
        fields_a = ["OBJECTID"] + compare_fields

    # --- Build dictionary from Layer B ---
    features_b = {}
    with arcpy.da.SearchCursor(layer_b, fields_b) as b_cursor:
        for row_b in b_cursor:
            b_oid = row_b[0]
            if "SHAPE" in compare_fields:
                attr_values = tuple(row_b[1:-1])
                geom_b = row_b[-1]
            else:
                attr_values = tuple(row_b[1:])
                geom_b = None
            features_b[attr_values] = (b_oid, geom_b)

    # --- Compare each feature from Layer A ---
    with arcpy.da.SearchCursor(layer_a, fields_a) as a_cursor:
        for row_a in a_cursor:
            a_oid = row_a[0]
            if "SHAPE" in compare_fields:
                attr_values = tuple(row_a[1:-1])
                geom_a = row_a[-1]
            else:
                attr_values = tuple(row_a[1:])
                geom_a = None

            if attr_values not in features_b:
                mismatches.append(("Missing in Layer B", a_oid, None, attr_values))
            else:
                b_oid, geom_b = features_b[attr_values]
                if "SHAPE" in compare_fields:
                    if not geom_a.equals(geom_b):
                        mismatches.append(("Geometry Mismatch", a_oid, b_oid, attr_values))

    # --- Print Summary ---
    arcpy.AddMessage(f"Total mismatches found: {len(mismatches)}")

    # --- Create the Output Table ---
    arcpy.management.CreateTable("in_memory", "temp_mismatch_table")

    # Add necessary fields
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "Type", "TEXT")
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "OBJECTID_A", "LONG")
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "OBJECTID_B", "LONG")

    # Add user-selected attribute fields (except SHAPE)
    out_fields = [f for f in compare_fields if f != "SHAPE"]
    for field in out_fields:
        arcpy.management.AddField("in_memory\\temp_mismatch_table", field, "TEXT")

    # Insert mismatches
    insert_fields = ["Type", "OBJECTID_A", "OBJECTID_B"] + out_fields
    with arcpy.da.InsertCursor("in_memory\\temp_mismatch_table", insert_fields) as ic:
        for mismatch_type, a_oid, b_oid, attr_values in mismatches:
            row_to_insert = [mismatch_type, a_oid, b_oid if b_oid else None] + [str(v) for v in attr_values]
            ic.insertRow(row_to_insert)

    # Save results to the output table
    arcpy.management.CopyRows("in_memory\\temp_mismatch_table", output_table)
    arcpy.management.Delete("in_memory\\temp_mismatch_table")  # Clean up

    arcpy.AddMessage(f"Results saved to: {output_table}")


# --- Script Tool Execution ---
if __name__ == "__main__":
    # Get parameters from the tool
    layer_a = arcpy.GetParameterAsText(0)
    layer_b = arcpy.GetParameterAsText(1)
    compare_fields = arcpy.GetParameterAsText(2).split(";")  # Multi-value selection
    output_table = arcpy.GetParameterAsText(3)  # This is now required

    compare_layers(layer_a, layer_b, compare_fields, output_table)
import arcpy

def compare_layers(layer_a, layer_b, compare_fields, output_table):
    """
    Compare two layers and check if each feature in layer_a has a matching feature in layer_b.
    
    :param layer_a: Path or layer name for Feature Class A
    :param layer_b: Path or layer name for Feature Class B
    :param compare_fields: List of fields to compare. If "SHAPE" is in the list, geometry is compared.
    :param output_table: Path to an output table to store mismatches (Mandatory)
    :return: None (Results are written to the output table)
    """
    mismatches = []

    # --- Ensure the output table is always created ---
    if not output_table:
        output_table = "in_memory\\Mismatches"
        arcpy.AddMessage("No output table specified. Results will be saved to: in_memory\\Mismatches")

    # --- Prepare the field lists ---
    if "SHAPE" in compare_fields:
        fields_b = ["OBJECTID"] + [f for f in compare_fields if f != "SHAPE"] + ["SHAPE@"]
        fields_a = ["OBJECTID"] + [f for f in compare_fields if f != "SHAPE"] + ["SHAPE@"]
    else:
        fields_b = ["OBJECTID"] + compare_fields
        fields_a = ["OBJECTID"] + compare_fields

    # --- Build dictionary from Layer B ---
    features_b = {}
    with arcpy.da.SearchCursor(layer_b, fields_b) as b_cursor:
        for row_b in b_cursor:
            b_oid = row_b[0]
            if "SHAPE" in compare_fields:
                attr_values = tuple(row_b[1:-1])
                geom_b = row_b[-1]
            else:
                attr_values = tuple(row_b[1:])
                geom_b = None
            features_b[attr_values] = (b_oid, geom_b)

    # --- Compare each feature from Layer A ---
    with arcpy.da.SearchCursor(layer_a, fields_a) as a_cursor:
        for row_a in a_cursor:
            a_oid = row_a[0]
            if "SHAPE" in compare_fields:
                attr_values = tuple(row_a[1:-1])
                geom_a = row_a[-1]
            else:
                attr_values = tuple(row_a[1:])
                geom_a = None

            if attr_values not in features_b:
                mismatches.append(("Missing in Layer B", a_oid, None, attr_values))
            else:
                b_oid, geom_b = features_b[attr_values]
                if "SHAPE" in compare_fields:
                    if not geom_a.equals(geom_b):
                        mismatches.append(("Geometry Mismatch", a_oid, b_oid, attr_values))

    # --- Print Summary ---
    arcpy.AddMessage(f"Total mismatches found: {len(mismatches)}")

    # --- Create the Output Table ---
    arcpy.management.CreateTable("in_memory", "temp_mismatch_table")

    # Add necessary fields
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "Type", "TEXT")
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "OBJECTID_A", "LONG")
    arcpy.management.AddField("in_memory\\temp_mismatch_table", "OBJECTID_B", "LONG")

    # Add user-selected attribute fields (except SHAPE)
    out_fields = [f for f in compare_fields if f != "SHAPE"]
    for field in out_fields:
        arcpy.management.AddField("in_memory\\temp_mismatch_table", field, "TEXT")

    # Insert mismatches
    insert_fields = ["Type", "OBJECTID_A", "OBJECTID_B"] + out_fields
    with arcpy.da.InsertCursor("in_memory\\temp_mismatch_table", insert_fields) as ic:
        for mismatch_type, a_oid, b_oid, attr_values in mismatches:
            row_to_insert = [mismatch_type, a_oid, b_oid if b_oid else None] + [str(v) for v in attr_values]
            ic.insertRow(row_to_insert)

    # Save results to the output table
    arcpy.management.CopyRows("in_memory\\temp_mismatch_table", output_table)
    arcpy.management.Delete("in_memory\\temp_mismatch_table")  # Clean up

    arcpy.AddMessage(f"Results saved to: {output_table}")


# --- Script Tool Execution ---
if __name__ == "__main__":
    # Get parameters from the tool
    layer_a = arcpy.GetParameterAsText(0)
    layer_b = arcpy.GetParameterAsText(1)
    compare_fields = arcpy.GetParameterAsText(2).split(";")  # Multi-value selection
    output_table = arcpy.GetParameterAsText(3)  # This is now required

    compare_layers(layer_a, layer_b, compare_fields, output_table)
