import os
import itertools
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField, SelectLayerByLocation as SelectByLocation,FindIdentical

from arcpy.edit import FlipLine

from Utils.Generalization import check_self_loops
from Utils.Helpers import count_features
from Utils.Lines import update_is_duplicated_field
from Utils.Configs import LOOPS_TEMPLATE
from FindLoops import find_loops

arcpy.env.overwriteOutput = True

def optimize_loop(selected_loop, levelling_lines,points_data,output_optimized_loop) -> None:


    feature_count = count_features(selected_loop)
    if feature_count != 1:
        AddError("You must select exactly one feature.")
        raise arcpy.ExecuteError

    #update_is_duplicated_field(levelling_lines)


    participating_lines = SelectByLocation(in_layer=levelling_lines,overlap_type="SHARE_A_LINE_SEGMENT_WITH",select_features=selected_loop,
                                            search_distance=None,selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")
    



    with SearchCursor(participating_lines, ["UniqueRowNum","UniqueLineNum","IsDuplicated"]) as cursor:
        non_duplicated_row_nums = []
        duplicated_line_nums = []
        for row in cursor:
            if row[2] == 0:
                non_duplicated_row_nums.append(row[0])
            else:
                duplicated_line_nums.append(row[1])

    if not duplicated_line_nums:
        AddMessage("No duplicated lines found in the selected loop.")
        return
    else:
        AddMessage(f"Found {len(set(duplicated_line_nums))} duplicated lines in the selected loop.")

    SelectByAttribute(in_layer_or_view=levelling_lines, selection_type="CLEAR_SELECTION")
    sr = Describe(participating_lines).spatialReference


    # Separate the path and name from the sorted_segments_output parameter
    out_path, out_name = os.path.split(output_optimized_loop)

    output_loops = CreateFeatureclass(out_path=out_path, out_name=out_name, geometry_type="POLYGON",template=LOOPS_TEMPLATE, spatial_reference=sr)
        
    AddField(in_table=output_loops,field_name="ParticipatingLines",field_type="TEXT",field_length=2500)

    options = []

    for line_number in set(duplicated_line_nums):
        temp_rows_list = []
        with SearchCursor(levelling_lines, ["UniqueRowNum"], where_clause=f"UniqueLineNum = {line_number}") as cursor:
            for row in cursor:
                temp_rows_list.append(row[0])
        options.append(temp_rows_list)
        AddMessage("Added options for line number: " + str(line_number) + " with rows: " + str(temp_rows_list))
    combinations = list(itertools.product(*options))

    num_of_options = 1
    for sublist in options:
        num_of_options *= len(sublist)

    AddMessage(f"Total number of combinations: {num_of_options}")



    

    participating_lines =[]
    for p in combinations:
        AddMessage(f"Processing combination: {p}")
        
        temp_list = non_duplicated_row_nums + list(p)
        temp_list_str = ",".join(str(temp_line) for temp_line in temp_list)
        participating_lines.append(temp_list_str)

        current_loop = SelectByAttribute(in_layer_or_view=levelling_lines,selection_type="NEW_SELECTION",where_clause=f"UniqueRowNum IN ({temp_list_str})",invert_where_clause=None)
        temp_loop = find_loops(current_loop, points_data, 'independent', r"memory\temp_loop")

        field_mapping = fr'MeasHeightDiff "MeasHeightDiff" true true false 8 Double 0 0,First,#,{temp_loop},MeasHeightDiff,-1,-1;'+\
                        fr'MeasDist "MeasDist" true true false 4 Long 0 0,First,#,{temp_loop},MeasDist,-1,-1;'+\
                        fr'MiscBF "MiscBF" true true false 8 Double 0 0,First,#,{temp_loop},MiscBF,-1,-1;;'+\
                        fr'SetUpsNum "SetUpsNum" true true false 4 Long 0 0,First,#,{temp_loop},SetUpsNum,-1,-1;'+\
                        fr'GravCorrection "GravCorrection" true true false 8 Double 0 0,First,#,{temp_loop},GravCorrection,-1,-1;'+\
                        fr'CorrHeightDiff "CorrHeightDiff" true true false 8 Double 0 0,First,#,{temp_loop},CorrHeightDiff,-1,-1;'+\
                        fr'NumOfLines "NumOfLines" true true false 4 Long 0 0,First,#,{temp_loop},NumOfLines,-1,-1;'+\
                        fr'NumOfSegments "NumOfSegments" true true false 4 Long 0 0,First,#,{temp_loop},NumOfSegments,-1,-1'
        
        Append(inputs=r"memory\temp_loop",target=output_loops,schema_type="NO_TEST",field_mapping=field_mapping,
               subtype="",expression="",match_fields=None,update_geometry="NOT_UPDATE_GEOMETRY",feature_service_mode="USE_FEATURE_SERVICE_MODE")

    with UpdateCursor(output_loops, ["ParticipatingLines"]) as cursor:
        for row in cursor:
            row[0] = participating_lines.pop(0)
            cursor.updateRow(row)
        


if __name__ == '__main__':

    # Get input parameters from user
    selected_loop = arcpy.GetParameterAsText(0) 
    levelling_lines = arcpy.GetParameterAsText(1) 
    points_data = arcpy.GetParameterAsText(2)
    output_optimized_loop = arcpy.GetParameterAsText(3) 


    optimize_loop(selected_loop, levelling_lines,points_data,output_optimized_loop)

    