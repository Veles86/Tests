import arcpy
from arcpy import GetParameterAsText, env, Delete_management as Delete
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField
from FindLoops import find_loops

env.overwriteOutput = True

def merge_loops(loops_input, lines_input, points_input, merged_loop_output) -> None:

    dissolved_loop = Dissolve(in_features=loops_input, out_feature_class=r"memory\dissolved_loop",
                        dissolve_field=None,statistics_fields=None, multi_part="SINGLE_PART", unsplit_lines="DISSOLVE_LINES", concatenation_separator="" )

    selected_lines = SelectByLocation(in_layer=lines_input,overlap_type="SHARE_A_LINE_SEGMENT_WITH",select_features=dissolved_loop,search_distance=None,
                                      selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")

    find_loops(selected_lines, points_input, 'independent', merged_loop_output)

    Delete(dissolved_loop)
    SelectByAttribute(lines_input, "CLEAR_SELECTION")

if __name__ == '__main__':

    # Get input parameters from user
    loops_input = GetParameterAsText(0)  
    lines_input = GetParameterAsText(1) 
    points_input = GetParameterAsText(2) 
    merged_loop_output = GetParameterAsText(3)  

    merge_loops(loops_input, lines_input, points_input, merged_loop_output)
    #gdb_path = fr"E:\git\LevellingThesis\Data\LoopsTests.gdb"
    #loops_input = gdb_path + fr"\LoopsDar1"
    #lines_input = gdb_path + fr"\Teum1998SegmentsFixed_generalized_Dar1_Only"
    #points_output = gdb_path+fr"\Teum1998SegmentsFixed_coordinates" 
    #merged_loop_output = gdb_path+fr"\MergedLoopsOutput" 
    