import arcpy
from arcpy import GetParameterAsText, env, Delete_management as Delete
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField, FeatureToPolygon
from FindLoops import find_loops

env.overwriteOutput = True


def find_max_loop(lines_input, points_input, max_loop_output) -> None:


    polygons_from_lines = FeatureToPolygon(in_features=lines_input,out_feature_class=r"memory\polygons_from_lines",cluster_tolerance=None,attributes="NO_ATTRIBUTES",label_features=None)

    dissolved_loop = Dissolve(in_features=polygons_from_lines, out_feature_class=r"memory\dissolved_loop",
                    dissolve_field=None,statistics_fields=None, multi_part="SINGLE_PART", unsplit_lines="DISSOLVE_LINES", concatenation_separator="" )

    selected_lines = SelectByLocation(in_layer=lines_input,overlap_type="SHARE_A_LINE_SEGMENT_WITH",select_features=dissolved_loop,search_distance=None,
                                      selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")

    find_loops(selected_lines, points_input, 'independent', max_loop_output)

    Delete(dissolved_loop)
    Delete(polygons_from_lines)
    SelectByAttribute(lines_input, "CLEAR_SELECTION")




if __name__ == '__main__':

    # Get input parameters from user
    lines_input = GetParameterAsText(0) 
    points_input = GetParameterAsText(1) 
    max_loop_output = GetParameterAsText(2)  

    find_max_loop(lines_input, points_input, max_loop_output)
    #gdb_path = fr"E:\git\LevellingThesis\Data\LoopsTests.gdb"
    #loops_input = gdb_path + fr"\LoopsDar1"
    #lines_input = gdb_path + fr"\Teum1998SegmentsFixed_generalized_Dar1_Only"
    #points_output = gdb_path+fr"\Teum1998SegmentsFixed_coordinates" 
    #merged_loop_output = gdb_path+fr"\MergedLoopsOutput" 