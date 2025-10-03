import os
import networkx as nx
import arcpy
from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField
from typing import Literal

from arcpy.analysis import Union
from arcpy.management import DeleteIdentical
from Utils.Helpers import start_editing, stop_editing,get_feature_layer_workspace,count_features

from Utils.Loops import create_graph, extract_data_from_cycles, create_polygons_from_loops, create_multigraph, num_of_independent_loops
from Utils.Points import extract_coordinates_from_external_source
from Utils.Calculations import calculate_missing_heights
from Utils.Configs import SPATIAL_REFERENCE, CONTROL_POINTS_FULL_DATA

arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True



def test_missing_h(levelling_lines,points_source_data) -> None:
    calculate_missing_heights(levelling_lines, points_source_data)


def find_loops(levelling_lines_input, source_points_input, loops_type:Literal['independent','all'], loops_output) -> None:
    # Main script execution

    # get coordinates data fow the participating lines
    points_data = extract_coordinates_from_external_source(levelling_lines_input, source_points_input)
    
    G = create_graph(levelling_lines_input, points_data)

    num_of_ind_loops = num_of_independent_loops(G)
    #TODO add check with selection for overlaping polygons after networkx usage

    polygons_from_cycles = "memory\polygons_from_cycles"
    #sr = SpatialReference(6991)
    polygons_data = extract_data_from_cycles(G)
    create_polygons_from_loops(polygons_data,polygons_from_cycles)
    #CopyFeatures(in_features=polygons_from_cycles,out_feature_class=loops_output+"_temp") #TODO delete later
    polygons_from_union = "memory\polygons_from_union"
    Union(in_features=polygons_from_cycles,out_feature_class=polygons_from_union,join_attributes="ALL",cluster_tolerance=None,gaps="GAPS")
    DeleteIdentical(in_dataset=polygons_from_union,fields="Shape",xy_tolerance=None,z_tolerance=0)

    n_features_in_polygons = count_features(polygons_from_cycles)
    #SelectByLocation(in_layer=polygons_from_union,overlap_type="ARE_IDENTICAL_TO",select_features=polygons_from_cycles,search_distance=None,selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")
    identical_polygons = SelectByLocation(in_layer=polygons_from_union,overlap_type="ARE_IDENTICAL_TO",select_features=polygons_from_cycles,
                            search_distance=None,
                            selection_type="NEW_SELECTION",
                            invert_spatial_relationship="NOT_INVERT"
                        )
    num_of_identical_features = count_features(identical_polygons)
    SelectByAttribute(in_layer_or_view=identical_polygons,selection_type="CLEAR_SELECTION")
    
    if num_of_identical_features != n_features_in_polygons:

        polygons_data_new = []
        undetected_polygons = [] #TODO
        with SearchCursor(polygons_from_union,["SHAPE@"]) as cursor:
            for row in cursor:
                
                temp_polygon = row[0]
                temp_lines = SelectByLocation(in_layer=levelling_lines_input,overlap_type="SHARE_A_LINE_SEGMENT_WITH",select_features=temp_polygon,search_distance=None,selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")
                temp_points = extract_coordinates_from_external_source(temp_lines,source_points_input)
                temp_graph = create_graph(temp_lines, temp_points)
                temp_polygons_data = extract_data_from_cycles(temp_graph)

                if temp_polygons_data:
                    polygons_data_new.append(temp_polygons_data[0])
                else:
                    undetected_polygons.append(temp_polygon)
            
        if undetected_polygons:
            temp_lines = SelectByLocation(in_layer=levelling_lines_input,overlap_type="SHARE_A_LINE_SEGMENT_WITH",select_features=undetected_polygons,search_distance=None,selection_type="NEW_SELECTION",invert_spatial_relationship="NOT_INVERT")
            temp_points = extract_coordinates_from_external_source(temp_lines,source_points_input)
            temp_graph = create_graph(temp_lines, temp_points)  
            temp_polygons_data = extract_data_from_cycles(temp_graph)
            if temp_polygons_data:
                polygons_data_new.extend(temp_polygons_data)
        create_polygons_from_loops(polygons_data_new,polygons_from_cycles)
    

    CopyFeatures(in_features=polygons_from_cycles,out_feature_class=loops_output)

    if num_of_ind_loops == count_features(polygons_from_cycles):
        AddMessage(f"✅{num_of_ind_loops} independent loops were found")
    else:
        AddMessage(f"⚠️{count_features(polygons_from_cycles)}/{num_of_ind_loops} independent loops were found")


    Delete(polygons_from_cycles)
    Delete(polygons_from_union)

    """
    cycles = nx.cycle_basis(G)

    # Print the cycles using the 'full_name'
    AddMessage("Independent cycles in the graph:")

    # Print the cycles
    AddMessage("Independent cycles in the graph:")
    for cycle in cycles:
        full_names = [G.nodes[node]['Name'] for node in cycle]
        cycle_coords = [G.nodes[node]['Shape'] for node in cycle]
        AddMessage(full_names)
        AddMessage(cycle_coords)

        # Print the cycles
    AddMessage("Independent cycles in the graph:")
    for cycle in cycles:
        dH_values = []
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            dH_values.append(G[u][v]['MeasHeightDiff'])
        AddMessage(dH_values)


    # Print the cycles with start node, end node, and attributes
    AddMessage("Independent cycles in the graph:")
    for cycle in cycles:
        height_diff = 0
        dist = 0
        AddMessage(f"len {len(cycle)}")
        AddMessage(f" range len {range(len(cycle))}")
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            edge_data = G[u][v]
            start_point = edge_data.get('StartPointCode', 'N/A')
            end_point = edge_data.get('EndPointCode', 'N/A')
            dH_value = edge_data.get('MeasHeightDiff', 'N/A')
            dist_value = edge_data.get('MeasDist', 'N/A')
            if u == end_point and v == start_point:
                height_diff = height_diff - dH_value
            elif u == start_point and v == end_point:
                height_diff = height_diff + dH_value
            dist = dist + dist_value
            AddMessage(f"Edge from {u} to {v}: start_point={start_point}, end_point={end_point}, dH={dH_value}")
        AddMessage(f"Cycle {cycle}: height_diff={round(height_diff,5)}, dist={round(dist/1000,1)}")

"""

'''
    cycles = nx.cycle_basis(G)

    # Print the cycles
    AddMessage("Independent cycles in the graph:")
    for cycle in cycles:
        AddMessage(cycle)
        '''



"""    
def create_polygons_from_cycles(G, output_fc):
    # Create a new feature class for storing polygons
    spatial_reference = arcpy.SpatialReference(4326)  # Adjust spatial reference as needed
    arcpy.CreateFeatureclass_management(
        out_path=arcpy.Describe(output_fc).path,
        out_name=arcpy.Describe(output_fc).name,
        geometry_type="POLYGON",
        spatial_reference=spatial_reference
    )
    
    # Add fields for storing attributes
    arcpy.AddField_management(output_fc, "HeightDiff", "DOUBLE")
    arcpy.AddField_management(output_fc, "Dist_km", "DOUBLE")
    
    # Insert polygons with attributes into the feature class
    with arcpy.da.InsertCursor(output_fc, ["SHAPE@", "HeightDiff", "Dist_km"]) as cursor:
        # Find independent cycles
        cycles = nx.cycle_basis(G)
        
        for cycle in cycles:
            height_diff = 0
            dist = 0
            coords = []  # Store coordinates for polygon geometry
            
            for i in range(len(cycle)):
                u = cycle[i]
                v = cycle[(i + 1) % len(cycle)]
                edge_data = G[u][v]
                start_point = edge_data.get('start_point', 'N/A')
                end_point = edge_data.get('end_point', 'N/A')
                dH_value = edge_data.get('dH', 0)
                dist_value = edge_data.get('dist', 0)
                u_coords = G.nodes[u]['coords']  # Ensure nodes have 'coords' attribute
                
                # Summarize attributes
                if u == end_point and v == start_point:
                    height_diff -= dH_value
                elif u == start_point and v == end_point:
                    height_diff += dH_value
                dist += dist_value
                
                # Add node coordinates
                coords.append(u_coords)
            
            # Close the polygon by adding the first point again
            coords.append(coords[0])
            
            # Create a polygon geometry
            polygon = arcpy.Polygon(arcpy.Array([arcpy.Point(*coord) for coord in coords]), spatial_reference)
            
            # Insert the polygon and attributes
            cursor.insertRow([polygon, round(height_diff, 5), round(dist / 1000, 1)])

    print(f"Polygons saved to {output_fc}")
"""


if __name__ == '__main__':

    # Get input parameters from user
    levelling_lines_input = GetParameterAsText(0)  # Input feature class of selected segments or lines
    source_points_input = GetParameterAsText(1)  # Input feature class with coordinates data
    loops_type = GetParameterAsText(2)  # Input feature class of loops type: independent or dependent
    loops_output = GetParameterAsText(3)  # Output feature class for calculated loops


    #gdb_path = fr"E:\git\LevellingThesis\Data\Testing4.gdb"
    #levelling_lines_input = gdb_path + fr"\SegmentsAreaB_generalized"
    #source_points_input = CONTROL_POINTS_FULL_DATA
    #loops_type = 'independent'  # Input feature class of loops type: independent or dependent
    #loops_output = gdb_path+fr"\AreaB_generalized_polygons"  # Output feature class for calculated loops

    find_loops(levelling_lines_input, source_points_input, loops_type, loops_output)
    #test_missing_h(levelling_lines_input,source_points_input)


    #gdb_path = fr"E:\git\LevellingThesis\Data\Testings.gdb"
    #levelling_segments_input = gdb_path + fr"\recent_filtered_3"
    #source_points_input = gdb_path + fr"\ControlPoints_full"
    #points_classification_output = gdb_path + fr"\even_yehuda_points_classification"
    #sorted_segments_output = gdb_path + fr"\even_yehuda_sorted_segments"
    #generalized_lines_output = gdb_path + fr"\even_yehuda_generalized"
  







