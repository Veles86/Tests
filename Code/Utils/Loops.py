import os
import networkx as nx
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, ListFields, SpatialReference, env,  Exists, Array, Point, Polygon
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from typing import Literal
from arcpy.management import CreateFeatureclass, AddField
from Utils.Configs import SPATIAL_REFERENCE, LOOPS_TEMPLATE



def num_of_independent_loops(G: nx.Graph) -> int:
    """
    Returns the number of independent loops in the graph G.
    """
    C = nx.number_connected_components(G)
    E = G.number_of_edges()
    N = G.number_of_nodes()

    num_of_independent_loops = E - N + C

    return num_of_independent_loops

def restore_loops_from_polygons(polygons: FeatureLayer, lines: FeatureLayer) -> None:
    pass
#TODO: check if number of polygons after union is larger than the number of polygons before union
#if yes, delete duplicates from union, select union polygons that are identical to the original polygons
#for remaining polygons, the data should be restored.

def CodeToXY(G, Code):
    """
    Returns (x, y) coordinate for a given code.
    You can implement this to query G.nodes[code]['SHAPE@XY'] or some other source
    of coordinate lookup. The exact implementation depends on your data.
    """

    return G.nodes[Code]['SHAPE@XY']
    

def get_cycle_coords_with_midpoints(G, cycle):
    """
    Given a cycle (list of node IDs) in graph G,
    return a list of XY coordinates including any midpoint codes found on edges.
    """
    coords_in_cycle = []
    n = len(cycle)
    for i in range(n):
        current_node = cycle[i]
        next_node = cycle[(i + 1) % n]  # to loop back for the last edge to first node

        # 1) Append the current node's coordinate
        #    Make sure you match the exact key you used for node geometry.
        current_node_xy = G.nodes[current_node]['SHAPE@XY']  # or 'Shape'
        coords_in_cycle.append(current_node_xy)

        # 2) Retrieve edge data (for an undirected graph, either [current][next] or [next][current] works the same)
        edge_data = G[current_node][next_node]

        # 3) If there's a MidPointCodes attribute, parse it
        if 'MidPointCodes' in edge_data and edge_data['MidPointCodes']:
            midpoint_codes_str = edge_data['MidPointCodes']
            midpoint_codes_list = [c.strip() for c in midpoint_codes_str.split(',')]

            # 4) Check the direction of the edge vs. the cycle path
            start_code = edge_data.get("StartPointCode")
            end_code = edge_data.get("EndPointCode")

            if start_code is not None and end_code is not None:
                # If cycle is going from Start to End, keep order
                if current_node == start_code and next_node == end_code:
                    midpoint_codes = midpoint_codes_list

                # If cycle is going from End to Start, reverse order
                elif current_node == end_code and next_node == start_code:
                    midpoint_codes = midpoint_codes_list[::-1]
                else:
                    # If there's some mismatch, either skip or default to forward
                    midpoint_codes = midpoint_codes_list

            # 5) Convert each code to its XY coordinate and append in order
            for m_code in midpoint_codes:
                midpoint_xy = CodeToXY(G,int(m_code))
                coords_in_cycle.append(midpoint_xy)

    # Optionally, if you need a closed loop for polygon geometry,
    # you might want to append coords_in_cycle[0] again at the end.
    return coords_in_cycle

def extract_data_from_cycles(G) -> list:
    """
    Extracts data from cycles in the graph and returns a list of dictionaries.

    Parameters:
        G (networkx.Graph): The input graph.

    Returns:
        list: A list of dictionaries, one per cycle, containing the extracted data.
    """


    generalized_lines = any("MidPointCodes" in data for _, _, data in G.edges(data=True))

    cycles = nx.cycle_basis(G)

    polygons_data = []  # List of dictionaries, one per cycle

    for cycle in cycles:
        # 1) Get all node coordinates for this cycle

        if generalized_lines:
            coords_in_cycle = get_cycle_coords_with_midpoints(G, cycle)
        else:
            coords_in_cycle = [G.nodes[node]['Shape'] for node in cycle]
        
        # 2) Initialize aggregate values
        meas_height_diff = 0
        meas_dist = 0
        misc_bf = 0
        g_corr = 0
        corr_height_diff = 0
        setups_num = 0

        if generalized_lines:
            segments_num = 0

        # 3) We'll store data for each edge in this list
        edges_data_list = []
        
        # 4) Iterate over each edge in the cycle
        #    (the last edge connects the last node back to the first)
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            
            edge_data = G[u][v]  # Access the edge attributes
            start_point = edge_data.get('StartPointCode', 'N/A')
            end_point   = edge_data.get('EndPointCode', 'N/A')
            
            dH_value    = edge_data.get('MeasHeightDiff', 'N/A')
            dist_value  = edge_data.get('MeasDist', 'N/A')
            bf_value    = edge_data.get('MiscBF', 'N/A')
            grav_corr_value = edge_data.get('GravCorrection', 'N/A')
            corr_dH_value = edge_data.get('CorrHeightDiff', 'N/A')
            setups_num_value = edge_data.get('SetUpsNum', 'N/A')

            if generalized_lines:
                mid_point_codes = edge_data['MidPointCodes']
                segments_num_value = edge_data.get('NumOfSegments', 'N/A')
            
            # Update aggregate height difference based on direction
            if u == end_point and v == start_point:
                meas_height_diff -= dH_value
                g_corr -= grav_corr_value
                corr_height_diff -= corr_dH_value
            elif u == start_point and v == end_point:
                meas_height_diff += dH_value
                g_corr += grav_corr_value
                corr_height_diff += corr_dH_value
            
            # Update aggregates
            meas_dist    += dist_value
            misc_bf += bf_value
            setups_num += setups_num_value

            if generalized_lines:
                segments_num += segments_num_value


            # Collect edge-level info in a dictionary
            edge_dict = {
                'u': u,
                'v': v,
                'StartPointCode': start_point,
                'EndPointCode': end_point,
                'MeasHeightDiff': dH_value,
                'MeasDist': dist_value,
                'MiscBF': bf_value,
                'SetUpsNum': setups_num_value,
                'GravCorrection': grav_corr_value,
                'CorrHeightDiff': corr_dH_value
            }

            if generalized_lines:
                edge_dict['MidPointCodes'] = mid_point_codes
                edge_dict['SegmentsNum'] = segments_num_value


            edges_data_list.append(edge_dict)
        
        # 5) Build a dictionary for the entire cycle/polygon
        polygon_dict = {
            'coords': coords_in_cycle,
            'MeasHeightDiff': meas_height_diff,
            'MeasDist': meas_dist,
            'MiscBF': misc_bf,
            'SetUpsNum': setups_num,
            'GravCorrection': g_corr,
            'CorrHeightDiff': corr_height_diff,
            #'NumOfEdges': len(cycle),
            #'NumOfEdges': segments_num,
            'EdgesData': edges_data_list
        }


        if generalized_lines:
            polygon_dict['NumOfEdges'] = segments_num
            polygon_dict['NumOfLines'] = len(cycle)
        else:
            polygon_dict['NumOfEdges'] = len(cycle)
        
        # 6) Append this polygon’s data to the overall list
        polygons_data.append(polygon_dict)

    # 'polygons_data' is now a list of dictionaries. Each dictionary represents one cycle.
    # You can iterate over 'polygons_data' to inspect or export the results.


    return polygons_data




def create_polygons_from_loops(polygons_data, out_fc):
    """
    Creates and populates a polygon feature class from the polygons_data list of dictionaries.
    
    :param polygons_data: List of dictionaries, each with:
        {
          'coords':       [(x1, y1), (x2, y2), ...],
          'MeasHeightDiff':  float,
          'MeasDist':         float,
          'MiscBF':      float,
          'SetUpsNum': int,
          'GravCorrection': float,
          'CorrHeightDiff': float,
          'NumOfEdges': int,
          'EdgesData':   [...]
        }
    :param out_fc: Path to the output feature class (or shapefile).
    :param spatial_ref: An arcpy.SpatialReference object (optional).
    """

    
    ## Extract the folder and the feature class name
    #out_folder = env.workspace  # or use something like os.path.dirname(out_fc)
    #fc_name = out_fc                  # if you're passing just a name, otherwise parse differently
    

    out_folder, fc_name = os.path.split(out_fc)

    CreateFeatureclass(
        out_path=out_folder,
        out_name=fc_name,
        geometry_type="POLYGON",
        spatial_reference=SPATIAL_REFERENCE
    )
    # Add fields for the attributes
    AddField(out_fc, "MeasHeightDiff", "DOUBLE")
    AddField(out_fc, "MeasDist", "DOUBLE")
    AddField(out_fc, "MiscBF", "DOUBLE")
    AddField(out_fc, "SetUpsNum", "SHORT")
    AddField(out_fc, "GravCorrection", "DOUBLE")
    AddField(out_fc, "CorrHeightDiff", "DOUBLE")
    AddField(out_fc, "NumOfSegments", "SHORT")



    # --- 2) Insert polygons ---
    # Define which fields you'll write to: SHAPE@ plus your attribute fields
    fields = ["SHAPE@", "MeasHeightDiff", "MeasDist", "MiscBF", "SetUpsNum", "GravCorrection", "CorrHeightDiff", "NumOfSegments"]

    # Check if any dictionary in polygons_data contains the key "NumOfLines"
    generalized_lines = any("NumOfLines" in polygon for polygon in polygons_data)

    # If "NumOfLines" exists, add the field to the feature class
    if generalized_lines:
        AddField(out_fc, "NumOfLines", "SHORT")
        fields.append("NumOfLines")
    
    with InsertCursor(out_fc, fields) as cur:
        for polygon_dict in polygons_data:
            coords = polygon_dict["coords"]
            meas_dH = polygon_dict["MeasHeightDiff"]
            meas_dist    = polygon_dict["MeasDist"]
            misc_bf     = polygon_dict["MiscBF"]
            setups_num = polygon_dict["SetUpsNum"]
            grav_corr = polygon_dict["GravCorrection"]
            corr_dH = polygon_dict["CorrHeightDiff"]
            num_segments   = polygon_dict["NumOfEdges"]

            if generalized_lines:
                num_lines = polygon_dict["NumOfLines"]

            
            # 2A) Build an arcpy.Array of arcpy.Point objects
            ring = Array()
            
            # Add each (x, y) as a Point
            for (x, y) in coords:
                ring.add(Point(x, y))
            
            # (Optional) Ensure the polygon is "closed" by re-adding the first point at the end.
            # ArcPy is usually fine if the polygon is "not closed," but sometimes you may want this:
            if coords and coords[0] != coords[-1]:
                ring.add(Point(*coords[0]))
            
            # 2B) Create an arcpy.Polygon geometry
            polygon_geom = Polygon(ring, SPATIAL_REFERENCE)
            
            # 2C) Insert the new row with geometry + attributes
            row = (polygon_geom, meas_dH, meas_dist, misc_bf, setups_num, grav_corr, corr_dH, num_segments)
            
            if generalized_lines:
                row += (num_lines,)
                
            cur.insertRow(row)
    
    #print(f"Polygons inserted into {out_fc} successfully!")

def create_graph(levelling_lines, points_data, lines_data_type:Literal['regular_segments', 'generalized_lines']='regular_segments') -> nx.Graph:
    """
    Create a networkx undirected graph from levelling_lines and points_data

    Parameters
    ----------
    levelling_lines : FeatureLayer
        Feature layer of levelling lines
    points_data : FeatureLayer
        Feature layer of control points
    lines_data_type : Literal['regular_segments', 'generalized_lines'], optional
        Type of lines data. The default is 'regular_segments'.
        defines the fields to be used in the graph

    Returns
    -------
    G : nx.Graph
        Networkx graph
    """


    # Initialize an empty undirected graph
    G = nx.Graph()
    
    # Add nodes from points_data with attributes (including geometry)
    fields = [f.name for f in ListFields(points_data)] + ["SHAPE@XY"]

    with SearchCursor(points_data, fields) as cursor:
        for row in cursor:
            node_attributes = {field: value for field, value in zip(fields, row)}
            point_code = node_attributes["Code"]  # Node identifier
            G.add_node(point_code, **node_attributes)  # Store all attributes, including geometry

    # Add edges from levelling_lines with attributes
    fields = [f.name for f in ListFields(levelling_lines)]  # List all fields
    with SearchCursor(levelling_lines, fields) as cursor:
        for row in cursor:
            edge_attributes = {field: value for field, value in zip(fields, row)}
            startPointCode, endPointCode = edge_attributes["StartPointCode"], edge_attributes["EndPointCode"]  # Extract edge nodes
            G.add_edge(startPointCode, endPointCode, **edge_attributes)  # Store all attributes
    
    return G


def create_multigraph(levelling_lines, points_data) -> nx.MultiGraph:
    """
    Create a networkx undirected multigraph from levelling_lines and points_data

    Parameters
    ----------
    levelling_lines : FeatureLayer
        Feature layer of levelling lines
    points_data : FeatureLayer
        Feature layer of control points

    Returns
    -------
    G : nx.MultiGraph
        Networkx multigraph
    """

    # Initialize an empty undirected multigraph
    G = nx.MultiGraph()
    
    # Add nodes from points_data with attributes
    with SearchCursor(points_data, ["code", "full_name", "H", "R"]) as cursor:
        for row in cursor:
            point_code = row[0]
            G.add_node(point_code, point_name=row[1], height=row[2], rank=row[2])

    with SearchCursor(levelling_lines, ["codeA", "codeB", "dH", "dist", "n", "diff"]) as cursor:
        for row in cursor:
            codeA, codeB = row[0], row[1]
            G.add_edge(codeA, codeB,  start_point=codeA, end_point=codeB, dH=round(row[2],5), dist=round(row[3],0), n=round(row[4],0), diff=round(row[5],2))

    return G