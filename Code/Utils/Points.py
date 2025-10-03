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
from arcpy.management import Project, AddXY,Append, CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount,Dissolve, Sort, DeleteField, AlterField
from arcpy import CopyFeatures_management as CopyFeatures
from arcpy import Delete_management as Delete, Describe, ListFields

from Utils.Configs import CONTROL_POINTS_FULL_DATA
from Utils.Helpers import get_feature_layer_workspace,start_editing, stop_editing

def get_point_height(point_code, points_source_data = CONTROL_POINTS_FULL_DATA) -> float:
    """
    Get the height of a point by its code from the control points source data.

    Parameters
    ----------
    point_code : int
        The code of the point.
    points_source_data : str, optional
        The path to the control points source data. The default is CONTROL_POINTS_FULL_DATA.

    Returns
    -------
    float
        The height of the point, or None if the point is not found.
    """
    
    where_clause = f"Code = {point_code}"
    
    with SearchCursor(points_source_data, ["H"], where_clause=where_clause) as cursor:
        for row in cursor:
            return row[0]
    
    return None

    
def get_point_names(codes_list) -> list:

    names_in_search_order = []
    codes_in_search_order = []
    where_clause = f"Code IN ({','.join(str(code) for code in codes_list)})"
    with SearchCursor(CONTROL_POINTS_FULL_DATA, ["Code", "Name"], where_clause=where_clause) as cursor:
        for row in cursor:
            code = row[0]
            codes_in_search_order.append(code)
            full_name = row[1]
            if full_name:
                names_in_search_order.append(full_name)
            else:
                names_in_search_order.append(str(code))

    # Ensure the names are in the same order as the codes
    code_to_name = dict(zip(codes_in_search_order, names_in_search_order))


    if len(codes_list) >1:
        ordered_names = [code_to_name[code] for code in codes_list]

        return ordered_names
    else:
        return names_in_search_order

def calculate_geographic_coordinates(control_points_fc) -> None:
    """
    Calculate geographic coordinates for control points

    Parameters
    ----------
    control_points_fc : FeatureLayer
        Control points feature layer.
    """

    gg_dictionary = project_IG0512_to_IGD0512_dict(control_points_fc)

    #wc = get_feature_layer_workspace(control_points_fc)
    #editor = start_editing(wc)


    # Update the control points feature class with the calculated geographic coordinates
    with UpdateCursor(control_points_fc, ["Code", "Longitude", "Latitude"]) as cursor:
        for row in cursor:
            code = row[0]
            if code in gg_dictionary:
                row[1] = gg_dictionary[code]["POINT_X"]
                row[2] = gg_dictionary[code]["POINT_Y"]
                cursor.updateRow(row)

    #stop_editing(editor)


def add_data_to_points_dictionary(main_points_dict:dict, additional_data:dict) -> None:
    """
    Add a data to the main points dictionary from an additional dictionary

    Parameters
    ----------
    main_points_dict : dict
        Main points dictionary to which the data will be added.
    additional_data : dict
        Additional data to be added to the main points dictionary.

    Returns
    -------
    None
        The main points dictionary is updated in place. 

    """ 

    for key, value in additional_data.items():
        if isinstance(value, dict) and key in main_points_dict and isinstance(main_points_dict[key], dict):
            # Recursively merge nested dictionaries
            additional_data(main_points_dict[key], value)
        else:
            # Overwrite or add the key-value pair
            additional_data[key] = value




def project_IG0512_to_IGD0512_dict(input_points:Literal['FeatureLayer'], output_fc = None) -> dict:
    """
    Project the control points from IG0512 to IGD0512

    Parameters
    ----------
    input_points : Literal['FeatureLayer']
        Control points to be projected. 
    add_xy : bool, optional
        Add X and Y coordinates to the output. The default is True.
    """
    
    if not output_fc:
        dont_save_output = True
        output_fc = os.path.join("memory", "gg_coordinates")
    else:
        dont_save_output = False

    arcpy.env.preserveGlobalIds = False

    projected_points = Project(in_dataset=input_points,out_dataset=output_fc,
        out_coor_system='GEOGCS["IGD05(2012)",DATUM["Israeli_Geodetic_Datum_2005(2012)",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]',
        transform_method="IGD05(2012)_To_IG05(2012)_Intermediate_CRS",
        in_coor_system='PROJCS["Israeli_Grid_05-12",GEOGCS["IG05(2012)_Intermediate_CRS",DATUM["IG05(2012)_Intermediate_Datum",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["False_Easting",219529.584],PARAMETER["False_Northing",626907.39],PARAMETER["Central_Meridian",35.20451694444444],PARAMETER["Scale_Factor",1.0000067],PARAMETER["Latitude_Of_Origin",31.73439361111111],UNIT["Meter",1.0]],VERTCS["Israeli Vertical Datum",VDATUM["IGSN 1971"],PARAMETER["Vertical_Shift",0.0],PARAMETER["Direction",1.0],UNIT["Meter",1.0]]',
        preserve_shape="NO_PRESERVE_SHAPE",max_deviation=None,vertical="NO_VERTICAL"
    )
    
    #wc = get_feature_layer_workspace(projected_points)
    #editor = start_editing(wc)
    AddXY(projected_points)

    #AlterField(in_table=projected_points,field="POINT_X",new_field_name="Longitude",new_field_alias="lambda")
    #AlterField(in_table=projected_points,field="POINT_Y",new_field_name="Latitude",new_field_alias="phi")

    #AlterField(in_table=projected_points,field="POINT_X",new_field_name="lambda",new_field_alias="Longitude")
    #AlterField(in_table=projected_points,field="POINT_Y",new_field_name="phi",new_field_alias="Latitude")
    
    #stop_editing(editor)
    #basic_point_dict = create_points_dictionary(input_points)

    #additional_dict_fields = ["longitude", "latitude"]
    additional_dict_fields = ["POINT_X", "POINT_Y"]
    additional_dict = create_points_dictionary(projected_points, additional_dict_fields)
    # Rename POINT_X and POINT_Y
    #additional_dict["Phi"] = additional_dict.pop("POINT_X")
    #additional_dict["Lambda"] = additional_dict.pop("POINT_Y")


    #add_data_to_points_dictionary(basic_point_dict, additional_dict)

    if dont_save_output:
        Delete(output_fc)
    
    del projected_points

    
    return additional_dict







def create_points_dictionary(control_points_fc, fields:list = None) -> dict:
    """
    Create a dictionary of points from a feature class.

    Parameters
    ----------
    control_points_fc : FeatureLayer
        Feature layer of control points.
    fields : list, optional
        List of fields to include in the dictionary. If None, basic fields will be included.
        The default is None.

    Returns
    -------
    dict
        Dictionary of points with code as key and other fields as nested dictionary.

    """

    if fields is None:
        # Fields to include in the dictionary
        fields = [field.name for field in arcpy.ListFields(control_points_fc)]
        #fields = ["code", "type", "full_name", "Y", "X", "H", "R"]
        
    elif "Code" not in fields:
        fields = ["Code"] + fields

    # Create the dictionary with 'code' as the key and other fields as nested dictionary
    points_dict = {}
    with SearchCursor(control_points_fc, fields) as cursor:
        for row in cursor:
            code = row[fields.index("Code")]  # Find the index of the 'code' field
            # Create a dictionary dynamically for each point's attributes
            point_data = {fields[i]: row[i] for i in range(len(fields))}
            # Add to the main dictionary with 'Code' as the key
            points_dict[code] = point_data

    return points_dict




def extract_coordinates_from_external_source(levelling_lines, control_points_source, print_messages = False) -> FeatureLayer:
    """"
    Extract coordinates of the participating point in levelling_lines from control_points_source

    Parameters
    ----------
    levelling_lines : FeatureLayer
        Feature layer of levelling segments or generalized levelling lines
    control_points_source : FeatureLayer
        Feature layer of control points in which the participating point will be searched
    print_messages : bool, optional
        If True, print messages, by default False

    Returns
    -------
    temp_layer : FeatureLayer
        Feature layer of control points participating in the given levelling lines
    """

    fields = ListFields(levelling_lines)
    field_names = [f.name for f in fields]

    if "MidPointCodes" in field_names:
        generalized_lines = True
    else:
        generalized_lines = False

    # Create a set to store unique codes
    unique_codes = set()

    # Read codes from LevellingLines
    if print_messages:
     AddMessage("Reading codes from levelling lines...")


    if generalized_lines:
        with SearchCursor(levelling_lines, ["StartPointCode", "EndPointCode", "MidPointCodes"]) as cursor:
            for row in cursor:
                unique_codes.add(row[0])
                unique_codes.add(row[1])

                if row[2]:
                    string_of_codes = row[2]
                    codes_list = [code.strip() for code in string_of_codes.split(",")]
                    for code in codes_list:
                        unique_codes.add(int(code))

    else:

        with SearchCursor(levelling_lines, ["StartPointCode", "EndPointCode"]) as cursor:
            for row in cursor:
                unique_codes.add(row[0])
                unique_codes.add(row[1])

    # Convert set to a list and format for SQL IN clause
    code_list = ", ".join(map(str, unique_codes))

    # Create a query to filter ControlPoints
    query = f"Code IN ({code_list})"

    # Create a feature layer with the query

    if print_messages:
        AddMessage("Filtering control points...")

    filtered_layer = SelectByAttribute(control_points_source, "NEW_SELECTION", query)

    # Generate a random string of characters
    random_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=8))

    # Define the memory path with the random suffix
    memory_path = fr"memory\extracted_coordinates_{random_suffix}"

    extracted_coordinates = CopyFeatures(filtered_layer, memory_path)
    SelectByAttribute(in_layer_or_view=filtered_layer,selection_type="CLEAR_SELECTION")


    return extracted_coordinates



def extract_coordinates_from_lines_data(levelling_segments) -> FeatureLayer:
    """"
    Extract coordinates of the participating point in levelling_segments from the lines table without additional source
    """
    pass

def get_distance_between_points(x1, y1, x2, y2) -> float:
    """
    Calculate the distance between two points.

    Parameters
    ----------
    x1 : float
        X-coordinate of the first point.
    y1 : float
        Y-coordinate of the first point.
    x2 : float
        X-coordinate of the second point.
    y2 : float
        Y-coordinate of the second point.

    Returns
    -------
    float
        The distance between the two points.
    """
    if (x1 and y1) and (x2 and y2):
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
    
    return None
