import os
import arcpy
import networkx as nx
import numpy as np
import math
import random
import string
import pandas as pd
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, ListFields, Describe
from arcpy.da import SearchCursor, Editor, UpdateCursor, InsertCursor
from arcpy import MakeFeatureLayer_management as MakeFeatureLayer
from arcpy.edit import FlipLine
from arcpy.management import Append, CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount,Dissolve, Sort, DeleteField
from arcpy import CopyFeatures_management as CopyFeatures
from arcpy import Delete_management as Delete, Describe


#arcpy.env.preserveGlobalIds = False

def feature_class_to_dataframe(feature_class, field_list, where_clause=None):
    """
    Converts a feature class or a subset of it to a Pandas DataFrame.

    Args:
        feature_class (str): The path to the feature class or layer name.
        field_list (list): A list of field names to include in the DataFrame.
        where_clause (str, optional): An SQL WHERE clause to filter rows. 
                                      Defaults to None, which means all features
                                      will be included.

    Returns:
        pandas.DataFrame: A DataFrame containing the data from the feature class.
                          Returns an empty DataFrame if an error occurs or no
                          features are found.
    """
    try:
        # Use SearchCursor to read the data into a list of tuples
        with SearchCursor(feature_class, field_list, where_clause) as cursor:
            # Use a list comprehension for efficient conversion
            # Ensure cursor is valid and field_list contains valid fields
            if cursor is None:
                raise ValueError("Cursor could not be initialized. Check the feature class and field list.")
            data_list = [row for row in cursor]

        # Create the Pandas DataFrame from the list of rows and field names
        if not data_list:
            print(f"Warning: No features found in '{feature_class}' with the given where clause.")
            return pd.DataFrame(columns=field_list)
        
        df = pd.DataFrame(data_list, columns=field_list)
        return df

    except arcpy.ExecuteError:
        # Catch errors from arcpy (e.g., feature class not found, invalid field)
        print(f"ArcPy Error: {arcpy.GetMessages(2)}")
        return pd.DataFrame(columns=field_list) # Return empty dataframe on error
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return pd.DataFrame(columns=field_list) # Return empty dataframe on error
    

def bool_str_to_bool(bool_str: str) -> bool:
    """
    Convert a string representation of a boolean to a boolean.

    Parameters
    ----------
    bool_str : str
        The string representation of a boolean.

    Returns
    -------
    bool
        The boolean value.
    """

    return bool_str.lower() == "true"


def count_features(feature_class) -> int:
    """
    returns the number of features in the given feature class

    Parameters
    ----------
    feature_class : FeatureLayer
        Feature layer of line segments

    Returns
    -------
    num_of_features : int
        Number of features in the feature class
    """
    
    return int(GetCount(feature_class)[0])


def get_feature_layer_workspace(feature_layer: FeatureLayer) -> str:
    
    path  =  Describe(feature_layer).path

    if path and not path.lower().endswith(".gdb"):
        path = os.path.dirname(path)


    return path

def start_editing(workspace: str) -> Editor:
    """
    Start editing in a geodatabase workspace using the arcpy.da.Editor class.

    Parameters:
    - workspace (str): The path to the geodatabase workspace.
    """

    editor: Editor = Editor(workspace=workspace, multiuser_mode=True)
    editor.startEditing(with_undo=False, multiuser_mode=True)
    editor.startOperation()
    return editor


def stop_editing(editor: Editor, save: bool = True) -> None:
    """
    Stop editing in a geodatabase workspace using the arcpy.da.Editor class.

    Parameters:
        editor (Editor): The current active editor object.
        save (bool): True to save the edits (This is the default), False to discard changes.
    """
    editor.stopOperation()
    editor.stopEditing(save_changes = save)








def mmyy_int_to_yyyymm_str(mmyy_int: int) -> str:
    """
    Convert an integer date in the format MMYY to a string date in the format YYYYMM.

    Parameters
    ----------
    mmyy_int : int
        The integer date in the format MMYY.

    Returns
    -------
    str
        The string date in the format YYYYMM.
    """

    year = mmyy_int % 100
    month = mmyy_int // 100

    if year < 25:
        year += 2000
    else:
        year += 1900

    return f"{year}/{str(month).zfill(2)}"



