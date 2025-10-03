import os
import itertools
import arcpy
from arcpy import GetParameterAsText
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField, SelectLayerByLocation as SelectByLocation,FindIdentical

from arcpy.edit import FlipLine
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Configs import GENERALIZED_LINES_TEMPLATE, MATLAB_FUNCTIONS_DIR
from Utils.Generalization import check_self_loops
from Utils.Helpers import count_features
from Utils.Configs import LOOPS_TEMPLATE

import matlab.engine
import pandas as pd
import numpy as np

arcpy.env.overwriteOutput = True



def feature_data_to_list(feature_layer, fields, where=None) -> list:
    
    output_data = []
    with SearchCursor(feature_layer, fields,where_clause=where) as cursor:
        for row in cursor:
            output_data.append(list(row))


    return output_data






# Helper function to perform the conversion
def matlab_table_to_dataframe(matlab_table, engine):
    """Converts a MATLAB table object to a Pandas DataFrame."""
    # Put the table in MATLAB's workspace to work with it
    engine.workspace['temp_table_from_py'] = matlab_table
    
    # Get column headers and data using table2cell
    headers = engine.eval("temp_table_from_py.Properties.VariableNames")
    data_cell = engine.eval("table2cell(temp_table_from_py)")
    
    # Clean up the MATLAB workspace
    engine.eval("clear temp_table_from_py", nargout=0)
    
    # Create the DataFrame
    return pd.DataFrame(data_cell, columns=headers)


def adjust_measurements(levelling_lines, known_points) -> None:


    # Set Pandas to display dataframes nicely in the console
    pd.set_option('display.width', 100)
    pd.set_option('display.max_columns', 10)

    AddMessage("Starting MATLAB engine...")
    # It's good practice to use a try/finally block to ensure MATLAB is closed
    try:
        eng = matlab.engine.start_matlab()
        AddMessage("MATLAB engine started successfully.")
        # Add the folder containing your .m file to the MATLAB path
        eng.addpath(MATLAB_FUNCTIONS_DIR, nargout=0)


        # 2. Convert Python lists to MATLAB data types
        obs_data_ml = matlab.double(levelling_lines)
        known_data_ml = matlab.double(known_points)

        # 3. Call the MATLAB function
        AddMessage("\nCalling 'adjustLevelingNetwork' in MATLAB...")
        # The MATLAB struct is returned as a Python dictionary
        results = eng.adjustLevellingNetwork(obs_data_ml, known_data_ml)
        AddMessage("Adjustment complete. Processing results in Python.")

        # 4. Process the results
        
        # Simple values can be accessed directly.

        dof = results['DOF']
        adj_type = results['AdjustmentType']
        AddMessage(f"\nAdjustment Type: {adj_type}, Degrees of Freedom: {dof:.0f}")

        # --- Convert MATLAB Tables to Pandas DataFrames ---
        


        # Convert the 'AdjustedHeights' table
        #df_heights = matlab_table_to_dataframe(results['AdjustedHeights'], eng)
        df_heights = results['AdjustedHeights']
        # Convert the 'Statistics' table
        #df_stats = matlab_table_to_dataframe(results['Statistics'], eng)
        df_stats = results['Statistics']
        # 5. Display the DataFrames in Python
        AddMessage("\n--- Adjusted Heights (as Pandas DataFrame) ---")
        AddMessage(df_heights)
        
        AddMessage("\n--- Statistical Analysis (as Pandas DataFrame) ---")
        AddMessage(df_stats)

        # Now you can work with the dataframes as you would with any other Pandas object
        # For example, find the observation with the largest residual
        max_residual_obs = df_stats.loc[df_stats['Residual'].abs().idxmax()]
        AddMessage("\n--- Example: Observation with max absolute residual ---")
        AddMessage(max_residual_obs)


    finally:
        # Ensure the MATLAB engine is stopped
        if 'eng' in locals():
            AddMessage("\nStopping MATLAB engine.")
            eng.quit()







     







if __name__ == '__main__':

    # Get input parameters from user
    levelling_lines = GetParameterAsText(0) 
    point_data = GetParameterAsText(1) 
    sql_filter = GetParameterAsText(2)
    output_results_folder = GetParameterAsText(3)
    output_calculated_points = GetParameterAsText(4)
    output_calculated_lines = GetParameterAsText(5)


    known_points = []
    if sql_filter and point_data:
        with SearchCursor(point_data, ["Code", "H"],where_clause= sql_filter) as cursor:
            for row in cursor:
                if row[1]:
                    known_points.append([row[0], row[1]])
                else:
                    AddWarning(f"Point {row[0]} has no height value. Skipping.")


    fields = ["OBJECTID", "StartPointCode", "EndPointCode", "CorrHeightDiff", "MeasDist"]
    #meas_data = feature_data_to_list(levelling_lines, fields)


    meas_data = []
    with SearchCursor(levelling_lines, fields) as cursor:
        for row in cursor:
            meas_data.append([row[0], row[1], row[2], row[3], 1/row[4]])

    adjust_measurements(meas_data, known_points)

    '''
        # Save df_results as a CSV file
    #output_csv_path = os.path.join(os.path.dirname(points_source_data), "calculated_heights_results.csv")
    output_csv_path = r"E:\OneDrive\Desktop\Code\temp\calculated_heights_results.csv"
    df_results.to_csv(output_csv_path, index=False)
    AddMessage(f"Results saved to {output_csv_path}.")

    '''