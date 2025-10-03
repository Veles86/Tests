import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

from arcpy.edit import FlipLine
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Configs import GENERALIZED_LINES_TEMPLATE, MATLAB_FUNCTIONS_DIR
from Utils.Generalization import check_self_loops
import matlab.engine

arcpy.env.overwriteOutput = True


def optimize_loop(py_levelling_data):



# --- 1. Setup and Configuration ---


    # Set the path to the folder containing your 'calculateLoopCombinations.m' file.
    # IMPORTANT: You MUST change this path to match the location on your system.
    matlab_script_folder = MATLAB_FUNCTIONS_DIR

    # Check if the path exists before proceeding
    if not os.path.isdir(matlab_script_folder):
        arcpy.AddError(f"Error: The specified MATLAB script folder does not exist: {matlab_script_folder}")
        # Exit the script if the path is invalid
        exit()

    arcpy.AddMessage("Starting MATLAB engine...")
    # It's good practice to wrap the engine in a try/finally block
    # to ensure it quits even if an error occurs.
    eng = None
    try:
        # Start a new MATLAB session
        eng = matlab.engine.start_matlab()
        arcpy.AddMessage("MATLAB engine started successfully.")

        # Add the folder containing your .m file to the MATLAB path
        eng.addpath(matlab_script_folder, nargout=0)
        arcpy.AddMessage(f"Added '{matlab_script_folder}' to MATLAB path.")



        # --- 2. Prepare Data ---
        # In a real-world script, you would get this data from a feature class table
        # using an arcpy.da.SearchCursor.
        # For this example, we use the same data from the MATLAB example.
        # Data format: [LineID, StartPt, EndPt, Distance, dH]
        # Duplicate line, but measured 102 -> 101

        arcpy.AddMessage("\nInput levelling data (from Python):")
        for row in py_levelling_data:
            arcpy.AddMessage(str(row))

        # Convert the Python list of lists into a MATLAB matrix (matlab.double)
        matlab_levelling_data = matlab.double(py_levelling_data)


        # --- 3. Call the MATLAB Function ---
        arcpy.AddMessage("\nCalling 'calculateLoopCombinations' function in MATLAB...")
        # The 'nargout=1' argument specifies that we expect one output from the function.
        results = eng.calculateLoopCombinations(matlab_levelling_data, nargout=1)
        arcpy.AddMessage("Calculation complete. Processing results from MATLAB.")


        # --- 4. Process and Display Results ---
        # The 'results' variable is a MATLAB struct array, which the engine
        # converts into a Python list of dictionaries.
        arcpy.AddMessage(f"\nFound {len(results)} possible loop combinations.")
        best_corr_misclosure = float('inf')
        for i, combo in enumerate(results):
            # The 'combination' field is a matlab.double matrix. We convert it to a
            # Python list of lists for easier handling and display.
            combo_lines = [list(row) for row in combo['combination']]

            # Print results to the ArcGIS geoprocessing messages window
            #arcpy.AddMessage(f"\n--- Combination {i+1} ---")
            #arcpy.AddMessage("  Lines Used ([ID, Start, End, Dist, dH, SetUps, MiscBF,CorrectedHeightDiff]):")
            #for line in combo_lines:
                #arcpy.AddMessage(f"    {line}")
            
            # Access the other fields directly
            total_dist = combo['totalDistance']
            misclosure = combo['misclosure']
            total_setups = combo['totalSetups']
            total_bf_misclosure = combo['totalBFMisclosure']
            total_corrected_height_diff = combo['misclosureWithGravity']

            
            #arcpy.AddMessage(f"  Total Distance: {total_dist:.0f}")
            #arcpy.AddMessage(f"  Misclosure: {misclosure:.5f}")
            #arcpy.AddMessage(f"  Total Setups: {total_setups:.0f}")
            #arcpy.AddMessage(f"  Misclosure BF: {total_bf_misclosure:.2f}")
            #arcpy.AddMessage(f"  Misclosure with Gravity: {total_corrected_height_diff:.5f}")

            best_corr_misclosure = min(best_corr_misclosure, total_corrected_height_diff)

    except Exception as e:
        arcpy.AddError(f"An error occurred: {e}")

    finally:
    # --- 5. Clean Up ---
        if eng:
            arcpy.AddMessage("\nQuitting MATLAB engine.")
            eng.quit()
        arcpy.AddMessage(f"Best corrected misclosure found: {best_corr_misclosure:.5f}")



        







if __name__ == '__main__':

    py_levelling_data = [
    [1,409385,565896,817,3.16529,12,-1.44,3.165288],
    [2,409385,565896,4620,3.15752,107,1.09,3.157519],
    [3,409385,414208,28486,-6.66405,701,9.69,-6.664264],
    [4,565896,414469,6125,-6.34492,145,5.65,-6.344961],
    [6,414235,566072,5546,-24.79981,164,1.99,-24.799812],
    [7,414235,414203,3136,-26.84697,103,-5.89,-26.846916],
    [8,566072,414212,3990,15.67358,99,-0.05,15.673575],
    [9,566072,414212,2712,15.67161,81,1.5,15.671613],
    [10,414208,414212,2623,-2.01137,73,-5.61,-2.011383],
    [12,414203,409366,10977,23.97645,354,5.76,23.976542],
    [13,414203,409366,22146,23.9778,743,10.5,23.977916],
    [14,409395,414469,6816,-6.94792,257,6.03,-6.947949],
    [15,409395,414469,9260,-6.9478,291,6.08,-6.947819],
    [16,409395,409366,1582,-6.17606,44,0.2,-6.176065]]

    optimize_loop(py_levelling_data)
