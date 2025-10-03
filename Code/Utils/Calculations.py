import os
import arcpy
import networkx as nx
import numpy as np
import math
import matlab.engine
import random
import string
from typing import Literal
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, ListFields
from arcpy.da import SearchCursor, Editor, UpdateCursor, InsertCursor
from arcpy import MakeFeatureLayer_management as MakeFeatureLayer
from arcpy.edit import FlipLine
from arcpy.management import Append, CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute, SelectLayerByLocation as SelectByLocation,GetCount,Dissolve, Sort, DeleteField
from arcpy import CopyFeatures_management as CopyFeatures
from arcpy import Delete_management as Delete, Describe

from Utils.Generalization import classify_points, sort_segments, generalize_sorted_segments
from Utils.Points import create_points_dictionary, extract_coordinates_from_external_source, get_distance_between_points
from Utils.Helpers import count_features, feature_class_to_dataframe

from Utils.Configs import HEIGHT_CODES_DICTIONARY, MATLAB_FUNCTIONS_DIR

def allowed_misclosure_mm(rank, dist_km) -> float:
    """ 
    Calculate the allowed misclosure based on the criteria and distance in kilometers.
    """

    criteria = 0


    if rank == 1:
        criteria = 3
    elif rank == 2:
        criteria = 5
    elif rank == 3:
        criteria = 10
    elif rank == 4:
        criteria = 20
    elif rank == 5:
        criteria = 30
    elif rank == 6:
        criteria = 60


    return criteria *math.sqrt(dist_km)


def get_max_rank(achived_misclosure_mm,dist_km):
    """
    Calculate the maximum rank based on the misclosure in millimeters.
    """


    if achived_misclosure_mm < allowed_misclosure_mm(1, dist_km):
        return HEIGHT_CODES_DICTIONARY[1]  # Rank H1
    elif achived_misclosure_mm < allowed_misclosure_mm(2, dist_km):
        return HEIGHT_CODES_DICTIONARY[2]  # Rank H2
    elif achived_misclosure_mm < allowed_misclosure_mm(3, dist_km):
        return HEIGHT_CODES_DICTIONARY[3]  # Rank H3
    elif achived_misclosure_mm < allowed_misclosure_mm(4, dist_km):
        return HEIGHT_CODES_DICTIONARY[4]  # Rank H4
    elif achived_misclosure_mm < allowed_misclosure_mm(5, dist_km):
        return HEIGHT_CODES_DICTIONARY[5]  # Rank H5
    elif achived_misclosure_mm < allowed_misclosure_mm(6, dist_km):
        return HEIGHT_CODES_DICTIONARY[6]  # Rank H6
    else:
        return 0


def leveling_adjustment_constrained_station_with_accuracies(
    in_lines_fc,         # Lines with pointA, pointB, dH, num_of_set_ups
    out_points_fc,       # Output FC for adjusted heights
    station_fc=None,     # (optional) FC with approximate heights: StationID, ApproxHeight, SHAPE@XY
    fixed_station_id=None,  # Station to fix (if any). If None => free-net
    fixed_height=100.0,     # The known height for the fixed station
    k=1.0               # Error-model constant
):
    """
    Leveling adjustment in 1D, optionally constraining one station to a known height.
    Also computes residuals, a posteriori variance factor, and station std. dev. (accuracies).

    Args:
        in_lines_fc : Feature class of leveling lines
            Fields required: pointA, pointB, dH, num_of_set_ups
            Meaning: dH = H_B - H_A (observed)
        out_points_fc : Output FC for station results
        station_fc : Optional FC with StationID, ApproxHeight, SHAPE@XY
        fixed_station_id : The ID of the station to fix. If None => free-net
        fixed_height : The height to which we fix that station
        k : used for weighting (e.g. sigma = k * sqrt(num_of_set_ups))

    Output:
        A new point FC with:
         - StationID
         - AdjustedHeight
         - StdDev (1-sigma)
    """
    arcpy.env.overwriteOutput = True

    # ---------------------------------------------------------------------
    # 1. Collect station IDs
    # ---------------------------------------------------------------------
    station_ids = set()
    fields_in = ["A", "B", "dH", "n"]
    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            station_ids.add(row[0])
            station_ids.add(row[1])

    station_ids = sorted(list(station_ids))
    num_stations = len(station_ids)

    # map station -> index
    station_index = {sid: i for i, sid in enumerate(station_ids)}

    # If we have a fixed station ID, check it's valid
    if fixed_station_id is not None and fixed_station_id not in station_index:
        raise ValueError(f"Fixed station '{fixed_station_id}' not in station IDs from {in_lines_fc}.")

    # ---------------------------------------------------------------------
    # 2. Approx heights
    # ---------------------------------------------------------------------
    station_heights = {}
    if station_fc:
        with arcpy.da.SearchCursor(station_fc, ["full_name", "H"]) as scur:
            for sid, approxH in scur:
                station_heights[sid] = float(approxH)

    for sid in station_ids:
        if sid not in station_heights:
            station_heights[sid] = 0.0

    # ---------------------------------------------------------------------
    # 3. Build A, l, W
    # ---------------------------------------------------------------------
    A_list = []
    l_list = []
    w_list = []
    num_obs = 0

    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            pointA, pointB, dh_obs, n_setup = row
            dh_obs = float(dh_obs)
            n_setup = float(n_setup)

            if pointA not in station_index or pointB not in station_index:
                continue

            iA = station_index[pointA]
            iB = station_index[pointB]

            # approximate
            H_A = station_heights[pointA]
            H_B = station_heights[pointB]
            dh_approx = H_B - H_A

            misclosure = dh_obs - dh_approx

            rowA = np.zeros(num_stations)
            rowA[iA] = -1
            rowA[iB] = +1

            A_list.append(rowA)
            l_list.append(misclosure)

            # weighting
            # sigma = k * sqrt(num_of_set_ups)
            if n_setup <= 0:
                sigma = k
            else:
                sigma = k * np.sqrt(n_setup)
            var = sigma**2
            w_val = 1.0 / var if var > 0 else 1e12
            w_list.append(w_val)

            num_obs += 1

    A = np.array(A_list)  # shape (num_obs, num_stations)
    l = np.array(l_list)  # shape (num_obs,)
    W = np.diag(w_list)   # shape (num_obs, num_obs)

    # Normal eqns
    N = A.T @ W @ A  # (num_stations, num_stations)
    t = A.T @ W @ l  # (num_stations,)

    # ---------------------------------------------------------------------
    # 4. Solve: either constrained or free-net
    # ---------------------------------------------------------------------
    if fixed_station_id is None:
        # Free-net -> rank deficient => pseudoinverse
        # x = N^+ t
        # Then Qx = sigma0^2 * N^+
        # We'll get Qx after computing sigma0^2
        # We'll also compute residuals v.
        N_pinv = np.linalg.pinv(N)
        x_solution = N_pinv @ t

        # Later, Q_x = sigma0^2 * N_pinv
        # We'll do that after we get sigma0^2.

        # We'll store N_pinv for the covariance formula
        N_inverse_like = N_pinv
        used_augmented = False

    else:
        # Minimal constraint with Lagrange multipliers
        idx_fixed = station_index[fixed_station_id]

        # c = row vector
        c = np.zeros(num_stations)
        c[idx_fixed] = 1.0
        d = fixed_height

        # Augmented system
        # S = [[N,   c^T],
        #      [c,   0  ]]
        # rhs = [t, d]
        S = np.block([
            [N, c[:, None]],             # c as column
            [c[None, :], np.array([[0.]])]
        ])
        rhs = np.concatenate([t, [d]])

        # Solve
        X_lambda = np.linalg.solve(S, rhs)
        x_solution = X_lambda[:num_stations]
        # lambda_ = X_lambda[num_stations]  # not used directly

        # For the covariance, we need the top-left block of S^-1
        # We'll invert S (which is (n+1)x(n+1)) and then extract
        # S_inv[:n_stations, :n_stations].
        # We'll do that after we have residuals => we get sigma0^2 => multiply by that block
        used_augmented = True
        S_big = S  # store for later
        rhs_big = rhs

        # We'll find S_inv after we compute sigma0^2
        # Q_x = sigma0^2 * [S_inv]_xx

        N_inverse_like = None  # We'll compute it after we compute S_inv

    # Update station heights from x_solution
    station_heights_adj = {}
    for i, sid in enumerate(station_ids):
        station_heights_adj[sid] = x_solution[i]

    # ---------------------------------------------------------------------
    # 5. Compute residuals, sigma0^2
    # ---------------------------------------------------------------------
    # v = l - A*x
    # shape: (num_obs,)
    v = l - A @ x_solution
    vTv = v.T @ W @ v  # scalar

    # Degrees of freedom
    # #obs - #params + #constraints
    # #params = num_stations
    # #constraints = 0 if free-net, else 1 if we fixed one station
    if fixed_station_id is None:
        num_constraints = 0
    else:
        num_constraints = 1

    dof = num_obs - num_stations + num_constraints
    if dof <= 0:
        # In a small or heavily constrained network, dof can be 0 or negative
        # We'll clamp to 1 for safety or raise an error
        dof = 1

    sigma0_sq = vTv / dof
    sigma0 = np.sqrt(sigma0_sq)

    # ---------------------------------------------------------------------
    # 6. Compute Q_x => standard deviations
    # ---------------------------------------------------------------------
    if fixed_station_id is None:
        # free-net: Q_x = sigma0^2 * N_pinv
        Q_x = sigma0_sq * N_inverse_like
    else:
        # constrained: Q_x = sigma0^2 * [S^-1]_xx
        # We must invert S_big
        S_inv = np.linalg.inv(S_big)
        # top-left block is (num_stations x num_stations)
        Q_x_block = S_inv[:num_stations, :num_stations]
        Q_x = sigma0_sq * Q_x_block

    station_stddevs = {}
    for i, sid in enumerate(station_ids):
        # standard deviation = sqrt(diagonal)
        var_i = Q_x[i, i]
        if var_i < 0:
            # numerical issues can cause tiny negative
            var_i = 0
        station_stddevs[sid] = np.sqrt(var_i)

    # ---------------------------------------------------------------------
    # 7. Write Output FC
    # ---------------------------------------------------------------------
    if arcpy.Exists(out_points_fc):
        arcpy.Delete_management(out_points_fc)

    out_path, out_name = os.path.split(out_points_fc)
    sr = Describe(in_lines_fc).spatialReference
    arcpy.management.CreateFeatureclass(
        out_path,
        out_name,
        "POINT",
        spatial_reference=sr
    )

    arcpy.AddField_management(out_points_fc, "StationID", "TEXT")
    arcpy.AddField_management(out_points_fc, "AdjustedHeight", "DOUBLE")
    arcpy.AddField_management(out_points_fc, "StdDev", "DOUBLE")

    # If we have geometry in station_fc, read it
    station_geometry = {}
    if station_fc:
        field_names = [f.name for f in arcpy.ListFields(station_fc)]
        if "SHAPE@XY" in field_names:
            with arcpy.da.SearchCursor(station_fc, ["StationID", "SHAPE@XY"]) as scur:
                for sid, xy in scur:
                    station_geometry[sid] = xy

    # Insert final
    with arcpy.da.InsertCursor(out_points_fc, ["SHAPE@XY", "StationID", "AdjustedHeight", "StdDev"]) as icur:
        for sid in station_ids:
            xy = station_geometry.get(sid, (0.0, 0.0))
            h_adj = station_heights_adj[sid]
            std_dev = station_stddevs[sid]
            icur.insertRow((xy, sid, h_adj, std_dev))

    # ---------------------------------------------------------------------
    # 8. Print Summary
    # ---------------------------------------------------------------------
    if fixed_station_id is None:
        print("Free-net leveling solution (no station fixed).")
    else:
        print(f"Station '{fixed_station_id}' fixed to {fixed_height}.")

    print(f"Number of observations: {num_obs}")
    print(f"Number of stations (params): {num_stations}")
    print(f"Constraints: {num_constraints}")
    print(f"Degrees of freedom: {dof}")
    print(f"Residual sum of squares v^T W v = {vTv:.5f}")
    print(f"Estimated sigma0 = {sigma0:.5f}, sigma0^2 = {sigma0_sq:.5f}")
    print("Results written to:", out_points_fc)

def calculate_missing_heights_matlab(levelling_lines, points_source_data, known_points_query:str = None) -> None:
    """
    Calculate missing heights in the levelling lines using the provided points data and optionally known points query.
    It selects subset from the input lines that have missing heights at least in one of the edges
    (By default missing height is defined as NULL, 0 or missing rank value) and then it generalize the lines data


    """

    if not known_points_query:
        known_points_query = "HRank in (1,2,3,4,5)"
    
    known_points_df = feature_class_to_dataframe(points_source_data, ["Code", "H"], where_clause=known_points_query)

    unknown_points_selection = SelectByAttribute(in_layer_or_view=points_source_data, selection_type="NEW_SELECTION", where_clause=known_points_query, invert_where_clause="INVERT")
    
    unknown_point_codes = []
    with SearchCursor(unknown_points_selection, ["Code"]) as cursor:
        for row in cursor:
            unknown_point_codes.append(row[0])
    unknown_point_codes_str = ",".join(map(str, unknown_point_codes))


    line_selection = SelectByAttribute(levelling_lines, where_clause=f"StartPointCode IN ({unknown_point_codes_str}) OR EndPointCode IN ({unknown_point_codes_str})", selection_type="NEW_SELECTION")
    
    fields = ["OBJECTID", "SegmentCode", "StartPointCode", "EndPointCode", "MeasHeightDiff","MeasDist", "SegmentDirection"]
    lines_for_calculation_df = feature_class_to_dataframe(line_selection, fields)


    # Save df_results as a CSV file
    #output_csv_path = os.path.join(os.path.dirname(points_source_data), "calculated_heights_results.csv")
    output_csv_path = r"E:\OneDrive\Desktop\Code\temp\known_points.csv"
    known_points_df.to_csv(output_csv_path, index=False)
    AddMessage(f"known_points saved to {output_csv_path}.")
    output_csv_path = r"E:\OneDrive\Desktop\Code\temp\lines_for_calculation_df.csv"
    lines_for_calculation_df.to_csv(output_csv_path, index=False)
    AddMessage(f"lines_for_calculation_df saved to {output_csv_path}.")

    lines_for_calculation_df = lines_for_calculation_df.astype({col: float for col in lines_for_calculation_df.select_dtypes(include='int').columns})
    AddMessage("Starting MATLAB engine...")
    # It's good practice to use a try/finally block to ensure MATLAB is closed
    try:
        eng = matlab.engine.start_matlab()
        AddMessage("MATLAB engine started successfully.")
        # Add the folder containing your .m file to the MATLAB path
        eng.addpath(MATLAB_FUNCTIONS_DIR, nargout=0)

        # 3. Call the MATLAB function
        AddMessage("\nCalling MATLAB function...")
        # The MATLAB struct is returned as a Python dictionary
        df_results = eng.calculateApproximatedHeights(lines_for_calculation_df, known_points_df)
        AddMessage("Calculation completed. Processing results in Python.")

        # 4. Process the results
        
        if df_results is not None and not df_results.empty:
            with UpdateCursor(points_source_data, ["Code", "H", "HRank"]) as cursor:
                for row in cursor:
                    code = row[0]
                    if code in df_results['Code'].values:
                        row[1] = df_results.loc[df_results['Code'] == code, 'H'].values[0]
                        row[2] = 31 # Calculated unregistered height
                        cursor.updateRow(row)    
            AddMessage("Filling heights completed.")
        
        # Save df_results as a CSV file
        #output_csv_path = os.path.join(os.path.dirname(points_source_data), "calculated_heights_results.csv")
        output_csv_path = r"E:\OneDrive\Desktop\Code\temp\calculated_heights_results.csv"
        df_results.to_csv(output_csv_path, index=False)
        AddMessage(f"Results saved to {output_csv_path}.")


    finally:
        # Ensure the MATLAB engine is stopped
        if 'eng' in locals():
            AddMessage("\nStopping MATLAB engine.")
            eng.quit()




    lines_with_missing_heights = CopyFeatures(line_selection, fr"memory\lines_with_missing_heights")
    SelectByAttribute(in_layer_or_view=points_source_data,selection_type="CLEAR_SELECTION")
    filtered_coordinates = extract_coordinates_from_external_source(lines_with_missing_heights, points_source_data)

    junction_points_codes, regular_points_codes, dead_end_points_codes = classify_points(lines_with_missing_heights, filtered_coordinates, classified_points_output = None, print_messages = False)


    known_points = []

    # Creating list of known points from the filtered coordinates according to the known_points_query
    with SearchCursor(filtered_coordinates, ['Code'],where_clause=known_points_query) as cursor:
        for row in cursor:
            known_points.append(row[0])



    junction_points_codes = list(set(known_points + junction_points_codes))
    regular_points_codes = [item for item in regular_points_codes if item not in known_points]

    sorted_segments = sort_segments(lines_with_missing_heights, filtered_coordinates,  junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output=None, print_messages = False )
   
    generalize_sorted_segments(sorted_segments, generalized_lines_output=None, print_messages = False)
 
    results = parametric_height_adjustment(sorted_segments, filtered_coordinates, known_points)

    if results:
        with UpdateCursor(points_source_data, ["Code", "H", "HRank"]) as cursor:
            for row in cursor:
                code = row[0]
                if code in results:
                    row[1] = results[code]
                    row[2] = 31 # Calculated unregistered height
                    cursor.updateRow(row)

    del lines_with_missing_heights, filtered_coordinates
    Delete(fr"memory\lines_with_missing_heights")
def calculate_missing_heights(levelling_lines, points_source_data, known_points_query:str = None) -> None:
    """
    Calculate missing heights in the levelling lines using the provided points data and optionally known points query.
    It selects subset from the input lines that have missing heights at least in one of the edges
    (By default missing height is defined as NULL, 0 or missing rank value) and then it generalize the lines data


    """

    if not known_points_query:
        known_points_query = "HRank in (1,2,3,4,5)"
    

    unknown_points_selection = SelectByAttribute(in_layer_or_view=points_source_data, selection_type="NEW_SELECTION", where_clause=known_points_query, invert_where_clause="INVERT")
    
    selected_codes = []
    with SearchCursor(unknown_points_selection, ["Code"]) as cursor:
        for row in cursor:
            selected_codes.append(row[0])
    selected_codes_str = ",".join(map(str, selected_codes))


    line_selection = SelectByAttribute(levelling_lines, where_clause=f"StartPointCode IN ({selected_codes_str}) OR EndPointCode IN ({selected_codes_str})", selection_type="NEW_SELECTION")
    lines_with_missing_heights = CopyFeatures(line_selection, fr"memory\lines_with_missing_heights")
    SelectByAttribute(in_layer_or_view=points_source_data,selection_type="CLEAR_SELECTION")
    filtered_coordinates = extract_coordinates_from_external_source(lines_with_missing_heights, points_source_data)

    junction_points_codes, regular_points_codes, dead_end_points_codes = classify_points(lines_with_missing_heights, filtered_coordinates, classified_points_output = None, print_messages = False)


    known_points = []

    # Creating list of known points from the filtered coordinates according to the known_points_query
    with SearchCursor(filtered_coordinates, ['Code'],where_clause=known_points_query) as cursor:
        for row in cursor:
            known_points.append(row[0])



    junction_points_codes = list(set(known_points + junction_points_codes))
    regular_points_codes = [item for item in regular_points_codes if item not in known_points]

    sorted_segments = sort_segments(lines_with_missing_heights, filtered_coordinates,  junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output=None, print_messages = False )
   
    generalize_sorted_segments(sorted_segments, generalized_lines_output=None, print_messages = False)
 
    results = parametric_height_adjustment(sorted_segments, filtered_coordinates, known_points)

    if results:
        with UpdateCursor(points_source_data, ["Code", "H", "HRank"]) as cursor:
            for row in cursor:
                code = row[0]
                if code in results:
                    row[1] = results[code]
                    row[2] = 31 # Calculated unregistered height
                    cursor.updateRow(row)

    del lines_with_missing_heights, filtered_coordinates
    Delete(fr"memory\lines_with_missing_heights")







def parametric_height_adjustment(levelling_lines:FeatureLayer, points_data:FeatureLayer,known_points:list) -> dict:

    # Retrieve unique point codes from codeA and codeB columns
    codes = set()
    with SearchCursor(levelling_lines, ["StartPointCode", "EndPointCode"]) as cursor:
        for row in cursor:
            codes.update([row[0], row[1]])

    # Create a sorted list of unique codes to define matrix columns

    known_points = sorted(list(set(known_points)))
    unique_codes = [item for item in sorted(codes) if item not in known_points]
    code_index = {code: idx for idx, code in enumerate(unique_codes)}

    all_points_dictionary = create_points_dictionary(points_data)
    known_points_dictionary = {k: v for k, v in all_points_dictionary.items() if k not in unique_codes}

    # Initialize incident matrix A and vector L
    num_segments = count_features(levelling_lines)
    A = np.zeros((num_segments, len(unique_codes)))
    Lb = np.zeros(num_segments)
    L0 = np.zeros(num_segments)

    fields = ["StartPointCode", "EndPointCode", "MeasHeightDiff"]
    with SearchCursor(levelling_lines, fields) as cursor:
        for i, row in enumerate(cursor):
            codeA, codeB, dH = row
            if codeA in known_points:
                L0[i] = L0[i] + known_points_dictionary[codeA].get('H', 0)
            else:
                A[i, code_index[codeA]] = -1  # Start point

            if codeB in known_points:
                L0[i] = L0[i] - known_points_dictionary[codeB].get('H', 0)
            else:
                A[i, code_index[codeB]] = 1   # End point

            Lb[i] = dH                     # dH value


    L = Lb - L0

    # Solve the system of linear equations
    A_T = A.T
    ATA = A_T@A
    ATL = A_T@L
    ATA_inv = np.linalg.inv(ATA)
    X = ATA_inv@ATL

    # Create a dictionary of point codes and their corresponding heights
    point_heights = {}
    for i, code in enumerate(unique_codes):
        point_heights[code] = X[i]

    #np.savetxt('matA.csv', A)
    #np.savetxt('L0.csv', L0)
    #np.savetxt('Lb.csv', Lb)
    #np.savetxt('L.csv', L)
    #np.savetxt('X.csv', X)


    return point_heights