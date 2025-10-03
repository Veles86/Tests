import os
import arcpy
from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField

#from Utils.Helpers import extract_coordinates_from_external_source,start_editing, stop_editing, reverse_line_segments, get_first_segment_with_point, add_segment, delete_segment,get_last_segment_end_point_code
#from Utils.Helpers import get_segments_end_point_names,get_segments_end_point_codes, check_self_loops, create_points_dictionary, count_features, mmyy_int_to_yyyymm_str, fill_midpoints_and_dates, fill_dH_and_dist_differences

from Utils.Points import extract_coordinates_from_external_source, create_points_dictionary, add_data_to_points_dictionary, project_IG0512_to_IGD0512_dict
arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True


import numpy as np

import arcpy
import numpy as np

import arcpy
import numpy as np

import arcpy
import numpy as np

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




def leveling_adjustment_constrained_station(
    in_lines_fc,         # Feature class of lines with: pointA, pointB, dH, num_of_set_ups
    out_points_fc,       # Output feature class for adjusted station heights
    station_fc=None,     # (Optional) station FC with approximate heights: StationID, ApproxHeight, SHAPE@XY
    fixed_station_id=None,  # The station ID we want to fix
    fixed_height=100.0,     # The height value for that station
    k=1.0               # Error-model constant for weighting
):
    """
    Perform a 1D leveling adjustment, constraining one known station to a known height
    via Lagrange multipliers.

    Parameters:
    -----------
    in_lines_fc : str
        Line feature class with fields:
          - pointA (start station ID)
          - pointB (end station ID)
          - dH (observed difference: H_B - H_A)
          - num_of_set_ups (used to derive weight)
    out_points_fc : str
        Output point feature class to store station IDs and adjusted heights.
    station_fc : str or None
        Optional feature class with fields:
          - StationID
          - ApproxHeight
          - SHAPE@XY (optional, for geometry)
        If None, we initialize all station heights to 0.0.
    fixed_station_id : str or None
        The ID of the station we want to fix. If None, no constraint is applied
        (this would revert to a free-net solution).
    fixed_height : float
        The height to which we fix the specified station.
    k : float
        Constant in the error model. Example usage: sigma = k * sqrt(num_of_set_ups).

    Result:
    -------
    A minimal-constrained 1D leveling solution, anchored at 'fixed_station_id = fixed_height'.
    """
    arcpy.env.overwriteOutput = True

    # ---------------------------------------------------------------------
    # 1. Collect unique station IDs from the line FC
    # ---------------------------------------------------------------------
    station_ids = set()
    fields_in = ["A", "B", "dH", "n"]
    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            station_ids.add(row[0])
            station_ids.add(row[1])

    station_ids = sorted(list(station_ids))
    num_stations = len(station_ids)

    # Create index mapping: station -> parameter index
    station_index = {sid: i for i, sid in enumerate(station_ids)}

    if fixed_station_id is not None and fixed_station_id not in station_index:
        raise ValueError(f"Fixed station ID '{fixed_station_id}' is not in the list of stations.")

    # ---------------------------------------------------------------------
    # 2. Approximate Heights
    # ---------------------------------------------------------------------
    station_heights = {}
    if station_fc:
        # Read approximate heights if available
        # Expecting fields: StationID, ApproxHeight
        with arcpy.da.SearchCursor(station_fc, ["StationID", "ApproxHeight"]) as scur:
            for sid, approxH in scur:
                station_heights[sid] = float(approxH)

    # Initialize any missing station to 0.0
    for sid in station_ids:
        if sid not in station_heights:
            station_heights[sid] = 0.0

    # ---------------------------------------------------------------------
    # 3. Build A, l, W from line FC
    # ---------------------------------------------------------------------
    A_rows = []
    l_vals = []
    w_vals = []

    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            pointA, pointB, dh_obs, num_setups = row
            dh_obs = float(dh_obs)
            num_setups = float(num_setups)

            if (pointA not in station_index) or (pointB not in station_index):
                continue

            iA = station_index[pointA]
            iB = station_index[pointB]

            # Approx difference
            H_A_approx = station_heights[pointA]
            H_B_approx = station_heights[pointB]
            dh_approx = H_B_approx - H_A_approx

            misclosure = dh_obs - dh_approx

            # Partial derivatives: dH_B=+1, dH_A=-1
            rowA = np.zeros(num_stations)
            rowA[iA] = -1.0
            rowA[iB] = +1.0

            A_rows.append(rowA)
            l_vals.append(misclosure)

            # Weighting
            # Example: sigma = k * sqrt(num_setups)
            if num_setups <= 0:
                # fallback if zero or negative
                sigma = k
            else:
                sigma = k * np.sqrt(num_setups)
            var = sigma**2
            w_val = 1.0/var if var > 0 else 1e12
            w_vals.append(w_val)

    # Convert to NumPy arrays
    A = np.array(A_rows)       # shape: (num_obs, num_stations)
    l = np.array(l_vals)       # shape: (num_obs,)
    W = np.diag(w_vals)        # shape: (num_obs, num_obs)

    # Normal equations
    N = A.T @ W @ A  # (num_stations, num_stations)
    t = A.T @ W @ l  # (num_stations,)

    # ---------------------------------------------------------------------
    # 4. Handle the constraint: H_fixedStation = fixed_height
    # ---------------------------------------------------------------------
    if fixed_station_id is None:
        # If no station is fixed, do a free-net solution with pseudoinverse
        # (Rank-deficient by 1)
        N_pinv = np.linalg.pinv(N)
        x_solution = N_pinv @ t
    else:
        # Minimal constraint using Lagrange multipliers
        idx_fixed = station_index[fixed_station_id]

        # C x = d  => single constraint row
        # c = [0..0, 1, 0..0], same dimension as x
        c = np.zeros(num_stations)
        c[idx_fixed] = 1.0
        d = fixed_height

        # Augment the system
        # S = [[N,   c^T],
        #      [c,    0 ]]
        # RHS = [t, d]
        S = np.block([
            [N, c[:, None]],
            [c[None, :], np.array([[0.]])]
        ])
        rhs = np.concatenate([t, np.array([d])])

        # Solve for [x, lambda]
        X_lambda = np.linalg.solve(S, rhs)
        x_solution = X_lambda[:num_stations]

    # Update station heights from x_solution
    for i, sid in enumerate(station_ids):
        station_heights[sid] = x_solution[i]

    # ---------------------------------------------------------------------
    # 5. Write Output Feature Class
    # ---------------------------------------------------------------------
    if arcpy.Exists(out_points_fc):
        arcpy.Delete_management(out_points_fc)

    # Create a new point FC
    out_path, out_name = os.path.split(out_points_fc)
    sr = Describe(in_lines_fc).spatialReference
    arcpy.management.CreateFeatureclass(
        out_path,
        out_name,
        "POINT",
        spatial_reference=sr
    )

    # Add fields
    arcpy.AddField_management(out_points_fc, "StationID", "TEXT")
    arcpy.AddField_management(out_points_fc, "AdjustedHeight", "DOUBLE")

    # If we have geometry in station_fc, read it:
    station_geometry = {}
    if station_fc:
        # We'll attempt to read SHAPE@XY
        field_names = [f.name for f in arcpy.ListFields(station_fc)]
        if "SHAPE@XY" in field_names:
            with arcpy.da.SearchCursor(station_fc, ["StationID", "SHAPE@XY"]) as scur:
                for sid, xy in scur:
                    station_geometry[sid] = xy

    # Insert final results
    with arcpy.da.InsertCursor(out_points_fc, ["SHAPE@XY", "StationID", "AdjustedHeight"]) as icur:
        for sid in station_ids:
            h_adj = station_heights[sid]
            xy = station_geometry.get(sid, (0.0, 0.0))
            icur.insertRow((xy, sid, h_adj))

    # ---------------------------------------------------------------------
    # 6. Report
    # ---------------------------------------------------------------------
    if fixed_station_id is None:
        print("Free-net leveling adjustment (no station fixed).")
    else:
        print(f"Station '{fixed_station_id}' fixed to {fixed_height}.")
    print("Results written to:", out_points_fc)









def free_net_leveling_lines(
    in_lines_fc,       # Feature class of lines with fields: pointA, pointB, dH, num_of_set_ups
    out_points_fc,     # Output feature class of stations with adjusted heights
    station_fc=None,   # (Optional) existing station FC with approximate heights
    k=1.0              # Error-model constant (example usage in weighting)
):
    """
    Perform a free-net leveling adjustment using line features with fields:
      - pointA (start station ID)
      - pointB (end station ID)
      - dH (observed difference: H_B - H_A)
      - num_of_set_ups (used to derive weight, e.g. W = 1 / (k^2 * num_of_set_ups))

    If station_fc is provided, we look up approximate heights from there.
    Otherwise, we initialize all station heights to 0.

    Args:
        in_lines_fc: Input line feature class of leveling observations.
        out_points_fc: Output feature class where we'll write the adjusted station heights.
        station_fc: (Optional) Feature class with fields [StationID, ApproxHeight].
        k: A constant for the error model (example: sigma = k * sqrt(num_of_set_ups)).
    """

    arcpy.env.overwriteOutput = True

    # ---------------------------------------------------------------------
    # 1. Collect station IDs from lines
    # ---------------------------------------------------------------------
    station_ids = set()  # we will gather unique station IDs (pointA, pointB)
    
    fields_in = ["A", "B", "dH", "n"]
    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            station_ids.add(row[0])  # pointA
            station_ids.add(row[1])  # pointB

    station_ids = sorted(list(station_ids))

    # Create an index to map station ID -> row in the parameter vector
    station_index = {stn_id: i for i, stn_id in enumerate(station_ids)}
    num_stations = len(station_ids)

    # ---------------------------------------------------------------------
    # 2. Approximate Heights
    # ---------------------------------------------------------------------
    # We'll store the approximate heights in a dictionary: station_heights[stn_id] = approx_height
    station_heights = {}

    if station_fc:
        # If we have an existing station feature class with approximate heights:
        #  Fields: StationID, ApproxHeight
        existing_fields = [f.name for f in arcpy.ListFields(station_fc)]
        if "StationID" not in existing_fields or "ApproxHeight" not in existing_fields:
            raise ValueError("station_fc must have fields 'StationID' and 'ApproxHeight'.")

        # read them into a dict
        with arcpy.da.SearchCursor(station_fc, ["StationID", "ApproxHeight"]) as scur:
            for sid, approxH in scur:
                station_heights[sid] = float(approxH)

    # For any station missing in station_heights, initialize with 0.0
    for sid in station_ids:
        if sid not in station_heights:
            station_heights[sid] = 0.0

    # ---------------------------------------------------------------------
    # 3. Build A, l, W from line FC
    # ---------------------------------------------------------------------
    A_rows = []
    l_vals = []
    w_vals = []

    with arcpy.da.SearchCursor(in_lines_fc, fields_in) as cursor:
        for row in cursor:
            pointA = row[0]
            pointB = row[1]
            dH_obs = float(row[2])           # measured difference: H_B - H_A
            num_setups = float(row[3])       # # of setups to derive weighting

            # skip if not recognized
            if (pointA not in station_index) or (pointB not in station_index):
                continue

            # approximate differences
            H_A_approx = station_heights[pointA]
            H_B_approx = station_heights[pointB]
            dH_approx  = H_B_approx - H_A_approx

            # misclosure
            misclosure = dH_obs - dH_approx

            # partial derivatives for (H_B - H_A)
            # d/dH_A = -1
            # d/dH_B = +1
            rowA = np.zeros(num_stations)
            rowA[station_index[pointA]] = -1.0
            rowA[station_index[pointB]] = +1.0

            A_rows.append(rowA)
            l_vals.append(misclosure)

            # Weight
            # example: sigma = k * sqrt(num_of_setups)
            # => variance = (k^2) * num_of_setups
            # => weight = 1 / variance
            sigma = k * np.sqrt(num_setups)
            variance = sigma**2
            w = 1.0 / variance if variance != 0.0 else 1e12  # safeguard
            w_vals.append(w)

    # Convert to NumPy
    A = np.array(A_rows)          # shape: (num_obs, num_stations)
    l = np.array(l_vals)          # shape: (num_obs,)
    W = np.diag(w_vals)           # shape: (num_obs, num_obs)

    # ---------------------------------------------------------------------
    # 4. Normal Equations & Free-Net Solve (pseudoinverse)
    # ---------------------------------------------------------------------
    # N = A^T W A,  t = A^T W l
    N = A.T @ W @ A  # (num_stations, num_stations)
    t = A.T @ W @ l  # (num_stations, )

    # Pseudoinverse to handle rank deficiency (1 for leveling)
    N_pinv = np.linalg.pinv(N)
    x_corrections = N_pinv @ t

    # Update station heights
    for i, sid in enumerate(station_ids):
        station_heights[sid] += x_corrections[i]

    # ---------------------------------------------------------------------
    # 5. Write Output Feature Class
    # ---------------------------------------------------------------------
    # We'll create a point feature class with a field 'StationID' and 'AdjustedHeight'.
    # If you have actual XY coords, you can incorporate them. 
    # For demonstration, we place all points at (0,0) or we can retrieve from station_fc if geometry is available.
    if arcpy.Exists(out_points_fc):
        arcpy.Delete_management(out_points_fc)

    # You might store these as points with dummy XY or if you want real geometry,
    # read geometry from station_fc or from something else.
    out_path, out_name = os.path.split(out_points_fc)
    sr = Describe(in_lines_fc).spatialReference
    arcpy.management.CreateFeatureclass(
        out_path,
        out_name,
        "POINT",
        spatial_reference=sr
    )

    # Add fields
    arcpy.AddField_management(out_points_fc, "StationID", "TEXT")
    arcpy.AddField_management(out_points_fc, "AdjustedHeight", "DOUBLE")

    # If we want XY from station_fc, let's read them into a dictionary
    station_geometry = {}
    if station_fc:
        with arcpy.da.SearchCursor(station_fc, ["StationID", "SHAPE@XY"]) as scur:
            for sid, xy in scur:
                station_geometry[sid] = xy

    # Insert adjusted stations
    with arcpy.da.InsertCursor(out_points_fc, ["SHAPE@XY", "StationID", "AdjustedHeight"]) as icur:
        for sid in station_ids:
            h_adj = station_heights[sid]
            if sid in station_geometry:
                xy = station_geometry[sid]
            else:
                xy = (0.0, 0.0)  # fallback if unknown
            icur.insertRow((xy, sid, h_adj))

    print("Free-net leveling adjustment completed.")
    print("Adjusted station heights written to:", out_points_fc)






def main_function(p1,p2,p3) -> None:
    # Main script execution
    pass



if __name__ == '__main__':

    # Get input parameters from user
    #p1 = arcpy.GetParameterAsText(0) 
    #p2 = arcpy.GetParameterAsText(1) 
    #p3 = arcpy.GetParameterAsText(2) 

  
    gdb_path = fr"E:\git\LevellingThesis\Data\Testings2.gdb"
    points_input = gdb_path + fr"\rosh_haayin_coordinates"
    lines_input = gdb_path + fr"\rosh_haayin_lines"
    output = gdb_path + fr"\rosh_haayin_freenet_output_wo_constrains"


    leveling_adjustment_constrained_station_with_accuracies(
        lines_input,         # Feature class of lines with: pointA, pointB, dH, num_of_set_ups
        output,       # Output feature class for adjusted station heights
        None,     # (Optional) station FC with approximate heights: StationID, ApproxHeight, SHAPE@XY
        fixed_station_id=None,  # The station ID we want to fix
        fixed_height=46.027,     # The height value for that station
        k=1.0               # Error-model constant for weighting
    )
    print("Done")  



    leveling_adjustment_constrained_station(
        lines_input,         # Feature class of lines with: pointA, pointB, dH, num_of_set_ups
        output,       # Output feature class for adjusted station heights
        station_fc=None,     # (Optional) station FC with approximate heights: StationID, ApproxHeight, SHAPE@XY
        fixed_station_id='1206b',  # The station ID we want to fix
        fixed_height=46.027,     # The height value for that station
        k=1.0               # Error-model constant for weighting
    )
    print("Done")
    free_net_leveling_lines(
        lines_input,       # Feature class of lines with fields: pointA, pointB, dH, num_of_set_ups
        output,     # Output feature class of stations with adjusted heights
        station_fc=None,   # (Optional) existing station FC with approximate heights
        k=1.0              # Error-model constant (example usage in weighting)
    )
    print("Done")

    points_dict = create_points_dictionary(points_input)
    #main_function(p1,p2,p3)
    print("Done")

