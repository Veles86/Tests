from typing import Literal
import math
from arcpy import AddError
from math import sin, radians, sqrt
from arcpy.da import SearchCursor, Editor, UpdateCursor, InsertCursor

from Utils.Points import create_points_dictionary
from Utils.Helpers import get_feature_layer_workspace,start_editing, stop_editing

def calculate_gravimetric_corrections(levelling_lines_fc, control_points_fc, method:Literal["Normal","Orthometric"]="Normal") -> None:
    """
    Calculate gravimetric corrections for levelling lines feature class
    
    Parameters
    ----------
    levelling_lines_fc : FeatureLayer
        Levelling lines feature layer.
    control_points_fc : FeatureLayer
        Control points feature layer.
    method : Literal["normal","orthometric"]
        The method to use. The default is "normal".

    Returns
    -------
    None
        The levelling lines feature class is updated in place with the calculated gravimetric corrections.

    """
    fields = ["Code", "H", "NormG"]
    points_dict = create_points_dictionary(control_points_fc,fields)

    #wc = get_feature_layer_workspace(levelling_lines_fc)
    #editor = start_editing(wc)
    with UpdateCursor(levelling_lines_fc, ["GlobalID", "StartPointCode", "EndPointCode", "MeasHeightDiff", "GravCorrection", "CorrHeightDiff"]) as cursor:
        for row in cursor:
            from_code = row[1]
            to_code = row[2]
            dH = row[3]
            grav_correction = get_gravimetric_correction(points_dict[from_code]["H"], dH, points_dict[from_code]["NormG"], points_dict[to_code]["NormG"], method=method)
            row[5] = points_dict[row[1]]["H"] - points_dict[row[2]]["H"] #?
            row[4] = grav_correction
            row[5] = dH + grav_correction

            cursor.updateRow(row)

    #stop_editing(editor)





def calculate_theoretic_normal_gravity(control_points_fc, method:Literal['IGF80','classic'] = 'IGF80') -> None:
    """
    Calculate theoretic normal gravity for control points feature class
    
    Parameters
    ----------
    control_points_fc : FeatureLayer
        Control points feature layer.
    method : Literal['IGF80','classic']
        method to use. The default is 'IGF80 (International Gravity Formula 1980)

    Returns
    -------
    None
        The control points feature class is updated in place with the calculated theoretic normal gravity.

    """


    #wc = get_feature_layer_workspace(control_points_fc)
    #editor = start_editing(wc)
    with UpdateCursor(control_points_fc, ["Code", "Latitude","NormG"]) as cursor:
        for row in cursor:
            code = row[0]
            phi = row[1]
            row[2] = get_theoretic_normal_gravity(phi, method)
            cursor.updateRow(row)

    #stop_editing(editor)


def get_theoretic_normal_gravity(phi, method:Literal['IGF80','classic'] = 'IGF80') -> float:
    """
    Calculate the theoretic normal gravity based on GRS80 ellipsoid with classic formula or improved IGF80 formula (that uses Somigliana equation)
    
    Parameters
    ----------
    phi : float
        Latitude in decimal degrees.

    method : Literal['IGF80','classic']
        method to use. The default is 'IGF80 (International Gravity Formula 1980)

    Returns
    -------
    float
        Theoretic normal gravity in mGal at the given latitude

    """

    phi_rad = radians(phi)

    
    equator_normal_g = 978032.67715 #normal gravity at the equator in mGal
    
    if method == 'IGF80':
        grav_phi = equator_normal_g * ((1 + 0.001931851353 * (sin(phi_rad) ** 2)) / sqrt(1 - 0.0066943800229 * (sin(phi_rad) ** 2)))
    else:
        grav_phi = equator_normal_g * (1 + 0.0053024 * (sin(phi_rad) ** 2) + 0.0000058 * (sin(2 * phi_rad) ** 2))
    
    return grav_phi

def helmert_gravity_correction(H, g) -> float:
    """
    Calculate the helmert gravity correction
    
    Parameters
    ----------
    H : float
        The height of the point in meters above the geoid.
    g : float
        The gravity value at the point in mGal.

    Returns
    -------
    float
        The helmert gravity correction in mGal

    """

    if H >= 0:
        return g + 0.0424*H
    else:
        return g + 0.1543*H

def get_gravimetric_correction(height_A, dH_AB, gA, gB, method:Literal["Normal","Orthometric"]="Normal") -> float:
    """
    Calculate the gravimetric_correction for levelling segment AB
    
    Parameters
    ----------
    height_A : float
        The height of point A in meters above the sea level.
    dH_AB : float
        The height difference between point A and B in meters.
    gA : float
        The gravity value at point A in mGal.
    gB : float
        The gravity value at point B in mGal.
    method : Literal["normal","orthometric"]
        The method to use. The default is "normal".

    Returns
    -------
    float
        The gravimetric correction in mGal

    """

    if method == "Normal":
        equator_normal_g = 978032.67715 #normal gravity at the equator in mGal
        correction = (1/equator_normal_g) * (height_A * (gA - gB) + ((gA + gB)/2 -gB) * dH_AB)
    else:
        helmert_gA = helmert_gravity_correction(height_A, gA)
        helmert_gB = helmert_gravity_correction(height_A + dH_AB, gB)
        correction = (1/helmert_gB) * (height_A * (helmert_gA - helmert_gB) + ((gA + gB)/2 -helmert_gB) * dH_AB)

    return correction


