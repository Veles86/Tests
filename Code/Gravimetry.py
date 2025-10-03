from typing import Literal
from arcpy import AddError
from math import sin, radians, sqrt

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

def calculate_gravimetric_correction(height_A, dH_AB, gA, gB, method:Literal["normal","orthometric"]="normal") -> float:
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

    if method == "normal":
        equator_normal_g = 978032.67715 #normal gravity at the equator in mGal
        correction = (1/equator_normal_g) * (height_A * (gA - gB) + ((gA + gB)/2 -gB) * dH_AB)
    else:
        helmert_gA = helmert_gravity_correction(height_A, gA)
        helmert_gB = helmert_gravity_correction(height_A + dH_AB, gB)
        correction = (1/helmert_gB) * (height_A * (helmert_gA - helmert_gB) + ((gA + gB)/2 -helmert_gB) * dH_AB)

    return correction

if __name__ == '__main__':

    phi1 = 32.5
    phi2 = 32.6

    grav1 = get_theoretic_normal_gravity(phi1)
    grav2 = get_theoretic_normal_gravity(phi2)

    height_A = 100
    dH_AB = 10.123
    corr = calculate_gravimetric_correction(height_A, dH_AB, grav1, grav2)
    print(f'Gravimetric correction: {corr} m')

