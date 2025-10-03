import os
import arcpy
from arcpy import AddMessage, AddError, GetParameterAsText
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy import Delete_management as Delete, Describe
from arcpy import CopyFeatures_management as CopyFeatures
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute
from arcpy import MakeFeatureLayer_management as MakeFeatureLayer
from Utils.Points import extract_coordinates_from_external_source, calculate_geographic_coordinates
from Utils.Helpers import bool_str_to_bool, count_features
from Utils.Gravimetry import calculate_theoretic_normal_gravity, calculate_gravimetric_corrections
from Utils.Calculations import calculate_missing_heights

def filter_control_points(levelling_lines_fc, control_points_fc, output_fc, calc_gg_coordinates, calc_norm_g, fill_missing_heights, known_points_query, calc_grav_corr):
    """
    Creates a filtered feature class of control points with only those that are participating in the levelling lines

    Parameters
    ----------
    levelling_lines_fc : FeatureLayer
        Feature layer of levelling lines
    control_points_fc : FeatureLayer
        Feature layer of control points
    output_fc : str
        Path to the output feature class
    calc_gg_coordinates : bool
        True to calculate gg coordinates in the control points feature class, False otherwise
    calc_norm_g : bool
        True to calculate norm g in the control points feature class, False otherwise
    fill_missing_heights : bool
        True to fill missing heights in the control points feature class, False otherwise
    known_points_query : str
        Query for known points
    calc_grav_corr : bool
        True to calculate gravity correction in the lines feature class, False otherwise
    """

    filtered_points_fc = extract_coordinates_from_external_source(levelling_lines_fc, control_points_fc, True)


    if calc_gg_coordinates:
        AddMessage("Calculating gg coordinates...")
        calculate_geographic_coordinates(filtered_points_fc)
        AddMessage("Geographic coordinates have been calculated")

        if calc_norm_g:
            AddMessage("Calculating theoretic normal gravity...")
            calculate_theoretic_normal_gravity(filtered_points_fc)
            AddMessage("Theoretic normal gravity has been calculated")



    if fill_missing_heights:
        AddMessage("Filling missing heights...")
        if not known_points_query:
            AddMessage("Empty query was recieved, default query will be used ('HRank in (1,2,3,4,5)')")
            known_points_query = "HRank in (1,2,3,4,5)"
        calculate_missing_heights(levelling_lines_fc, filtered_points_fc, known_points_query)


        AddMessage("Missing heights have been filled")

        if calc_grav_corr and calc_gg_coordinates and calc_norm_g:
            AddMessage("Calculating gravity corrections...")
            SelectByAttribute(levelling_lines_fc,"CLEAR_SELECTION")
            calculate_gravimetric_corrections(levelling_lines_fc, filtered_points_fc)
            AddMessage("Gravity corrections have been calculated")




    # Copy the filtered features to a new feature class
    AddMessage(f"Creating output feature class: {output_fc}")
    CopyFeatures(filtered_points_fc, output_fc)

    # Clean up
    Delete(filtered_points_fc)

    AddMessage(f"Filtered control points have been created")

if __name__ == '__main__':

    # Get input parameters
    levelling_lines_fc = GetParameterAsText(0)
    control_points_fc = GetParameterAsText(1)
    output_fc = GetParameterAsText(2)
    calc_gg_coordinates = bool_str_to_bool(GetParameterAsText(3))
    calc_norm_g = bool_str_to_bool(GetParameterAsText(4))
    fill_missing_heights = bool_str_to_bool(GetParameterAsText(5))
    known_points_query = GetParameterAsText(6)
    calc_grav_corr = bool_str_to_bool(GetParameterAsText(7))

    #levelling_lines_fc = "E:\git\LevellingThesis\Data\Testing4.gdb" + "\\" +  "SegmentsAreaB"
    #control_points_fc = "E:\git\LevellingThesis\Data\LevellingDB_new.gdb" + "\\" + "ControlPoints"
    #output_fc = "E:\git\LevellingThesis\Data\Testing4.gdb" + "\\" +  "SegmentsAreaB_coordinates"
    #calc_gg_coordinates = True
    #calc_norm_g = True
    #fill_missing_heights = True
    #known_points_query = "HRank in (1,2,3,4,5)"
    #calc_grav_corr = True



    filter_control_points(levelling_lines_fc, control_points_fc, output_fc, calc_gg_coordinates, calc_norm_g, fill_missing_heights, known_points_query, calc_grav_corr)
