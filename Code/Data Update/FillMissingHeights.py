
import arcpy
from arcpy import GetParameterAsText
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Calculations import calculate_missing_heights,calculate_missing_heights_matlab





if __name__ == "__main__":

    points_fc = GetParameterAsText(0)
    lines_fc = GetParameterAsText(1)
    known_points_query = GetParameterAsText(2)



    #lines_fc = r"E:\git\LevellingThesis\Project\temp.gdb\Segments1"
    #points_fc = r"E:\git\LevellingThesis\Project\temp.gdb\ControlPointsLevellingOnly1"
    #known_points_query = "HRank in (1,2,3,4)"
    calculate_missing_heights_matlab(lines_fc, points_fc, known_points_query)
    