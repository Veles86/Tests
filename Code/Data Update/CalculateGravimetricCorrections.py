
import arcpy
import os
import sys
sys.path.append(os.path.abspath(".."))
from Utils.Gravimetry import calculate_gravimetric_corrections




if __name__ == "__main__":

    points_fc = arcpy.GetParameterAsText(0)
    lines_fc = arcpy.GetParameterAsText(1)
    method = arcpy.GetParameterAsText(2)  # "Normal" or "Orthometric"
    if method not in ["Normal", "Orthometric"]:
        method = "Normal"

    #lines_fc = "E:\git\LevellingThesis\Data\Testing5.gdb" + "\\" +  "Izun_lines"
    #points_fc = "E:\git\LevellingThesis\Data\Testing5.gdb" + "\\" +  "Izun_coordinates2"
    calculate_gravimetric_corrections(lines_fc, points_fc,method=method)
