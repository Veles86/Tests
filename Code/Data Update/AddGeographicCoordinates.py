
from arcpy import GetParameterAsText


import os
import sys
sys.path.append(os.path.abspath(".."))

from Utils.Points import calculate_geographic_coordinates

if __name__ == "__main__":

    points_fc = GetParameterAsText(0)
    calculate_geographic_coordinates(points_fc)
