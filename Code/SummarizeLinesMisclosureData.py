import os
import math
import networkx as nx
import arcpy
from arcpy import GetParameterAsText, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, AddMessage, AddError, AddWarning, Exists, Array, Point, Polygon, SpatialReference,env
from arcgis.features import FeatureLayer
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,SelectLayerByLocation as SelectByLocation,GetCount, Dissolve, Sort, Append, DeleteField
from typing import Literal

from arcpy.analysis import Union
from arcpy.management import DeleteIdentical
from Utils.Helpers import start_editing, stop_editing,get_feature_layer_workspace,count_features

from Utils.Loops import create_graph, extract_data_from_cycles, create_polygons_from_loops, create_multigraph, num_of_independent_loops
from Utils.Points import extract_coordinates_from_external_source
from Utils.Calculations import get_max_rank
from Utils.Configs import SPATIAL_REFERENCE, CONTROL_POINTS_FULL_DATA

arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True





def summarize_loops_data(loops_input):

    with UpdateCursor(loops_input, ['MeasHeightDiff', 'MeasDist', 'MiscBF','CorrHeightDiff','AbsMeasHeightDiffMM','AbsCorrHeightDiffMM','MaxRankByMeasHeightDiff','MaxRankByCorrHeightDiff']) as cursor:
        for row in cursor:
            meas_dist_km = row[1] / 1000.0 
            row[4] = abs(row[0]) * 1000.0
            row[5] = abs(row[3]) * 1000.0

            row[6] = get_max_rank(row[4], meas_dist_km)
            row[7] = get_max_rank(row[5], meas_dist_km)
            cursor.updateRow(row)

if __name__ == '__main__':

    # Get input parameters from user
    loops_input = GetParameterAsText(0)  # Input feature class of selected segments or lines


    summarize_loops_data(loops_input)
    AddMessage("Loops data summarized successfully.")

  







