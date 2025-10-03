from arcpy import SpatialReference

SPATIAL_REFERENCE = SpatialReference(6991)
ROOT_DATA_LOCATION = fr"E:\git\LevellingThesis\Data\LevellingDB.gdb"
GENERALIZED_LINES_TEMPLATE = ROOT_DATA_LOCATION + fr"\Lines"
LOOPS_TEMPLATE = ROOT_DATA_LOCATION + fr"\Loops"
SEGMENTS_DATA = ROOT_DATA_LOCATION + fr"\Segments"
# CONTROL_POINTS_FULL_DATA = ROOT_DATA_LOCATION + fr"\ControlPoints"
CONTROL_POINTS_FULL_DATA = ROOT_DATA_LOCATION + fr"\ControlPointsLevellingOnly"
DUPLICATED_CONTROL_POINTS_DATA = ROOT_DATA_LOCATION + fr"\DuplicatedControlPoints"

ATTRIBUTE_RULES_PATH = r"E:\git\LevellingThesis\Scheme\Attribute Rules"
MATLAB_FUNCTIONS_DIR = r"E:\git\LevellingThesis\Code\MATLAB\functions"


HEIGHT_CODES_DICTIONARY = {
    1: 21, # H1
    2: 22, # H2
    3: 23, # H3
    4: 24, # H4
    5: 25, # H5
    6: 26  # H6
}
