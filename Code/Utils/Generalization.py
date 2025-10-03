import os
import arcpy

from arcgis.features import FeatureLayer
from arcpy import AddMessage, AddError, AddWarning, Delete_management as Delete, Describe, CopyFeatures_management as CopyFeatures, CreateTable_management as CreateTable, ListFields
from arcpy.da import SearchCursor, InsertCursor, UpdateCursor
from arcpy.management import CreateFeatureclass, AddField, SelectLayerByAttribute as SelectByAttribute,GetCount, Dissolve, Sort, Append, DeleteField,EnableEditorTracking, SelectLayerByLocation as SelectByLocation, MakeTableView
from arcpy.analysis import Statistics, Frequency

from collections import Counter
from collections import defaultdict


from Utils.Points import create_points_dictionary, create_points_dictionary, get_distance_between_points,extract_coordinates_from_external_source, get_point_names

from Utils.Lines import reverse_segments_new, reverse_line_segments_without_copy,reverse_line_segments_from_tuple, add_segment_from_tuple, get_single_segment_fc, fill_dH_and_dist_differences, reverse_line_segments, get_first_segment_with_point, add_segment, delete_segment,get_last_segment_end_point_code, get_last_segment,reverse_generalized_lines
from Utils.Helpers import count_features, mmyy_int_to_yyyymm_str,start_editing, stop_editing, get_feature_layer_workspace

from Utils.Configs import GENERALIZED_LINES_TEMPLATE, MATLAB_FUNCTIONS_DIR

arcpy.env.overwriteOutput = True
#arcpy.env.preserveGlobalIds = True


def check_proper_division(sorted_segments, points_classification) -> None:

    junctions_and_dead_ends = SelectByAttribute(in_layer_or_view=points_classification,selection_type="NEW_SELECTION",where_clause="Type <> 'regular'")
    temp_dissolved = Dissolve(in_features=sorted_segments, out_feature_class=r"memory\temp_dissolved", dissolve_field="line_num", statistics_fields="StartPointCode FIRST;EndPointCode LAST", multi_part="MULTI_PART",unsplit_lines="DISSOLVE_LINES",concatenation_separator=",")

    num_of_junctions = count_features(junctions_and_dead_ends)

    intersected_junctions = SelectByLocation(in_layer=points_classification,overlap_type="BOUNDARY_TOUCHES",select_features=temp_dissolved,search_distance=None,selection_type="SUBSET_SELECTION")

    num_of_intersected_junctions = count_features(intersected_junctions)

    if num_of_junctions != num_of_intersected_junctions:
        AddWarning(f"Found {num_of_junctions - num_of_intersected_junctions} from {num_of_junctions} junction points that are not connected to any segment. Please check the input data.")
    else:
        AddMessage(f"✅ All {num_of_junctions} junction points are connected to at least one segment.")


    regular_points = SelectByAttribute(in_layer_or_view=points_classification, selection_type="NEW_SELECTION", where_clause="Type = 'regular'")
    num_of_regular_points = count_features(regular_points)
    regular_points_at_edges = SelectByLocation(in_layer=regular_points,overlap_type="BOUNDARY_TOUCHES",select_features=temp_dissolved,search_distance=None,selection_type="SUBSET_SELECTION")
    num_of_regular_points_at_edges = count_features(regular_points_at_edges)

    if num_of_regular_points_at_edges != 0:
        AddWarning(f"Found {num_of_regular_points_at_edges} from {num_of_regular_points} regular points that are touching generalized lines' edges. Please check the input data.")
    else:
        AddMessage(f"✅ All {num_of_regular_points} regular points are not touching generalized lines' edges.")

    Delete(temp_dissolved)
    SelectByAttribute(in_layer_or_view=points_classification, selection_type="CLEAR_SELECTION")


def check_lines_connectivity(sorted_segments, output_table):
    """
    Checks connectivity of segments grouped by line_num and writes the result to an output table.
    Each row will contain line_num and disconnection_count.
    """

    # Step 1: Group segments by line_num
    lines_dict = defaultdict(list)
    with SearchCursor(sorted_segments, ["line_num", "StartPointCode", "EndPointCode"]) as cursor:
        for line_num, start_code, end_code in cursor:
            lines_dict[line_num].append((start_code, end_code))

    # Step 2: Analyze connectivity
    disconnection_counts = {}
    for line_num, segments in lines_dict.items():
        if len(segments) <= 1:
            disconnection_counts[line_num] = 0
            continue

        count = 0
        for i in range(1, len(segments)):
            prev_end = segments[i - 1][1]
            curr_start = segments[i][0]
            if prev_end != curr_start:
                count += 1
        disconnection_counts[line_num] = count

    total_disconnections = sum(disconnection_counts.values())
    if total_disconnections == 0:
        AddMessage("✅ All lines are connected without disconnections.")
    else:
        AddWarning(f"Found {total_disconnections} disconnections in the lines.")

    # Step 3: Create output table
           
    out_path, out_name = os.path.split(output_table)
    CreateTable(out_path=out_path, out_name=out_name)

    AddField(output_table, "line_num", "LONG")
    AddField(output_table, "disconnection_count", "LONG")

    # Step 4: Insert results into table
    with InsertCursor(output_table, ["line_num", "disconnection_count"]) as cursor:
        for line_num, count in disconnection_counts.items():
            cursor.insertRow((line_num, count))





def fill_midpoints_and_adj2025(generalized_lines, dissolved_lines, sorted_segments) -> None:


    # Step 1: Read and process data from the source feature class
    adj2025_dict = {}
    mid_points_dict = {}

    with SearchCursor(dissolved_lines, ["line_num", "CONCATENATE_Adj2025","CONCATENATE_StartPointCode","COUNT_line_num"]) as cursor:
        for row in cursor:
            key = row[0]
            adj2025_str = row[1]
            mid_points_str = row[2]


            if adj2025_str:

                #digits = [int(d) for d in adj2025_str.split(',') if d.strip().isdigit()]
                codes = [int(d.strip()) for d in adj2025_str.split(',')]
                codes = list(set(codes))  # Remove duplicates

                #count = Counter(digits)
                #sorted_items = sorted(count.items())  # sort by digit
                #output = ", ".join(f"{v}x{k}" for k, v in sorted_items)

                if len(codes) == 1:
                    output = codes[0]
                else:
                    output = 17 # If there are multiple codes, we set it to 17 as defined in the domain
                adj2025_dict[key] = output

            if mid_points_str and row[3] > 1:  # Only process if there are multiple segments
                output = ",".join(mid_points_str.split(",")[1:])
                mid_points_dict[key] = output

    # Step 2: Write transformed values to the target feature class
    with UpdateCursor(generalized_lines, ["UniqueRowNum", "Adj2025","MidPointCodes","MidPointNames"]) as cursor:
        for row in cursor:
            key = row[0]
            updated_flag = False
            if key in adj2025_dict:
                row[1] = adj2025_dict[key]
                updated_flag = True

            if key in mid_points_dict:
                temp_str = mid_points_dict[key]
                row[2] = temp_str
                try:
                    point_codes = [int(x.strip()) for x in temp_str.split(",")]
                except:
                    AddWarning(f"Failed to convert midpoints string '{temp_str}' to point codes for line {key}. Unique Row Num: {row[0]}")
                    # For rare cases where the midpoints string contains code in scientific notation e.g. "3e+04" instead of "30000"
                    point_codes = []

                    is_first_row = True
                    with SearchCursor(sorted_segments, ["StartPointCode","StartPointName"], f"line_num = {key}") as search_cursor:
                        for sorted_row in search_cursor:
                            if is_first_row:
                                is_first_row = False
                            else:
                                point_codes.append(sorted_row[0])

                    

                
                point_names = get_point_names(point_codes)
                if row[2]:
                    row[3] = " to ".join(point_names)
                else:
                    row[3] = None
                updated_flag = True



            if updated_flag:
                cursor.updateRow(row)
                









def fill_precise_levelling_source_summary(generalized_lines, dissolved_lines) -> None:

    # Field names
    source_id_field = "line_num"         # join field in source
    target_id_field = "UniqueRowNum"         # join field in target
    code_field = "CONCATENATE_PreciseLevellingSource"            # field in source with values like "1,1,2,0,2,1"
    converted_field = "PreciseLevellingSourceSummarized"       # field in target to store "1x0, 3x1, 2x2"

    # Step 1: Read and process data from the source feature class
    code_dict = {}

    with SearchCursor(dissolved_lines, [source_id_field, code_field]) as cursor:
        for row in cursor:
            key = row[0]
            code_str = row[1]
            if code_str:
                digits = [int(d) for d in code_str.split(',') if d.strip().isdigit()]
                count = Counter(digits)
                sorted_items = sorted(count.items())  # sort by digit
                output = ", ".join(f"{v}x{k}" for k, v in sorted_items)
                code_dict[key] = output

    # Step 2: Write transformed values to the target feature class
    with UpdateCursor(generalized_lines, [target_id_field, converted_field]) as cursor:
        for row in cursor:
            key = row[0]
            if key in code_dict:
                row[1] = code_dict[key]
                cursor.updateRow(row)

def  fill_missing_data_for_generalized_lines(generalized_lines) -> None:



    fields = ["UniqueLineNum", "LineCode","LineName","LineDirection","StartPointCode","EndPointCode","StartPointName","EndPointName", "numOfSegments","NumOfPoints"]

    # filling basic data for the generalized lines 
    with UpdateCursor(generalized_lines, fields) as update_cursor:
        for row in update_cursor:
            start_point_code = row[4]
            end_point_code = row[5]
            start_point_name = row[6]
            end_point_name = row[7]
            num_of_segments = row[8]

            if start_point_name and end_point_name:
                row[2] = f"{start_point_name} to {end_point_name}"
            
            if start_point_code and end_point_code:
                min_code = min(start_point_code, end_point_code)
                max_code = max(start_point_code, end_point_code)
                row[1] = f"{min_code}_{max_code}"
                if min_code == start_point_code:
                    row[3] = 1
                else:
                    row[3] = -1

            if num_of_segments:
                row[9] = num_of_segments + 1

            update_cursor.updateRow(row)

    # numerating the unique line codes
    codes_list = []
    line_code_to_unique_num = {}
    unique_num = 1

    with UpdateCursor(generalized_lines, ["LineCode", "UniqueLineNum"]) as cursor:
        for row in cursor:
            line_code = row[0]
            codes_list.append(line_code)
            if line_code not in line_code_to_unique_num:
                line_code_to_unique_num[line_code] = unique_num
                unique_num += 1
            row[1] = line_code_to_unique_num[line_code]
            cursor.updateRow(row)


    code_counts = Counter(codes_list)

    # marking duplicated lines
    with UpdateCursor(generalized_lines, ["LineCode","IsDuplicated"]) as cursor:
        for row in cursor:
            row[1] = 1 if code_counts[row[0]] > 1 else 0
            cursor.updateRow(row)



def fill_midpoints(generalized_lines, sorted_segments) -> None:

    """
    fillind the data for MidPointCodes and MidPointNames in the generalized lines feature class
    according to the data in the sorted_segments feature class
    
    """
    with UpdateCursor(generalized_lines, ["UniqueRowNum","MidPointCodes", "MidPointNames","NumOfPoints"]) as update_cursor:
        for generalized_row in update_cursor:
            line_number = generalized_row[0]
            mid_points_codes = []
            mid_points_names = []
            num_of_segments = count_features(SelectByAttribute(sorted_segments, where_clause=f"line_num = {line_number}", selection_type="NEW_SELECTION"))
            is_first_row = True
            with SearchCursor(sorted_segments, ["StartPointCode","StartPointName"], f"line_num = {line_number}") as search_cursor:
                for sorted_row in search_cursor:
                    if num_of_segments == 1:
                        pass  # No midpoints for single-segment lines
                    else:
                        if is_first_row:
                            is_first_row = False
                        else:
                            mid_points_codes.append(sorted_row[0])
                            mid_points_names.append(sorted_row[1])


            if mid_points_codes:
                generalized_row[1] = ", ".join(str(temp_code) for temp_code in mid_points_codes)

            if mid_points_names:
                generalized_row[2] = " to ".join(mid_points_names)

            generalized_row[3] = 2 + len(mid_points_codes)

            update_cursor.updateRow(generalized_row)





def fill_midpoints_and_dates(generalized_lines, sorted_segments) -> None:

    """
    fillind the data for MidPointCodes, MidPointNames, MeasDateStart, MeasDateEnd fields in the generalized lines feature class
    according to the data in the sorted_segments feature class
    
    """
    with UpdateCursor(generalized_lines, ["UniqueRowNum","MidPointCodes", "MidPointNames", "MeasDateStart", "MeasDateEnd","NumOfPoints"]) as update_cursor:
        for generalized_row in update_cursor:
            line_number = generalized_row[0]
            mid_points_codes = []
            mid_points_names = []
            dates = []
            num_of_segments = count_features(SelectByAttribute(sorted_segments, where_clause=f"line_num = {line_number}", selection_type="NEW_SELECTION"))
            is_first_row = True
            with SearchCursor(sorted_segments, ["StartPointCode","StartPointName","MeasDate"], f"line_num = {line_number}") as search_cursor:
                for sorted_row in search_cursor:
                    if num_of_segments == 1:
                        #temp_date = mmyy_int_to_yyyymm_str(sorted_row[2])
                        dates.append(sorted_row[2])
                    else:
                        if is_first_row:
                            is_first_row = False
                        else:
                            mid_points_codes.append(sorted_row[0])
                            mid_points_names.append(sorted_row[1])

                            #temp_date = mmyy_int_to_yyyymm_str(sorted_row[2])
                            dates.append(sorted_row[2])

            if mid_points_codes:
                generalized_row[1] = ", ".join(str(temp_code) for temp_code in mid_points_codes)

            if mid_points_names:
                generalized_row[2] = " to ".join(mid_points_names)

            generalized_row[5] = 2 + len(mid_points_codes)
            dates = sorted(dates)
            generalized_row[3] = dates[0]
            generalized_row[4] = dates[-1]


            update_cursor.updateRow(generalized_row)





def break_self_loop(sorted_segments, point_coordinates, line_num,start_point) -> int:
    """
    Recieving a set of segments that creating a self-loop and breaking it at one point so the would be loop turns to 2 lines
    The breakage point will be set according to 2 criterias: the participating point's ranks and their distance from the start point:
    If one (and only one) of the participating points (without the starting point) has rank 1 or 2, the breakage point will be that point
    otherwise, the breakage point will be the farthest point from the starting point


    Parameters
    ----------
    sorted_segments : FeatureLayer
        Feature layer of sorted segments
    point_coordinates : FeatureLayer
        Feature layer of control points
    line_num : int
        Line number of the generalized line
    start_point : int
        Starting point of the generalized line

    Returns
    -------
    breakage_point_code : int
        Code of the breakage point

    
    """
    breakage_point_code = None

    lines_segments = SelectByAttribute(sorted_segments, where_clause=f"line_num = {line_num}", selection_type="NEW_SELECTION")
    num_of_segments = count_features(lines_segments)
    participating_points = extract_coordinates_from_external_source(lines_segments, point_coordinates)
    points_dictionary = create_points_dictionary(participating_points)

    high_ranked_points = []
    farthest_point = None
    max_distance = 0
    start_point_data = points_dictionary[start_point]
    for point_code in points_dictionary:
        if point_code != start_point:
            point_data = points_dictionary[point_code]
            if point_data["HRank"] in [1, 2, 11, 12]:
                high_ranked_points.append(point_code)
            temp_dist = get_distance_between_points(start_point_data["X"], start_point_data["Y"], point_data["X"], point_data["Y"])
            if temp_dist > max_distance:
                max_distance = temp_dist
                farthest_point = point_code



    if len(high_ranked_points) ==1:
        breakage_point_code = high_ranked_points[0]
    else:
        breakage_point_code = farthest_point


   
    break_objectID = 0
    break_num = 0

    with SearchCursor(lines_segments, ["OBJECTID", "UniqueRowNum", "StartPointCode","EndPointCode"]) as cursor:
        for row in cursor:
            if row[3] == breakage_point_code:
                break_objectID = row[0]
                break_num = row[1]
                break

    
    increase_line_num = False
    SelectByAttribute(in_layer_or_view=sorted_segments,selection_type="CLEAR_SELECTION")
    with UpdateCursor(sorted_segments, ["OBJECTID", "UniqueRowNum","line_num"]) as cursor:
        for row in cursor:

            if increase_line_num:
                row[2] = row[2]+1
                cursor.updateRow(row)
            elif row[0] == break_objectID and row[1] == break_num:
                increase_line_num = True

    return breakage_point_code

def check_self_loops_new(sorted_segments, point_coordinates) -> None:
    """
    Checking if the sorted segments data will have self-loops and breaking them at one point so the would be loop turns to 2 lines

    Parameters
    ----------
    sorted_segments : FeatureLayer
        Feature layer of sorted segments
    point_coordinates : FeatureLayer
        Feature layer of control points

    
    """
    num_of_lines = 0
    num_of_self_loops = 0
    points_dictionary = create_points_dictionary(point_coordinates)
    self_looping_points = []
    breakage_points = []
    self_looping_lines = []




    temp_dissolve = Dissolve(in_features=sorted_segments,out_feature_class=r"memory\temp_dissolve",statistics_fields="StartPointCode FIRST;EndPointCode LAST",dissolve_field="line_num",concatenation_separator="")
    num_of_lines = count_features(temp_dissolve)
    #self_looping_lines_selection = SelectByAttribute(in_layer_or_view=temp_dissolve,selection_type="NEW_SELECTION",where_clause="FIRST_StartPointCode = LAST_EndPointCode",invert_where_clause=None)
    #num_of_self_loops = count_features(self_looping_lines_selection)

    num_of_fixed_lines = 0
    with SearchCursor(temp_dissolve, ["line_num", "FIRST_StartPointCode","LAST_EndPointCode"]) as cursor:
        for row in cursor:
            if row[1] == row[2]:
                self_looping_lines.append(row[0]+num_of_fixed_lines)
                start_point = row[1]
                self_looping_points.append(start_point)
                #AddMessage(f"line {row[0]+num_of_fixed_lines} is self-looping")
                temp_point = break_self_loop(sorted_segments, point_coordinates, row[0]+num_of_fixed_lines ,start_point)
                breakage_points.append(temp_point)
                num_of_fixed_lines += 1

    '''
    AddMessage(f"The initial data has {num_of_lines} lines and {num_of_self_loops} self-looping lines")
    if num_of_self_loops > 0:
        AddMessage(f"Breaking self-loops...")
        num_of_fixed_lines = 0
        with SearchCursor(self_looping_lines_selection, ["line_num", "FIRST_StartPointCode"]) as cursor:
            for row in cursor:
                self_looping_lines.append(row[0]+num_of_fixed_lines)
                start_point = row[1]
                self_looping_points.append(start_point)
                AddMessage(f"line {row[0]+num_of_fixed_lines} is self-looping")
                temp_point = break_self_loop(sorted_segments, point_coordinates, row[0]+num_of_fixed_lines ,start_point)
                breakage_points.append(temp_point)
                num_of_fixed_lines += 1

    ''' 

    if num_of_fixed_lines > 0:
        AddMessage(f"Found {num_of_fixed_lines} self-looping lines from the initial {num_of_lines}. Breaking them at one point each.")
        for i in range(num_of_fixed_lines):
            #current_line = self_looping_lines[i] + num_of_fixed_lines
            self_loop_point_code = self_looping_points[i]
            self_loop_point_data = points_dictionary[self_loop_point_code]
            breake_point_code = breakage_points[i]
            breake_point_data = points_dictionary[breake_point_code]
            AddMessage(f"   The self loop around point {self_loop_point_data['Name']} ({self_loop_point_code}) was broken at point {breake_point_data['Name']} ({breake_point_code})")
        AddMessage(f"✅ Successfully broke {num_of_fixed_lines} self-looping lines. Final sorted data has {num_of_lines+num_of_fixed_lines} lines")
    else:
        AddMessage(f"✅ No self-looping lines found in the initial data. Final sorted data has {num_of_lines} lines")
    #SelectByAttribute(in_layer_or_view=sorted_segments,selection_type="CLEAR_SELECTION")
    Delete(temp_dissolve)
    del points_dictionary   




def check_self_loops(sorted_segments, point_coordinates) -> None:
    """
    Checking if the sorted segments data will have self-loops and breaking them at one point so the would be loop turns to 2 lines

    Parameters
    ----------
    sorted_segments : FeatureLayer
        Feature layer of sorted segments
    point_coordinates : FeatureLayer
        Feature layer of control points

    
    """
    num_of_lines = 0
    num_of_self_loops = 0
    points_dictionary = create_points_dictionary(point_coordinates)
    self_looping_points = []
    breakage_points = []
    last_segment = get_last_segment(sorted_segments)

    with SearchCursor(last_segment, ["line_num"]) as cursor:
        row = next(cursor, None)
        if row:
            num_of_lines = row[0]

    for i in range(num_of_lines):
        lines_segments = SelectByAttribute(sorted_segments, where_clause=f"line_num = {i+1}", selection_type="NEW_SELECTION")
        num_of_segments = count_features(lines_segments)
        if num_of_segments > 1:
            start_point = 0
            end_point = 0
            ind = 0
            with SearchCursor(lines_segments, ["StartPointCode", "EndPointCode"]) as cursor:
                for row in cursor:
                    ind += 1
                    if ind == 1:
                        start_point = row[0]
                    elif ind == num_of_segments:
                        end_point = row[1]
                    
            if start_point == end_point:
                num_of_self_loops += 1
                self_looping_points.append(start_point)
                AddMessage(f"line {i+1} is self-looping")
                temp_point = break_self_loop(sorted_segments, point_coordinates, i+1 ,start_point)
                breakage_points.append(temp_point)
                num_of_lines += 1

    if num_of_self_loops>0:
        AddMessage(f"Found {num_of_self_loops} self-looping set of segments:")
        for i in range(num_of_self_loops):
            self_loop_point_code = self_looping_points[i]
            self_loop_point_data = points_dictionary[self_loop_point_code]
            breake_point_code = breakage_points[i]
            breake_point_data = points_dictionary[breake_point_code]
            AddMessage(f"   The self loop around point {self_loop_point_data['Name']} ({self_loop_point_code}) was broken at point {breake_point_data['Name']} ({breake_point_code})")

    SelectByAttribute(in_layer_or_view=sorted_segments,selection_type="CLEAR_SELECTION")
    Delete(last_segment)
    del points_dictionary



def classify_points(line_segments, point_coordinates, classified_points_output = None, print_messages = True) -> tuple:
    """
    Classify points into junction points, regular points, and dead-end points using codeA and codeB, and save to output feature class.
    
    Parameters:
    - line_segments (FeatureLayer): Input feature class of line segments
    - point_coordinates (FeatureLayer): Input feature class of control points
    - classified_points_output (str): Output feature class for classified points (optional)
    
    Returns:
    - tuple: Tuple containing the codes of junction points, regular points, and dead-end points
    """ 
    point_counts = {}

    if print_messages:
        AddMessage("Classifying points...")

    # Count occurrences of points in StartPointCode and EndPointCode
    with SearchCursor(line_segments, ['StartPointCode', 'EndPointCode']) as cursor:
        for row in cursor:
            point_a = row[0]
            point_b = row[1]

            # Increment counts for both points
            point_counts[point_a] = point_counts.get(point_a, 0) + 1
            point_counts[point_b] = point_counts.get(point_b, 0) + 1

    # Get the spatial reference from the input segments feature class
    sr = Describe(line_segments).spatialReference

    
    # Check if the field SegmentCode exists in the line_segments feature class
    field_names = [field.name for field in ListFields(line_segments)]
    if "SegmentCode" in field_names:
        # Setting the counting of dead end points from duplicated segments to 1
        unique_segments = Statistics(in_table=line_segments,out_table=r"memory\unique_segments",statistics_fields="StartPointCode FIRST;EndPointCode FIRST",case_field="SegmentCode")
    elif "LineCode" in field_names:
        unique_segments = Statistics(in_table=line_segments,out_table=r"memory\unique_segments",statistics_fields="StartPointCode FIRST;EndPointCode FIRST",case_field="LineCode")
    else:
        AddError("The fields 'LineCode' or 'SegmentCode' does not exist in the line_segments feature class.")
        raise ValueError("The field 'SegmentCode' or 'LineCode' is missing.")
    
    all_checked_codes = []
    with SearchCursor(unique_segments, ['FIRST_StartPointCode', 'FIRST_EndPointCode']) as cursor:
        for start_code, end_code in cursor:
            all_checked_codes.append(start_code)
            all_checked_codes.append(end_code)
    
    code_counts = Counter(all_checked_codes)
    updated_dead_ends = [code for code, count in code_counts.items() if count == 1]

    #for key, count in code_counts.items():
    #    AddMessage(f"{key}: {count}")

    for code in updated_dead_ends:
        point_counts[code] = 1  # Set count to 1 for dead-end points from duplicated segments


    Delete(unique_segments)


    if classified_points_output:
        # Separate the path and name from the classified_points_output parameter
        out_path, out_name = os.path.split(classified_points_output)

        if not out_path:
            out_path = arcpy.env.workspace  # Use the current workspace if no path is provided


    else:
        out_path = "memory"
        out_name = "classified_points"

    # Create the output points feature class
    CreateFeatureclass(out_path=out_path, out_name=out_name, geometry_type="POINT", spatial_reference=sr)

    

    # Add fields for code, name, num, and type
    AddField(os.path.join(out_path, out_name), "Code", "LONG")
    AddField(os.path.join(out_path, out_name), "Name", "TEXT")
    AddField(os.path.join(out_path, out_name), "Count", "LONG")  # Number of appearances
    AddField(os.path.join(out_path, out_name), "Type", "TEXT")  # Point type (junction, regular, dead-end)


    # Cache point coordinates and names in a dictionary for efficient lookup

    point_cache = create_points_dictionary(point_coordinates, ["SHAPE@XY", "Code", "Name"])
    """
    point_cache = {}
    with SearchCursor(point_coordinates, ['SHAPE@XY',  "Code", "Name"]) as cursor:
        for row in cursor:
            xy = row[0]
            code = row[1]
            full_name = row[2]
            point_cache[code] = (xy, full_name)

    """
    junction_points_codes = []
    regular_points_codes = []
    dead_end_points_codes = []

    #ws = get_feature_layer_workspace(out_path)
    #editor = start_editing(ws)
    with InsertCursor(os.path.join(out_path, out_name), ['SHAPE@XY', 'Code', 'Name', 'Count', 'Type']) as insert_cursor:
        for point_code, point_count in point_counts.items():
            #xy, name = point_cache[point_code]
            xy, name = point_cache[point_code]["SHAPE@XY"], point_cache[point_code]["Name"]        
            if point_count == 1:
                point_type = "dead-end"
                dead_end_points_codes.append(point_code)
            elif point_count == 2:
                point_type = "regular"
                regular_points_codes.append(point_code)
            elif point_count >= 3:
                point_type = "junction"
                junction_points_codes.append(point_code)
            else:
                point_type = "unknown"
                AddError(f"Point {point_code} was not found in the segments data")
                    
            # Insert the classified point with XY geometry
            insert_cursor.insertRow([xy, point_code, name, point_count, point_type])
                
    #stop_editing(editor)

    if print_messages:
        AddMessage(f"Found {len(junction_points_codes)} junction points, {len(regular_points_codes)} regular points and {len(dead_end_points_codes)} dead-end points.")

    if not classified_points_output:
        Delete(os.path.join(out_path, out_name))
        del out_path, out_name
    return junction_points_codes,regular_points_codes,dead_end_points_codes

def sort_segments_matlab(segments, point_coordinates, junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output=None, print_messages=True) -> None:
    """
    Sort the segments between junction points and create output feature class.

    Parameters:
    - segments (FeatureLayer): Input feature class of line segments
    - point_coordinates (list): List of point coordinates (used for self-loop check)
    - junction_points_codes (list): List of junction point codes
    - regular_points_codes (list): List of regular point codes
    - dead_end_points_codes (list): List of dead-end point codes
    - sorted_segments_output (str): Optional path to output feature class
    - print_messages (bool): Whether to print progress messages

    Returns:
    - FeatureLayer: Output feature class of sorted segments
    """

    if print_messages:
        AddMessage("Sorting segments...")

    junctions_and_dead_ends_codes = junction_points_codes + dead_end_points_codes



    # Output feature class
    if sorted_segments_output:
        out_path, out_name = os.path.split(sorted_segments_output)
        if not out_path:
            out_path = arcpy.env.workspace
    else:
        out_path = "memory"
        out_name = "sorted_segments"
        sorted_segments_output = os.path.join(out_path, out_name)

    # Create a copy of the segments feature class to avoid modifying the original data

    segments_copy = CopyFeatures(segments, r"memory\segments_copy")
    AddField(segments_copy, "line_num", "LONG")
    AddField(segments_copy, "new_order", "LONG")
    AddField(segments_copy, "reverse_flag", "SHORT")


    # === Step 1: Extract segment endpoints ===
    #segments_list = []
    segment_endpoints = []
    with SearchCursor(segments_copy, ["OID@", "StartPointCode", "EndPointCode"]) as cursor:
        for row in cursor:
            #segments_list.append(row[0])
            segment_endpoints.append([row[1], row[2]])

    # === Step 2: Call MATLAB ===
    import matlab.engine
    eng = matlab.engine.start_matlab()

    try:
        AddMessage("Calling MATLAB function to sort segments...")
        # Add the path to MATLAB
        eng.addpath(os.path.abspath(MATLAB_FUNCTIONS_DIR), nargout=0)
        segment_endpoints_ml = matlab.double(segment_endpoints)
        junctions_ml = matlab.double(junctions_and_dead_ends_codes)
        regulars_ml = matlab.double(regular_points_codes)
        '''
        vec1 = matlab.double(junctions_and_dead_ends_codes)
        vec2 = matlab.double(regular_points_codes)

        

        sz1, rows1, cols1 = eng.check_vector_size(vec1, nargout=3)
        sz2, rows2, cols2 = eng.check_vector_size(vec2, nargout=3)

        AddMessage(f"vec1 -> size: {sz1}, rows: {rows1}, cols:, {cols1}")
        AddMessage(f"vec2 -> size: {sz2}, rows: {rows2}, cols:, {cols2}")
        '''

        sorted_info = eng.sortSegments(segment_endpoints_ml, junctions_ml, regulars_ml)

    except Exception as e:
        AddError(f"MATLAB error: {e}")
        raise
    finally:
        # Stop MATLAB engine
        eng.quit()
        AddMessage("MATLAB function call completed.")



    # === Step 3: Process sorted segments ===
    #sorted_indices = list(sorted_info["sorted_index"][0])
    #segment_indices = list(sorted_info["segment_index"][0])
    sorted_indices = list(sorted_info["sorted_index"][0])
    #sorted_indices = [x - 1 for x in sorted_info["sorted_index"][0]] # Convert to 0-based indexes list
    segment_indices = [x - 1 for x in sorted_info["segment_index"][0]]  # Convert to 0-based indexes list

    line_numbers    = list(sorted_info["line_number"][0])
    reverse_flags   = list(sorted_info["reverse_flag"][0])


    unmarked_segments = []
    with UpdateCursor(segments_copy, ["OBJECTID", "line_num", "new_order", "reverse_flag"]) as cursor:
        idx = 0
        for row in cursor:
            if idx == segment_indices[idx-len(unmarked_segments)]:
                row[1] = line_numbers[idx-len(unmarked_segments)]
                row[2] = sorted_indices[idx-len(unmarked_segments)]
                row[3] = reverse_flags[idx-len(unmarked_segments)]
                cursor.updateRow(row)
            else:
                unmarked_segments.append(row[0])  # Store unmarked segment OIDs
            idx += 1

    # dealing with cases where some segments were not marked with line numbers due to the self looping of the segments in the input
    # when all points in the loop classified as regular
    if unmarked_segments:
        unmarked_segments_str = ", ".join(map(str, unmarked_segments))
        AddMessage(f"The segments with the following OIDs were not marked with line numbers in the sorting process: {unmarked_segments_str}")
        AddMessage(f"Each segment will be marked as a separate line")
        last_line_num = max(line_numbers)
        last_new_order = max(sorted_indices)
        with UpdateCursor(segments_copy, ["OBJECTID", "line_num","new_order"],where_clause=f"OBJECTID IN ({unmarked_segments_str})") as cursor:
            for row in cursor:
                row[1] = last_line_num + 1
                row[2] = last_new_order + 1
                last_line_num += 1
                last_new_order += 1
                cursor.updateRow(row)



        '''
        participating_segment_codes = []
        with SearchCursor(segments_copy, ["OBJECTID", "SegmentCode"], where_clause=f"OBJECTID IN ({unmarked_segments_str})") as cursor:
            for row in cursor:
                participating_segment_codes.append(row[1])
        participating_segment_codes = list(set(participating_segment_codes))  # Remove duplicates
        segment_code_to_line_num = {}
        for segment_code in participating_segment_codes:
            last_line_num += 1
            segment_code_to_line_num[segment_code] = last_line_num
        with UpdateCursor(segments_copy, ["OBJECTID", "SegmentCode", "line_num"],where_clause=f"OBJECTID IN ({unmarked_segments_str})") as cursor:
            for row in cursor:
                if row[1] in segment_code_to_line_num:
                    row[2] = segment_code_to_line_num[row[1]]
                    cursor.updateRow(row)
        '''
    # === Step 4: Sort segments ===

    sorted_segments = Sort(in_dataset=segments_copy, out_dataset=sorted_segments_output,sort_field="new_order ASCENDING",spatial_sort_method="UR")
    sorted_segments_to_reverse = SelectByAttribute(in_layer_or_view=sorted_segments, where_clause="reverse_flag = 1", selection_type="NEW_SELECTION")
    reverse_segments_new(sorted_segments_to_reverse)
    SelectByAttribute(in_layer_or_view=sorted_segments, selection_type="CLEAR_SELECTION")

    # Get count of original segments
    count_original_segments = count_features(segments)
    # Get count of sorted segments
    count_sorted_segments = count_features(sorted_segments)

    if count_sorted_segments != count_original_segments and print_messages:
        AddWarning(f"{count_original_segments - count_sorted_segments} unsorted segments remain.")

    AddMessage(f"checking for self-loops in the sorted segments...")
    check_self_loops_new(sorted_segments, point_coordinates)


    Delete(segments_copy)
    DeleteField(in_table=sorted_segments,drop_field="new_order;reverse_flag;ORIG_FID",method="DELETE_FIELDS")

    AddMessage("Sorting process completed.")
    #AddMessage(f"The {count_sorted_segments} segments have been sorted into {int(max(line_numbers))} lines.")
    return sorted_segments


def sort_segments(segments, point_coordinates,junction_points_codes, regular_points_codes, dead_end_points_codes, sorted_segments_output = None,print_messages=True) -> FeatureLayer:
    """
    Sort the segments between junction points and create output feature class.

    Parameters:
    - segments (FeatureLayer): Input feature class of line segments
    - junction_points_codes (list): List of junction point codes
    - regular_points_codes (list): List of regular point codes
    - dead_end_points_codes (list): List of dead-end point codes

    Returns:
    - FeatureLayer: Output feature class of sorted segments with additional field "line_num" to represent the number of the generalized line
    """
    if print_messages:
        AddMessage("Sorting segments...")

    # Combine the lists of junction points and dead-end points since in both cases the search should be stopped
    junctions_and_dead_ends_codes = junction_points_codes + dead_end_points_codes
    del junction_points_codes, dead_end_points_codes

    # Create a copy of the segments feature class to avoid modifying the original data
    # arcpy.env.preserveGlobalIds = True
    segments_input_copy = CopyFeatures(segments, r"memory\segments_input_copy")

    # Get the spatial reference from the input segments feature class
    sr = Describe(segments).spatialReference

    if sorted_segments_output:
        # Separate the path and name from the sorted_segments_output parameter
        out_path, out_name = os.path.split(sorted_segments_output)

        if not out_path:
            out_path = arcpy.env.workspace  # Use the current workspace if no path is provided


    else:
        out_path = "memory"
        out_name = "sorted_segments"

    ## Separate the path and name from the classified_points_output parameter
    #out_path, out_name = os.path.split(sorted_segments_output)

    #if not out_path:
    #    out_path = arcpy.env.workspace  # Use the current workspace if no path is provided


    sorted_segments = CreateFeatureclass(out_path=out_path, out_name=out_name,geometry_type="POLYLINE",template=segments,spatial_reference=sr)

    # Add fields for storing line number
    AddField(sorted_segments, "line_num", "SHORT")

    line_number = 0
    for temp_code in junctions_and_dead_ends_codes:
        #_#_# AddMessage(f"temp code: {temp_code}")
        point_finished = False
        while not point_finished:
            last_point_code = temp_code
            line_finished = False
            line_number += 1
            #_#_# AddMessage(f"starting line number: {line_number}")

            while not line_finished:
                temp_segment = get_first_segment_with_point(segments_input_copy, last_point_code)
                
                if temp_segment:
                    prev_point_code = last_point_code
                    last_point_code = get_last_segment_end_point_code(temp_segment)
                    #_#_# AddMessage(f"Pulled segment with points {get_segments_end_point_names(temp_segment)}")
                    #_#_# AddMessage(f"Pulled segment with points {get_segments_end_point_codes(temp_segment)}")
                    #_#_# AddMessage(f"prev point code: {prev_point_code}, last point code: {last_point_code}")

                    if prev_point_code == last_point_code:
                        #_#_# AddMessage(f"points are equal, need to reverse")
                        reversed_temp_segment = reverse_line_segments(temp_segment)
                        add_segment(sorted_segments, reversed_temp_segment,line_number)
                        last_point_code = get_last_segment_end_point_code(reversed_temp_segment)
                        #_#_# AddMessage(f"reversed segment points {get_segments_end_point_names(reversed_temp_segment)}")
                        #_#_# AddMessage(f"reversed segment codes {get_segments_end_point_codes(reversed_temp_segment)}")
                        #_#_# AddMessage(f"prev point code: {prev_point_code}, last point code: {last_point_code}")
                        Delete(reversed_temp_segment)
                    else:
                        add_segment(sorted_segments, temp_segment,line_number)

                    delete_segment(segments_input_copy, temp_segment)

                    if last_point_code in junctions_and_dead_ends_codes:
                        line_finished = True
                        if print_messages:
                            AddMessage(f"Line {line_number} finished")

                    elif last_point_code in regular_points_codes:
                        # Normally nothing should happen here
                        #_#_# AddMessage(f"Line {line_number} continues")
                        pass  # This is a placeholder indicating that nothing needs to be done
                    else:
                        #_#_# AddMessage(f"Point {last_point_code} was not found in the segments data")
                        pass

                else:
                    #_#_# AddMessage(f"Point {temp_code} not found in the remaining segments data")
                    line_finished = True
                    point_finished = True 
                    line_number -= 1


    if (int(GetCount(segments_input_copy)[0]) != 0):
        if print_messages:
            AddWarning(f"There are still {int(GetCount(segments_input_copy)[0])} segments left in the segments data")
    Delete(segments_input_copy)
    del junctions_and_dead_ends_codes,regular_points_codes
    # arcpy.env.preserveGlobalIds = False
    check_self_loops(sorted_segments, point_coordinates)
    return sorted_segments


def generalize_sorted_segments(sorted_segments, generalized_lines_output=None, print_messages=True) -> FeatureLayer:



    if print_messages:
        AddMessage("Generalizing sorted segments...")

    dissolved_lines = Dissolve(in_features=sorted_segments,out_feature_class=r"memory\dissolved_lines",dissolve_field="line_num",
                    #statistics_fields="A FIRST;B LAST;dH SUM;dist SUM;n SUM;diff SUM;codeA FIRST;typeA FIRST;nameA FIRST;groupA FIRST;rankA FIRST;codeB LAST;typeB LAST;nameB LAST;groupB LAST;rankB LAST;Y1 FIRST;X1 FIRST;H1 FIRST;Y2 LAST;X2 LAST;H2 LAST;line_num COUNT",

                    statistics_fields="StartPointCode FIRST;StartPointName FIRST;EndPointCode LAST;EndPointName LAST;MeasHeightDiff SUM;MeasDist SUM;SetUpsNum SUM;MiscBF SUM;GravCorrection SUM;CorrHeightDiff SUM;line_num COUNT;PreciseLevellingSource CONCATENATE;Adj2025 CONCATENATE; StartPointCode CONCATENATE;MeasDate MIN;MeasDate MAX",

                    multi_part="MULTI_PART",unsplit_lines="DISSOLVE_LINES",concatenation_separator=",")

    sorted_dissolved_lines = Sort(in_dataset=dissolved_lines,out_dataset=r"memory\sorted_dissolved_lines",sort_field="line_num ASCENDING",spatial_sort_method="UR")
    

    AddMessage("Dissolved lines created.")
    Delete(dissolved_lines)

    # Get the spatial reference from the input segments feature class
    sr = Describe(sorted_segments).spatialReference

    if generalized_lines_output:
    # Separate the path and name from the sorted_segments_output parameter
        out_path, out_name = os.path.split(generalized_lines_output)

        if not out_path:
            out_path = arcpy.env.workspace  # Use the current workspace if no path is provided


    else:
        out_path = "memory"
        out_name = "generalized_lines"
    
    # Separate the path and name from the classified_points_output parameter
    #out_path, out_name = os.path.split(generalized_lines_output)
    #if not out_path:
    #    out_path = arcpy.env.workspace  # Use the current workspace if no path is provided

    # Create the output points feature class
    generalized_lines = CreateFeatureclass(out_path=out_path, out_name=out_name, geometry_type="POLYLINE",template=GENERALIZED_LINES_TEMPLATE, spatial_reference=sr)

    #DeleteField(in_table=generalized_lines,drop_field="date;calc_file;text_file;group_;Original;month;year;line_num",method="DELETE_FIELDS")
    #DeleteField(in_table=generalized_lines,drop_field="MeasYear;MeasMonth;MeasDate;MeasFile;OriginalMeasRow;SubTextFile;line_num",method="DELETE_FIELDS")
    
    #AddField(in_table=generalized_lines,field_name="num_of_segments",field_type="SHORT")

    field_mapping = fr'StartPointCode "קוד נקודת התחלה" true true false 4 Long 0 0,First,#,{sorted_dissolved_lines},FIRST_StartPointCode,-1,-1;'+\
                    fr'EndPointCode "קוד נקודת סיום" true true false 4 Long 0 0,First,#,{sorted_dissolved_lines},LAST_EndPointCode,-1,-1;'+\
                    fr'StartPointName "שם נקודת התחלה" true true false 255 Text 0 0,First,#,{sorted_dissolved_lines},FIRST_StartPointName,0,7;'+\
	                fr'EndPointName "שם נקודת סיום" true true false 255 Text 0 0,First,#,{sorted_dissolved_lines},LAST_EndPointName,0,7;'+\
                    fr'MeasHeightDiff "הפרש גובה מדוד" true true false 8 Double 0 0,First,#,{sorted_dissolved_lines},SUM_MeasHeightDiff,-1,-1;'+\
                    fr'MeasDist "אורך מדידה" true true false 2 Short 0 0,First,#,{sorted_dissolved_lines},SUM_MeasDist,-1,-1;'+\
                    fr'SetUpsNum "מספר תחנות" true true false 2 Short 0 0,First,#,{sorted_dissolved_lines},SUM_SetUpsNum,-1,-1;'+\
                    fr'MiscBF "אי סגירה בין מהלכים במ[double_quote]מ" true true false 8 Double 0 0,First,#,{sorted_dissolved_lines},SUM_MiscBF,-1,-1;'+\
                    fr'GravCorrection "תיקון גרבימטרי" true true false 8 Double 0 0,First,#,{sorted_dissolved_lines},SUM_GravCorrection,-1,-1;'+\
                    fr'CorrHeightDiff "הפרש גובה מתוקן" true true false 8 Double 0 0,First,#,{sorted_dissolved_lines},SUM_CorrHeightDiff,-1,-1;'+\
                    fr'UniqueRowNum "מספר שורה ייחודי" true true false 4 Long 0 0,First,#,{sorted_dissolved_lines},line_num,-1,-1;'+\
                    fr'MeasDateStart "תאריך תחילת מדידה" true true false 255 Text 0 0,First,#,{sorted_dissolved_lines},MIN_MeasDate,0,7;'+\
                    fr'MeasDateEnd "תאריך סיום מדידה" true true false 255 Text 0 0,First,#,{sorted_dissolved_lines},MAX_MeasDate,0,7;'+\
                    fr'NumOfSegments "מספר קטעים" true true false 2 Short 0 0,First,#,{sorted_dissolved_lines},COUNT_line_num,-1,-1'
    

                    


    
    Append(inputs=sorted_dissolved_lines,target=generalized_lines,schema_type="NO_TEST",
        field_mapping=field_mapping,  subtype="",expression="",match_fields=None,update_geometry="NOT_UPDATE_GEOMETRY")
    


    AddMessage("Lines were appended to the generalized lines feature class.")
    #fill_precise_levelling_source_summary(generalized_lines, sorted_dissolved_lines)
    fill_midpoints_and_adj2025(generalized_lines, sorted_dissolved_lines, sorted_segments)
    AddMessage("Midpoints and adj2025 classification were filled in the generalized lines feature class.")

    Delete(in_data=sorted_dissolved_lines)

    



    EnableEditorTracking(in_dataset=generalized_lines,creator_field="created_user",creation_date_field="created_date",
                        last_editor_field="last_edited_user",last_edit_date_field="last_edited_date",add_fields="NO_ADD_FIELDS",record_dates_in="UTC")



    #fill_midpoints(generalized_lines, sorted_segments)

    fill_missing_data_for_generalized_lines(generalized_lines)
    AddMessage("Missing data for generalized lines was filled.")



    return generalized_lines