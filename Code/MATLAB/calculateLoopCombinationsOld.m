function [resultsCell] = calculateLoopCombinations(levellingData)
    % calculateLoopCombinations: Analyzes levelling loop data with duplicate lines.
    %
    % Syntax:
    %   resultsCell = calculateLoopCombinations(levellingData)
    %
    % Input:
    %   levellingData - An N-by-5 matrix where each row represents a
    %                   measured levelling line. The data should represent
    %                   a single, closed loop.
    %                   Columns:
    %                   1: Line Number (ID) - Can be any numeric identifier.
    %                   2: Start Point Code (numeric)
    %                   3: End Point Code (numeric)
    %                   4: Distance
    %                   5: Measured Height Difference (dH from Start to End)
    %
    % Output:
    %   resultsCell - A cell array where each cell contains a struct with
    %                 information for one possible loop combination. This format
    %                 is compatible with the MATLAB Engine for Python.
    %                 Each struct has the fields:
    %                   .combination     - A K-by-5 matrix of the specific lines.
    %                   .misclosure      - The calculated misclosure for the loop.
    %                   .totalDistance   - The total distance of the loop.
    %
    % Description:
    %   This function takes levelling data for a single closed loop that may
    %   contain multiple measurements for the same line segment. It first
    %   identifies all unique line segments and groups their corresponding
    %   measurements. It then traces the loop's path to establish a
    %   consistent direction for calculation. Finally, it generates all
    %   possible unique loop combinations and, for each one, calculates the
    %   total misclosure and the total distance.

    %% --- Input Validation ---
    if nargin < 1 || isempty(levellingData)
        error('Input levellingData is required.');
    end
    if size(levellingData, 2) ~= 5
        error('Input levellingData must have 5 columns.');
    end

    %% --- Step 1: Group Measurements by Unique Segment ---
    % A segment is defined by its two endpoints, regardless of order.
    % We use a containers.Map to group all measurements for each segment.
    segmentMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    allPoints = unique([levellingData(:, 2); levellingData(:, 3)]);

    for i = 1:size(levellingData, 1)
        line = levellingData(i, :);
        p1 = line(2);
        p2 = line(3);

        % Create a canonical key for the segment (e.g., '101-105') by
        % ordering the point codes. This ensures (P1, P2) and (P2, P1)
        % map to the same segment.
        segmentKey = sprintf('%d-%d', min(p1, p2), max(p1, p2));

        if isKey(segmentMap, segmentKey)
            % Append this measurement to the existing list for this segment
            segmentMap(segmentKey) = [segmentMap(segmentKey); {line}];
        else
            % Create a new list for this new segment
            segmentMap(segmentKey) = {line};
        end
    end

    %% --- Step 2: Find the Ordered Path of the Loop ---
    % We build an adjacency list to represent the network as a graph and
    % traverse it to find the ordered sequence of points in the loop.
    adj = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for i = 1:numel(allPoints)
        adj(allPoints(i)) = [];
    end

    segmentKeys = keys(segmentMap);
    for i = 1:numel(segmentKeys)
        points = sscanf(segmentKeys{i}, '%d-%d');
        p1 = points(1);
        p2 = points(2);
        adj(p1) = [adj(p1), p2];
        adj(p2) = [adj(p2), p1];
    end

    % Traverse the graph to get the ordered loop path.
    % This assumes a single, simple loop structure where each node has
    % a degree of 2.
    startNode = allPoints(1);
    path = [startNode];
    prevNode = -1; % Placeholder for the previous node in the traversal
    currNode = startNode;

    for i = 1:(numel(allPoints) - 1)
        neighbors = adj(currNode);
        % Find the next node in the path that is not the one we just came from.
        if numel(neighbors) == 1 && i < (numel(allPoints) - 1)
             error('Loop is not closed or data is invalid. Point %d is a dead end.', currNode);
        end

        if neighbors(1) ~= prevNode
            nextNode = neighbors(1);
        else
            nextNode = neighbors(2);
        end
        
        path = [path, nextNode];
        prevNode = currNode;
        currNode = nextNode;
    end
    
    % Final check to ensure the discovered path is a closed loop
    neighborsOfLast = adj(path(end));
    if ~any(neighborsOfLast == startNode)
        error('Could not determine a single closed loop from the data provided.');
    end


    %% --- Step 3: Prepare for Combination Generation ---
    % Get the measurement groups for each segment in the ordered path
    numSegments = numel(path);
    orderedSegmentGroups = cell(1, numSegments);
    groupSizes = zeros(1, numSegments);

    for i = 1:numSegments
        p1 = path(i);
        % The last segment connects the last point back to the start point
        if i < numSegments
            p2 = path(i+1);
        else
            p2 = path(1);
        end

        segmentKey = sprintf('%d-%d', min(p1, p2), max(p1, p2));
        if ~isKey(segmentMap, segmentKey)
            error('Internal logic error: Path contains a segment with no measurements: %d-%d', p1, p2);
        end
        orderedSegmentGroups{i} = segmentMap(segmentKey);
        groupSizes(i) = numel(orderedSegmentGroups{i});
    end

    %% --- Step 4: Generate Combinations and Calculate Results ---
    % We use ndgrid to efficiently generate all combinations of indices.
    % The total number of combinations is the product of the number of
    % measurements for each segment.
    numCombinations = prod(groupSizes);
    
    % Create cell array of vectors, e.g., {1:2, 1:1, 1:3} for ndgrid
    gridVectors = arrayfun(@(n) 1:n, groupSizes, 'UniformOutput', false);
    
    % Generate index grids. Each grid corresponds to a segment.
    indexGrids = cell(1, numSegments);
    [indexGrids{:}] = ndgrid(gridVectors{:});

    % Flatten the index grids to get a list of index combinations
    combinationIndices = zeros(numCombinations, numSegments);
    for i = 1:numSegments
        combinationIndices(:, i) = indexGrids{i}(:);
    end

    % Pre-allocate the results struct for efficiency
    results = struct('combination', cell(numCombinations, 1), ...
                     'misclosure', cell(numCombinations, 1), ...
                     'totalDistance', cell(numCombinations, 1));

    % Iterate over each combination to perform the calculations
    for i = 1:numCombinations
        currentCombinationLines = cell(1, numSegments);
        totalDistance = 0;
        misclosure = 0;

        % Build the combination and calculate sums
        for j = 1:numSegments
            % Get the specific measurement for this segment and combination
            group = orderedSegmentGroups{j};
            lineIndex = combinationIndices(i, j);
            lineData = group{lineIndex};
            
            currentCombinationLines{j} = lineData;
            
            % Sum the distance (always positive)
            totalDistance = totalDistance + lineData(4);
            
            % Sum the height difference, carefully considering direction.
            % 'path_p1' is the start of the segment in our traced path.
            path_p1 = path(j); 
            
            if lineData(2) == path_p1 
                % Measurement direction matches the path traversal direction
                misclosure = misclosure + lineData(5);
            else 
                % Measurement direction is opposite to the path direction
                misclosure = misclosure - lineData(5);
            end
        end
        
        % Store the results for this combination
        results(i).combination = vertcat(currentCombinationLines{:});
        results(i).misclosure = misclosure;
        results(i).totalDistance = totalDistance;
    end
    
    %% --- Step 5: Package for Python Engine ---
    % The MATLAB engine for Python cannot return a non-scalar struct array.
    % We convert the struct array into a cell array. The engine will convert
    % this cell array into a Python list of dictionaries.
    if numCombinations > 0
        resultsCell = num2cell(results);
    else
        resultsCell = {}; % Return an empty cell if no combinations found
    end
end
