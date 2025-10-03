function [resultsCell] = calculateLoopCombinations(levellingData)
    % calculateLoopCombinations: Analyzes levelling loop data with duplicate lines.
    %
    % Syntax:
    %   resultsCell = calculateLoopCombinations(levellingData)
    %
    % Input:
    %   levellingData - An N-by-9 matrix where each row represents a
    %                   measured levelling line. The data should represent
    %                   a single, closed loop.
    %                   Columns:
    %                   1: Line Number (ID)
    %                   2: Start Point Code (numeric)
    %                   3: End Point Code (numeric)
    %                   4: Measured Height Difference (dH from Start to End)
    %                   5: Distance
    %                   6: Number of Setups
    %                   7: Back/Forward Misclosure per line
    %                   8: Gravimetric Correction per line
    %                   9: Corrected Height Difference (Measured + Correction)
    %
    % Output:
    %   resultsCell - A cell array where each cell contains a struct with
    %                 information for one possible loop combination. This format
    %                 is compatible with the MATLAB Engine for Python.
    %                 Each struct has the fields:
    %                   .codes                - A string with all the line ids, separated by comma
    %                   .combination          - A K-by-8 matrix of the specific lines.
    %                   .misclosure           - The calculated geodetic misclosure.
    %                   .totalDistance        - The total distance of the loop.
    %                   .totalSetups          - The total number of instrument setups.
    %                   .totalBFMisclosure    - The sum of back/forward misclosures.
    %                   .totalGravCorrection  - The sum of the gravity correction
    %                   .misclosureWithGravity - The misclosure after applying gravity corrections.
    %
    % Description:
    %   This function takes levelling data for a single closed loop that may
    %   contain multiple measurements for the same line segment. It first
    %   identifies all unique line segments and groups their corresponding
    %   measurements. It then traces the loop's path to establish a
    %   consistent direction for calculation. Finally, it generates all
    %   possible unique loop combinations and, for each one, calculates the
    %   total misclosure, distance, and other relevant metrics.

    %% --- Input Validation ---
    if nargin < 1 || isempty(levellingData)
        error('Input levellingData is required.');
    end
    if size(levellingData, 2) ~= 9
        error('Input levellingData must have 9 columns.');
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
    startNode = allPoints(1);
    path = [startNode];
    prevNode = -1; 
    currNode = startNode;

    for i = 1:(numel(allPoints) - 1)
        neighbors = adj(currNode);
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
    
    neighborsOfLast = adj(path(end));
    if ~any(neighborsOfLast == startNode)
        error('Could not determine a single closed loop from the data provided.');
    end

    %% --- Step 3: Prepare for Combination Generation ---
    numSegments = numel(path);
    orderedSegmentGroups = cell(1, numSegments);
    groupSizes = zeros(1, numSegments);

    for i = 1:numSegments
        p1 = path(i);
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
    numCombinations = prod(groupSizes);
    gridVectors = arrayfun(@(n) 1:n, groupSizes, 'UniformOutput', false);
    
    indexGrids = cell(1, numSegments);
    [indexGrids{:}] = ndgrid(gridVectors{:});

    combinationIndices = zeros(numCombinations, numSegments);
    for i = 1:numSegments
        combinationIndices(:, i) = indexGrids{i}(:);
    end

    % Pre-allocate the results struct with the new fields
    results = struct('codes', cell(numCombinations, 1), ...
                     'combination', cell(numCombinations, 1), ...
                     'misclosure', cell(numCombinations, 1), ...
                     'totalDistance', cell(numCombinations, 1), ...
                     'totalSetups', cell(numCombinations, 1), ...
                     'totalBFMisclosure', cell(numCombinations, 1), ...
                     'totalGravCorrection', cell(numCombinations, 1), ...
                     'misclosureWithGravity', cell(numCombinations, 1));

    % Iterate over each combination to perform the calculations
    for i = 1:numCombinations
        currentCombinationLines = cell(1, numSegments);
        currentCodesString = "";

        totalDistance = 0;
        misclosure = 0;
        totalSetups = 0;
        totalBFMisclosure = 0;
        totalGravCorrection = 0;
        misclosureWithGravity = 0;

        % Build the combination and calculate sums
        for j = 1:numSegments
            group = orderedSegmentGroups{j};
            lineIndex = combinationIndices(i, j);
            lineData = group{lineIndex};
            
            currentCombinationLines{j} = lineData;
            currentCodesString = currentCodesString + string(lineData(1)) +",";
            % Sum direction-independent fields
            totalDistance = totalDistance + lineData(5);
            totalSetups = totalSetups + lineData(6);
            totalBFMisclosure = totalBFMisclosure + lineData(7);
            
            % Sum direction-dependent fields
            path_p1 = path(j); 
            
            
            if lineData(2) == path_p1 
                % Measurement direction matches the path traversal direction
                misclosure = misclosure + lineData(4);
                misclosureWithGravity = misclosureWithGravity + lineData(9);
                totalGravCorrection = totalGravCorrection + lineData(8);
            else 
                % Measurement direction is opposite to the path direction
                misclosure = misclosure - lineData(4);
                misclosureWithGravity = misclosureWithGravity - lineData(9);
                totalGravCorrection = totalGravCorrection + lineData(8);
            end
        end
        
        % Store the results for this combination
        results(i).codes = extractBefore(currentCodesString, length(char(currentCodesString)));
        results(i).combination = vertcat(currentCombinationLines{:});
        results(i).misclosure = misclosure;
        results(i).totalDistance = totalDistance;
        results(i).totalSetups = totalSetups;
        results(i).totalBFMisclosure = totalBFMisclosure;
        results(i).totalGravCorrection = totalGravCorrection;
        results(i).misclosureWithGravity = misclosureWithGravity;
    end
    
    %% --- Step 5: Package for Python Engine ---
    if numCombinations > 0
        resultsCell = num2cell(results);
    else
        resultsCell = {}; % Return an empty cell if no combinations found
    end
end
