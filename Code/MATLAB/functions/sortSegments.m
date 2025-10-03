function sortedSegments = sortSegments(segmentEnds, junctionCodes, regularCodes)
% segmentEnds: Nx2 double (each row is [start_code, end_code])
% junctionCodes: vector of junction + dead-end point codes
% regularCodes: vector of regular point codes
%
% Input validation:
% Ensure segmentEnds is a non-empty Nx2 matrix
% Ensure junctionCodes and regularCodes are non-empty vectors
% Ensure all codes in segmentEnds are present in either junctionCodes or regularCodes
%
% Output:
%   sortedInfo - struct with fields:
%       sorted_index: vector of indices of segments after the sorting
%                       process
%       segment_index: vector of indices of segments
%       line_number: vector of line numbers assigned
%       reverse_flag: 1 if the segment was reversed, 0 otherwise
%    The output is organized in the same order as the input

% Validate inputs
if isempty(segmentEnds) || size(segmentEnds, 2) ~= 2
    error('segmentEnds must be a non-empty Nx2 matrix.');
end
if isempty(junctionCodes) || isempty(regularCodes)
    error('junctionCodes and regularCodes must be non-empty vectors.');
end
%if any(~ismember(segmentEnds(:), [junctionCodes; regularCodes]))
%    error('All codes in segmentEnds must be present in junctionCodes or regularCodes.');
%end

% Initialize
numSegments = size(segmentEnds, 1);
used = false(numSegments, 1);
segment_index = [];
line_number = [];
reverse_flag = [];

line_num = 0;

for i = 1:length(junctionCodes)
    start_code = junctionCodes(i);
    point_finished = false;

    while ~point_finished
        current_code = start_code;
        line_num = line_num + 1;
        line_done = false;

        while ~line_done
            % Find the next unused segment connected to current_code
            idx = find(~used & (segmentEnds(:,1) == current_code | segmentEnds(:,2) == current_code), 1);

            if isempty(idx)
                point_finished = true;
                line_num = line_num - 1;  % Rollback this line_num
                break;
            end

            used(idx) = true;
            seg = segmentEnds(idx,:);

            if seg(1) == current_code
                next_code = seg(2);
                rev = 0;
            else
                next_code = seg(1);
                rev = 1;
            end

            segment_index(end+1) = idx;
            line_number(end+1) = line_num;
            reverse_flag(end+1) = rev;

            if ismember(next_code, junctionCodes)
                line_done = true;
            elseif ismember(next_code, regularCodes)
                current_code = next_code;
            else
                line_done = true;
                point_finished = true;
            end
        end
    end
end

% build the scalar struct with arrays (compatible with matlab.engine in Python)
sortedSegments = struct( ...
    'sorted_index', (1:length(segment_index)),...
    'segment_index', segment_index, ...
    'line_number', line_number, ...
    'reverse_flag', reverse_flag ...
);

% Sort in original order of the input segments and return
[~, sortedIndices] = sort(sortedSegments.segment_index);
sortedSegments.line_number = sortedSegments.line_number(sortedIndices);
sortedSegments.sorted_index = sortedSegments.sorted_index(sortedIndices);
sortedSegments.reverse_flag = sortedSegments.reverse_flag(sortedIndices);
sortedSegments.segment_index = sortedSegments.segment_index(sortedIndices);
