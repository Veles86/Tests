function results = analyzeLevellingLoopGraph(line_data)
% ANALYZELEVELLINGLOOPGRAPH Finds valid closed loop paths using DFS and computes misclosure.
% INPUT:
%   line_data - N×5 table: LineNumber, StartPoint, EndPoint, Distance, MeasuredHeightDiff
% OUTPUT:
%   results - struct array with fields: LineSequence, TotalDistance, HeightMisclosure

    % Extract data
    n = height(line_data);
    lines = table2struct(line_data);

    % Build directed graph: edges from StartPoint to EndPoint with line indices
    G = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    for i = 1:n
        key = lines(i).StartPoint;
        if ~isKey(G, key)
            G(key) = [];
        end
        G(key) = [G(key); i]; % Store index of line
    end

    results = [];
    result_idx = 1;

    % Try starting from each line
    for i = 1:n
        start_code = lines(i).StartPoint;
        visited = false(1, n);
        path = [];
        dfs(i, lines(i).StartPoint, lines(i).EndPoint, 1);
    end

    function dfs(curr_idx, start_code, current_code, depth)
        % DFS to find closed loop paths
        visited(curr_idx) = true;
        path(depth) = curr_idx;

        if current_code == start_code && depth > 1
            % Valid loop found
            dist = 0;
            dh = 0;
            for k = 1:depth
                dist = dist + lines(path(k)).Distance;
                dh = dh + lines(path(k)).MeasuredHeightDiff;
            end
            results(result_idx).LineSequence = [lines(path).LineNumber]; %#ok<AGROW>
            results(result_idx).TotalDistance = dist;
            results(result_idx).HeightMisclosure = dh;
            result_idx = result_idx + 1;
        else
            % Continue DFS
            if isKey(G, current_code)
                next_lines = G(current_code);
                for ni = 1:length(next_lines)
                    ni_idx = next_lines(ni);
                    if ~visited(ni_idx)
                        dfs(ni_idx, start_code, lines(ni_idx).EndPoint, depth + 1);
                    end
                end
            end
        end

        visited(curr_idx) = false;
    end
end
