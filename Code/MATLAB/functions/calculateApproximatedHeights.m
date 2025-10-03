function results = calculateApproximatedHeights(segmentTable, pointTable)
% Performs geodetic height network adjustment on subnetworks.
% Input:
%   segmentTable - table of segments (OBJECTID, SegmentCode, StartPointCode, EndPointCode, MeasHeightDiff, MeasDist, SegmentDirection)
%   pointTable   - table of known points (Code, H)
% Output:
%   results      - table with computed heights, standard deviations and sigma0

    % Step 1: Average repeated segments by SegmentCode
    codes = unique(segmentTable.SegmentCode);
    meanSegments = table();

    for i = 1:length(codes)
        code = codes{i};
        inds = strcmp(segmentTable.SegmentCode, code);
        sub = segmentTable(inds, :);

        dh_corrected = sub.MeasHeightDiff .* sub.SegmentDirection;
        k = 2;
        Dist_km = sub.MeasDist / 1000;
        weights = 1 ./ (k^2 * Dist_km);

        dH_avg = sum(dh_corrected .* weights) / sum(weights);
        std_eff = sqrt(1 / sum(weights));

        main = sub(find(sub.SegmentDirection == 1, 1), :);
        if isempty(main)
            main = sub(1, :);
        end

        meanSegments = [meanSegments; ...
            table({code}, string(main.StartPointCode), string(main.EndPointCode), ...
                  dH_avg, std_eff, main.MeasDist, main.OBJECTID, ...
                'VariableNames', {'SegmentCode','StartPointCode','EndPointCode','MeasHeightDiff','std','MeasDist','OBJECTID'})];
    end

    % Step 2: Remove anchor points and split graph
    G = graph(meanSegments.StartPointCode, meanSegments.EndPointCode);
    breakPoints = string(pointTable.Code);
    G_reduced = rmnode(G, breakPoints);

    % Step 3: Find connected components
    bins = conncomp(G_reduced);
    uniqueBins = unique(bins);
    results = table();

    for b = uniqueBins
        node_names_reduced = G_reduced.Nodes.Name;
        nodes = node_names_reduced(bins == b);
        subG = subgraph(G, nodes);

        all_neighbors = [];
        for i = 1:length(nodes)
            current_neighbors = neighbors(G, nodes(i));
            all_neighbors = [all_neighbors; current_neighbors];
        end

        neighborNames = unique(string(all_neighbors));
        nodes = string(nodes);
        neighborNames = string(neighborNames);
        breakPoints = string(breakPoints);
        junction_neighbors = intersect(neighborNames, breakPoints);
        allNodes = unique([nodes; junction_neighbors]);

        % Extract matching segments
        idx = ismember(meanSegments.StartPointCode, allNodes) & ...
              ismember(meanSegments.EndPointCode, allNodes);
        segs = meanSegments(idx, :);

        known_in_sub = intersect(allNodes, string(pointTable.Code));
        unknown_in_sub = setdiff(allNodes, known_in_sub);

        if length(known_in_sub) < 1
            warning("No known height point in subnetwork — skipping...");
            continue;
        end

        n_unknown = length(unknown_in_sub);
        m = height(segs);
        A = zeros(m, n_unknown);
        l = zeros(m, 1);
        P = zeros(m, m);

        for i = 1:m
            from = segs.StartPointCode{i};
            to = segs.EndPointCode{i};
            l(i) = segs.MeasHeightDiff(i);
            P(i, i) = 1 / segs.std(i)^2;

            if ismember(from, unknown_in_sub)
                A(i, find(strcmp(unknown_in_sub, string(from)))) = -1;
            end
            if ismember(to, unknown_in_sub)
                A(i, find(strcmp(unknown_in_sub, string(to)))) = 1;
            end
        end

        for i = 1:m
            from = segs.StartPointCode{i};
            to = segs.EndPointCode{i};
            deltaH = 0;
            if ismember(from, known_in_sub)
                deltaH = deltaH - pointTable.H(strcmp(string(pointTable.Code), from));
            end
            if ismember(to, known_in_sub)
                deltaH = deltaH + pointTable.H(strcmp(string(pointTable.Code), to));
            end
            l(i) = l(i) - deltaH;
        end

        x_hat = (A' * P * A) \ (A' * P * l);
        v = A * x_hat - l;
        sigma0 = sqrt(v' * P * v / (m - n_unknown));
        Cov = sigma0^2 * inv(A' * P * A);
        sigmas = sqrt(diag(Cov));

        if ~isempty(x_hat)
            results = [results; table(unknown_in_sub, x_hat, sigmas, repmat(sigma0, n_unknown, 1), ...
                'VariableNames', {'Code', 'CalculatedHeight', 'StdDev', 'Sigma0'})];
        else
            warning("Adjustment failed — skipping...");
        end
    end
end