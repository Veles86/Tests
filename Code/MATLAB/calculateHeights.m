% קוד MATLAB מלא ומעודכן לביצוע תיאום רשתות גובה בתתי-רשתות
% כולל מיצוע קטעים חוזרים, חישוב מטריציוני, ושמירת סטיית תקן

% שלב 1: קריאת הנתונים
segments = readtable('Segments.csv');        % כולל: SegmentCode, FROM, TO, dH, MeasDist, std, SegmentDirection
known = readtable('Points.csv');       % כולל: Point, H, Grade

% שלב 2: מיצוע קטעים חוזרים לפי SegmentCode
codes = unique(segments.SegmentCode);
meanSegments = table();

for i = 1:length(codes)
    code = codes{i};
    inds = strcmp(segments.SegmentCode, code);
    sub = segments(inds, :);

    % התאמת כיוון הפרשי גובה לפי SegmentDirection
    dh_corrected = sub.MeasHeightDiff .* sub.SegmentDirection;

    k = 2;
    Dist_km = sub.MeasDist / 1000;
    weights = 1 ./ (k^2 * Dist_km);



    dH_avg = sum(dh_corrected .* weights) / sum(weights);
    std_eff = sqrt(1 / sum(weights));

    % שמירת כיוון הקטע לפי ההגדרה הראשית (SegmentDirection == 1)
    main = sub(find(sub.SegmentDirection == 1, 1), :);
    if isempty(main)
        main = sub(1, :); % אם אין מופע בכיוון חיובי
    end

    meanSegments = [meanSegments; table({code}, string(main.StartPointCode), string(main.EndPointCode), dH_avg, std_eff, main.MeasDist, main.UniqueRowNum, ...
        'VariableNames', {'SegmentCode','StartPointCode','EndPointCode','MeasHeightDiff','std','MeasDist','UniqueRowNum'})];
end

% שלב 3: בניית גרף קשרים והסרת נקודות עוגן (GRADE <= 3) לפירוק לרשתות
G = graph(meanSegments.StartPointCode, meanSegments.EndPointCode);
known(~(known.HRank <= 3),:)=[];
breakPoints = known.Code(known.HRank <= 4);
G_reduced = rmnode(G, breakPoints);

% שלב 4: מציאת רכיבים קשירים
bins = conncomp(G_reduced);
uniqueBins = unique(bins);
results = table();

for b = uniqueBins
    node_names_reduced = G_reduced.Nodes.Name;  % this gives the node names in order
    nodes = node_names_reduced(bins == b);      % now safe to index

    nodes = G_reduced.Nodes.Name(bins == b);
    subG = subgraph(G, nodes);

    all_neighbors = [];

    for i = 1:length(nodes)
        current_neighbors = neighbors(G, nodes(i));
        all_neighbors = [all_neighbors; current_neighbors];
    end
    
    % הסר כפילויות והמר למחרוזות
    neighborNames = unique(string(all_neighbors));
    nodes = string(nodes);
    neighborNames = string(neighborNames);
    breakPoints = string(breakPoints);
    junction_neighbors = intersect(neighborNames, breakPoints);
    allNodes = unique([nodes; junction_neighbors]);



    % שליפת קטעים מתאימים
    idx = ismember(meanSegments.StartPointCode, allNodes) & ismember(meanSegments.EndPointCode, allNodes);
    segs = meanSegments(idx, :);
    known = known(~isnan(known.H), :);
    % נקודות ידועות ולא ידועות
    known_in_sub = intersect(allNodes, string(known.Code));
    unknown_in_sub = setdiff(allNodes, string(known_in_sub));
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

    % תיקון l עבור נקודות ידועות
    for i = 1:m
        from = segs.StartPointCode{i};
        to = segs.EndPointCode{i};
        deltaH = 0;
        if ismember(from, known_in_sub)
            deltaH = deltaH - known.H(strcmp(string(known.Code), from));
        end
        if ismember(to, known_in_sub)
            deltaH = deltaH + known.H(strcmp(string(known.Code), to));
        end
        l(i) = l(i) - deltaH;
    end

    % פתרון Least Squares
    x_hat = (A' * P * A) \ (A' * P * l);
    v = A * x_hat - l;
    sigma0 = sqrt(v' * P * v / (m - n_unknown));
    Cov = sigma0^2 * inv(A' * P * A);
    sigmas = sqrt(diag(Cov));
    if ~isempty(x_hat)
        results = [results; table(unknown_in_sub, x_hat, sigmas, repmat(sigma0, n_unknown, 1), ...
         'VariableNames', {'Point', 'ComputedHeight', 'StdDev', 'Sigma0'})];
    else
        warning("No known height point in subnetwork — skipping...");
    end
end

% תצוגה או שמירה לקובץ
%disp(results);
writetable(results, 'computed_heights.csv');