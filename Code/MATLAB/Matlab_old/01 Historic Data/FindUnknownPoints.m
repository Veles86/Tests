clc
clear all
close all



opts = detectImportOptions('coordinates\all_from_DB_english.csv');
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinates = readtable('coordinates\all_from_DB_english.csv',opts);
coordinates = [coordinates;readtable('coordinates\all_from_DB_hebrew.csv',opts)];
coordinates = [coordinates;readtable('coordinates\missing_points.csv',opts)];
coordinates = [coordinates;readtable('coordinates\temporary_missing_points.csv',opts)];
coordinates = [coordinates;readtable('coordinates\old_BM_not_in_DB.csv',opts)];
coordinates = [coordinates;readtable('coordinates\pickets.csv',opts)];
coordinates = [coordinates;readtable('coordinates\eccentric_points.csv',opts)];


[fileName,path] = uigetfile('*.csv','Select Merged CSV file');


opts = detectImportOptions(strcat(path,fileName),'Delimiter',',');
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(9) = {'char'};
opts.VariableTypes(10) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
opts.VariableTypes(17) = {'char'};
opts.VariableTypes(22) = {'char'};
inputFile = readtable(strcat(path,fileName),opts);
fileName = fileName(1:end-4);


%% CHOOSING DATA FOR WANTING POINTS - SELECTING FROM CSV OR ENTERING MANUALLY

% [fileName2,path2] = uigetfile('*.csv','Select list of unknown points you want to test',path);
% opts = detectImportOptions(strcat(path2,fileName2));
% opts.VariableTypes(1) = {'char'};
% wantedPoints = readtable(strcat(path2,fileName2),opts);


wanted_points_names = {'4068MPI'};

wantedPoints = table(wanted_points_names');
wantedPoints.Properties.VariableNames = {'pointName'};



%%
% the radius for searching similar points
RADIUS_XY = 1000;

if sum(mod(inputFile.dH,1))==0
    REGULAR_DH_MODE = 0;
else
    REGULAR_DH_MODE = 1;
end

n = height(wantedPoints);

approximatedCoordinates = table('Size',[n 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',["pointName","calcY","calcX","calcH","AB_dist","diff_dH","num_similar_pnts"]);

allSimilarPoints = [];
counter=1;
for i = 1: n
    tempPoint = wantedPoints.pointName{i};
    [connectedFromA,connectedFromB] = GetConnectedPoints(inputFile,tempPoint);
    [approxY,approxX,approxH,dist_between_AB,diff_dH] = CalculateApproximatedCoordinates(connectedFromA,connectedFromB);

    approximatedCoordinates.pointName(i) = tempPoint;
    approximatedCoordinates.calcY(i) = approxY;
    approximatedCoordinates.calcX(i) = approxX;
    approximatedCoordinates.calcH(i) = approxH;
    approximatedCoordinates.AB_dist(i) = dist_between_AB;
    approximatedCoordinates.diff_dH(i) = diff_dH;

    [tempSimilarPoints] = FindSimilarPoints(coordinates,approxY,approxX,approxH,RADIUS_XY);
    if width(tempSimilarPoints)~=1 % if not NaN
        h = height(tempSimilarPoints);
        approximatedCoordinates.num_similar_pnts(i) = h;

        tempIndexesTable = table('Size',[h 2],'VariableTypes',{'double','double'},'VariableNames',["index","pointNum"]);
        tempIndexesTable.index = [counter:counter+h-1]';
        tempIndexesTable.pointNum = repmat(i,h,1);

        tempApproxTable = repmat(approximatedCoordinates(i,:),h,1);

        tempTable = [tempIndexesTable,tempApproxTable,tempSimilarPoints];

        allSimilarPoints = [allSimilarPoints;tempTable];

        counter = counter + h;
    else
        approximatedCoordinates.num_similar_pnts(i) = 0;
    end


end


writetable(approximatedCoordinates,strcat('output_finding_unknown_points\',fileName,'_approximated_coordinates.csv'));
writetable(allSimilarPoints,strcat('output_finding_unknown_points\',fileName,'_similar_points.csv'));


function [connectedFromA,connectedFromB] = GetConnectedPoints(inputFile,pointName)


connectedIndexesA = (string(inputFile.B)==pointName);
connectedFromA = inputFile(connectedIndexesA,:);

connectedIndexesB = (string(inputFile.A)==pointName);
connectedFromB = inputFile(connectedIndexesB,:);

[~,ia,~] = unique(string(connectedFromA.A),'rows','stable');
connectedFromA = connectedFromA(ia,:);

[~,ia,~] = unique(string(connectedFromB.B),'rows','stable');
connectedFromB = connectedFromB(ia,:);

end

function [swappedRows] = SwapRows(dataRows)

    swappedRows = dataRows;
    swappedRows(:,["A","B"]) = dataRows(:,["B","A"]);
    swappedRows.dH = -dataRows.dH;

    swappedRows(:,["codeA","typeA","nameA","groupA","rankA"]) = dataRows(:,["codeB","typeB","nameB","groupB","rankB"]);
    swappedRows(:,["codeB","typeB","nameB","groupB","rankB"]) = dataRows(:,["codeA","typeA","nameA","groupA","rankA"]);
    swappedRows(:,["Y1","X1","H1"]) = dataRows(:,["Y2","X2","H2"]);
    swappedRows(:,["Y2","X2","H2"]) = dataRows(:,["Y1","X1","H1"]);


end


function [Y,X,H,totalDist,diff_dH] = CalculateApproximatedCoordinates(connectedFromA,connectedFromB)

if (isnan(sum(mod(connectedFromA.dH,1))) && sum(mod(connectedFromB.dH,1))==0) ||...
        (isnan(sum(mod(connectedFromB.dH,1))) && sum(mod(connectedFromA.dH,1))==0) ||...
        (sum(mod(connectedFromA.dH,1))+sum(mod(connectedFromB.dH,1)))==0
    REGULAR_DH_MODE = 0;
else
    REGULAR_DH_MODE = 1;
end


if height(connectedFromB)~=0
dataRows = [connectedFromA;SwapRows(connectedFromB)];
else
dataRows = connectedFromA;
end

[~,ia,~] = unique(string(dataRows.A),'rows','stable');
dataRows = dataRows(ia,:);

h = height(dataRows);

if h>=2
    num_of_couples = (h*h-h)/2;
    tempY = zeros(num_of_couples,1);
    tempX = zeros(num_of_couples,1);
    tempH = zeros(num_of_couples,1);
    tempTotalDist = zeros(num_of_couples,1);
    tempDiff_dH = zeros(num_of_couples,1);
    counter = 0;

    for i = 1: (h-1)
        for j = (i+1) : h
            counter = counter + 1;
            tempTotalDist(counter) = dataRows.dist(i) + dataRows.dist(j);
            relativeDist = dataRows.dist(i)/tempTotalDist(counter);

            deltaY = dataRows.Y1(j) - dataRows.Y1(i);
            deltaX = dataRows.X1(j) - dataRows.X1(i);

            tempY(counter) = round(dataRows.Y1(i) + deltaY*relativeDist);
            tempX(counter) = round(dataRows.X1(i) + deltaX*relativeDist);

            if REGULAR_DH_MODE
                tempH(counter) = round(nanmean([dataRows.H1(i)+dataRows.dH(i),dataRows.H1(j)+dataRows.dH(j)]),3);
                deltaH_db = dataRows.H1(j) - dataRows.H1(i);
                deltaH_file = dataRows.dH(i) - dataRows.dH(j);
                tempDiff_dH(counter) = round(abs(deltaH_db-deltaH_file),3);

            else
                deltaH_db = dataRows.H1(j) - dataRows.H1(i);
                deltaH_file = dataRows.dH(i) - dataRows.dH(j);

                if (abs(deltaH_file/1000-(deltaH_db))<abs(deltaH_file/10000-(deltaH_db))) && (abs(deltaH_file/10000-(deltaH_db))<abs(deltaH_file/100000-(deltaH_db)))
                    factor = 1000;
                elseif (abs(deltaH_file/10000-(deltaH_db))<abs(deltaH_file/100000-(deltaH_db))) && (abs(deltaH_file/10000-(deltaH_db))<abs(deltaH_file/1000-(deltaH_db)))
                    factor = 10000;
                else

                    factor = 100000;
                end

                tempH(counter) = round(nanmean([dataRows.H1(i)+dataRows.dH(i)/factor,dataRows.H1(j)+dataRows.dH(j)/factor]),3);
                tempDiff_dH(counter) = round(abs(deltaH_db-deltaH_file/factor),3);

            end

        end
    end


    Y = round(nanmean(tempY));
    X = round(nanmean(tempX));
    H = round(nanmean(tempH),3);
    totalDist = round(nanmean(tempTotalDist));
    diff_dH = round(nanmean(tempDiff_dH),3);
else

    Y = NaN;
    X = NaN;
    H = NaN;
    totalDist = NaN;
    diff_dH = NaN;


end



end

function [similarPoints] = FindSimilarPoints(coordinates,Y,X,H,deltaDist)

if isnan(Y) || isnan(X)
    similarPoints = NaN;

else
    rows = sqrt((coordinates.Y-Y).^2+(coordinates.X-X).^2)<deltaDist;
    foundPoints = coordinates(rows,:);

    n = height(foundPoints);

    if n==0
        similarPoints = NaN;
    else

        deltas = table('Size',[n 2],'VariableTypes',{'double','double'},'VariableNames',["dist","dH"]);

        deltas.dist = sqrt((foundPoints.Y-Y).^2+(foundPoints.X-X).^2);
        deltas.dH = abs(foundPoints.H - H);


        similarPoints = [foundPoints,deltas];
        similarPoints = sortrows(similarPoints,["dH","dist"]);


    end
end






end

% files = {'dar1.csv'};
% for u = 1:length(files)
%
% fileName = files{u};
% opts = detectImportOptions(strcat('data\',fileName));
% opts.VariableTypes(2) = {'char'};
% opts.VariableTypes(3) = {'char'};
% opts.VariableTypes(13) = {'double'};
% opts.VariableTypes(14) = {'double'};
% inputFile = readtable(strcat('data\',fileName),opts);
%
% end
%
%
% Y = 169738;
% X = 600870;
% dist = 100;
% H = 144.29;
% dH = 0.1;
%
%
% rows = sqrt((coordinates.Y-Y).^2+(coordinates.X-X).^2)<dist;
%
% foundPointsByXY = coordinates(rows,:);
%
% rows = (coordinates.H>(H-dH)) & (coordinates.H<(H+dH));
% foundPointsByH = coordinates(rows,:);