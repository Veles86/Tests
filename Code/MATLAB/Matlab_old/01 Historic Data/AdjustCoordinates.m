clc
clear all
close all

[fileName,path] = uigetfile('*.csv','Select CSV file','MultiSelect','on');

if ~iscell(fileName)
    fileName = {fileName};
    k=1;
else
    k = length(fileName);
end



% REGULAR_DH_MODE = 0; % 1 if dh values in files have decimal point, 0 if without decimal point

FILTER_COORDINATES = 1; % 1 if want to remove less significant groups (like polygon and private) from the coordinates table

% Parameters for splitting tables between good and bad
% default values 0.2 and 2000
DH_THRESHOLD = 0.2;
DIST_THRESHOLD = 2000;



opts = detectImportOptions('coordinates\all_from_DB_english.csv');
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinates = readtable('coordinates\all_from_DB_english.csv',opts);
coordinates = [coordinates;readtable('coordinates\all_from_DB_hebrew.csv',opts)];
coordinates = [coordinates;readtable('coordinates\missing_points.csv',opts)];
coordinates = [coordinates;readtable('coordinates\temporary_missing_points.csv',opts)];
% coordinates = [coordinates;readtable('coordinates\old_BM_not_in_DB.csv',opts)];
coordinates = [coordinates;readtable('coordinates\pickets.csv',opts)];
coordinates = [coordinates;readtable('coordinates\eccentric_points.csv',opts)];

if FILTER_COORDINATES

coordinates(coordinates.group==8,:)=[];
coordinates(coordinates.group==11,:)=[];

end

%%


for u = 1:k

tempFileName = fileName{u};
opts = detectImportOptions(strcat(path,tempFileName));
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(9) = {'char'};
opts.VariableTypes(10) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
inputFile = readtable(strcat(path,tempFileName),opts);

tempDH = inputFile.dH;
if sum(mod(tempDH,1))==0
    REGULAR_DH_MODE = 0;
else
    REGULAR_DH_MODE = 1;
end


tempFileName = tempFileName(1:end-4);
n = height(inputFile);
%%
output_point_names = table('Size',[n 3],'VariableTypes',{'double','string','string'},'VariableNames',["num","pntA","pntB"]);
output_point_idn = table('Size',[n 10],'VariableTypes',{'double','string','string','double','double','double','string','string','double','double'},'VariableNames',["codeA","typeA","nameA","groupA","rankA","codeB","typeB","nameB","groupB","rankB"]);
output_coordinates = table('Size',[n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["Y1","X1","H1","Y2","X2","H2"]);
output_idn_quality = table('Size',[n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["dH_from_DB","dH_from_file","dH_diff","dist_from_DB","dist_from_file","dist_diff"]);

fullOutput = [];

for i = 1 : n
    %     temp = inputFile(i,:);
    pntA = string(inputFile{i,"A"});
    pntB = string(inputFile{i,"B"});
    dH_file = inputFile{i,"dH"};
    dist_file = inputFile{i,"dist"};

    rows = strcmpi(string(coordinates.full_name),(pntA));
    optionsA = coordinates(rows, :);
    optionsA = sortrows(optionsA,["R","H"]);

    rows = strcmpi(string(coordinates.full_name),(pntB));
    optionsB = coordinates(rows, :);
    optionsB = sortrows(optionsB,["R","H"]);

    output_point_names(i,:) = table(inputFile{i,"num"},pntA,pntB);

    if ~strcmpi(pntA,pntB)
        allOptions = FindAllOptions(optionsA,optionsB,dH_file,dist_file,REGULAR_DH_MODE);
        bestOption = FindBestOption(allOptions,DH_THRESHOLD,DIST_THRESHOLD);
        output_point_idn(i,:) = [bestOption(:,1:10)];
        output_coordinates(i,:) = [bestOption(:,11:16)];
        output_idn_quality(i,:) = bestOption(:,17:22);
        m = height(allOptions);
        tempTable = table(ones(m,1)*inputFile{i,"num"},repmat(pntA,m,1),repmat(pntB,m,1),'VariableNames',["num","pntA","pntB"]);
        allOptions = [tempTable,allOptions];

        fullOutput = [fullOutput;allOptions];


    end

end

tempTable = table((1:height(fullOutput))','VariableNames',"i");
fullOptionsTable = [tempTable,fullOutput];
bestOptionsTable = [output_point_names,output_point_idn,output_coordinates,output_idn_quality];

writetable(fullOptionsTable,strcat('output_adjusting_coordinates\',tempFileName,'_full.csv'));
writetable(bestOptionsTable,strcat('output_adjusting_coordinates\',tempFileName,'_best.csv'));


mergedTable = [inputFile,output_point_idn,output_coordinates,output_idn_quality];
writetable(mergedTable,strcat('output_adjusting_coordinates\',tempFileName,'_merged_all.csv'));



[mergedTableGood, mergedTableBad, mergedTableUnknown] = SplitMergedTable(mergedTable,DH_THRESHOLD,DIST_THRESHOLD);

writetable(mergedTableUnknown,strcat('output_adjusting_coordinates\',tempFileName,'_merged_unknown.csv'));

mergedTableBad = AddWKTField(mergedTableBad);
% if ~REGULAR_DH_MODE
%    mergedTableBad.dH = mergedTableBad.dH_from_file;
%    mergedTableBad.diff = mergedTableBad.diff/100;
% end
writetable(mergedTableBad,strcat('output_adjusting_coordinates\',tempFileName,'_merged_bad.csv'));

mergedTableGood.A = mergedTableGood.nameA;
mergedTableGood.B = mergedTableGood.nameB;
if ~REGULAR_DH_MODE
   mergedTableGood.dH = mergedTableGood.dH_from_file;
   mergedTableGood.diff = mergedTableGood.diff/100;
end

mergedTableGood = AddWKTField(mergedTableGood);
writetable(mergedTableGood,strcat('output_adjusting_coordinates\',tempFileName,'_merged_good.csv'));

remainedInput = [mergedTableBad(:,1:14);mergedTableUnknown(:,1:14)];
writetable(remainedInput,strcat('output_adjusting_coordinates\',tempFileName,'_remained_input.csv'));


coordinatesGood = GetCoordinates(mergedTableGood);
coordinatesBad = GetCoordinates(mergedTableBad);
coordinatesUnknown = GetCoordinates(mergedTableUnknown);

writetable(coordinatesGood,strcat('output_adjusting_coordinates\',tempFileName,'_coordinates_good.csv'));
writetable(coordinatesBad,strcat('output_adjusting_coordinates\',tempFileName,'_coordinates_bad.csv'));
writetable(coordinatesUnknown,strcat('output_adjusting_coordinates\',tempFileName,'_coordinates_unknown.csv'));

rowsA = isnan(mergedTable.codeA);
rowsB = isnan(mergedTable.codeB);

pointsNotFound = sort(unique(string([mergedTable.A(rowsA);mergedTable.B(rowsB)])));
notFoundPoints = table(pointsNotFound,'VariableNames',"pointName");
writetable(notFoundPoints,strcat('output_adjusting_coordinates\',tempFileName,'_not_found_points.csv'));


end




function [allOptions] = FindAllOptions(optionsA,optionsB,dH_file,dist_file,REGULAR_DH_MODE)

m = size(optionsA,1);
n = size(optionsB,1);

if (m~=0 && n~=0)
    temp1 = table('Size',[m*n 10],'VariableTypes',{'double','string','string','double','double','double','string','string','double','double'},'VariableNames',["codeA","typeA","nameA","groupA","rankA","codeB","typeB","nameB","groupB","rankB"]);
    temp2 = table('Size',[m*n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["Y1","X1","H1","Y2","X2","H2"]);
    temp3 = table('Size',[m*n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["dH_from_DB","dH_from_file","dH_diff","dist_from_DB","dist_from_file","dist_diff"]);

    counter = 1;
    for i = 1 : m
        for j = 1 : n

            temp1.codeA(counter) = optionsA.code(i);
            temp1.typeA(counter) = optionsA.type(i);
            temp1.nameA(counter) = optionsA.full_name(i);
            temp1.groupA(counter) = optionsA.group(i);
            temp1.rankA(counter) = optionsA.R(i);

            y1 = optionsA.Y(i);
            x1 = optionsA.X(i);
            h1 = optionsA.H(i);

            temp2.Y1(counter) = y1;
            temp2.X1(counter) = x1;
            temp2.H1(counter) = h1;


            temp1.codeB(counter) = optionsB.code(j);
            temp1.typeB(counter) = optionsB.type(j);
            temp1.nameB(counter) = optionsB.full_name(j);
            temp1.groupB(counter) = optionsB.group(j);
            temp1.rankB(counter) = optionsB.R(j);

            y2 = optionsB.Y(j);
            x2 = optionsB.X(j);
            h2 = optionsB.H(j);

            temp2.Y2(counter) = y2;
            temp2.X2(counter) = x2;
            temp2.H2(counter) = h2;

            if (~isnan(h1) && ~isnan(h2))
                temp3.dH_from_DB(counter) = h2-h1;

                if (REGULAR_DH_MODE)
                    temp3.dH_from_file(counter) = dH_file;
                else

                    if (abs(dH_file/1000-(h2-h1))<abs(dH_file/10000-(h2-h1))) && (abs(dH_file/10000-(h2-h1))<abs(dH_file/100000-(h2-h1)))
                        temp3.dH_from_file(counter) = dH_file/1000;
                    elseif (abs(dH_file/10000-(h2-h1))<abs(dH_file/100000-(h2-h1))) && (abs(dH_file/10000-(h2-h1))<abs(dH_file/1000-(h2-h1)))
                        temp3.dH_from_file(counter) = dH_file/10000;
                    else
                        temp3.dH_from_file(counter) = dH_file/100000;

                    end
                end

                temp3.dH_diff(counter) = abs(temp3.dH_from_DB(counter)-temp3.dH_from_file(counter));

            else
                temp3.dH_from_DB(counter) = NaN;

                if (REGULAR_DH_MODE)
                    temp3.dH_from_file(counter) = dH_file;
                else
                    temp3.dH_from_file(counter) = dH_file/100000;
                end
                temp3.dH_diff(counter) = NaN;

            end

            temp3.dist_from_DB(counter) = round(sqrt((x2-x1)^2+(y2-y1)^2));
            temp3.dist_from_file(counter) = dist_file;
            temp3.dist_diff(counter) = abs(temp3.dist_from_file(counter)-temp3.dist_from_DB(counter));
            counter = counter + 1;

        end
    end


elseif (m==0 && n==0)

    temp1 = table(NaN,"NO_DATA","NO_DATE",NaN,NaN,NaN,"NO_DATA","NO_DATA",NaN,NaN,'VariableNames',["codeA","typeA","nameA","groupA","rankA","codeB","typeB","nameB","groupB","rankB"]);
    temp2 = table(NaN,NaN,NaN,NaN,NaN,NaN,'VariableNames',["Y1","X1","H1","Y2","X2","H2"]);
    temp3 = table(NaN,dH_file/100000,NaN,NaN,dist_file,NaN,'VariableNames',["dH_from_DB","dH_from_file","dH_diff","dist_from_DB","dist_from_file","dist_diff"]);


elseif (m==0)

    temp1 = table('Size',[n 10],'VariableTypes',{'double','string','string','double','double','double','string','string','double','double'},'VariableNames',["codeA","typeA","nameA","groupA","rankA","codeB","typeB","nameB","groupB","rankB"]);
    temp2 = table('Size',[n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["Y1","X1","H1","Y2","X2","H2"]);
    temp3 = table('Size',[n 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["dH_from_DB","dH_from_file","dH_diff","dist_from_DB","dist_from_file","dist_diff"]);


    for i = 1 :n
        temp1.codeA(i) = NaN;
        temp1.typeA(i) = "NO_DATA";
        temp1.nameA(i) = "NO_DATA";
        temp1.rankA(i) = NaN;
        temp1.groupA(i) = NaN;

        temp1.codeB(i) = optionsB.code(i);
        temp1.typeB(i) = optionsB.type(i);
        temp1.nameB(i) = optionsB.full_name(i);
        temp1.rankB(i) = optionsB.R(i);
        temp1.groupB(i) = optionsB.group(i);

        temp2.Y1(i) = NaN;
        temp2.X1(i) = NaN;
        temp2.H1(i) = NaN;

        temp2.Y2(i) = optionsB.Y(i);
        temp2.X2(i) = optionsB.X(i);
        temp2.H2(i) = optionsB.H(i);

        temp3.dH_from_DB(i) = NaN;

        if REGULAR_DH_MODE
            temp3.dH_from_file(i) = dH_file;
        else
            temp3.dH_from_file(i) = dH_file/100000;
        end

        temp3.dH_diff(i) = NaN;
        temp3.dist_from_DB(i) = NaN;
        temp3.dist_from_file(i) = dist_file;
        temp3.dist_diff(i) = NaN;
    end


else % n==0

    temp1 = table('Size',[m 10],'VariableTypes',{'double','string','string','double','double','double','string','string','double','double'},'VariableNames',["codeA","typeA","nameA","groupA","rankA","codeB","typeB","nameB","groupB","rankB"]);
    temp2 = table('Size',[m 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["Y1","X1","H1","Y2","X2","H2"]);
    temp3 = table('Size',[m 6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',["dH_from_DB","dH_from_file","dH_diff","dist_from_DB","dist_from_file","dist_diff"]);

    for i = 1 :m
        
        temp1.codeA(i) = optionsA.code(i);
        temp1.typeA(i) = optionsA.type(i);
        temp1.nameA(i) = optionsA.full_name(i);
        temp1.rankA(i) = optionsA.R(i);
        temp1.groupA(i) = optionsA.group(i);
        
        temp1.codeB(i) = NaN;
        temp1.typeB(i) = "NO_DATA";
        temp1.nameB(i) = "NO_DATA";
        temp1.rankB(i) = NaN;
        temp1.groupB(i) = NaN;

        temp2.Y1(i) = optionsA.Y(i);
        temp2.X1(i) = optionsA.X(i);
        temp2.H1(i) = optionsA.H(i);

        temp2.Y2(i) = NaN;
        temp2.X2(i) = NaN;
        temp2.H2(i) = NaN;

        temp3.dH_from_DB(i) = NaN;

        if REGULAR_DH_MODE
            temp3.dH_from_file(i) = dH_file;
        else
            temp3.dH_from_file(i) = dH_file/100000;
        end

        temp3.dH_diff(i) = NaN;
        temp3.dist_from_DB(i) = NaN;
        temp3.dist_from_file(i) = dist_file;
        temp3.dist_diff(i) = NaN;
    end

end

allOptions = [temp1,temp2,temp3];


end



function [bestOption] = FindBestOption(allOptions,DH_THRESHOLD,DIST_THRESHOLD)

if height(allOptions)==1
    bestOption=allOptions;
else

    allOptions = sortrows(allOptions,["dH_diff","dist_diff"]);

    rows = (~isnan(allOptions.dH_diff));
    by_dH_diff = allOptions(rows, :);
    if height(by_dH_diff)>0
        best_by_dH_diff = by_dH_diff(1,:);
    end

    rows = (isnan(allOptions.dH_diff));
    by_dist_diff = allOptions(rows, :);
    if height (by_dist_diff)>0
        best_by_dist_diff = by_dist_diff(1,:);
    end

    if height(by_dH_diff)==0
        bestOption = best_by_dist_diff;
    elseif height (by_dist_diff)==0
        bestOption = best_by_dH_diff;
    elseif (best_by_dH_diff.dH_diff < DH_THRESHOLD)
        bestOption = best_by_dH_diff;
    elseif (best_by_dist_diff.dist_diff<DIST_THRESHOLD)
        bestOption = best_by_dist_diff;
    else
        bestOption = best_by_dH_diff;
    end


end





end


function [goodTable, badTable, unknownTable] = SplitMergedTable(mergedTable,DH_THRESHOLD,DIST_THRESHOLD)

rows_unknown = isnan(mergedTable.dist_diff);
unknownTable = mergedTable(rows_unknown,:);
mergedTable = mergedTable(~rows_unknown,:);
rows_good = (mergedTable.dH_diff<DH_THRESHOLD & mergedTable.dist_diff<DIST_THRESHOLD) | ((isnan(mergedTable.dH_diff) & mergedTable.dist_diff<DIST_THRESHOLD));
goodTable = mergedTable(rows_good,:);
badTable = mergedTable(~rows_good,:);

end

function [coordinatesTable] = GetCoordinates(InputMergedTable)

n = height(InputMergedTable);

coordinatesTable = table('Size',[2*n 8],'VariableTypes',{'double','string','string','double','double','double','double','double'},'VariableNames',["code","type","name","group","rank","Y","X","H"]);
coordinatesTable(1:n,:) = [InputMergedTable(:,15:19),InputMergedTable(:,25:27)];
coordinatesTable(n+1:end,:) = [InputMergedTable(:,20:24),InputMergedTable(:,28:30)];

rows = isnan(coordinatesTable.code);
coordinatesTable(rows,:) = [];

rows = isnan(coordinatesTable.rank);
coordinatesTable.rank(rows) = 9999;
rows = isnan(coordinatesTable.group);
coordinatesTable.group(rows) = 9999;
rows = isnan(coordinatesTable.H);
coordinatesTable.H(rows) = 9999;
coordinatesTable = unique(coordinatesTable,'rows');


rows = coordinatesTable.rank==9999;
coordinatesTable.rank(rows) = NaN;
rows = coordinatesTable.group==9999;
coordinatesTable.group(rows) = NaN;
rows = coordinatesTable.H==9999;
coordinatesTable.H(rows) = NaN;

coordinatesTable = sortrows(coordinatesTable,1);

end

function [outputMergedTable] = AddWKTField(InputMergedTable)

n = height(InputMergedTable);

WKT_column = repmat("",n,1);

for i=1:n

WKT_column(i) = strcat("LINESTRING(",string(InputMergedTable.Y1(i))," ",string(InputMergedTable.X1(i)),":",string(InputMergedTable.Y2(i))," ",string(InputMergedTable.X2(i)),")");

end

outputMergedTable = [InputMergedTable,table(WKT_column,'VariableNames',["WKT"])];



end

