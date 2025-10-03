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
n = height(inputFile);




rowsA = string(inputFile.typeA)=='duplicate';
rowsB = string(inputFile.typeB)=='duplicate';

codesA = inputFile.codeA(rowsA);
codesB = inputFile.codeB(rowsB);

[namesA,groupsA] = GetMainNames(codesA,coordinates);
[namesB,groupsB] = GetMainNames(codesB,coordinates);

not_found_counter = sum(namesA=="NOT FOUND") + sum(namesB=="NOT FOUND");
found_too_much_counter = sum(namesA=="FOUND MORE THAN 1") + sum(namesB=="FOUND MORE THAN 1");
%%
if not_found_counter~=0
    warning('There are points for which the main name was not found');
end

if found_too_much_counter~=0
    warning('There are points for which more than 1 main name was found');
end
%%
inputFile.A(rowsA) = cellstr(namesA);
inputFile.B(rowsB) = cellstr(namesB);

inputFile.typeA(rowsA) = cellstr("main");
inputFile.typeB(rowsB) = cellstr("main");

inputFile.nameA(rowsA) = cellstr(namesA);
inputFile.nameB(rowsB) = cellstr(namesB);

inputFile.groupA(rowsA) = groupsA;
inputFile.groupB(rowsB) = groupsB;

writetable(inputFile,strcat('output_swapping_duplicate_names/',fileName,'_swapped.csv'));

updatedCoordinatesTable = GetCoordinates(inputFile);

writetable(updatedCoordinatesTable,strcat('output_swapping_duplicate_names/',fileName,'_swapped_coordinates.csv'));

function [mainNames,mainGroups] = GetMainNames(codesList,coordinates)

k = length(codesList);
mainNames = repmat("",k,1);
mainGroups=zeros(k,1);

for i = 1:k

    tempCode = codesList(i);
    allPointNames = coordinates(coordinates.code==tempCode,:);
    tempRows = string(allPointNames.type)=='main';
    tempName = allPointNames.full_name(tempRows);
    if size(tempName,1)==1
    mainNames(i)=tempName;
    mainGroups(i) = allPointNames.group(tempRows);
    elseif size(tempName,1)==0
        mainNames(i)='NOT FOUND';
        mainGroups(i) = NaN;
        
    else
        mainNames(i)='FOUND MORE THAN 1';
        mainGroups(i) = NaN;
        
    end


end

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
