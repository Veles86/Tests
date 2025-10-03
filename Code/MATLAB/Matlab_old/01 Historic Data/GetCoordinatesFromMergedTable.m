clc
clear all
close all


[fileName,path] = uigetfile({'*.csv';'*.xlsx'},'Select Merged file');

temp = strsplit(fileName,'.');

if temp{end}=="csv"

opts = detectImportOptions(strcat(path,fileName),'Delimiter',',');

else

opts = detectImportOptions(strcat(path,fileName));
end

opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(9) = {'char'};
opts.VariableTypes(10) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
opts.VariableTypes(17) = {'char'};
opts.VariableTypes(22) = {'char'};
inputFile = readtable(strcat(path,fileName),opts);
fileName = temp{1};
n = height(inputFile);


coordinatesTable = GetCoordinates(inputFile);
writetable(coordinatesTable,strcat(path,fileName,'_coordinates.csv'));




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



