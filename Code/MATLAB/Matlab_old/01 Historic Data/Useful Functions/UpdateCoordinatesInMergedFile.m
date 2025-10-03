clc
clear
close all

% The script updates coordinates in merged CSV file


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


opts = detectImportOptions('..\coordinates\all_from_DB_english.csv');
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinatesUpdated = readtable('..\coordinates\backup_before_update_2024-03-08\updated_points-all_till_now.csv',opts);

%%


originalInputFile = inputFile;

[tf_A,idx_upd_1] = ismember(inputFile.codeA,coordinatesUpdated.code);
[tf_B,idx_upd_2] = ismember(inputFile.codeB,coordinatesUpdated.code);

indicesA = find(idx_upd_1);
indicesB = find(idx_upd_2);


for i = 1 : length(indicesA)
inputFile(indicesA(i),25:27)=coordinatesUpdated(idx_upd_1(indicesA(i)),7:9);
inputFile(indicesA(i),19)=coordinatesUpdated(idx_upd_1(indicesA(i)),10);

end


for i = 1 : length(indicesB)
inputFile(indicesB(i),28:30)=coordinatesUpdated(idx_upd_2(indicesB(i)),7:9);
inputFile(indicesB(i),24)=coordinatesUpdated(idx_upd_2(indicesB(i)),10);
end
%%
inputFile.dH_from_DB = round(inputFile.H2 - inputFile.H1,3);
inputFile.dH_diff = round(abs(inputFile.dH_from_DB-inputFile.dH_from_file),5);

inputFile.dist_from_DB = round(sqrt((inputFile.X2-inputFile.X1).^2+(inputFile.Y2-inputFile.Y1).^2));
inputFile.dist_diff = abs(inputFile.dist_from_file-inputFile.dist_from_DB);

inputFile.WKT = strcat("LINESTRING(",string(inputFile.Y1)," ",string(inputFile.X1),":",string(inputFile.Y2)," ",string(inputFile.X2),")");



coordinates = GetCoordinatesFromMergedTable(inputFile);


writetable(inputFile,strcat(path,fileName,'_updated.csv'));
writetable(coordinates,strcat(path,fileName,'_updated_coordinates.csv'));


% Remove NaN values from both files at the end



function [coordinatesTable] = GetCoordinatesFromMergedTable(InputMergedTable)

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