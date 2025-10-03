clc
clear
close all

% The script updates coordinates in the DB file

opts = detectImportOptions('..\coordinates\all_from_DB_english.csv');
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinatesEN = readtable('..\coordinates\all_from_DB_english.csv',opts);
coordinatesHE = readtable('..\coordinates\all_from_DB_hebrew.csv',opts);
coordinatesUpdated = readtable('E:\OneDrive\Research\Matlab\01 Historic Data\coordinates\backup_before_update_2024-03-08\updated_points.csv',opts);


[tf_EN,idx_upd_1] = ismember(coordinatesEN.code,coordinatesUpdated.code);
[tf_HE,idx_upd_2] = ismember(coordinatesHE.code,coordinatesUpdated.code);

indicesEN = find(idx_upd_1);
indicesHE = find(idx_upd_2);


for i = 1 : length(indicesEN)
coordinatesEN(indicesEN(i),7:11)=coordinatesUpdated(idx_upd_1(indicesEN(i)),7:11);

% A = coordinatesEN(indicesEN(i),7:11);
% B = coordinatesUpdated(idx_upd_1(indicesEN(i)),7:11);


end


for i = 1 : length(indicesHE)
coordinatesHE(indicesHE(i),7:11)=coordinatesUpdated(idx_upd_2(indicesHE(i)),7:11);


end

writetable(coordinatesEN,strcat('..\coordinates\','updated_EN.csv'));
writetable(coordinatesHE,strcat('..\coordinates\','updated_HE.csv'));

% after updating coordinates, replace NaN values with empty char and copy the notes column for hebrew points
