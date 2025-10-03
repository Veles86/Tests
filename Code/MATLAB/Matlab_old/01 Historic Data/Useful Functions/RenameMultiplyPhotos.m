clc
clear all
close all


[fileName,path] = uigetfile('*.csv','Select CSV file with lists of names');

pointsTable = readtable(strcat(path,fileName));

h = height(pointsTable);

for i=1:h

    movefile(strcat(path,pointsTable.old_name{i}),strcat(path,pointsTable.new_name{i}));


end