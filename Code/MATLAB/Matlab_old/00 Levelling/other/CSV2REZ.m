clc
clear
close all


input_rez_table = readtable('south2.csv');

% defining output folder's name
formatOut = 'yyyy-mm-dd_HH-MM-SS';
outputFolder = strcat('output_',datestr(now,formatOut));
mkdir(outputFolder);

WriteREZFile( outputFolder,'output.rez',input_rez_table );