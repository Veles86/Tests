clc
clear all

addpath('functions\')



% coordinates = ReadCoordinatesFiles('E:\OneDrive\Research\matlab\coordinates\');
% 
% coordinates = coordinates(string(coordinates.type)~='duplicate',:);
% coordinates = coordinates(coordinates.group~=8,:);
% coordinates = coordinates(coordinates.group~=11,:);
% coordinates = coordinates(coordinates.group~=13,:);
opts = detectImportOptions('input\test\all_coordinates.csv');
opts.VariableTypes = {'double','char','char','double','double','double','double','double'};
coordinates = readtable('input\test\all_coordinates.csv',opts);

%%
tic
[phi,lambda] = IG12toIGD12(coordinates.Y,coordinates.X);
toc
%%
temp_table = table(phi,lambda,'VariableNames',["phi","lambda"]);

output_coordinates = [coordinates,temp_table];

writetable(output_coordinates,strcat('E:\OneDrive\Research\matlab\coordinates\Gravimetry\','GG_coordinates.csv'));