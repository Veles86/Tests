clc
clear all

addpath('functions\')

[fileNameXYH,path] = uigetfile('*.csv','Select coordinates file',path);
fileNameXYH = fileNameXYH(1:end-4);
inputFileCoordinates = readtable(strcat(path,fileNameXYH));

% coordinates = ReadCoordinatesFiles('E:\OneDrive\Research\matlab\coordinates\');
% 
% coordinates = coordinates(string(coordinates.type)~='duplicate',:);
% coordinates = coordinates(coordinates.group~=8,:);
% coordinates = coordinates(coordinates.group~=11,:);
% coordinates = coordinates(coordinates.group~=13,:);
% opts = detectImportOptions('input\test\all_coordinates.csv');
% opts.VariableTypes = {'double','char','char','double','double','double','double','double'};
% coordinates = readtable('input\test\all_coordinates.csv',opts);

%%
tic
[phi,lambda] = IG12toIGD12(inputFileCoordinates.Y,inputFileCoordinates.X);
toc
%%
temp_table = table(phi,lambda,'VariableNames',["phi","lambda"]);

output_coordinates = [inputFileCoordinates,temp_table];

writetable(output_coordinates,strcat('output_phi_lambda_coordinates_list\',fileNameXYH,'_GG.csv'));