clc
clear all
close all

addpath('functions\')

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

[fileNameXYH,pathXYH] = uigetfile('*.csv','Select coordinates file',path);

inputCoordinates = readtable(strcat(pathXYH,fileNameXYH));

[fileNameG,pathG] = uigetfile('*.csv','Select gravity file',path);

inputG = readtable(strcat(pathG,fileNameG));
%%

full_coordinates = [inputCoordinates,inputG(:,"grav"),table(zeros(height(inputCoordinates),1),'VariableNames',"gG")];

full_coordinates.gG = CalcGGag(full_coordinates.H,full_coordinates.grav);

temp_table = table(zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),'VariableNames',["gA","gB","gG_A","gG_B","OC","dHO"]);
full_merged_file = [inputFile,temp_table];

[~,idxA] = ismember(full_merged_file.codeA,full_coordinates.code);
[~,idxB] = ismember(full_merged_file.codeB,full_coordinates.code);

full_merged_file.gA = full_coordinates{idxA,"grav"};
full_merged_file.gB = full_coordinates{idxB,"grav"};
full_merged_file.gG_A = full_coordinates{idxA,"gG"};
full_merged_file.gG_B = full_coordinates{idxB,"gG"};

full_merged_file.OC = (1./full_merged_file.gG_B).*(full_merged_file.H1.*(full_merged_file.gG_A-full_merged_file.gG_B)+((full_merged_file.gG_A+full_merged_file.gG_B)/2-full_merged_file.gG_B).*full_merged_file.dH);
full_merged_file.dHO = full_merged_file.dH - full_merged_file.OC;

writetable(full_merged_file,strcat('output_filling_grav_data/',fileName,'_wth_grav.csv'));
writetable(full_coordinates,strcat('output_filling_grav_data/',fileName,'_wth_grav_coordinates.csv'));

full_merged_file.WKT = strrep(full_merged_file.WKT,':',',');

writetable(full_merged_file,strcat('output_filling_grav_data/',fileName,'_wth_grav_WKT.csv'),'Delimiter',';');


function [gGag] = CalcGGag(H,g)

gGag = zeros(length(H),1);
posH = (H>0);
negH = (H<=0);


gGag(posH) = g(posH)+0.0424.*H(posH);
gGag(negH) = g(negH)+0.1543.*H(negH);




end