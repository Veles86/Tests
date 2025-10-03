clc
clear all
close all



[fileName,path] = uigetfile('*.csv','Select Merged CSV file','MultiSelect', 'on');



if ~iscell(fileName)
    fileName = {fileName};
    n = 1;
else
    n = length(fileName);
end

for i = 1:n
opts = detectImportOptions(strcat(path,fileName{i}),'Delimiter',',');
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(9) = {'char'};
opts.VariableTypes(10) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
opts.VariableTypes(17) = {'char'};
opts.VariableTypes(22) = {'char'};
inputFile = readtable(strcat(path,fileName{i}),opts);
tempFileName = fileName{i};
tempFileName = tempFileName(1:end-4);


inputFile.WKT = strrep(inputFile.WKT,':',',');

writetable(inputFile,strcat('output_WKT/',tempFileName,'_WKT.csv'),'Delimiter',';');

end