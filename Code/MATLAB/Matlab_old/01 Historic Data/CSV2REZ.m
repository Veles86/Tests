clc
clear
close all

[fileName,path] = uigetfile('*.csv','Select CSV file');

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

inputFile = inputFile(:,2:9);
inputFile.Properties.VariableNames = {'Point_A','Point_B','Height_Difference','Distance','Num_Across','Difference_Between_BF','Date_Measured','Name_File'};

% 
% % defining output folder's name
% formatOut = 'yyyy-mm-dd_HH-MM-SS';
% outputFolder = strcat('output_',datestr(now,formatOut));
% mkdir(outputFolder);


WriteREZFile('output_rez',strcat(fileName,'.rez'),inputFile );
%%
[fileNameXYH,path] = uigetfile('*.csv','Select coordinates file',path);

inputFileCoordinates = readtable(strcat(path,fileNameXYH));
h = height(inputFileCoordinates);


outputCoordinates = [inputFileCoordinates(:,[3 6 7 8]),table(repmat("",h,1)),inputFileCoordinates(:,5)];
outputCoordinates.Properties.VariableNames = {'POINT_NAME','X','Y','HEIGHT','gGal','RANK'};



writetable(outputCoordinates,strcat('output_rez\',fileName,'_coordinates.csv'));



function [ number, letter ] = SplitPointName( fullPointName )

if ~isa(fullPointName,'char')
    fullPointName = char(fullPointName);
end

ascii = double(fullPointName);
digits = (ascii > 47 & ascii < 58);

if (sum(digits)==length(digits))
    if length(digits)<=4
        number = string(fullPointName);
        letter = string('');
    else
        number = string(fullPointName(1:4));
        letter = string(fullPointName(5:end));
    end
elseif sum(digits)==0
    number = string('');
    letter = string(fullPointName);
else
    index = find(digits==0,1);
    if index<=5
        number = string((fullPointName(1:index-1)));
        letter = string(fullPointName(index:end));
    else
        number = string((fullPointName(1:4)));
        letter = string(fullPointName(5:end));
    end
end

end


function [  ] = WriteREZFile( folder,output_file_name,rez_table )

fileID = fopen(fullfile(folder,output_file_name),'w');

n = size(rez_table,1);
%formatSpec = '%4s%-4s %4s%-4s %11.5f %5g. %3g %6.2f  %04g %+24s %+24s\n';

 formatSpec = '%4s%-4s %4s%-4s %11.5f %5g. %3g %6.2f  %04g %+24s\n';
%fprintf(fileID,'Point_A   Point_B   Height_Difference    Distance    Num_Across Difference_Between_BF Date_Measured       Name_File     Source\n');

fprintf(fileID,'Point_A   Point_B   Height_Difference    Distance    Num_Across Difference_Between_BF Date_Measured       Name_File\n');
for i = 1:n
    tempRow = rez_table(i,:);
    
    
    [fromNum,fromLetter] = SplitPointName(string(tempRow.Point_A));
    [toNum,toLetter] = SplitPointName(string(tempRow.Point_B));
    
    
%     fprintf(fileID,formatSpec,fromNum,fromLetter,toNum,toLetter,tempRow.Height_Difference,tempRow.Distance,...
%         tempRow.Num_Across,tempRow.Difference_Between_BF,tempRow.Date_Measured,string(tempRow.Name_File),string(tempRow.Source));
    
    fprintf(fileID,formatSpec,fromNum,fromLetter,toNum,toLetter,tempRow.Height_Difference,tempRow.Distance,...
        tempRow.Num_Across,tempRow.Difference_Between_BF,tempRow.Date_Measured,string(tempRow.Name_File));
    
    

end

fclose(fileID);
end