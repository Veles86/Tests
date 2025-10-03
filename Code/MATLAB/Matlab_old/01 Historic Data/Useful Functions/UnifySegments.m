clc
clear
close all




[fileName,path] = uigetfile('*.rez','Select REZ file','MultiSelect','on');

if ~iscell(fileName)
    fileName = {fileName};
    u=1;
else
    u = length(fileName);
end

for k = 1:u

%Insert names of pickets that are different than "PIK" names

pointsToUnify={'117J','120J','1624bb'}';
% pointsToUnify = {'10K','100GPIK'}';
opts = detectImportOptions(strcat(path,fileName{k}),'FileType','text');
opts.VariableTypes(1) = {'char'};
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(8) = {'char'};

inputFile = readtable(strcat(path,fileName{k}),opts);

rows = contains(string(inputFile.Point_A),'PIK');
pikNames = (inputFile.Point_A(rows));
rows = contains(string(inputFile.Point_B),'PIK');
pikNames = [pikNames;inputFile.Point_B(rows)];
pointsToUnify = [pointsToUnify;pikNames];
pointsToUnify = unique(pointsToUnify);

tempName = fileName{k};
tempName = tempName(1:end-4);
n = length(pointsToUnify);

for i=1:n
    h = height(inputFile);
    tempPoint = pointsToUnify{i};
    
    indexesA = string(inputFile.Point_A)==tempPoint;
    indexesB = string(inputFile.Point_B)==tempPoint;
    indexes_all = indexesA + indexesB;
    num_of_appearances = sum(indexes_all);
    dataRows = [inputFile(indexesB,:);SwapRows(inputFile(indexesA,:))];

    if num_of_appearances==0
        warning(strcat('No appearance found for point: ',tempPoint));
    elseif num_of_appearances==1
        warning(strcat('Only one appearance found for point: ',tempPoint));
    elseif num_of_appearances==2
        newRow = Unify2Rows(dataRows(1,:),dataRows(2,:));
        idx = find(indexes_all);
        inputFile(idx(1),:) = newRow;
        inputFile(idx(2),:) = [];
    else
        for a = 1: (num_of_appearances -1)
            for b = (a+1) :num_of_appearances
                newRow = Unify2Rows(dataRows(a,:),dataRows(b,:));
                h = h+1;
                inputFile(h,:) = newRow;
            end
        end
        idx = find(indexes_all);
        inputFile(idx,:) = [];

    end



end

WriteREZFile(path,strcat(tempName,'_wo_pickets','.rez'),inputFile );

end

function [swappedRows] = SwapRows(dataRows)

    swappedRows = dataRows;
    swappedRows.Point_A = dataRows.Point_B;
    swappedRows.Point_B = dataRows.Point_A;
    swappedRows.Height_Difference = -dataRows.Height_Difference;


end

function [date] = FindLateDate(date1,date2)

if date1==date2
    date = date1;
elseif mod(date1,100)>mod(date2,100)
    date = date1;
elseif mod(date1,100)<mod(date2,100)
    date = date2;
elseif round(date1/100,0)>round(date2/100,0)
    date = date1;
else
    date = date2;
end


end

function [unifiedRow] = Unify2Rows(dataRow1,dataRow2)

if string(dataRow1.Point_B)==string(dataRow2.Point_B)
    dataRow2 = SwapRows(dataRow2);
elseif string(dataRow1.Point_A)==string(dataRow2.Point_B)
    temp = dataRow2;
    dataRow2 = dataRow1;
    dataRow1 = temp;
elseif string(dataRow1.Point_A)==string(dataRow2.Point_A)
    dataRow1 = SwapRows(dataRow1);   
end

unifiedRow = dataRow1;
unifiedRow.Point_B = dataRow2.Point_B;
unifiedRow.Height_Difference = dataRow1.Height_Difference + dataRow2.Height_Difference;
unifiedRow.Distance = dataRow1.Distance + dataRow2.Distance;
unifiedRow.Num_Across = dataRow1.Num_Across + dataRow2.Num_Across;
unifiedRow.Difference_Between_BF = dataRow1.Difference_Between_BF + dataRow2.Difference_Between_BF;
unifiedRow.Date_Measured = FindLateDate(dataRow1.Date_Measured,dataRow2.Date_Measured);

end


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