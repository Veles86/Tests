clc
clear
close all




[fileName,path] = uigetfile('*.xlsx','Select excel file');

opts = detectImportOptions(strcat(path,fileName));
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(9) = {'char'};
opts.VariableTypes(10) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
opts.VariableTypes(17) = {'char'};
opts.VariableTypes(22) = {'char'};
inputFile = readtable(strcat(path,fileName),opts);
fileName = fileName(1:end-5);
n = height(inputFile);
%%
counter = 1;

while counter<=height(inputFile)

    tempRow = inputFile(counter,:);

    rows = (inputFile.dH_abs == tempRow.dH_abs) & (inputFile.dist == tempRow.dist) & (inputFile.n == tempRow.n) & (inputFile.diff == tempRow.diff) & (inputFile.date == tempRow.date);
    
    if sum(rows)>1
        
        indexes = find(rows);
        if indexes(1)~=counter
            error('something wrong');
        end

        for i = 2 : length(indexes)
          
            inputFile{counter,end-7:end} = inputFile{counter,end-7:end} + inputFile{indexes(i),end-7:end};

            
        end
        rows(counter) = 0;
        inputFile(rows,:) = [];

    end

    counter = counter + 1;
end

writetable(inputFile,strcat(path,fileName,'_merged_duplicates.xlsx'));
