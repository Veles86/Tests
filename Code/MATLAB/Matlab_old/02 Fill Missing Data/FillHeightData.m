clc
clear
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
%%

RANK_LIMIT = 5;
knownSegments = (inputFile.rankA<RANK_LIMIT) & (inputFile.rankB<RANK_LIMIT);

unknownData = inputFile;
knownData = inputFile(knownSegments,:);
unknownData(knownSegments,:) = [];


knownHeights = GetCoordinatesFromMergedTable(unknownData);
knownHeights(~(knownHeights.rank<RANK_LIMIT),:)=[];

unknownData{~(unknownData.rankA<RANK_LIMIT),"H1"}=NaN;
unknownData{~(unknownData.rankB<RANK_LIMIT),"H2"}=NaN;
unknownData{~(unknownData.rankA<RANK_LIMIT),"rankA"}=99; %code for manually removed ranks
unknownData{~(unknownData.rankB<RANK_LIMIT),"rankB"}=99; %code for manually removed ranks


%%

pointsA = string(unknownData.A);
pointsB = string(unknownData.B);

allPoints = [pointsA;pointsB];

distinctPoints = unique(allPoints);

[tf,idx] = ismember(allPoints,distinctPoints);


counts = hist(idx,length(distinctPoints))';


edgePoints = distinctPoints(counts==1);
midPoints = distinctPoints(counts==2);
junctionPoints = distinctPoints (counts>2);

%%
M = ismember(unknownData.A,edgePoints);
N = ismember(unknownData.B,edgePoints);
A = isnan(unknownData.H1);
B = isnan(unknownData.H2);

%%

finished = false;
i=1;

while ~finished

    if i>height(knownHeights)
        finished = true;
    else
        tempCode = knownHeights{i,"code"};
        tempH = knownHeights{i,"H"};

        knownA = (unknownData.codeA==tempCode) & (isnan(unknownData.H2));
        knownB = (unknownData.codeB==tempCode) & (isnan(unknownData.H1));

        indA = find(knownA);
        indB = find(knownB);

        newKnownCodes=zeros(length(indA)+length(indB),1);
        newKnownNames=repmat("",length(indA)+length(indB),1);
        newKnownHeights=zeros(length(indA)+length(indB),1);

        for j = 1 : length(indA)
            newKnownCodes(j) = unknownData.codeB(indA(j));
            newKnownNames(j) = string(unknownData.B(indA(j)));
            newKnownHeights(j) = round(tempH + unknownData.dH(indA(j)),3);
        end


        for k = 1 :length(indB)
            newKnownCodes(length(indA)+k) = unknownData.codeA(indB(k));
            newKnownNames(length(indA)+k) = string(unknownData.A(indB(k)));
            newKnownHeights(length(indA)+k) = round(tempH - unknownData.dH(indB(k)),3);
        end

        %averages duplicates, if they exist
        if length(newKnownCodes)~=length(unique(newKnownCodes))
            uniqueCodes = unique(newKnownCodes);
            [~,idx] = ismember(newKnownCodes,uniqueCodes);
            counts = hist(idx,length(uniqueCodes))';
            duplicates = uniqueCodes(counts>1);

            for l = 1:length(duplicates)
                tempNewCode = duplicates(l);
                tempAvgH = round(mean(newKnownHeights(newKnownCodes==tempNewCode)),3);
                
                newKnownCodes(newKnownCodes==tempNewCode) = [];
                newKnownHeights(newKnownCodes==tempNewCode) = [];
                newKnownCodes(length(newKnownCodes)+1) = tempNewCode;
                newKnownHeights(length(newKnownCodes)+1) = tempAvgH;
            end

        end
        
        totalInd = zeros(height(unknownData),1);
        for p = 1:length(newKnownCodes)
            indAA = unknownData.codeA==newKnownCodes(p);
            indBB = unknownData.codeB==newKnownCodes(p);
            unknownData.H1(indAA) = newKnownHeights(p);
            unknownData.H2(indBB) = newKnownHeights(p);
            unknownData.rankA(indAA) = 99;
            unknownData.rankB(indBB) = 99;

            totalInd = totalInd + indAA + indBB;
        end

       tempCoordinates = GetCoordinatesFromMergedTable(unknownData(logical(totalInd),:));
       [~,tempIdx] = ismember(tempCoordinates.code,newKnownCodes);
       knownHeights = [knownHeights;tempCoordinates(logical(tempIdx),:)];

    end


    i = i+1;
end

unknownData.dH_from_DB = round(unknownData.H2-unknownData.H1,3);
unknownData.dH_diff = round(abs(unknownData.dH_from_file-unknownData.dH_from_DB),5);

fixedData = [knownData;unknownData];
fixedData = sortrows(fixedData,"num");
writetable(fixedData,strcat('output_filling_height_data/',fileName,'_filled_height.csv'));
writetable(GetCoordinatesFromMergedTable(fixedData),strcat('output_filling_height_data/',fileName,'_filled_height_coordinates.csv'));


fixedData.WKT = strrep(fixedData.WKT,':',',');

writetable(fixedData,strcat('output_filling_height_data/',fileName,'_filled_height_WKT.csv'),'Delimiter',';');

%%



% function [coordinatesTable] = GetCoordinatesFromMergedTable(InputMergedTable)
% 
% n = height(InputMergedTable);
% 
% coordinatesTable = table('Size',[2*n 8],'VariableTypes',{'double','string','string','double','double','double','double','double'},'VariableNames',["code","type","name","group","rank","Y","X","H"]);
% coordinatesTable(1:n,:) = [InputMergedTable(:,15:19),InputMergedTable(:,25:27)];
% coordinatesTable(n+1:end,:) = [InputMergedTable(:,20:24),InputMergedTable(:,28:30)];
% 
% rows = isnan(coordinatesTable.code);
% coordinatesTable(rows,:) = [];
% 
% rows = isnan(coordinatesTable.rank);
% coordinatesTable.rank(rows) = 9999;
% rows = isnan(coordinatesTable.group);
% coordinatesTable.group(rows) = 9999;
% rows = isnan(coordinatesTable.H);
% coordinatesTable.H(rows) = 9999;
% coordinatesTable = unique(coordinatesTable,'rows');
% 
% 
% rows = coordinatesTable.rank==9999;
% coordinatesTable.rank(rows) = NaN;
% rows = coordinatesTable.group==9999;
% coordinatesTable.group(rows) = NaN;
% rows = coordinatesTable.H==9999;
% coordinatesTable.H(rows) = NaN;
% 
% coordinatesTable = sortrows(coordinatesTable,1);
% 
% end