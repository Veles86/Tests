clc
clear
close all


input_rez_table = readtable('lines_all_sources.csv');

original_input = input_rez_table;

fromPntList = string(input_rez_table.Point_A);
toPntList = string(input_rez_table.Point_B);

allPntsInstances = sort([fromPntList;toPntList]);

uniquePnts = unique(allPntsInstances);
instances = zeros(length(uniquePnts),1);
for i=1:length(uniquePnts)
    instances(i) = sum(strcmp(allPntsInstances,uniquePnts(i)));
    
end

instancesT = table(instances,uniquePnts);
junctionPoints = instancesT((instancesT.instances>2),:);
junctionPoints = string(junctionPoints{:,2});
deadEndPoints = instancesT((instancesT.instances==1),:);
deadEndPoints = string(deadEndPoints{:,2});


% defining output folder's name
formatOut = 'yyyy-mm-dd_HH-MM-SS';
outputFolder = strcat('output_',datestr(now,formatOut));
mkdir(outputFolder);

fileNum = 0;

for i=1:length(junctionPoints)
    startPoint =  junctionPoints(i);
    
    startIndices = find(input_rez_table.Point_A==startPoint);
    startLines = input_rez_table(startIndices,:);
    endIndices = find(input_rez_table.Point_B==startPoint);
    endLines = input_rez_table(endIndices,:);
    
    startLines = sortrows([startLines; InvertLines(endLines)],'Point_B');
    
    input_rez_table = DeleteRows(input_rez_table,[startIndices; endIndices]);
    
    skipIndex = 0;
    for j = 1 : size(startLines,1)
        if ~(j==skipIndex)
            tempFile = startLines(j,:);
            nextPoint = string(tempFile.Point_B(end));
            
            %         if (sum(junctionPoints==nextPoint) || sum(deadEndPoints==nextPoint)) %if the next point is a junction point or a dead end point
            %             fileNum = fileNum + 1;
            %             output_path = sprintf('output\\%03d_%s-%s.rez',fileNum,startPoint,nextPoint);
            %             WriteREZFile( output_path,tempFile );
            %         else
            
            closedLoop = false;
            
            while ~(sum(junctionPoints==nextPoint) || sum(deadEndPoints==nextPoint) || isempty(input_rez_table))
                index = find(input_rez_table.Point_A==nextPoint);
                tempLine = input_rez_table(index,:);
                
                
                
                if isempty(index)
                    index = find(input_rez_table.Point_B==nextPoint);
                    tempLine = InvertLines(input_rez_table(index,:));
                    if isempty(index)  %then it's closed loop which starts and ends at the same point
                        index = find(startLines.Point_B==nextPoint);
                        tempLine = InvertLines(startLines(index,:));
                        closedLoop = true;
                    end
                end
                
                if ~closedLoop
                    input_rez_table = DeleteRows(input_rez_table,index);
                else
                    skipIndex = index;
                    
                end
                
                tempFile = [tempFile;tempLine];
                nextPoint = string(tempFile.Point_B(end));
            end
            fileNum = fileNum + 1;
            
            fileName = sprintf('%04d_%s-%s.rez',fileNum,startPoint,nextPoint);
            WriteREZFile( outputFolder,fileName,tempFile );
            
        end
    end
    
end


function [ reversedTable ] = InvertLines( inputTable )

reversedTable = inputTable;

reversedTable.Height_Difference = -inputTable.Height_Difference;
reversedTable.Point_A = inputTable.Point_B;
reversedTable.Point_B = inputTable.Point_A;

end

function [ newTable ] = DeleteRows( inputTable, indices )

newTable = inputTable;
indices = sort(indices,'descend');
for i=1:size(indices,1)
    
    newTable(indices(i),:) = [];
    
end

end

