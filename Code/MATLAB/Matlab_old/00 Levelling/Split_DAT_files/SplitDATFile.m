clc
clear
close all

% edgePoints = string({'23TM','3J'})';

%for joining files - join them manually and then execute "splitting" with
%only 2 points - start and end ant it will calculate the final results
% edgePoints = string({'21TM','22J'})';
edgePoints = string({'3J','J03b'})';

numOfLevellingLines = length(edgePoints)-1;
outputFiles = cell(numOfLevellingLines,1);
indices = cell(length(edgePoints),1);
linesNames = strings(numOfLevellingLines,1);
[fileName,path] = uigetfile({'*.DAT';'*.*'}, 'Select Input Files','MultiSelect', 'off');

% defining output folder's name
formatOut = 'yyyy-mm-dd_HH-MM-SS';
outputFolder = strcat(fileName(1:end-4),'_splitted_',datestr(now,formatOut));
mkdir(outputFolder);

inputData=readtable(fullfile(path,fileName),'Delimiter','|','HeaderLines',0,'ReadVariableNames',false);
inputData = inputData(:,1:6);
inputData.Var4(strcmp(inputData.Var4,'')) = {'                      '};
inputData.Var5(strcmp(inputData.Var5,'')) = {'                      '};
inputData.Var6(strcmp(inputData.Var6,'')) = {'                      '};


endLine = find(contains(inputData.Var3,'End-Line'));

if length(endLine==1)
    
    inputData(endLine,:) = [];
end
startLocation = find(contains(inputData.Var3,'Start-Line'));

if length(startLocation)~=1
    
    error('more than 1 start line');
    
else
    
    header = inputData(1:startLocation,:);

    tempColumn = string(inputData.Var3(startLocation+1:end));
    
    for j=1:length(edgePoints)
        tempIndices = find(contains(tempColumn,strcat({' '},edgePoints(j),{' '})));
        if isempty(tempIndices)
            error('point not found')
        else
            indices{j} = tempIndices+startLocation;
        end
    end
    
    for i=1:numOfLevellingLines
        fromPointIndices = indices{i};
        toPointIndices = indices{i+1};
        
        if (length(fromPointIndices)==2) && (length(toPointIndices)==3)
            startLine = indices{i}(1);
            endLine = indices{i+1}(2);
            outputFiles{i} = [header;inputData(startLine:endLine,:)];
        elseif (((length(fromPointIndices)==3) && (length(toPointIndices)==3)) ||...
                ((length(fromPointIndices)==3) && (length(toPointIndices)==4)))
            startLine = indices{i}(2);
            endLine = indices{i+1}(2);
            outputFiles{i} = [header;inputData(startLine:endLine,:)];
        elseif (length(fromPointIndices)==2) && (length(toPointIndices)==4)
            outputFiles{i} = inputData(1:end-2,:);
        else
            error('wrong data')
        end
        
        linesNames(i) = strcat(edgePoints(i),'-',edgePoints(i+1),'.dat');
        %         outputFiles{i}{1,3} = {strcat('TO',{'  '},linesNames(i),{'                '})};
        
        outputFiles{i}{1,3} = {sprintf('TO  % -27s',linesNames(i))};
    end
    
    for i=1:numOfLevellingLines
        tempColumn1 = string(outputFiles{i}.Var4);
        tempColumn2 = string(outputFiles{i}.Var5);
        
        Rb_locations = find(contains(tempColumn1,'Rb'));
        Rf_locations = find(contains(tempColumn1,'Rf'));
        
        if length(Rb_locations)~=length(Rf_locations)
            error('wrong data')
        else
            
            Rb = (split(tempColumn1(Rb_locations)));
            
            if size(Rb,2)==1
                Rb = Rb';
            end
            
            Rb = str2double(Rb(:,2));
            
            
            Rf = (split(tempColumn1(Rf_locations)));
            if size(Rf,2)==1
                Rf = Rf';
            end
            
            Rf = str2double(Rf(:,2));
            
            
            HD_Rb = split(tempColumn2(Rb_locations));
            
            if size(HD_Rb,2)==1
                HD_Rb = HD_Rb';
            end
            
            HD_Rb = str2double(HD_Rb(:,2));
            
            
            HD_Rf = split(tempColumn2(Rf_locations));
            
            if size(HD_Rf,2)==1
                HD_Rf = HD_Rf';
            end
                      
            HD_Rf = str2double(HD_Rf(:,2));
            
            
            sh = sum(Rb)-sum(Rf);
            
            n = length(Rb_locations);
            Db = sum(HD_Rb);
            Df = sum(HD_Rf);
            
            lastLine = outputFiles{i}(end,:);
            
            lastIndex = split(lastLine{1,2});
            lastIndex = str2double(lastIndex(2));
            
            lastPoint = split(lastLine{1,3});
            lastPoint = lastPoint(2);
            
            lineNum = split(lastLine{1,3});
            lineNum = lineNum(4);
            
            lastZ = split(lastLine{1,6});
            lastZ = (lastZ(2));
            
           
            
            newLine1 = {'For M5',...
                sprintf('Adr% 6s',num2str(lastIndex+1)),...
                sprintf('KD1% 9s% 19s',lastPoint,lineNum),...
                sprintf('Sh% 15s m   ',num2str(sh)),...
                '                      ',...
                '                      '};
            
            
            
            newLine2 = {'For M5',...
                sprintf('Adr% 6s',num2str(lastIndex+2)),...
                sprintf('KD2% 9s% 9s% 10s',lastPoint,num2str(n),lineNum),...
                sprintf('Db% 15s m   ',num2str(round(Db,2))),...
                sprintf('Df% 15s m   ',num2str(round(Df,2))),...
                sprintf('Z% 16s m   ',lastZ)};
            
            outputFiles{i} = [outputFiles{i};newLine1;newLine2];
            
            writetable(outputFiles{i},char(strcat(outputFolder,'\',linesNames(i))),'Delimiter','|','WriteVariableNames',false);
       
        end
        
        
    end
    
    
end

    
    

