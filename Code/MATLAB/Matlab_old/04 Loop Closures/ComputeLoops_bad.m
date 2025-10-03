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
inputFile = inputFile(inputFile.year>1995,:);
% duplicates = FindDuplicates(inputFile);
fileName = fileName(1:end-4);

%%
G = graph(inputFile.A, inputFile.B, inputFile.dist);
%%
% tic
% [allCycles,allEdgeCycles] = allcycles(G,"MinCycleLength",50,"MaxNumCycles",1000000);
% toc
%%
[basisCycles,basisEdgeCycles] = cyclebasis(G);
n = height(basisCycles);

[cycles,edgecycles] = cyclebasis(G);
tiledlayout flow
for k = 1:length(cycles)
    nexttile
    highlight(plot(G),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
    title("Cycle " + k)
end

%%
loopData = cell(n,1);
for i = 1 : n
loopData{i} = GetLoopData(string(basisCycles{i}),inputFile);
end
%%
loopSummary = SummarizeLoops(loopData);

function [loopSummary] = SummarizeLoops(loopData)

k = height(loopData);

loopSummary = table('Size',[k 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',["number","num_of_points","length_km","misc_mm_wo_grav","misc_mm_w_grav"]);


for i = 1 : k

temp = loopData{i};
loopSummary.number(i) = i;
loopSummary.num_of_points(i) = height(temp);
loopSummary.length_km(i) = round(sum(temp.dist)/1000,1);
loopSummary.misc_mm_wo_grav(i) = round(sum(temp.dH)*1000,2);
loopSummary.misc_mm_w_grav(i) = 0;




end

end


% function [loopDataWithG] = AddGravData(loopData,gravData)
% 
% 
% 
% 
% 
% 
% 
% 
% end




function [duplicates] = FindDuplicates(inputFile)

duplicates = [];
num_of_duplicates = 0;

index = 1;

while index<=height(inputFile)

    temp1 = fullData((fullData.A==loop(i)) & ((fullData.B==loop(i+1))),:);
    temp2 = fullData((fullData.A==loop(i+1)) & ((fullData.B==loop(i))),:);

    temp = [temp1;InvertRows(temp2)];


end
end

function cycles = find_loops(graph)
  edge_list = cyclebasis(graph);
  cycles = [];
  for i = 1:length(edge_list)
    new_cycle = edge_list(i);
    for j = 1:length(cycles)
      if setxor(new_cycle, cycles(j)) == []
        break;
      end
    end
    if j == length(cycles)
      cycles = [cycles new_cycle];
    end
  end
end


