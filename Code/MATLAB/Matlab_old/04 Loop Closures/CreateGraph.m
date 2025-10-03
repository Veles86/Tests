clc
clear all
close all

% % open the file
% fileID = fopen('input.rez');
% 
% % read the data
% data = textscan(fileID, '%s %s %f %f %d %d');
% 
% % close the file
% fclose(fileID);
% 
% % Create a table from the cell array
% T = table(data{1},data{2},data{3},data{4},data{5},data{6},...
%     'VariableNames',{'A','B','dH','dist','n','diff'});

[fileName,path] = uigetfile('*.csv','Select Merged CSV file');

T = readtable(strcat(path,fileName));

% % create the graph
% G = digraph(T.A, T.B, T.dH);

G = graph(T.A, T.B, T.dist);

% plot the graph
figure;
plot(G);

[cycles,edgecycles] = allcycles(G);

p = plot(G,'EdgeLabel',1:numedges(G));
figure
tiledlayout flow
for k = 1:length(cycles)
    nexttile
    highlight(plot(G),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
    title("Cycle " + k)
end


[cycles,edgecycles] = cyclebasis(G);
tiledlayout flow
for k = 1:length(cycles)
    nexttile
    highlight(plot(G),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
    title("Cycle " + k)
end