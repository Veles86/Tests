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

% [fileNameXYH,path] = uigetfile('*.csv','Select coordinates file',path);
% inputFileCoordinates = readtable(strcat(path,fileNameXYH));

inputFileCoordinates = GetCoordinatesFromMergedTable(inputFile);
%%
[fileG,path] = uigetfile('*.csv','Select gravity data file',path);
inputGrav = readtable(strcat(path,fileG));

[~,ind] = ismember(inputFileCoordinates.code, inputGrav.code);
gTable = table(inputGrav.g(ind),'VariableNames',{'g'});
inputFileCoordinates = [inputFileCoordinates,gTable];
inputFileCoordinates.rank(inputFileCoordinates.rank==0)=NaN;
inputFileCoordinates.H(inputFileCoordinates.H==0)=NaN;
%%

EdgeTable = table([string(inputFile.A),string(inputFile.B)],inputFile.dist,ones(height(inputFile),1),'VariableNames',{'EndNodes','Weight','Direction'});
EdgeTable = [EdgeTable,inputFile];

% %ordering points in the order in which they appear in the edges table so
% %that the edges order will remain similar to the order from the input
% pointsOrder = unique([inputFile.codeA;inputFile.codeB],'stable');
% [~,ind] = ismember(pointsOrder, inputFileCoordinates.code);
% 
% NodeTable = table(inputFileCoordinates.name(ind),'VariableNames',{'Name'});
% NodeTable = [NodeTable,inputFileCoordinates(ind,:)];
% G = graph(EdgeTable,NodeTable);



NodeTable = table(inputFileCoordinates.name,'VariableNames',{'Name'});
NodeTable = [NodeTable,inputFileCoordinates];

G = graph(EdgeTable,NodeTable);




G = simplify(G);


% removes open edges
G=RemoveFreeEdges(G);
%marking rows in which matlab reversed the direction of the measured line
G=MarkDirection(G);


% freeEdges = 1;
% while freeEdges
%     x = degree(G)<2;
%     if sum(x)==0
%         freeEdges = false;
%     else
%         G=subgraph(G,(degree(G)>1));
%     end
% end

% edgesOrderInGraph = G.Edges.EndNodes;
% edgesOrderInTable = [G.Edges.A,G.Edges.B];
% idx = (string(edgesOrderInGraph(:,1))==string(edgesOrderInTable(:,1)) & string(edgesOrderInGraph(:,2))==string(edgesOrderInTable(:,2)));
% G.Edges.Direction(~idx)=-1;

% p = plot(G,'EdgeLabel',G.Edges.Weight,'XData',G.Nodes.Y,'YData',G.Nodes.X);

% p.XData = inputFileCoordinates.Y(ind);
% p.YData = inputFileCoordinates.X(ind);
%%

% [independentLoops,A,B,w] = FindIndependentLoops(G);
% independentGLoops = cell(length(independentLoops),1);
% removeIdx = zeros(length(independentLoops),1);
% additionalGLoops=[];
% additionalIdx=1;
% 
% for i = 1:length(independentLoops)
% 
% tempGLoop = subgraph(G,independentLoops{i});
% 
% if numnodes(tempGLoop)==numedges(tempGLoop)
% independentGLoops{i} = MarkDirection(tempGLoop);
% elseif  numnodes(tempGLoop)<numedges(tempGLoop)
% [tempIndependentLoops,~,~,~] = FindIndependentLoops(tempGLoop);
% 
% 
% 
% 
% for j=1:length(tempIndependentLoops)
% 
% tempGLoop2 = subgraph(tempGLoop,tempIndependentLoops{j});
% additionalGLoops{additionalIdx,1} = MarkDirection(tempGLoop2);
% additionalIdx = additionalIdx+1;
% 
% end
% 
% 
% removeIdx(i)=1;
% else
%     removeIdx(i)=1;
% end
% 
% 
% 
% end
% independentGLoops=independentGLoops(logical(~removeIdx));
% independentGLoops = [independentGLoops;additionalGLoops];

% figure;
% [T,pred] = minspantree(G);
% p = plot(G,'EdgeLabel',G.Edges.Weight,'XData',G.Nodes.Y,'YData',G.Nodes.X);
% highlight(p,T)
% 
% [~,ind] = ismember(T.Edges.num, G.Edges.num);
% S = rmedge(G,ind); %remained part besides the span tree
% % S = rmnode(S,find(degree(S)==0));
% 
% [~,ind] = ismember(T.Edges.num, G.Edges.num);
% S = rmedge(G,ind); %remained part besides the span tree
% % S = rmnode(S,find(degree(S)==0));
% 
% At = incidence(T)';
% At = At(:,1:end-1);
% Ac = incidence(S)';
% Ac = Ac(:,1:end-1);
% Af = [At;Ac];
% 
% Bt = -Ac/(At);
% Bf = [Bt,eye(size(Bt,1))];
% 
% edgesOrderInAf = [T.Edges.num;S.Edges.num];
% [~,ind] = ismember(G.Edges.num,edgesOrderInAf);
% 
% A = Af(ind,:);
% B = Bf(:,ind);
% 
% w = B*(G.Edges.dH.*G.Edges.Direction);
% n = length(w);

%%
% basisCycles = cell(n,1);
% basisCycles = independentLoops;
%     nProblems = 0;
% diff=[];
% 
% for i = 1 :n 
% 
% 
% % idx = logical(abs(B(i,:)));
% % points = unique([G.Edges.A(idx);G.Edges.B(idx)],'stable');
% % temp = subgraph(G,points);
% 
% temp = subgraph(G,basisCycles{i});
% if numnodes(temp)~=numedges(temp)
% 
% nProblems=nProblems+1;
% diff(nProblems)=numedges(temp)-numnodes(temp);
% 
% end
% edgesOrderInGraph = temp.Edges.EndNodes;
% edgesOrderInTable = [temp.Edges.A,temp.Edges.B];
% idx = (string(edgesOrderInGraph(:,1))==string(edgesOrderInTable(:,1)) & string(edgesOrderInGraph(:,2))==string(edgesOrderInTable(:,2)));
% temp.Edges.Direction(~idx)=-1;
% basisCycles{i} = temp;
% 
% 
% end
%
%%
[cycles,edgecycles] = cyclebasis(G);
% tiledlayout flow
% for k = 1:length(cycles)
%     nexttile
%     highlight(plot(G),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
%     title("Cycle " + k)
% end

independentGLoops2 = cell(length(cycles),1);
for i=1:length(cycles)
% independentGLoops2{i} = subgraph(G,cycles{i});
[~,ind] = ismember(cycles{i}, G.Nodes.Name);
independentGLoops2{i} = graph(G.Edges(edgecycles{i},:),G.Nodes(ind,:));
end
loopSummary2 = SummarizeLoops(independentGLoops2);

% tiledlayout flow
% for k = 1:length(cycles)
%     nexttile
%     highlight(plot(G),cycles{k},'Edges',edgecycles{k},'EdgeColor','r','NodeColor','r')
%     title("Cycle " + k)
% end

%%





% loopSummary = SummarizeLoops(independentGLoops);




function H = RemoveFreeEdges(G)
freeEdges = 1;
while freeEdges
    x = degree(G)<2;
    if sum(x)==0
        freeEdges = false;
    else
        G=subgraph(G,(degree(G)>1));
    end
end
H=G;
end

function H = MarkDirection(G)
edgesOrderInGraph = G.Edges.EndNodes;
edgesOrderInTable = [G.Edges.A,G.Edges.B];
idx = (string(edgesOrderInGraph(:,1))==string(edgesOrderInTable(:,1)) & string(edgesOrderInGraph(:,2))==string(edgesOrderInTable(:,2)));
G.Edges.Direction(~idx)=-1;
G.Edges.Direction(idx)=1;
H=G;
end


function [independentLoops,A,B,w] = FindIndependentLoops(G)

figure;
[T,~] = minspantree(G);
p = plot(G,'EdgeLabel',G.Edges.Weight,'XData',G.Nodes.Y,'YData',G.Nodes.X);
highlight(p,T)

[~,ind] = ismember(T.Edges.num, G.Edges.num);
S = rmedge(G,ind); %remained part besides the span tree
% S = rmnode(S,find(degree(S)==0));

[~,ind] = ismember(T.Edges.num, G.Edges.num);
S = rmedge(G,ind); %remained part besides the span tree
% S = rmnode(S,find(degree(S)==0));

At = incidence(T)';
At = At(:,1:end-1);
Ac = incidence(S)';
Ac = Ac(:,1:end-1);
Af = [At;Ac];

Bt = -Ac/(At);
Bf = [Bt,eye(size(Bt,1))];

edgesOrderInAf = [T.Edges.num;S.Edges.num];
[~,ind] = ismember(G.Edges.num,edgesOrderInAf);

A = Af(ind,:);
B = Bf(:,ind);

w = B*(G.Edges.dH.*G.Edges.Direction);

n = length(w);

independentLoops = cell(n,1);


for i = 1 :n 

idx = logical(abs(B(i,:)));
independentLoops{i} = unique([G.Edges.A(idx);G.Edges.B(idx)],'stable');

end
end

function [ Result ] = Factorial2( Value )
%Factorial2 - Calculates the value of n!
% Outputs the factorial value of the input number.
 if Value > 1
 Result = Factorial2(Value - 1) * Value;
 else
 Result = 1;
 end
end


function [loopSummary] = SummarizeLoops(loopData)

k = length(loopData);

loopSummary = table('Size',[k 5],'VariableTypes',{'double','double','double','double','double'},'VariableNames',["number","num_of_points","length_km","misc_mm_wo_grav","misc_mm_w_grav"]);


for i = 1 : k

temp = loopData{i};
[~, edgesOrder] = cyclebasis(temp);


edgesOrder = cell2mat(edgesOrder);


if temp.Edges{1,"Direction"} == 1
lastPoint = string(temp.Edges{1,"B"});
else
    lastPoint = string(temp.Edges{1,"A"});
end

% temp.Edges{1,"Direction"} = 1;

for j = 2 : height(temp.Edges)



        if (temp.Edges.A{edgesOrder(j)}==lastPoint)
            temp.Edges{edgesOrder(j),"Direction"} = 1;
            lastPoint = string(temp.Edges.B{edgesOrder(j)});
        else
            temp.Edges{edgesOrder(j),"Direction"} = -1;
            lastPoint = string(temp.Edges.A{edgesOrder(j)});
        end


end


loopSummary.number(i) = i;
loopSummary.num_of_points(i) = height(temp.Edges);
loopSummary.length_km(i) = round(sum(temp.Edges.dist)/1000,1);
loopSummary.misc_mm_wo_grav(i) = round(sum(temp.Edges.dH.*temp.Edges.Direction)*1000,2);
loopSummary.misc_mm_w_grav(i) = 0;




end

end

% function [orderedGLoop] = OrderGLoop(gLoop)
% 
% n = size(gLoop.Edges,1);
% visitedIndexes = zeros(n,1);
% 
% finished = false;
% 
% currentIndex=1;
% visitedIndexes(currentIndex)=1;
% 
% while ~finished
% 
% if sum(visitedIndexes)==n
%     finished = true;
% else
% 
% end
% 
% end
% 
% orderedGLoop = gLoop;
% 
% end



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




% function [duplicates] = FindDuplicates(inputFile)
% 
% duplicates = [];
% num_of_duplicates = 0;
% 
% index = 1;
% 
% while index<=height(inputFile)
% 
%     temp1 = fullData((fullData.A==loop(i)) & ((fullData.B==loop(i+1))),:);
%     temp2 = fullData((fullData.A==loop(i+1)) & ((fullData.B==loop(i))),:);
% 
%     temp = [temp1;InvertRows(temp2)];
% 
% 
% end
% end




