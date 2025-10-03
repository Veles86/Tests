clc
clear
close all

pointNames = string({'RFH95';'RFH96';'RFH97';'A';'B';'C';'D'});
knownH = [142.645,103.082,92.242,NaN,NaN,NaN,NaN]';
deltaH = [21.243,18.224,3.759,25.085,25.340,32.465,7.372,7.422,10.525,0.325]';
numOfStations = [15,13,12,17,18,22,13,12,11,12]';

A=[1 0 0 -1 0 0 0;
    0 -1 0 1 0 0 0;
    0 0 0 1 -1 0 0;
    1 0 0 0 -1 0 0;
    0 0 -1 0 1 0 0;
    1 0 0 0 0 -1 0;
    0 0 0 0 1 -1 0;
    0 0 0 0 0 1 -1;
    0 0 -1 0 0 0 1;
    0 1 0 0 0 0 -1];

numOfPoints = length(pointNames);
numOfLines = size(A,1);

nLoops = numOfLines - numOfPoints + 1;
nConditions = nLoops + sum(~isnan(knownH))-1;


