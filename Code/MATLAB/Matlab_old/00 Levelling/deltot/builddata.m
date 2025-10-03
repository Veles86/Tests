clc
clear
close all

% lines = readtable('deltot_table.xlsx');
% points = readtable('ORT-points.xlsx');

load('lines.mat')
load('points.mat')

indices = zeros(size(lines,1),2);
from_X = zeros(size(lines,1),1);
from_Y = zeros(size(lines,1),1);
from_H = zeros(size(lines,1),1);
to_X = zeros(size(lines,1),1);
to_Y = zeros(size(lines,1),1);
to_H = zeros(size(lines,1),1);
indices = zeros(size(lines,1),2);
status = zeros(size(lines,1),1);

for i=1:size(lines,1)
  fromP = string(lines.A(i));  
  toP = string(lines.B(i));  
  dist = lines.dist(i);  
    
  indexFrom = find(strcmp(string(points.name),fromP)==1);  
  indexTo =   find(strcmp(string(points.name),toP)==1);
  
  if size(indexFrom,1)==1
      indices(i,1)=indexFrom;
      from_X(i) = points.X(indexFrom);
      from_Y(i) = points.Y(indexFrom);
      from_H(i) = points.H(indexFrom);
  elseif size(indexFrom,1)>1
      indices(i,1)=-1;
  end
  
  if size(indexTo,1)==1
      indices(i,2)=indexTo;
      to_X(i) = points.X(indexTo);
      to_Y(i) = points.Y(indexTo);
      to_H(i) = points.H(indexTo);
  elseif size(indexTo,1)>1
      indices(i,2)=-1;
  end
  
  if size(indexFrom,1)==1 && size(indexTo,1)==1
     d1 = str2double(lines.dist{i});
     d2 = sqrt((to_X(i)-from_X(i))^2+(to_Y(i)-from_Y(i))^2);
      
     status(i) = abs(d1-d2); 
  end
    
end

output = table(lines.A,lines.B,lines.dH,lines.dist,lines.n,lines.disclosure,lines.month,lines.year,from_X,from_Y,to_X,to_Y,status);

