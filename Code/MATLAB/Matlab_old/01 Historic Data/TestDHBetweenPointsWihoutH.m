clc
clear all
close all

opts = detectImportOptions('coordinates\all_from_DB_english.csv');
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinates = readtable('coordinates\all_from_DB_english.csv',opts);
coordinates = [coordinates;readtable('coordinates\all_from_DB_hebrew.csv',opts)];
coordinates = [coordinates;readtable('coordinates\missing_points.csv',opts)];
coordinates = [coordinates;readtable('coordinates\temporary_missing_points.csv',opts)];
coordinates = [coordinates;readtable('coordinates\old_BM_not_in_DB.csv',opts)];
coordinates = [coordinates;readtable('coordinates\pickets.csv',opts)];
coordinates = [coordinates;readtable('coordinates\eccentric_points.csv',opts)];

%%
[fileName,path] = uigetfile('*.csv','Select Merged CSV file');

opts = detectImportOptions(strcat(path,fileName),'Delimiter',',');
opts.VariableTypes(2) = {'char'};
opts.VariableTypes(3) = {'char'};
opts.VariableTypes(13) = {'double'};
opts.VariableTypes(14) = {'double'};
opts.VariableTypes(17) = {'char'};
opts.VariableTypes(22) = {'char'};
inputFile = readtable(strcat(path,fileName),opts);
fileName = fileName(1:end-4);
n = height(inputFile);


allPoints = GetCoordinates(inputFile);
rows = isnan(allPoints.H);
knownHeights = allPoints(~rows,:);
% unknownHeights = allPoints(rows,:);


[~,ia] = unique(table(abs(inputFile.dH),inputFile(:,5:7)));
inputFile = inputFile(ia,:);
%%
rows_for_deletion = zeros(height(inputFile),1);
for i = 1: height(inputFile)
    if sum(string(knownHeights.name)==inputFile.A(i))>0 && sum(string(knownHeights.name)==inputFile.B(i))>0
        rows_for_deletion(i) = 1;
    end
end

inputFile(logical(rows_for_deletion),:)=[];
allPoints = GetCoordinates(inputFile);
rows = isnan(allPoints.H);
knownHeights = allPoints(~rows,:);
unknownHeights = allPoints(rows,:);

knownH_table = knownHeights(:,[3,8]);
rez_table = inputFile(:,2:7);
%%

[results, v,l_a,l_b] = Adjust(knownHeights(:,[3,8]), rez_table);










function [coordinatesTable] = GetCoordinates(InputMergedTable)

n = height(InputMergedTable);

coordinatesTable = table('Size',[2*n 8],'VariableTypes',{'double','string','string','double','double','double','double','double'},'VariableNames',["code","type","name","group","rank","Y","X","H"]);
coordinatesTable(1:n,:) = [InputMergedTable(:,15:19),InputMergedTable(:,25:27)];
coordinatesTable(n+1:end,:) = [InputMergedTable(:,20:24),InputMergedTable(:,28:30)];

rows = isnan(coordinatesTable.code);
coordinatesTable(rows,:) = [];

rows = isnan(coordinatesTable.rank);
coordinatesTable.rank(rows) = 9999;
rows = isnan(coordinatesTable.group);
coordinatesTable.group(rows) = 9999;
rows = isnan(coordinatesTable.H);
coordinatesTable.H(rows) = 9999;
coordinatesTable = unique(coordinatesTable,'rows');


rows = coordinatesTable.rank==9999;
coordinatesTable.rank(rows) = NaN;
rows = coordinatesTable.group==9999;
coordinatesTable.group(rows) = NaN;
rows = coordinatesTable.H==9999;
coordinatesTable.H(rows) = NaN;

coordinatesTable = sortrows(coordinatesTable,1);

end



function [output, v,l_a,l_b] = Adjust(knownH_table, deltaH_table)

% knownPoints = sort(unique(string(knownH_table{:,1})));
% variables = sort(unique(string([deltaH_table{:,1};deltaH_table{:,2}])));

knownPoints = (unique(string(knownH_table{:,1})));
variables = (unique(string([deltaH_table{:,1};deltaH_table{:,2}])));

for i=1:length(knownPoints)
   index = find(strcmpi(variables,knownPoints(i)));
   variables(index)=[]; 
end

%number of equations
n = size(deltaH_table,1);

%number of variables
u = length(variables);

%choose weigting parameter and build P
% sigma_apriori = 1;
WEIGHTING_PARAMETER = 4;

try
    
    if WEIGHTING_PARAMETER==1 %weighting by distance
        P = diag(1./((deltaH_table{:,4}).^2));
    elseif WEIGHTING_PARAMETER==2 %weighting by number of stations
        precision_per_station = 0.1; %mm
%       P = diag(1./(deltaH_table{:,5}));
        P = diag(1./(((precision_per_station/1000)*deltaH_table{:,5}).^2));
    elseif WEIGHTING_PARAMETER==3 %weigting by diff BF
        P = diag(1./((deltaH_table{:,6}./1000).^2));
    elseif WEIGHTING_PARAMETER==4 %custom weighting
        P = eye(n);
    else %without weighting
        P = eye(n);
    end
    
    
catch
    P = eye(n);
    
end
        






l_b = deltaH_table{:,3};

%build l_0
l_0 = zeros(n,1);

for i=1:size(knownH_table,1)
   fromPoint = strcmp(string(deltaH_table{:,1}),knownPoints(i));
   toPoint = strcmp(string(deltaH_table{:,2}),knownPoints(i));
   tempH = knownH_table{i,2};
   
   l_0 = l_0 + tempH.*(toPoint - fromPoint);
end


l = l_b - l_0;

%build A
A = zeros(n,u);
for i=1:length(variables)
    fromPoint = strcmp(string(deltaH_table{:,1}),variables(i));
    toPoint = strcmp(string(deltaH_table{:,2}),variables(i));
    
    A(:,i) = toPoint - fromPoint;
    
end

N_mat = A'*P*A;
U_mat = A'*P*l;

% X = N_mat\U_mat;
X=pinv(N_mat)*(U_mat);
l_a = A*X + l_0;

v = l_a - l_b;

sigma_apost = sqrt((v'*P*v)/(n-u));

% SigmaX = (sigma_apost^2)./N_mat;

SigmaX = (sigma_apost^2)*pinv(N_mat);

tolerance = 1.96; % 95%

output = table(variables,round(X,3),round(tolerance*sqrt(diag(SigmaX))*1000,2));

end