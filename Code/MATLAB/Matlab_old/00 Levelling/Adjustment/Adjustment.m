clc
clear
close all

% %select input data for known H
% [fileName,path] = uigetfile({'*.csv';'*.*'}, 'Select known data','MultiSelect', 'off');
% knownH_table = sortrows(readtable(fullfile(path,fileName)));
% 
% %select input data for dH
% [fileName,path] = uigetfile({'*.csv';'*.*'}, 'Select measured data','MultiSelect', 'off');
% deltaH_table = readtable(fullfile(path,fileName));


fileName_H = 'recent_input_H.csv';
fileName_dH = 'recent_input_dH.csv';
opts = detectImportOptions(fileName_dH);
opts.VariableTypes(1) = {'char'};
opts.VariableTypes(2) = {'char'};
deltaH_table = readtable(fileName_dH,opts);
opts = detectImportOptions(fileName_H);
opts.VariableTypes(1) = {'char'};
knownH_table = sortrows(readtable(fileName_H,opts));


knownPoints = sort(unique(string(knownH_table{:,1})));
variables = sort(unique(string([deltaH_table{:,1};deltaH_table{:,2}])));

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

X = N_mat\U_mat;

l_a = A*X + l_0;

v = l_a - l_b;

sigma_apost = sqrt((v'*P*v)/(n-u));

% SigmaX = (sigma_apost^2)./N_mat;

SigmaX = (sigma_apost^2)*inv(N_mat);

tolerance = 1.96; % 95%

output = table(variables,round(X,3),round(tolerance*sqrt(diag(SigmaX))*1000,2));