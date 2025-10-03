function [coordinates] = ReadCoordinatesFiles(path)

if path(end)~='\'
    path = strcat(path,'\');
end

opts = detectImportOptions(strcat(path,'all_from_DB_english.csv'));
opts.VariableTypes = {'double','char','char','char','double','double','double','double','double','double','char','char'};
coordinates = readtable(strcat(path,'all_from_DB_english.csv'),opts);
coordinates = [coordinates;readtable(strcat(path,'all_from_DB_hebrew.csv'),opts)];
coordinates = [coordinates;readtable(strcat(path,'missing_points.csv'),opts)];
coordinates = [coordinates;readtable(strcat(path,'temporary_missing_points.csv'),opts)];
% coordinates = [coordinates;readtable(strcat(path,'old_BM_not_in_DB.csv'),opts)];
coordinates = [coordinates;readtable(strcat(path,'pickets.csv'),opts)];
coordinates = [coordinates;readtable(strcat(path,'eccentric_points.csv'),opts)];


end

