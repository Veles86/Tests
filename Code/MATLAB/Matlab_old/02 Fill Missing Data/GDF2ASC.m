clc
clear
close all


[fileName,path] = uigetfile('*.gdf','Select gdf file');
% fileName = 'E:\OneDrive\Research\Gravimetry and Stability\ICGEM\EIGEN-6C4\EIGEN-6C4_gravity_earth.gdf';


FID = fopen(strcat(path,fileName));
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});
fileName = fileName(1:end-4);
%%
modelName = stringData(find(stringData=='modelname')+1);
functional = stringData(find(stringData=='functional')+1);
x_min = double(stringData(find(stringData=='longlimit_west')+1));
x_max = double(stringData(find(stringData=='longlimit_east')+1));
y_min = double(stringData(find(stringData=='latlimit_south')+1));
y_max = double(stringData(find(stringData=='latlimit_north')+1));
gridstep = double(stringData(find(stringData=='gridstep')+1));
nrows = double(stringData(find(stringData=='latitude_parallels')+1));
ncols = double(stringData(find(stringData=='longitude_parallels')+1));
gapvalue = double(stringData(find(stringData=='gapvalue')+1));

startInd = find(stringData=='end_of_head')+2;

height_data = [];
%%
if stringData(startInd+1)==stringData(startInd+5) % 4 columns of data
    main_data = double(stringData(startInd+3:4:end));
    height_data = double(stringData(startInd+2:4:end));

elseif stringData(startInd+1)==stringData(startInd+4) % 3 columns of data

    main_data = double(stringData(startInd+2:3:end));

else

    warning('unrecognized format');

end
%%
main_data = reshape(main_data,ncols,nrows)';

if ~isempty(height_data)
height_data = reshape(height_data,ncols,nrows)';

end
%%
WriteASC(main_data,height_data,x_min,y_min,gridstep,gapvalue,modelName,functional)

function [] = WriteASC(grid_data,height_data,xmin,ymin,cellsize,nodata_value,modelname,functional)

ncols = size(grid_data,2);
nrows = size(grid_data,1);

fileID = fopen(strcat('data_grav_rasters\',modelname,'_',functional,'.asc'),'w');
fprintf(fileID,'ncols         %d\n',ncols);
fprintf(fileID,'nrows         %d\n',nrows);
fprintf(fileID,'xllcorner     %.10f\n',xmin);
fprintf(fileID,'yllcorner     %.10f\n',ymin);
fprintf(fileID,'cellsize      %.7f\n',cellsize);
fprintf(fileID,'nodata_value  %d\n',nodata_value);



for i = 1:nrows
    for j = 1:ncols
        if j~=ncols
            fprintf(fileID,'%.3f ',grid_data(i,j));
        else
            fprintf(fileID,'%.3f\n',grid_data(i,j));
        end
    end
end

fclose(fileID);

if ~isempty(height_data)
    fileID = fopen(strcat('data_grav_rasters\',modelname,'_',functional,'_heights','.asc'),'w');

    for i = 1:nrows
        for j = 1:ncols
            if j~=ncols
                fprintf(fileID,'%.3f ',height_data(i,j));
            else
                fprintf(fileID,'%.3f\n',height_data(i,j));
            end
        end
    end


    fclose(fileID);

end


end








