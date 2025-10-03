clc
clear
close all


param_file = 'E:\OneDrive\Research\Matlab\02 Fill Missing Data\data_ILUM2\ILUM2.0.par';
grid_file = 'E:\OneDrive\Research\Matlab\02 Fill Missing Data\data_ILUM2\ILUM2.0_Jan2012-UND-DEV.DAT';


%%
FID = fopen(param_file);
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});

%%

x_min = double(stringData(find(stringData=='ymin')+1));
x_max = double(stringData(find(stringData=='ymax')+1));
y_min = double(stringData(find(stringData=='xmin')+1));
y_max = double(stringData(find(stringData=='xmax')+1));
gridstep = double(stringData(find(stringData=='deltaX')+1));
nrows = double(stringData(find(stringData=='lines')+1));
ncols = double(stringData(find(stringData=='columns')+1));
gapvalue = 9999;

FID = fopen(grid_file);
data = textscan(FID,'%s');
fclose(FID);
stringData = string(data{:});

undulation = double(stringData(3:4:end));
stdv_undulation = double(stringData(4:4:end));

undulation = reshape(undulation,ncols,nrows)';
stdv_undulation = reshape(stdv_undulation,ncols,nrows)';

%%
WriteASC(undulation,x_min,y_min,gridstep,gapvalue,'ILUM2_IL05_undulation.asc');
WriteASC(stdv_undulation,x_min,y_min,gridstep,gapvalue,'ILUM2_IL05_stdev.asc');

% E = [x_min;x_max];
% N = [y_min;y_max];
% 
% [ phi, lambda ] = IG12toIGD12(E,N );


function [] = WriteASC(grid_data,xmin,ymin,cellsize,nodata_value,filename)

ncols = size(grid_data,2);
nrows = size(grid_data,1);

fileID = fopen(strcat('data_grav_rasters\',filename),'w');

if mod(xmin,1)==0 && mod(ymin,1)==0 && mod(cellsize,1)==0
    fprintf(fileID,'ncols         %d\n',ncols);
    fprintf(fileID,'nrows         %d\n',nrows);
    fprintf(fileID,'xllcorner     %d\n',xmin);
    fprintf(fileID,'yllcorner     %d\n',ymin);
    fprintf(fileID,'cellsize      %d\n',cellsize);
    fprintf(fileID,'nodata_value  %d\n',nodata_value);

else

    fprintf(fileID,'ncols         %d\n',ncols);
    fprintf(fileID,'nrows         %d\n',nrows);
    fprintf(fileID,'xllcorner     %.10f\n',xmin);
    fprintf(fileID,'yllcorner     %.10f\n',ymin);
    fprintf(fileID,'cellsize      %.7f\n',cellsize);
    fprintf(fileID,'nodata_value  %d\n',nodata_value);

end


for i = 1:nrows
    for j = 1:ncols
        if j~=ncols
            fprintf(fileID,'%.4f ',grid_data(i,j));
        else
            fprintf(fileID,'%.4f\n',grid_data(i,j));
        end
    end
end

fclose(fileID);


end








