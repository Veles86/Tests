function [  ] = WriteREZFile( folder,output_file_name,rez_table )

fileID = fopen(fullfile(folder,output_file_name),'w');

n = size(rez_table,1);
formatSpec = '%4s%-4s %4s%-4s %11.5f %6g. %4g %7.2f  %04g %+24s %+24s\n';

% formatSpec = '%4s%-4s %4s%-4s %11.5f %6g. %4g %7.2f  %04g %+24s\n';
fprintf(fileID,'Point_A   Point_B   Height_Difference    Distance    Num_Across Difference_Between_BF Date_Measured       Name_File     Source\n');

% fprintf(fileID,'Point_A   Point_B   Height_Difference    Distance    Num_Across Difference_Between_BF Date_Measured       Name_File\n');
for i = 1:n
    tempRow = rez_table(i,:);
    
    
    [fromNum,fromLetter] = SplitPointName(string(tempRow.Point_A));
    [toNum,toLetter] = SplitPointName(string(tempRow.Point_B));
    
    
    fprintf(fileID,formatSpec,fromNum,fromLetter,toNum,toLetter,tempRow.Height_Difference,tempRow.Distance,...
        tempRow.Num_Across,tempRow.Difference_Between_BF,tempRow.Date_Measured,string(tempRow.Name_File),string(tempRow.Source));
    
%     fprintf(fileID,formatSpec,fromNum,fromLetter,toNum,toLetter,tempRow.Height_Difference,tempRow.Distance,...
%         tempRow.Num_Across,tempRow.Difference_Between_BF,tempRow.Date_Measured,string(tempRow.Name_File));
    
    

end

fclose(fileID);