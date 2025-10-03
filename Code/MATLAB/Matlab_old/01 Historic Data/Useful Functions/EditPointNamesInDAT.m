clc
clear all
close all

[fileName,path] = uigetfile('*.csv','Select csv file with changes to make','MultiSelect', 'off');
changes = readtable(strcat(path,fileName));
filename_changes = changes(string(changes.type)=='filename',:);
internal_changes = changes(string(changes.type)=='internal',:);
h1 = height(filename_changes);
h2 = height(internal_changes);

[fileName,path] = uigetfile('*.DAT','Select DAT files','MultiSelect', 'on',path);

mkdir(strcat(path,'output'));



if ~iscell(fileName)
    fileName = {fileName};
    n = 1;
else
    n = length(fileName);
end


%%
for i = 1:n

    newFileName = fileName{i};

    for j = 1 : h1
        newFileName = strrep(newFileName,filename_changes.from{j},filename_changes.to{j});
    end


    fid = fopen(strcat(path,fileName{i}),'r');
    f = fread(fid,'*char')';
    fclose(fid);

    for k = 1 : h2
        f(121:end) = strrep(f(121:end),char(strrep(internal_changes.from{k},'_',' ')),char(strrep(internal_changes.to{k},'_',' ')));
    end
    formatSpec = ' 1|TO  %-27s';
    f = char(regexprep(string(f),'\s1\|TO  \w*\-*\w*\.*\w*\s*',sprintf(formatSpec,newFileName)));
    fid = fopen(strcat(path,'output/',newFileName),'w');
    fprintf(fid,'%s',f);
    fclose(fid);

end