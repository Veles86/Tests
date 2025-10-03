clc
clear all
close all

changes = readtable('changes.csv');
h = height(changes);
fid = fopen('1546JAD-96J.DAT','r');

f = fread(fid,inf,'*char');

fclose(fid);

% f = strrep(f',char('  21J '),char('1521U '));

for i = 1: h
f = strrep(f,char(strrep(changes.from{i},'_',' ')),char(strrep(changes.to{i},'_',' ')));
end



fid = fopen('1546JAD-96J_out.DAT','w');
fprintf(fid,'%s',f);
fclose(fid);