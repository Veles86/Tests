clc
clear 
close all

load("segments.mat")
load("junctions.mat")
load("regulars.mat")

results = sortSegments(segments,junctions,regulars);


results.numeration = (1:length(results.segment_index))';

[~, sortedIndices] = sort(results.segment_index);
results.line_number = results.line_number(sortedIndices);
results.numeration = results.numeration(sortedIndices);
results.reverse_flag = results.reverse_flag(sortedIndices);
results.segment_index = results.segment_index(sortedIndices);

numeration = results.numeration;
reverse_flag = results.reverse_flag';
segment_indexes = results.segment_index';
line_num = results.line_number';