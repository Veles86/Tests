clc
clear
close all

deltotTable = readtable('deltot_table.xlsx');

pointA = '40F';
pointB = '5090';

results = SearchConnections(deltotTable, pointA,pointB);

function [ found ] = SearchConnections(table, p1, p2 )






end