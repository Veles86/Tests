
clc
clear 
close all

addpath 'E:\git\LevellingThesis\Code\MATLAB\functions'

load("input_adj2025_main_and_secondary.mat")

known_points = zeros(0,2);

%%
counter = 0;
max_counter = 500;
is_finished = false;
removed_lines = zeros(0,1);
final_redundancies = [];
while counter <=max_counter && ~is_finished

    counter = counter + 1;
    results = adjustLevellingNetwork(observations, known_points);

    line_idx = results.removeLine;
    
    if line_idx == 0
        is_finished = true;
        final_redundancies = [results.Statistics.Obs_ID,results.Statistics.Redundancy];
    else
        removed_lines(counter) = observations(line_idx,1);
        observations(line_idx, :) = [];
    end

end



