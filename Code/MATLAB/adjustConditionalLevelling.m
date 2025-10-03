function results = adjustConditionalLevelling(observations, known_points, conditions)
% adjustConditionalLeveling performs a full conditional least-squares adjustment.
% It handles both closed loops and open traverses between known control points.
%
% INPUT:
%   observations - [n x 5] matrix: [Obs_ID, From_ID, To_ID, dH, Weight]
%   known_points - [k x 2] matrix: [Point_ID, Height]. Can be empty.
%   conditions   - A cell array where each cell defines a condition. A condition
%                  is a struct with two fields:
%                  .type: 'loop' or 'traverse'
%                  .path: A vector of observation IDs. Use negative ID for
%                         reverse direction.
%                  For a 'traverse', the path MUST start and end at points
%                  that are defined in 'known_points'.
%
% OUTPUT:
%   results - A struct containing the adjustment results.

%% 1. Data Preparation
num_obs = size(observations, 1);
num_cond = length(conditions);
dof = num_cond; % Degrees of freedom = number of conditions

L = observations(:,4); % Vector of observations
P = diag(observations(:,5)); % Weight matrix
Q = inv(P); % Cofactor matrix

% Create maps for easy lookup
obs_map = containers.Map(observations(:,1), 1:num_obs);
obs_info = containers.Map('KeyType','double','ValueType','any');
for i=1:num_obs
    obs_info(observations(i,1)) = struct('from', observations(i,2), 'to', observations(i,3));
end
known_heights = containers.Map(known_points(:,1), known_points(:,2));

%% 2. Build Condition Matrix (B) and Misclosure Vector (w)
B = zeros(num_cond, num_obs);
w = zeros(num_cond, 1);

for i = 1:num_cond
    cond = conditions{i};
    path = cond.path;
    
    % Build the i-th row of the B matrix
    for j = 1:length(path)
        obs_id = abs(path(j));
        obs_idx = obs_map(obs_id);
        B(i, obs_idx) = sign(path(j));
    end
    
    % Calculate the i-th misclosure w(i)
    measured_sum = B(i, :) * L;
    
    if strcmpi(cond.type, 'loop')
        % For a closed loop, misclosure is the sum itself
        w(i) = -measured_sum;
        
    elseif strcmpi(cond.type, 'traverse')
        % For a traverse, misclosure is (H_known_end - H_known_start) - measured_sum
        % Find the start and end points of the traverse from the path
        start_obs_info = obs_info(abs(path(1)));
        end_obs_info = obs_info(abs(path(end)));
        
        % This logic assumes a simple, non-branching path for finding endpoints
        % A more robust implementation might require a full path reconstruction.
        % For this example, we trace the path to find the absolute start/end points.
        [start_point, end_point] = trace_path(path, obs_info);
        
        if ~isKey(known_heights, start_point) || ~isKey(known_heights, end_point)
            error('Traverse path must start and end at known points. Check condition %d.', i);
        end
        
        known_diff = known_heights(end_point) - known_heights(start_point);
        w(i) = known_diff - measured_sum;
        
    else
        error('Unknown condition type: %s', cond.type);
    end
end

%% 3. Solve the System
N_e = B * Q * B'; % Normal matrix for the correlates
k = N_e \ w; % Solve for correlates (Lagrange Multipliers)
v = Q * B' * k; % Solve for the residuals vector v

%% 4. Post-Adjustment Statistical Analysis
post_var_factor = (w' * k) / dof;
L_a = L + v; % Adjusted observations

%% 5. Package Results
results = struct();
results.AdjustedObservations = table(observations(:,1), L, v, L_a, 'VariableNames', {'Obs_ID', 'Measured_dH', 'Residual_v', 'Adjusted_dH'});
results.PostVarFactor = post_var_factor;
results.DOF = dof;
results.Misclosure_w = w;
results.Correlates_k = k;
end

function [start_node, end_node] = trace_path(path_vector, obs_info_map)
    % Helper function to trace a path of observations and find its absolute start and end nodes.
    % This is a simplified implementation.
    current_node = -1;
    
    % Find starting point
    first_obs_id = abs(path_vector(1));
    first_obs_dir = sign(path_vector(1));
    info = obs_info_map(first_obs_id);
    if first_obs_dir == 1
        current_node = info.from;
    else
        current_node = info.to;
    end
    start_node = current_node;

    % Trace the path
    for i = 1:length(path_vector)
        obs_id = abs(path_vector(i));
        direction = sign(path_vector(i));
        info = obs_info_map(obs_id);
        
        if direction == 1 && info.from == current_node
            current_node = info.to;
        elseif direction == -1 && info.to == current_node
            current_node = info.from;
        else
            % This indicates a broken path or a more complex topology than assumed.
            % For this example, we assume paths are continuous.
        end
    end
    end_node = current_node;
end