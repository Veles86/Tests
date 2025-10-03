function results = adjustLevellingNetwork(observations, known_points)
% adjustLevelingNetwork performs a least-squares adjustment of a 1D leveling network.
%
% Syntax:
%   results = adjustLevelingNetwork(observations, known_points)
%
% INPUT:
%   observations - An [n x 5] matrix with the observation data:
%                  Col 1: Unique Observation ID
%                  Col 2: Start Point ID (From)
%                  Col 3: End Point ID (To)
%                  Col 4: Measured Height Difference (dH)
%                  Col 5: Weight of the observation (p)
%
%   known_points - A [k x 2] matrix with the control point data:
%                  Col 1: Point ID
%                  Col 2: Height of the point
%                  If this matrix is empty ([]), a free-network adjustment is performed.
%
% OUTPUT:
%   results - A struct containing all adjustment results:
%             .AdjustedHeights    - Table of point IDs and their adjusted heights.
%             .PointStdDev        - Table of point IDs and their standard deviations.
%             .Statistics         - Table with detailed statistical analysis per observation
%                                   (Residuals, Redundancy, Normalized Residuals).
%             .PostVarFactor      - The a-posteriori variance factor (sigma-naught-squared).
%             .GlobalTestStat     - The Chi-Squared global test statistic.
%             .DOF                - Degrees of freedom.
%             .AdjustmentType     - The type of adjustment performed ('Constrained' or 'Free').
%             .PointCovMatrix     - The covariance matrix of the adjusted points.

%% 1. Prepare Data and Identify All Points
obs_ids = observations(:,1);
point_ids = unique([observations(:,2); observations(:,3); known_points(:,1)]);
num_points = length(point_ids);
num_obs = size(observations, 1);

% Create a map from point ID to matrix index (1, 2, 3...)
point_map = containers.Map(point_ids, 1:num_points);

%% 2. Build the Linear System (v = Ax - l)
% Assumption: Approximate heights (H0) are zero for all points.
A = zeros(num_obs, num_points);
l = zeros(num_obs, 1);
P = diag(observations(:,5)); % Diagonal weight matrix

for i = 1:num_obs
    from_idx = point_map(observations(i,2));
    to_idx = point_map(observations(i,3));
    
    A(i, from_idx) = -1;
    A(i, to_idx) = 1;
    
    % l = H0_j - H0_i - dH_measured. Since H0 = 0, l = dH_measured.
    l(i) = observations(i,4);
end

%% 3. Identify Adjustment Type and Solve
N = A' * P * A; % Normal equations matrix
t = A' * P * l;

if isempty(known_points)
    % ------------------
    % Case 1: Free Network
    % ------------------
    adjustment_type = 'Free';
    
    % N is singular. We use the Moore-Penrose pseudo-inverse to solve.
    % This solution minimizes the norm of the solution vector ||x||.
    x_hat = pinv(N) * t;
    
    % Covariance matrix for free network
    N_inv = pinv(N);
    
    % Degrees of freedom for a free network (1 datum defect)
    dof = num_obs - (num_points - 1);

else
    % ---------------------------------
    % Case 2: Constrained Network
    % ---------------------------------
    adjustment_type = 'Constrained';
    
    % Identify known and unknown points
    is_known = ismember(point_ids, known_points(:,1));
    known_indices = find(is_known);
    unknown_indices = find(~is_known);
    
    num_unknowns = length(unknown_indices);
    
    % Partition matrices into unknown and known parts
    A_u = A(:, unknown_indices);
    A_k = A(:, known_indices);
    
    % Create vector of known heights
    H_k = zeros(length(known_indices), 1);
    for i=1:length(known_indices)
        H_k(i) = known_points(known_points(:,1) == point_ids(known_indices(i)), 2);
    end
    
    % Adjust the 'l' vector
    l_adj = l - A_k * H_k;
    
    % Solve for unknown points only
    N_u = A_u' * P * A_u;
    t_u = A_u' * P * l_adj;
    
    x_u_hat = N_u \ t_u;
    
    % Assemble the full solution vector
    x_hat = zeros(num_points, 1);
    x_hat(unknown_indices) = x_u_hat;
    x_hat(known_indices) = H_k;
    
    % Assemble the full inverse normal matrix
    N_inv_u = inv(N_u);
    N_inv = zeros(num_points, num_points);
    N_inv(unknown_indices, unknown_indices) = N_inv_u;

    % Degrees of freedom for a constrained network
    dof = num_obs - num_unknowns;
end

%% 4. Post-Adjustment Statistical Analysis
% Residual vector
v = A * x_hat - l;

% A-posteriori variance factor
post_var_factor = (v' * P * v) / dof;
if post_var_factor < 0; post_var_factor = 0; end % Prevent small negative floating point results

% Covariance matrix of adjusted parameters
C_x = post_var_factor * N_inv;

% Standard deviations of adjusted points
point_std_dev = sqrt(abs(diag(C_x)));

% Covariance matrix of residuals (for data snooping)
Qv = inv(P) - A * N_inv * A';
Cv = post_var_factor * Qv;
std_dev_residuals = sqrt(abs(diag(Cv)));

% Redundancy numbers (internal reliability)
redundancy_numbers = diag(Qv * P);

% Normalized residuals (w-statistic) for Data Snooping
% Avoid division by zero for perfectly determined observations
norm_residuals = zeros(num_obs, 1);
non_zero_std_dev_idx = std_dev_residuals > 1e-9;
norm_residuals(non_zero_std_dev_idx) = v(non_zero_std_dev_idx) ./ std_dev_residuals(non_zero_std_dev_idx);

significance_level_baarda = 0.05;
probability_baarda = 1 - (significance_level_baarda/num_obs)/2;
critical_value_baarda = norminv(probability_baarda, 0, 1);


% Global Test Statistic (assuming a-priori variance is 1)
% Chi square
global_test_stat = dof * post_var_factor; 

significance_level_chi_square = 0.05;
critical_value_chi_square = chi2inv(1 - significance_level_chi_square, dof);


% Find the index of the maximum absolute normalized residual
[max_w, max_w_idx] = max(abs(norm_residuals));

if max_w > critical_value_baarda
    remove_line = max_w_idx;
else
    remove_line = 0;
end

%% 5. Package Results into Output Struct
results = struct();
results.AdjustmentType = adjustment_type;
results.AdjustedHeights = table(point_ids, x_hat, 'VariableNames', {'Point_ID', 'Adjusted_Height'});
results.PointStdDev = table(point_ids, point_std_dev, 'VariableNames', {'Point_ID', 'Std_Dev'});
results.Statistics = table(obs_ids, v, redundancy_numbers, norm_residuals, 'VariableNames', {'Obs_ID', 'Residual', 'Redundancy', 'Normalized_Residual_W'});
results.PostVarFactor = post_var_factor;
results.GlobalTestStat = global_test_stat;
results.DOF = dof;
results.PointCovMatrix = C_x;

results.removeLine = remove_line;

%fprintf('--- Adjustment Summary ---\n');
%fprintf('Adjustment Type: %s\n', adjustment_type);
%fprintf('Degrees of Freedom: %d\n', dof);
%fprintf('A-posteriori Variance Factor (Sigma Naught^2): %.6f\n', post_var_factor);
%disp(results.AdjustedHeights);
%fprintf('--- Statistical Analysis ---\n');
%disp(results.Statistics);

end