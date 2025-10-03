function results = adjustParametricLevelling(observations, known_points)
% adjustParametricLeveling performs a parametric least-squares adjustment.
% This function requires at least one known control point.
%
% INPUT:
%   observations - [n x 5] matrix: [Obs_ID, From_ID, To_ID, dH, Weight]
%   known_points - [k x 2] matrix: [Point_ID, Height]. CANNOT be empty.
%
% OUTPUT:
%   results - A struct containing the adjustment results.

%% 1. Input Validation and Data Preparation
if isempty(known_points)
    error('Parametric adjustment requires at least one known point. For free network, use a different function.');
end

obs_ids = observations(:,1);
point_ids = unique([observations(:,2); observations(:,3); known_points(:,1)]);
num_points = length(point_ids);
num_obs = size(observations, 1);
point_map = containers.Map(point_ids, 1:num_points);

%% 2. Build Linear System for Unknowns
% Identify known and unknown points
is_known = ismember(point_ids, known_points(:,1));
known_indices = find(is_known);
unknown_indices = find(~is_known);
num_unknowns = length(unknown_indices);

% Build full design matrix A and misclosure vector l
A = zeros(num_obs, num_points);
l = -observations(:,4); % Since H0=0, l = -dH
for i = 1:num_obs
    A(i, point_map(observations(i,2))) = -1;
    A(i, point_map(observations(i,3))) = 1;
end
P = diag(observations(:,5));

% Partition matrices for unknowns (u) and knowns (k)
A_u = A(:, unknown_indices);
A_k = A(:, known_indices);

% Build vector of known heights
H_k = zeros(length(known_indices), 1);
for i=1:length(known_indices)
    H_k(i) = known_points(known_points(:,1) == point_ids(known_indices(i)), 2);
end

% Adjust the 'l' vector based on known heights
l_adj = l - A_k * H_k;

%% 3. Solve the System
% Normal equations for the unknown parameters
N_u = A_u' * P * A_u;
t_u = A_u' * P * l_adj;
x_u_hat = N_u \ t_u; % Solve for corrections for unknowns

% Assemble the full vector of adjusted heights
x_hat = zeros(num_points, 1);
x_hat(unknown_indices) = x_u_hat;
x_hat(known_indices) = H_k;

%% 4. Post-Adjustment Statistical Analysis
v = A * x_hat - l; % Residuals
dof = num_obs - num_unknowns; % Degrees of freedom
post_var_factor = (v' * P * v) / dof;

% Covariance matrix for unknown parameters
C_x_u = post_var_factor * inv(N_u);
% Assemble full covariance matrix
C_x = zeros(num_points, num_points);
C_x(unknown_indices, unknown_indices) = C_x_u;
point_std_dev = sqrt(abs(diag(C_x)));

%% 5. Package Results
results = struct();
results.AdjustedHeights = table(point_ids, x_hat, 'VariableNames', {'Point_ID', 'Adjusted_Height'});
results.PointStdDev = table(point_ids, point_std_dev, 'VariableNames', {'Point_ID', 'Std_Dev'});
results.Residuals = table(obs_ids, v, 'VariableNames', {'Obs_ID', 'Residual'});
results.PostVarFactor = post_var_factor;
results.DOF = dof;

end