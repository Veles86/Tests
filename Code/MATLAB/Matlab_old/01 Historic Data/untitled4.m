clc
clear

% Open the file
fileID = fopen('input.rez');

% Read the data
data = textscan(fileID, '%s %s %f %f %d %d');

% Close the file
fclose(fileID);

% Create a table from the cell array
T = table(data{1},data{2},data{3},data{4},data{5},data{6},...
    'VariableNames',{'A','B','dH','dist','n','diff'});

% Create the graph
% G = digraph(T.A, T.B, T.dH);

G = digraph(T.A, T.B, T.dH);


% Create the adjacency matrix
adjacencyMatrix = full(adjacency(G));

% Call the DFS function
cycles = DFS(adjacencyMatrix);

% Display the cycles
for i = 1:length(cycles)
    disp(['Cycle ', num2str(i), ':']);
    disp(cycles{i});
end


function cycles = DFS(adjacencyMatrix)
    n = size(adjacencyMatrix, 1);
    color = zeros(1, n); % 0 = WHITE, 1 = GRAY, 2 = BLACK
    cycles = {};
    parent = zeros(1, n);
    
    for i = 1:n
        if color(i) == 0
            stack = i;
            while ~isempty(stack)
                v = stack(end);
                if color(v) == 0
                    color(v) = 1;
                    neighbors = find(adjacencyMatrix(v, :));
                    for w = neighbors
                        if color(w) == 0
                            parent(w) = v;
                            stack = [stack, w];
                        elseif color(w) == 1 && parent(v) ~= w
                            % A cycle is found
                            cycle = [w, v];
                            while parent(v) ~= w
                                v = parent(v);
                                cycle = [v, cycle];
                            end
                            cycles = [cycles; {cycle}];
                        end
                    end
                else
                    color(v) = 2;
                    stack(end) = [];
                end
            end
        end
    end
end
