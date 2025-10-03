function [sz, rows, cols] = checkVectorSize(vec)
    % Return size, number of rows and columns of a vector
    sz = size(vec);
    rows = size(vec, 1);
    cols = size(vec, 2);
end
