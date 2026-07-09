function grid = vector_to_symmetric_grid(vec,num_net)
    % Initialize a matrix of zeros
    grid = zeros(num_net, num_net);

    % Index to keep track of the position in the vector
    idx = 1;

    % Fill in the upper triangle of the matrix, including the diagonal
    for i = 1:num_net
        for j = i:num_net
            grid(i, j) = vec(idx);
            idx = idx + 1;
        end
    end

    % Since the matrix is symmetric, copy the upper triangle to the lower triangle
    grid = grid + triu(grid, 1)';
end
