function [i, j] = get_neighborhood4(n, m)
    % Return edges in 4-connected grid to feed sparse matrix:
    % i(1) <-> j(1) is the first horisontal edge 
    % ...

    % Horisontal edges
    hor_i = [1:(n * (m-1))]';
    hor_j = [(n+1):(n*m)]';

    % Vertical edges
    ver_i = [];
    ver_j = [];
    for column_idx = 1:m
        column = get_column(n, m, column_idx);
        ver_i = [ver_i; column(1:end-1)];
        ver_j = [ver_j; column(2:end)];
    end
    i = [hor_i; ver_i];
    j = [hor_j; ver_j];
end

function column = get_column(n, m, idx)
    start_idx =  (idx - 1) * n + 1;
    column = [start_idx:(start_idx + n - 1)]';
end
