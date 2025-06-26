function no_sign_changes = find_non_zeros(f, X, N)
    % Function to find intervals containing zeros of f
    % f: function handle
    % X: interval of type 'infsup' from INTLAB
    % N: number of subintervals

    % Subdivide the interval X into N subintervals
    x_points = linspace(inf(X), sup(X), N+1);

    % Evaluate function values at subinterval endpoints
    f_values = f(I_intval(x_points));

    % Identify sign changes between adjacent points
    no_sign_changes = (f_values(1:end-1) .* f_values(2:end) > 0);

    % Construct intervals where sign changes occur
    non_zero_intervals = arrayfun(@(i) I_infsup(x_points(i), x_points(i+1)), find(no_sign_changes), 'UniformOutput', false);

    % Convert cell array to matrix if needed
    non_zero_intervals = vertcat(non_zero_intervals{:});
end


% 


