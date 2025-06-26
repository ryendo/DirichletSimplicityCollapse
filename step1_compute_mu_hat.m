% This script computes a rigorous upper bound for the second Dirichlet eigenvalue (mu_hat_2)
% on collapsing triangles, as described in Section 4.1 of the paper.
% It utilizes the Rayleigh-Ritz method with INTLAB's interval arithmetic for rigorous enclosures.

tic
format compact short infsup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Symbolic Definitions and Matrix Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% Define the number of basis functions for the Rayleigh-Ritz method.
N = 17;  % default: N=17
% Define the number of subintervals for the 's' parameter in [0, 1].
N_s = 100; % default: N_s=100
% Set the height parameter 't' for the collapsing triangle.
t_val = tan(intval('pi')/60);

% Generate CSV file name for storing computed upper bounds.
csvFileName = sprintf('results/step1_1_results_N%d_Ns%d.csv', N, N_s);

% Initialize or load existing results from the CSV file.
if exist(csvFileName, 'file')
    fprintf('Loading existing results from %s\n', csvFileName);
    existingData = readtable(csvFileName);
    completedIndices = existingData.s_ind;
else
    fprintf('Creating new results file: %s\n', csvFileName);
    fid = fopen(csvFileName, 'w');
    fprintf(fid, 's_ind,s_inf,s_sup,upper_bounds\n');
    fclose(fid);
    completedIndices = [];
end

% Declare symbolic variables for matrix indices, parameters t, s, and pi.
syms idx1 idx2 iInd jInd t s a pi real
assume(idx1, 'integer');
assume(idx2, 'integer');
assume(iInd <= sym(N));
assume(jInd <= sym(N));
assume(t > 0);
assume(abs(s) < 1);

% Define paths for saving and loading symbolic matrices.
symFilenameBase = sprintf('upper_bound_matrix/K_M_sym_with_N_%d', N);
symMatFile = [symFilenameBase, '.mat'];
symMFile = [symFilenameBase, '.m'];

% Load or construct symbolic matrices (K_mat_sym, M_mat_sym) representing the integrals.
if exist(symMatFile, 'file')
    fprintf('Loading symbolic integrals from file: %s\n', symMatFile);
    load(symMatFile);
else
    fprintf('No existing symbolic integrals file. Building them now...\n');

    % Define the interval boundary points a_left and a_right for the 1D domain I^t.
    a_left_sym  = -t^(-sym(2)/sym(3))*(sym(1) + s);
    a_right_sym =  t^(-sym(2)/sym(3))*(sym(1) - s);
    x_split_sym = sym(0);

    syms x real

    % Define the coordinate transformation as per the problem formulation.    
    Y_sym = t^(sym(2)/sym(3))*x + s;  

    % Define the symbolic height function h(Y) and its derivative h'(Y) for the triangle.
    % h_p_sym for Y > s, h_m_sym for Y < s.
    h_p_sym       = @(Y) -t/(sym(1) - s)*(Y - sym(1));
    h_m_sym       = @(Y)  t/(sym(1) + s)*(Y + sym(1));
    h_prime_p_sym = @(Y) -t/(sym(1) - s);
    h_prime_m_sym = @(Y)  t/(sym(1) + s);

    % Define the symbolic potential V(t,x) for the Schrödinger operator.
    % V_plus_sym for x > 0, V_minus_sym for x < 0. (after coordinate transform)
    V_plus_sym = t^(sym(4)/sym(3))*( (pi^2)/(h_p_sym(Y_sym)^2) ...
             + ((sym(3) + sym(4)*pi^2)*(h_prime_p_sym(Y_sym)^2))/(sym(12)*h_p_sym(Y_sym)^2) ...
             - (pi^2/t^2) );

    V_minus_sym = t^(sym(4)/sym(3))*( (pi^2)/(h_m_sym(Y_sym)^2) ...
              + ((sym(3) + sym(4)*pi^2)*(h_prime_m_sym(Y_sym)^2))/(sym(12)*h_m_sym(Y_sym)^2) ...
              - (pi^2/t^2) );
    
    % Define the symbolic centers 'c_sym' for the basis functions 'phi_sym'.
    c_sym = symfun(...
    (a_left_sym + (a_right_sym - a_left_sym)/sym(2)) + ...
    (idx1 - sym(1)) * ((a_right_sym - a_left_sym) - (a_right_sym - a_left_sym)/sym(2)) / (sym(N) - sym(1)), ...
    idx1);
    
    % Define the symbolic Rayleigh-Ritz basis function phi_i(x).
    phi_sym = symfun((x - a_left_sym) * (a_right_sym - x) * exp(-sym(1)*(x - c_sym(idx1))^2), [x, idx1]);
    dphi_sym = symfun(diff(phi_sym(x,idx1), x), [x, idx1]);

    % Initialize symbolic stiffness (K) and mass (M) matrices.
    K_mat_sym = sym(zeros(N,N));
    M_mat_sym = sym(zeros(N,N));

    % Prepare symbolic expressions for integrals to build K and M matrices.
    phi_i_expr  = simplify(phi_sym(x, iInd));
    phi_j_expr  = simplify(phi_sym(x, jInd));
    dphi_i_expr = simplify(dphi_sym(x, iInd));
    dphi_j_expr = simplify(dphi_sym(x, jInd));

    % Define integrands for the stiffness matrix (K) for both positive and negative x.
    integrand_K_left  = simplify(dphi_i_expr .* dphi_j_expr + V_minus_sym.*phi_i_expr.*phi_j_expr);
    integrand_K_right = simplify(dphi_i_expr .* dphi_j_expr + V_plus_sym .*phi_i_expr.*phi_j_expr);
    % Define integrands for the mass matrix (M).
    integrand_M_left  = simplify(phi_i_expr.*phi_j_expr);
    integrand_M_right = simplify(phi_i_expr.*phi_j_expr);
    
    % Perform symbolic integration for K and M matrix elements.
    K_left  = int(integrand_K_left,  x, a_left_sym,  x_split_sym);
    K_right = int(integrand_K_right, x, x_split_sym, a_right_sym);
    M_left  = int(integrand_M_left,  x, a_left_sym, x_split_sym);
    M_right = int(integrand_M_right, x, x_split_sym, a_right_sym);

    % Sum the integrals over the split domain.
    valK = simplify(K_left + K_right);
    valM = simplify(M_left + M_right);

    % Populate the symmetric symbolic K and M matrices by substituting specific indices.
    for iInd_val = 1:N
        for jInd_val = iInd_val:N
            K_mat_sym(iInd_val,jInd_val) = subs(valK, [iInd, jInd], [iInd_val, jInd_val]);
            K_mat_sym(jInd_val,iInd_val) = K_mat_sym(iInd_val,jInd_val);
            M_mat_sym(iInd_val,jInd_val) = subs(valM, [iInd, jInd], [iInd_val, jInd_val]);
            M_mat_sym(jInd_val,iInd_val) = M_mat_sym(iInd_val,jInd_val);
        end
    end

    % Compute Taylor expansion of K and M matrices with respect to 's'.
    % This handles the dependency on 's' for later interval evaluation.
    TaylorK = taylor(K_mat_sym, s, 'ExpansionPoint', a, 'Order', 1);
    TaylorM = taylor(M_mat_sym, s, 'ExpansionPoint', a, 'Order', 1);

    fprintf('Saving symbolic integrals to file: %s\n', symMatFile);
    save(symMatFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.5) Generating .m Function for Interval Matrix Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, funcName] = fileparts(symMFile);
if ~exist(symMFile, 'file')
    fprintf('Creating .m file: %s\n', symMFile);
    fid = fopen(symMFile, 'w');
    fprintf(fid, 'function [K, M] = %s(a_val, s_val, t_val)\n', funcName);

    fprintf(fid, 'K = [\n');
    for ii = 1:N
        rowStr = '';
        for jj = 1:N
            exprStr = char(TaylorK(ii,jj));
            exprStr = regexprep(exprStr, '(\d+(\.\d*)?|\d*\.\d+)', 'I_intval(''$1'')');
            exprStr = strrep(exprStr, 'pi', 'I_pi');
            exprStr = regexprep(exprStr, '\<s\>', 's_val');
            exprStr = regexprep(exprStr, '\<t\>', 't_val');
            exprStr = regexprep(exprStr, '\<a\>', 'a_val');
            rowStr = [rowStr, ' ', exprStr];
            if jj < N
                rowStr = [rowStr, ','];
            end
        end
        fprintf(fid, '%s;\n', rowStr);
    end
    fprintf(fid, '];\n\n');

    fprintf(fid, 'M = [\n');
    for ii = 1:N
        rowStr = '';
        for jj = 1:N
            exprStr = char(TaylorM(ii,jj));
            exprStr = regexprep(exprStr, '(\d+(\.\d*)?|\d*\.\d+)', 'I_intval(''$1'')');
            exprStr = strrep(exprStr, 'pi', 'I_pi');
            exprStr = regexprep(exprStr, '\<s\>', 's_val');
            exprStr = regexprep(exprStr, '\<t\>', 't_val');
            exprStr = regexprep(exprStr, '\<a\>', 'a_val');
            rowStr = [rowStr, ' ', exprStr];
            if jj < N
                rowStr = [rowStr, ','];
            end
        end
        fprintf(fid, '%s;\n', rowStr);
    end
    fprintf(fid, '];\n\n');
    fprintf(fid, 'end\n');
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Interval Matrix Construction and Eigenvalue Problem Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting eigenvalue computation for N=%d basis functions over %d s-intervals.\n', N, N_s);

for i = 1:N_s
    % Skip computation for s-intervals that have already been completed.
    if ismember(i, completedIndices)
        fprintf('Skipping s-interval %d/%d (already computed).\n', i, N_s);
        continue;
    end

    % Define the current s-interval for computation.
    s_inf = inf((i-1)/intval(N_s));
    s_sup = sup(i/intval(N_s));
    a_val = intval(s_inf);
    s_val = infsup(s_inf, s_sup);

    fprintf('\nProcessing s-interval %d/%d: s in [%.4f, %.4f]\n', i, N_s, inf(s_val), sup(s_val));

    % Construct interval matrices K and M using the auto-generated function.
    tic
    % K_inf, M_inf are evaluated at a single point (s_inf) for error estimation base.
    [K_inf, M_inf] = feval(funcName, a_val, s_inf, t_val);
    % K_I, M_I are evaluated over the interval (s_val) to capture parameter uncertainty.    
    [K_I, M_I] = feval(funcName, a_val, s_val, t_val);
    fprintf('Matrix construction time for interval %d: %.2f seconds.\n', i, toc);    
    toc
    tic

    % Calculate eigenvalue bounds using interval arithmetic.
    % This part computes rigorous upper bounds for the eigenvalues.
    
    % Compute the inverse of the smallest eigenvalue of M_I (related to norm of M).
    % This provides a scaling factor for error propagation.    
    mu = 1/veigs(M_I, I_eye(N, N), 1); 
    mu = sqrt(min(mu));

    % Diagonalize K_inf and M_inf using eigenvectors of mid(M_inf) for stability.    
    [V,D]=eig(mid(M_inf));
    K_inf_ = V'*K_inf*V; K_inf_=hull(K_inf_,K_inf_');
    M_inf_ = V'*M_inf*V; M_inf_=hull(M_inf_,M_inf_');

    % Compute the initial rigorous upper bounds for the smallest 3 eigenvalues ('sm', 3).    
    upper_bounds_tmp = veigs(K_inf_, M_inf_, 'sm', 3);

    % Calculate the error introduced by the Taylor expansion of K and M over the s-interval.    
    Kdiff = norm(K_I - K_inf, 2);
    Mdiff = norm(M_I - M_inf, 2);
    error_s_upper = mu * Kdiff + mu^2 * Mdiff * Kdiff;

    % Apply the error term to the computed eigenvalue bounds to obtain the final rigorous upper bounds.
    upper_bounds = I_sup(upper_bounds_tmp + I_hull(-error_s_upper, error_s_upper) * ones(size(upper_bounds_tmp)));
    toc

    % Append the computed results (s-interval and upper bounds) to the CSV file.
    fid = fopen(csvFileName, 'a');
    upper_bounds_str = sprintf('%.15f,', upper_bounds);
    upper_bounds_str = upper_bounds_str(1:end-1);
    fprintf(fid, '%d,%.15f,%.15f,%s\n', i, s_inf, s_sup, upper_bounds_str);
    fclose(fid);

end

toc
fprintf('\nAll computations for step 1.1 are completed. Results are saved in %s\n', csvFileName);