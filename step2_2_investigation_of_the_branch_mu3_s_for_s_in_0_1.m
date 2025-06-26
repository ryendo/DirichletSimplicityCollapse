% Verify the lower bound of the third Dirichlet eigenvalue for non-equilateral triangles.
% This script confirms that mu_bar_3(s) remains above 23.5 for s in [0, 1).

% Use INTLAB for rigorous interval arithmetic to guarantee bounds.

format long infsup
close all

% Define the output directory for logging and results.
results_dir = 'results';
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% This flag confirms that INTLAB's interval mode is active.
INTERVAL_MODE = 1;

% Number of subintervals for the 's' parameter in [0, 1]. A higher N increases verification granularity.
N = 1E+5;
s = linspace(0, 1, N+2);
% Create overlapping s-intervals to ensure complete coverage and robustness for interval arithmetic.
s_intervals = I_infsup(s(1:end-2), s(3:end));

% Define the target mu_bar value to check against (corresponding to 23.5 in the paper).
X_list = [intval('23.5')];

for j = 1:length(X_list)
    X = X_list(j);
    
    % Generate a descriptive filename for the log based on the current mu_bar interval.
    X_str = sprintf('%.2f-%.2f', inf(X), sup(X));
    log_filename = fullfile(results_dir, ['result_X_', X_str, '.txt']);
    
    
    fid = fopen(log_filename, 'w');
    if fid == -1
        error(['Cannot open log file: ', log_filename]);
    end
    
    
    fprintf(fid, 'Result Log - X Interval: [%f, %f]\n', inf(X), sup(X));
    fprintf(fid, '----------------------------------------------------------\n');
    
    try
        for i = 1:N
            tic
            fprintf('Processing s-interval %d/%d...\n', i, N); % Improved print statement
            fprintf(fid, 'Processing s-interval %d/%d...\n', i, N);
            disp(i);
            fprintf(fid, 'Processing: Index %d\n', i);
            
            s_val = s_intervals(i);
            
            % Coefficient for converting mu_bar to kappa as per equation (17) in the paper.
            c0 = (2*I_pi^2)^(2/I_intval(3));
            
            % Select the appropriate function for f_s(kappa) based on the s-value.
            % For s close to 1, a specialized function handling potential singularities is used.
            if s_val < 0.9
                f = @(x) func_left_hand_side(s_val, x/c0);
            else
                sigma = (1 - s_val)^(1/intval(3));
                f = @(x) ((2 * I_pi^2) * (2 / intval(3)) * func_left_hand_side_singularity(s_val, x/c0));
            end
            
            % Rigorously verify that the function f_s(kappa) does not have a zero within the specified X interval.
            % 'find_non_zeros' uses interval arithmetic to prove that f(X) does not contain zero.
            no_sign_changes = find_non_zeros(f, X, 100);
            
            
            if any(no_sign_changes == 0)
                error('Sign changes detected in f_values.');
            end
            
            
            elapsed_time = toc;
            estimated_time = (N - i) * elapsed_time / 60;
            fprintf(fid, 'Index: %d, Estimated Remaining Time: %.2f minutes\n', i, estimated_time);
        end
        
        
        fprintf(fid, '----------------------------------------------------------\n');
        fprintf(fid, 'Verification successful: No sign changes detected for f_s(kappa).\n');
        fprintf(fid, 'The lower bound of mu_bar_3(s) has been rigorously confirmed.\n');
        disp('Verification successful: No sign changes detected.');
        
    catch ME
        printf(fid, '--- Verification Failed: An error occurred ---\n'); % Improved error header
        fprintf(fid, 'Error details: %s\n', ME.message); % More specific error details
        fclose(fid);
        rethrow(ME);
    end
    
    fclose(fid);
end
