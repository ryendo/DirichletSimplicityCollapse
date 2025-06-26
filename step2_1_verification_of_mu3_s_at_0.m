% This script rigorously verifies the third Dirichlet eigenvalue (kappa_3) at s=0.
% This is a critical step in confirming the lower bound of mu_bar_3(s) as per Section 4.2 of the paper.

% Use INTLAB for rigorous interval arithmetic to guarantee all computed enclosures.

format long infsup

% Flag indicating that INTLAB's interval mode is active.
INTERVAL_MODE=1;

% Define the search interval for kappa.
% This range is chosen to encompass the expected third eigenvalue.
x0 = infsup(0,4);

% Define the function f_s(kappa) from equation (19) in the paper.
% 'func_left_hand_side' is a helper function that implements the left-hand side of this implicit equation.
s = 0;
f = @(x)(func_left_hand_side(s,x));

% Configure options for verifynlssall, INTLAB's rigorous nonlinear system solver.
% These parameters control the precision and efficiency of the zero-finding process.
options = verifynlssallset( ...
    'Boxes', 2^10, ...         % Number of initial partition boxes for global search.
    'TolXAbs', 1e-14, ...      % Absolute tolerance for the width of the final solution boxes.
    'TolXRel', 1e-14, ...      % Relative tolerance for the width of the final solution boxes.
    'NIT', 5, ...              % Maximum number of global iterations.
    'ND', 10 ...               % Maximum number of local Newton iterations.
);

% Perform rigorous search for all zeros of f(kappa) within x0.
% X will contain validated enclosures of the zeros.
% XS will contain 'suspicious' boxes that could not be rigorously processed.
[X, XS] = verifynlssall(f, x0,options)

% Crucially, verify that no suspicious boxes remain.
% If XS_inf is not empty, it means not all zeros were rigorously enclosed, 
% indicating a potential issue or insufficient precision.
if length(XS) > 0
    fprintf('Error: Unresolved suspicious boxes remain after verification.\n');
    XS % Display the suspicious boxes for debugging.
    error('Verification failed: Not all zeros could be rigorously enclosed.');
else
    fprintf('Verification at s=0 completed successfully. No suspicious boxes were found.\n');
end

kappa_3_inf = intval(inf(X(3)));

% Define the conversion coefficient from kappa to mu_bar as per equation (17) in the paper.
% (2*pi^2)^(2/3)
conversion_coeff = (2 * I_pi^2)^(2/I_intval(3));

% Calculate the rigorous lower bound for mu_bar_3(0).
mu_bar_3_at_0_inf = kappa_3_inf * conversion_coeff;

% Display the rigorous lower bound of mu_bar_3(0).
fprintf('\nRigorous lower bound for mu_bar_3(0): %.4f\n', inf(mu_bar_3_at_0_inf));
