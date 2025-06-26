%The path of INTLAB toolbox and initialization.
addpath('Intlab_V12')

%The path of the codes for switch between verified computing and approximate computing.
addpath('mode_swith_interface')
addpath('verified_eig_estimation')
addpath('upper_bound_matrix')
addpath('functions')
addpath('results')

startintlab;

global INTERVAL_MODE;

INTERVAL_MODE=1;