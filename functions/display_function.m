INTERVAL_MODE = 1;
i = 100; N = 100;
kappa = 4;
% s = I_infsup(x,x+ep);
f = @(x)(func_left_hand_side(0,x));
% f = @(x)(partial_kappa_s(x, kappa));
% f = @(x)(test_func(x));
% f = @(x)(func_left_hand_side_singularity(0.01,x));
% f = @(x)(I_minus_Ai(x));
% f = @(x)(I_d_minus_Ai(x));
% f = @(x)(I_besselj(I_intval(3)/2,x));
% f = @(x)(x.^2-2);
x = linspace(0.1,5);
% f = @(x)(test_func(x));
plot(x,I_mid(f(x)),x,0*x)
x = linspace(0.001,5,500);
y = I_mid(f(x));  % Calculated y values
save('data.mat', 'x', 'y');
% % Initial setup
% clear; close all; clc;
% 
% INTERVAL_MODE = 1;
% ep = 1E-4;
% 
% % Slider initial values
% x_start = 0;
% x_end = 5;
% 
% % x values for plotting
% x = linspace(x_start, x_end, 500);
% 
% % Function definition
% func_left_hand_side = @(s, x) (func_left_hand_side(s,x)); % Replace with appropriate function
% I_mid = @(x) mid(x); % Example function to extract midpoint for simplicity
% 
% % Prepare the figure
% fig = figure('Name', 'Slider-controlled Plot', 'NumberTitle', 'off', 'Position', [100, 100, 600, 400]);
% ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.65]);
% 
% % Slider settings
% slider = uicontrol('Style', 'slider', ...
%                    'Min', 0, 'Max', 1, ...
%                    'Value', x_start, ...
%                    'Units', 'normalized', ...
%                    'Position', [0.1, 0.1, 0.8, 0.05], ...
%                    'Callback', @(src, event) updatePlot(src, ax, x, func_left_hand_side, I_mid, ep));
% 
% % Initial plot
% updatePlot(slider, ax, x, func_left_hand_side, I_mid, ep);
% 
% % Function to update the plot
% function updatePlot(slider, ax, x, func_left_hand_side, I_mid, ep)
%     s_val = get(slider, 'Value');
%     s = I_infsup(s_val, s_val + ep); % Consider INTERVAL_MODE
%     f = @(x)(func_left_hand_side(s, x));
%     y = I_mid(f(x));
% 
%     % Update plot
%     plot(ax, x, y, x, zeros(size(x)), 'r--'); % Add y=0 line
%     xlabel(ax, 's');
%     ylabel(ax, 'f(s)');
%     grid(ax, 'on');
% end