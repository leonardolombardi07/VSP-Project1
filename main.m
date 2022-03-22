%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleaning Workspace & Importing Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; % Closes all figures
clear; # Clears variables from workspace
clc; % Clears all text from command window

% Importing functions from folders with given relative path
addpath(genpath("./formulas"));
addpath(genpath("./plotting"));
addpath(genpath("./response_integration"));
addpath(genpath("./direct_integration"));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initializing Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.005; % Time step | Seconds
t = 0:dt:2; % Time vector (t(i), t(i + 1), t(i + 2), ...) | Seconds
F = get_force(t); % Force vector (F(t(i), F(t(i + 1)), ...) | Newtons
x0 = 0; % Initial position | Meters
xdot0 = 0; % Initial velocity | Meters per Second
M = 1; % Equivalent Mass | Kilograms
K = (600 + 400 + 100); % Spring Constant | Newtons per Meter

% Calculating Damping Constant C
N = 5; % Number of oscilations
A_ratio = 1/0.25; % Ratio between initial amplitude (A0) and amplitude in time NT (A(NT))
C = get_damping_coefficient(A_ratio, M, K, N); % Damping Constant | Kilograms per Second

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Integration Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Direct Integration

% Wilson Method
teta = 1.4;
[x_wilson, xdot_wilson, x2dot_wilson] = wilson_integration(t, F, x0, xdot0, M, K, C, teta);
plot_method(t, x_wilson, xdot_wilson, x2dot_wilson, "Integration with Method of Wilson");

% Newmark Method
alfa = 1/6; beta_ = 1/1.5;
[x3_newmark, xdot_newmark, x2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_);
plot_method(t, x3_newmark, xdot_newmark, x2dot_newmark, "Integration with Method of Newmark 1");

alfa = 1/6; beta_ = 1/2;
[x3_newmark, xdot_newmark, x2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_);
plot_method(t, x3_newmark, xdot_newmark, x2dot_newmark, "Integration with Method of Newmark 2");

alfa = 1/6; beta_ = 1/3;
[x3_newmark, xdot_newmark, x2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_);
plot_method(t, x3_newmark, xdot_newmark, x2dot_newmark, "Integration with Method of Newmark 3");

% Central Difference Method
[x_central, xdot_central, x2dot_central] = central_difference_integration(t, F, x0, xdot0, M, K, C);
plot_method(t, x_central, xdot_central, x2dot_central, "Integration with Method of Central Difference");

% Response Integration

% Constant Approximation Method
[x_const, xdot_const, x2dot_const] = constant_approximation_int(t, F, x0, xdot0, M, K, C);
plot_method(t, x_const, xdot_const, x2dot_const, "Response Integration by Constant Approximation");

% Linear Approximation Method
[x_linear, xdot_linear, x2dot_linear] = linear_approximation_int(t, F, x0, xdot0, M, K, C);
plot_method(t, x_linear, xdot_linear, x2dot_linear, "Response Integration by Linear Approximation");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical Solution using Laplace Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_analytical = get_inverse_laplace_transform(t);
xdot_analytical = diff(x_analytical) / dt;
x2dot_analytical = diff(xdot_analytical) / dt;
plot_method(t, x_analytical, xdot_analytical, x2dot_analytical, "Analytical solution (using Laplace transform and Heaviside step function");
