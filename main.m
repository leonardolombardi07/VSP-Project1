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
addpath(genpath("./analytical"));

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

teta = 1.4; alfa = 1/6; beta_a = 1/1.5; beta_b = 1/2; beta_c = 1/3;
[x_wilson, xdot_wilson, x2dot_wilson] = wilson_integration(t, F, x0, xdot0, M, K, C, teta);
[xa_newmark, xadot_newmark, xa2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_a);
[xb_newmark, xbdot_newmark, xb2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_b);
[xc_newmark, xcdot_newmark, xc2dot_newmark] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_c);
[x_central, xdot_central, x2dot_central] = central_difference_integration(t, F, x0, xdot0, M, K, C);
[x_const, xdot_const, x2dot_const] = constant_approximation_int(t, F, x0, xdot0, M, K, C);
[x_linear, xdot_linear, x2dot_linear] = linear_approximation_int(t, F, x0, xdot0, M, K, C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical Solution using Laplace Transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_analytical = solve_hardcoded(t);
xdot_analytical = diff(x_analytical) / dt;
x2dot_analytical = diff(xdot_analytical) / dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots and Comparisons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot all methods separately
plot_method(t, x_wilson, xdot_wilson, x2dot_wilson, "Integration with Method of Wilson");
plot_method(t, xa_newmark, xadot_newmark, xa2dot_newmark, "Integration with Method of Newmark 1");
plot_method(t, xb_newmark, xbdot_newmark, xb2dot_newmark, "Integration with Method of Newmark 2");
plot_method(t, xc_newmark, xcdot_newmark, xc2dot_newmark, "Integration with Method of Newmark 3");
plot_method(t, x_central, xdot_central, x2dot_central, "Integration with Method of Central Difference");
plot_method(t, x_const, xdot_const, x2dot_const, "Response Integration by Constant Approximation");
plot_method(t, x_linear, xdot_linear, x2dot_linear, "Response Integration by Linear Approximation");
plot_method(t, x_analytical, xdot_analytical, x2dot_analytical, "Analytical solution (using Laplace transform and Heaviside step function");

% Plot all methods together
plot_numerical_vs_analytical(t,
x_wilson, xdot_wilson, x2dot_wilson,
xa_newmark, xadot_newmark, xa2dot_newmark,
x_central, xdot_central, x2dot_central,
x_const, xdot_const, x2dot_const,
x_linear, xdot_linear, x2dot_linear,
x_analytical, xdot_analytical, x2dot_analytical);

% Plot Newmark comparison for different betas
plot_newmark_comparison(t, alfa,
xa_newmark, xadot_newmark, xa2dot_newmark, beta_a,
xb_newmark, xbdot_newmark, xb2dot_newmark, beta_b,
xc_newmark, xcdot_newmark, xc2dot_newmark, beta_c);

% Write text files containing data from displacement, velocity and acceleration
X = [t; x_wilson'; xa_newmark'; xb_newmark'; xc_newmark'; x_central'; x_const'; x_linear'; x_analytical'];
write_text_file(X, "Displacement");

XDOT = [t; xdot_wilson'; xadot_newmark'; xbdot_newmark'; xcdot_newmark'; xdot_central'; xdot_const'; xdot_linear'; [xdot_analytical; 0]'];
write_text_file(XDOT, "Velocity");

X2DOT = [t; x2dot_wilson'; xa2dot_newmark'; xb2dot_newmark'; xc2dot_newmark'; x2dot_central'; x2dot_const'; x2dot_linear'; [x2dot_analytical; 0; 0]'];
write_text_file(X2DOT, "Acceleration");
