close all;
clear;
clc;

% Importing functions
addpath(genpath("./response_integration"));
addpath(genpath("./direct_integration"));

dt = 0.005; # Time step | Seconds
t = 0:dt:2; # Time vector (t(i), t(i + 1), t(i + 2), ...) | Seconds
F = get_force(t); # Force vector (F(t(i), F(t(i + 1)), ...) | Newtons
x0 = 0; % Initial position | Meters
xdot0 = 0; % Initial velocity | Meters per Second
M = 1; % Equivalent Mass | Kilograms
K = (600 + 400 + 200); % Spring Constant | Newtons per Meter

% Calculating Damping Constant C
A0_over_An = 1/0.25; # Ratio between amplitudes
wn = sqrt(K / M); # Natural frequency
n = 5; # Number of oscilations
sigma = log(A0_over_An) / n;
zeta_ = sigma / sqrt((2 * pi)^2 + (sigma)^2);
C = 2 * M * wn * zeta_; % Damping Constant | Kilograms per Second

% Integration Methods

[x, xdot, x2dot] = central_difference_integration(t, F, x0, xdot0, M, K, C);
figure;
plot(t, x);
hold on;

teta = 1.0;
[x2, xdot, x2dot] = wilson_integration(t, F, x0, xdot0, M, K, C, teta);
plot(t, x2);
hold on;

alfa = 1/6; beta_ = 0.5;
[x3, xdot, x2dot] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_);
plot(t, x3);
hold on;

% Response Integration
[x4, xdot, x2dot] = constant_approximation_int(t, F, x0, xdot0, M, K, C);
plot(t, x4);
hold on;

[x5, xdot, x2dot] = linear_approximation_int(t, F, x0, xdot0, M, K, C);
plot(t, x5);
hold off;
