% t = 0:0.28:3.36; # Time Vector
% M = [2 0; 0 1]; % Mass Matrix
% K = [6 -2; -2 4]; % Stiffness Matrix
% DOF = length(K); # Number of Degrees of Freedom
% F = zeros(DOF, length(t)); # External Force
% U0 = zeros(DOF, 1); # Initial Positions
% Udot0 = zeros(DOF, 1); # Initial Velocities
% C = zeros(size(K)); % Damping Matrix

dt = 0.02;
t = 0:dt:0.2;
% F = get_force(t);
F = zeros(length(t), 1);
x0 = 1; % m
xdot0 = 0; % m / s
% M = 1; % Kg
M = 4; % Kg
% K = (600 + 400 + 200); % N/M
K = 2000;
% C = 3.0542; % Kg / s
C = 0; % Kg / s

[u1, udot, u2dot] = central_difference_integration(t, F, x0, xdot0, M, K, C);
disp(u1);
disp("\n");

teta = 1.0;
[u2, udot, u2dot] = wilson_integration(t, F, x0, xdot0, M, K, C, teta);
disp(u2);
disp("\n")

alfa = 1/6; beta_ = 0.5;
[u3, udot, u2dot] = newmark_integration(t, F, x0, xdot0, M, K, C, alfa, beta_);
% [u4, udot4, u2dot4] = newmark_int(t, F, x0, xdot0, M, K, C, alfa, beta_);

disp(u3);
