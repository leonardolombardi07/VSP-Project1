function [U, Udot, U2dot, t] = wilson_integration(t, F, U0, Udot0, M, K, C, teta)
    % Wilson-teta Integration Method
    %--------------------------------------------------------------------------
    % Integrates a N-DOF system wih a mass matrix "M", stiffness matrix "K" and
    % damping matrix "C" subjected to an external force F(t).
    % Returns the displacement, velocity and acceleration of the system with
    % respect to an inertial frame of reference.
    %
    % Input
    % ----------
    %       [t] :       Time Vector               [n, 1]
    %       [F] :       External Force Matrix     [DOF, n]
    %       [U0]:       Initial Position Vector   [DOF, 1]
    %       [Udot0]:    Initial Velocity Vector   [DOF, 1]
    %       [M]:        Mass Matrix               [DOF, DOF]
    %       [K]:        Stiffness Matrix          [DOF, DOF]
    %       [C]:        Damping Matrix            [DOF, DOF]
    %       [teta]:     Teta factor               scalar
    %
    % Output
    % ----------
    %       [U]:        Displacement Vector        [DOF, n]
    %       [Udot]:     Velocity Vector            [DOF, n]
    %       [U2dot]:    Acceleration Vector        [DOF, n]

    dt = t(2) - t(1); DOF = length(K); n = length(t);
    U = zeros(DOF, n); Udot = U; U2dot = U; Feff = U;

    % Initial Conditions
    U(:, 1) = U0;
    Udot(:, 1) = Udot0;
    U2dot(:, 1) = M \ (F(:, 1) - K * U(:, 1) - C * Udot(:, 1));

    % Integration Constants
    a0 = 6 / (teta * dt)^2; a1 = 3 / (teta * dt); a2 = 2 * a1;
    a3 = teta * dt / 2; a4 = a0 / teta; a5 = -a2 / teta;
    a6 = 1 - 3 / teta; a7 = dt / 2; a8 = dt^2/6;

    % Form Effective Stiffness Matrix
    Keff = K + a0 * M + a1 * C;

    for i = 1:(n - 1)
        % Calculating Effective Force
        Feff(:, i) = F(:, i) + teta * (F(:, i + 1) - F(:, i)) + M * (a0 * U(:, i) + a2 * Udot(:, i) + 2 * U2dot(:, i)) + C * (a1 * U(:, i) + 2 * Udot(:, i) + a3 * U2dot(:, i));
        % Solving for displacements at time (t+dt)
        U(:, i + 1) = Keff \ Feff(:, i);
        % Calculating displacements, velocities and accelerations at time t+dt
        U2dot(:, i + 1) = a4 * (U(:, i + 1) - U(:, i)) + a5 * (Udot(:, i)) + a6 * U2dot(:, i);
        Udot(:, i + 1) = Udot(:, i) + a7 * (U2dot(:, i + 1) + U2dot(:, i));
        U(:, i + 1) = U(:, i) + dt * Udot(:, i) + a8 * (U2dot(:, i + 1) + 2 * U2dot(:, i));
    end

end
