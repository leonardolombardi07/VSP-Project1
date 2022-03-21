function [x, xdot, x2dot] = wilson_integration(t, F, x0, xdot0, m, k, c, teta)
    % Newmark Integration Method
    %--------------------------------------------------------------------------
    % Integrates a 1-DOF system with an equivalent mass "m", spring stiffness "k" and damping
    % coeffiecient "c" subjected to an external force F(t).
    % Returns the displacement, velocity and acceleration of the system with
    % respect to an inertial frame of reference.
    %
    % Input
    % ----------
    %       [t] :       Time Vector             [n,1]
    %       [F] :       External Force          [n,1]
    %       [x0]:       Initial Position        scalar
    %       [xdot0]:    Initial Velocity        scalar
    %       [m]:        Equivalent Mass         scalar
    %       [k]:        System Stiffness        scalar
    %       [c]:        System Damping          scalar
    %       [teta]:     Integration Parameter   scalar
    %
    % Output
    % ----------
    %       [x]:        Displacement Response   [n,1]
    %       [xdot]:     Velocity                [n,1]
    %       [x2dot]:    Acceleration            [n,1]

    % Initializing Variables
    dt = t(2) - t(1);
    x = zeros(length(t), 1); xdot = x; x2dot = x;

    % Initial Conditions
    x(1) = x0; xdot(1) = xdot0;
    x2dot(1) = (F(1) - k * x0 - c * xdot0) / m;

    % Integration Constants
    a0 = 6 / (teta * dt)^2;
    a1 = 3 / (teta * dt);
    a2 = 2 * a1;
    a3 = teta * dt / 2;
    a4 = a0 / teta;
    a5 = -a2 / teta;
    a6 = 1 - (3 / teta);
    a7 = dt / 2;
    a8 = dt^2/6;

    % Form Effective Stiffness
    Keff = k + a0 * m + a1 * c;

    for i = 1:(length(t) - 1)
        Feff = F(i) + teta * (F(i + 1) - F(i)) + ...
            m * (a0 * x(i) + a2 * xdot(i) + 2 * x2dot(i) + ...
            c * (a1 * x(i) + 2 * xdot(i) + a3 * x2dot(i)));
        xi_plus_teta = Feff / Keff;

        % ALERT: DO NOT CHANGE THE ORDER OF THE VARIABLES BELOW
        % xdot(i + 1), for example, accesses the variable x2dot(i + 1),
        % which is defined in the previous line
        x2dot(i + 1) = a4 * (xi_plus_teta - x(i)) + a5 * xdot(i) + a6 * x2dot(i);
        xdot(i + 1) = xdot(i) + a7 * (x2dot(i + 1) + x2dot(i));
        x(i + 1) = x(i) + dt * xdot(i) + a8 * (x2dot(i + 1) + 2 * x2dot(i));
    end

end
