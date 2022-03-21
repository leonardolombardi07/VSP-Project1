function [x, xdot, x2dot] = central_difference_integration(t, F, x0, xdot0, m, k, c)
    % Central Difference Integration Method
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
    %
    % Output
    % ----------
    %       [x]:        Displacement Response   [n,1]
    %       [xdot]:     Velocity                [n,1]
    %       [x2dot]:    Acceleration            [n,1]

    % Initializing Variables
    dt = t(2) - t(1);
    x = zeros(length(t) + 1, 1);
    xdot = zeros(length(t), 1); x2dot = xdot;

    % Initial Conditions
    x2dot0 = (F(1) - c * xdot0 - k * x0) / m;
    x(1) = x0; xdot(1) = xdot0; x2dot(1) = x2dot0;

    % Integration Constants
    a0 = ((2 * m) / (dt^2)) - k;
    a1 = (c / (2 * dt)) - (m / (dt^2));
    K_circunflex = (m / (dt^2)) + (c / (2 * dt));

    for i = 1:(length(t))

        % We define the variable x_iminus1 because, in the first iteration
        % (when i=1), x(i-1) is x(0), which does not exist
        if i == 1
            x_iminus1 = x0 - dt * xdot0 + (((dt^2) * x2dot0) / 2);
        else
            x_iminus1 = x(i - 1);
        end

        F_circunflex = a0 * x(i) + a1 * x_iminus1 + F(i);

        % ALERT: DO NOT CHANGE THE ORDER OF THE VARIABLES BELOW
        % xdot(i + 1), for example, accesses the variable x(i + 1),
        % which is defined in the previous
        x(i + 1) = F_circunflex / K_circunflex;
        xdot(i) = (x(i + 1) - x_iminus1) / (2 * dt);
        x2dot(i) = (x(i + 1) - 2 * x(i) + x_iminus1) / (dt^2);
    end

    x = x(1:end - 1);
end
