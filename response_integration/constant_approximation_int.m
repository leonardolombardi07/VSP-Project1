function [x, xdot, x2dot] = constant_approximation_int(t, F, x0, xdot0, m, k, c)
    % Constant Approximation Integration Method
    %--------------------------------------------------------------------------
    % Integrates a 1-DOF system with an equivalent mass "m", spring stiffness "k" and damping
    % coeffiecient "c" subjected to a non-periodic external force F(t).
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
    x = zeros(length(t), 1); xdot = x; x2dot = x;

    % Initial Conditions
    x(1) = x0; xdot(1) = xdot0;
    x2dot(1) = (F(1) - k * x0 - c * xdot0) / m;

    % Integration Constants
    wn = sqrt(k / m);
    zeta_ = c / (2 * m * wn);
    wd = wn * sqrt(1 - zeta_^2);

    for i = 2:(length(t))
        dti = t(i) - t(i - 1);

        % Calculating x(i)
        first_factor = (F(i) / k) * (1 - ...
            exp(-zeta_ * wn * dti) * ...
            (cos(wd * dti) + (zeta_ * wn / wd) * sin(wd * dti)));

        second_factor = (exp(-zeta_ * wn * dti)) * ...
            (x(i - 1) * cos(wd * dti) + ...
            ((xdot(i - 1) + zeta_ * wn * x(i - 1)) / wd) * sin(wd * dti));

        x(i) = first_factor + second_factor;

        % Calculating xdot(i)
        first_factor = (F(i) * wd / k) * exp(-zeta_ * wn * dti)* ...
            (1 + (zeta_ * wn / wd)^2) * sin(wd * dti);

        second_factor = wd * exp(-zeta_ * wn * dti);

        third_factor = (-x(i - 1) * sin(wd * dti));

        fourth_factor = ((xdot(i - 1) + zeta_ * wn * x(i - 1)) / wd) * ...
            cos(wd * dti);

        fifth_factor = (-zeta_ * wn / wd) * ...
            (x(i - 1) * cos(wd * dti) + ...
            ((xdot(i - 1) + zeta_ * wn * x(i - 1)) / wd) * sin(wd * dti));

        xdot(i) = first_factor + second_factor * ...
            (third_factor + fourth_factor + fifth_factor);

        % Calculating x2dot(i)
        x2dot(i) = (F(i) - c * xdot(i) - k * x(i)) / m;
    end

end
