function [x, xdot, x2dot] = linear_approximation_int(t, F, x0, xdot0, m, k, c)
    % Linear Approximation Integration Method
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
        first_multiplier = ((F(i) - F(i - 1)) / (k * dti));
        f1 = dti; f2 = (2 * zeta_ / wn);
        f3 = exp(-zeta_ * wn * dti);
        f31 = ((2 * zeta_) / wn) * cos(wd * dti);
        f32 = ((wd^2 - zeta_^2 * wn^2) / (wn^2 * wd)) * sin(wd * dti);
        first_factor = first_multiplier * (f1 - f2 + f3 * (f31 - f32));

        second_multiplier = (F(i - 1) / k);
        s1 = 1; s2 = exp(-zeta_ * wn * dti);
        s21 = cos(wd * dti);
        s22 = (zeta_ * wn / wd) * sin(wd * dti);
        second_factor = second_multiplier * (s1 - s2 * (s21 + s22));

        third_multiplier = exp(-zeta_ * wn * dti);
        t1 = x(i - 1) * cos(wd * dti);
        t2 = ((xdot(i - 1) + zeta_ * wn * x(i - 1)) / wd) * sin(wd * dti);
        third_factor = third_multiplier * (t1 + t2);

        x(i) = first_factor + second_factor + third_factor;

        % Calculating xdot(i)
        first_multiplier = ((F(i) - F(i - 1)) / (k * dti));
        f1 = 1; f2 = exp(-zeta_ * wn * dti);
        f21 = cos(wd * dti);
        f22 = (zeta_ * wn / wd) * sin(wd * dti);
        first_factor = first_multiplier * (f1 - f2 * (f21 + f22));

        s1 = (F(i - 1) / k) * exp(-zeta_ * wn * dti) * ((wn^2) / wd) * sin(wd * dti);
        s2 = exp(-zeta_ * wn * dti);
        s21 = xdot(i - 1) * cos(wd * dti);
        s22 = (zeta_ * wn / wd) * (xdot(i - 1) + wn * x(i - 1) / zeta_) * sin(wd * dti);
        second_factor = s1 + s2 * (s21 - s22);

        xdot(i) = first_factor + second_factor;

    end

end
