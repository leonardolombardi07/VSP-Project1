function C = get_damping_coefficient(A_ratio, m, k, n)
    % Damping Coefficient Calculation
    %--------------------------------------------------------------------------
    % Calculates the damping coefficient from a system with the ratio between the
    % initial amplitude (A0) and a amplitude A(nT) at time n*T (where n is the number of
    % cycles and T is the period), equivalent mass "m", and spring stiffness "k"
    % Returns the damping coefficient
    %
    % Input
    % ----------
    %       [A_ratio]:  Ratio between amplitudes (A(t=0)/A(t=nT))     scalar
    %       [m]:        Equivalent Mass                               scalar
    %       [k]:        System Stiffness                              scalar
    %       [c]:        System Damping                                scalar
    %       [n]:        Number of cycles elapsed                      scalar
    %
    % Output
    % ----------
    %       [C]:        Damping Coefficient                           scalar

    wn = sqrt(k / m); # Natural frequency
    sigma = log(A_ratio) / n;
    zeta_ = sigma / sqrt((2 * pi)^2 + (sigma)^2);
    C = 2 * m * wn * zeta_; % Damping Constant | Kilograms per Second
end
