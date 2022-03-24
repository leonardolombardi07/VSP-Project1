function [x, xdot, x2dot] = solve_by_laplace(times, m, c, k)
    % Analytical Method using Laplace Transform
    %--------------------------------------------------------------------------
    % Solves a 1-DOF system with an equivalent mass "m", spring stiffness "k" and damping
    % coeffiecient "c" subjected to an external force F(t) using Laplace Transform.
    % Returns the displacement, velocity and acceleration of the system with
    % respect to an inertial frame of reference.
    %
    % Input
    % ----------
    %       (To avoid variable name clashes, we don't use the name "t". See below)
    %       [times]:   Time Vector              [n, 1]
    %
    % Output
    % ----------
    %       [x]:        Displacement Response   [n,1]
    %       [xdot]:     Velocity                [n,1]
    %       [x2dot]:    Acceleration            [n,1]

    % "syms" creates symbolic variables. They are special variables that can be manipulated
    % by MATLAB in order to find analytical expressions
    syms x t s;

    % We are solving a system on the form of: mx'' + cx' + kx = F(t)
    % The variable "F" represents the exciting force F(t)
    F = 200 + (100 - 500 * t) * heaviside(t - 0.2) + (500 * t - 300) * heaviside(t - 0.4);

    % Calculate the laplace transform of F(t)
    L_right = laplace(F, t, s); % L{F(t)}(s)

    % Calculate, manually, the laplace transform of mx'' + cx' + kx
    % considering x(0) = 0 and x'(0) = 0
    L_left = (m * s^2 + c * s + k); % L{mx'' + cx' + kx}(s)

    % The laplace transform of x(t) then becomes:
    L = L_right / L_left; % L{x(t)}(s)

    % To find x(t), calculate the inverse laplace transform of L{x(t)}(s)
    % ilaplace returns a symbolic variable representing x(t)
    x_symbolic = ilaplace(L);
    xdot_symbolic = diff(x_symbolic, t); % Taking the derivative to find velocity
    x2dot_symbolic = diff(xdot_symbolic, t); % Taking the derivative to find acceleration

    % Convert x(t) into numeric values (i.e, the vectors x, xdot and x2dot)
    % NOTE:
    % The symbolic variable "x_symbolic" is a function that has the symbolic variable
    % "t" as parameter. To substitute "t" for the numeric values in the vector "times",
    % we apply the function "subs": subs(symbolic_function, t,  times).
    % But this still returns symbolic variables. If, for example, symbolic_function = exp(t),
    % then subs(symbolic_function, t,  [2, 3]) = [exp(2), exp(3)]),
    % where exp(2) and exp(3) are still symbolic. To convert this symbolic variables to actual
    % numbers, apply double(subs((symbolic_function, t,  [2, 3]))) to get, in this example,
    % [7.3891   20.0855]
    x = double(subs(x_symbolic, t, times));
    xdot = double(subs(xdot_symbolic, t, times));
    x2dot = double(subs(x2dot_symbolic, t, times));
end
