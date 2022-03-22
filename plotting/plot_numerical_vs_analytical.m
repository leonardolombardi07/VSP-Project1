function _ = plot_numerical_vs_analytical(t,
    x_wilson, xdot_wilson, x2dot_wilson,
    x_newmark, xdot_newmark, x2dot_newmark,
    x_central, xdot_central, x2dot_central,
    x_const, xdot_const, x2dot_const,
    x_linear, xdot_linear, x2dot_linear,
    x_analytical, xdot_analytical, x2dot_analytical)

    xlabel_ = "Time [s]";
    title_ = "Comparison of Wilson, Newmark, Central Difference, Constant Approximation, Linear Approximation and Analytical methods";
    legend_ = {"Wilson", "Newmark", "Central", "Constant Approx", "Linear Approx", "Analytical"};

    figure_ = figure("name", title_);
    title(title_);

    subplot(3, 1, 1);
    plot(t, x_wilson,
    t, x_newmark,
    t, x_central,
    t, x_const,
    t, x_linear,
    t, x_analytical);
    xlabel(xlabel_);
    ylabel("Displacement [m]");
    legend(legend_);
    grid on;

    subplot(3, 1, 2);
    plot(t, xdot_wilson,
    t, xdot_newmark,
    t, xdot_central,
    t, xdot_const,
    t, xdot_linear,
    t(:, 1:length(xdot_analytical)), xdot_analytical);
    xlabel(xlabel_);
    ylabel("Velocity [m/s]");
    legend(legend_);
    grid on;

    subplot(3, 1, 3);
    plot(t, x_wilson,
    t, x2dot_newmark,
    t, x2dot_central,
    t, x2dot_const,
    t, x2dot_linear,
    t(:, 1:length(x2dot_analytical)), x2dot_analytical);
    xlabel(xlabel_);
    ylabel("Acceleration [m/s^2]");
    legend(legend_);
    grid on;

    % Comment/Uncomment line below to save plot as png file
    saveas(figure_, strcat("./generated_figures/", title_, ".png"));
end