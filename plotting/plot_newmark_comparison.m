function _ = plot_newmark_comparison(t,
    x_a, xdot_a, x2dot_a,
    x_b, xdot_b, x2dot_b,
    x_c, xdot_c, x2dot_c)
    xlabel_ = "Time [s]";

    figure_ = figure("name", title_);
    title(title_);

    subplot(3, 1, 1);
    plot(t, x_a, t, x_b, t, x_c);
    xlabel(xlabel_);
    ylabel("Displacement [m]");
    grid on;

    subplot(3, 1, 2);
    plot(t, xdot_a, t, xdot_b, t, xdot_c);
    xlabel(xlabel_);
    ylabel("Velocity [m/s]");
    grid on;

    subplot(3, 1, 3);
    plot(t, x2dot_a, t, x2dot_b, t, x2dot_c);
    xlabel(xlabel_);
    ylabel("Acceleration [m/s^2]");
    grid on;

    % Comment/Uncomment line below to save plot as png file
    % saveas(figure_, strcat("./generated_figures/", title_, ".png"));
end
