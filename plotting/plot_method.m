function _ = plot_method(t, x, xdot, x2dot, title_)
    xlabel_ = "Time [s]";

    figure_ = figure("name", title_);
    title(title_);

    subplot(3, 1, 1);
    plot(t, x);
    xlabel(xlabel_);
    ylabel("Displacement [m]");
    grid on;

    subplot(3, 1, 2);
    % Ensuring t and xdot have the same length
    plot(t(:, 1:length(xdot)), xdot);
    xlabel(xlabel_);
    ylabel("Velocity [m/s]");
    grid on;

    subplot(3, 1, 3);
    % Ensuring t and x2dot have the same length
    plot(t(:, 1:length(x2dot)), x2dot);
    xlabel(xlabel_);
    ylabel("Acceleration [m/s^2]");
    grid on;

    % Comment/Uncomment line below to save plot as png file
    saveas(figure_, strcat("./generated_figures/", title_, ".png"));
end
