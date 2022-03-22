function _ = plot_newmark_comparison(t, alfa,
    x_a, xdot_a, x2dot_a, beta_a,
    x_b, xdot_b, x2dot_b, beta_b,
    x_c, xdot_c, x2dot_c, beta_c)

    xlabel_ = "Time [s]";
    title_ = strcat("Comparison of Newmark Method for alfa = ", num2str(alfa), " and different betas");
    beta2display = @(beta_) strcat("Beta: ", num2str(beta_));
    legend_ = {beta2display(beta_a), beta2display(beta_b), beta2display(beta_c)};

    figure_ = figure("name", title_);
    title(title_);

    subplot(3, 1, 1);
    plot(t, x_a, t, x_b, t, x_c);
    xlabel(xlabel_);
    ylabel("Displacement [m]");
    legend(legend_);
    grid on;

    subplot(3, 1, 2);
    plot(t, xdot_a, t, xdot_b, t, xdot_c);
    xlabel(xlabel_);
    ylabel("Velocity [m/s]");
    legend(legend_);
    grid on;

    subplot(3, 1, 3);
    plot(t, x2dot_a, t, x2dot_b, t, x2dot_c);
    xlabel(xlabel_);
    ylabel("Acceleration [m/s^2]");
    legend(legend_);
    grid on;

    % Comment/Uncomment line below to save plot as png file
    saveas(figure_, strcat("./generated_figures/", title_, ".png"));
end
