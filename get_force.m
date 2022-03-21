function F = get_force(t)
    F = zeros(length(t), 1);

    for i = 1:length(t)

        if t(i) <= 0.2
            F(i) = 200;
        elseif t(i) > 0.4
            F(i) = 0;
        else
            F(i) = -500 * t(i) + 300;
        end

    end

end
