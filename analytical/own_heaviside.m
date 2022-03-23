function teta = own_heaviside(t)
    % This is a polyfill for the MATLAB heaviside function,
    % as it doesn't work well on Octave
    if t <= 0
        teta = 0;
    else
        teta = 1;
    end

end
