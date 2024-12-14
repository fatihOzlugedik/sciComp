% Explicit Euler Function
function [t, y] = explicit_euler(f, y0, t)
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    for i = 2:length(t)
        dt = t(i) - t(i-1);
        y(:,i) = y(:,i-1) + dt * f(y(:,i-1),t(i-1));
    end
end

