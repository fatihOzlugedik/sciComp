% Implicit Euler Function
function [t, y] = implicit_euler(f, df, y0, t, dt)
    y = zeros(size(t));
    y(1) = y0;
    for i = 2:length(t)
        G = @(yn) yn - y(i-1) - dt * f(yn);
        dG = @(yn) 1 - dt * df(yn);
        y(i) = newton_solver(y(i-1), G, dG);
    end
end
