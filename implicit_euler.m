% Implicit Euler Function
function [t, y] = implicit_euler(f, df, y0, t, dt)
    dim=length(y0);
    y = zeros(dim,length(t));
    y(:,1) = y0;

    for i = 2:length(t)
        G = @(yn) yn - y(i-1) - dt * f(yn,t(i)-t(i-1));
        dG = @(yn) eye(dim) - dt * df(yn,t(i)-t(i-1));

       %Checks if the non-linear equation solvable.  
        if abs(det(dG(y(:,i-1)))) < 2*eps
            t = t(1:i-1);
            y = y(:,1:i-1); 
            fprintf('Time step %d: The equation is not solvable due to singular Jacobian.\n', i);
            return;
         end
        y(:,i) = newton_solver(y(:,i-1), G, dG);
    end
end
