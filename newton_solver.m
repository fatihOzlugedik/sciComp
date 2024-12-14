% Newton Solver
function root = newton_solver(x0, G, dG)
    max_iter = 100;
    tol = 1e-8;
    iter = 0;
    x = x0;
    while iter < max_iter
        x_new = x - dG(x)\G(x);
        if norm(G(x_new)) < tol
            root = x_new;
            return;
        end
        x = x_new;
        iter = iter + 1;
    end
    error('Newton method did not converge');
end
