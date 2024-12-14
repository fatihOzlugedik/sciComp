% Main script: euler_analysis.m
% Script for solving Dahlquist's test equation using explicit and implicit Euler methods.

% Clear workspace and command window
clear; clc;

% Parameters
lambda = -7;
y0 = 1;
t_end = 5;
dt_values = [1/2, 1/4, 1/8, 1/16, 1/32];
exact_solution = @(t) exp(-7 * t);

% Time vector for exact solution
fine_t = linspace(0, t_end, 1000);
fine_y = exact_solution(fine_t);

% Results storage
errors_explicit = zeros(1, length(dt_values));
errors_implicit = zeros(1, length(dt_values));
reduction_factors_explicit = zeros(1, length(dt_values));
reduction_factors_implicit = zeros(1, length(dt_values));

% Explicit Euler Analysis
figure;
hold on;
for i = 1:length(dt_values)
    dt = dt_values(i);
    t = 0:dt:t_end;
    [t_exp, y_exp] = explicit_euler(@(y,t) lambda * y, y0, t);

    % Plot explicit Euler solution
    plot(t_exp, y_exp, 'DisplayName', sprintf('Explicit Euler \Delta t = %.3f', dt));

    % Compute error
    exact_y = exact_solution(t);
    errors_explicit(i) = sqrt(dt / t_end * sum((y_exp - exact_y).^2));

    % Compute error reduction factor
    if i > 1
        reduction_factors_explicit(i) = errors_explicit(i-1) / errors_explicit(i);
    end
end

% Plot exact solution
plot(fine_t, fine_y, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
legend;
xlim([0, t_end]); ylim([-1, 1]);
title('Explicit Euler and Exact Solution');
xlabel('Time (t)'); ylabel('x(t)');

% Implicit Euler Analysis
figure;
hold on;
for i = 1:length(dt_values)
    dt = dt_values(i);
    t = 0:dt:t_end;
    [t_imp, y_imp] = implicit_euler(@(y,t) lambda * y, @(y,t) lambda, y0, t, dt);

    % Plot implicit Euler solution
    plot(t_imp, y_imp, 'DisplayName', sprintf('Implicit Euler \Delta t = %.3f', dt));

    % Compute error
    exact_y = exact_solution(t);
    errors_implicit(i) = sqrt(dt / t_end * sum((y_imp - exact_y).^2));

    % Compute error reduction factor
    if i > 1
        reduction_factors_implicit(i) = errors_implicit(i-1) / errors_implicit(i);
    end
end

% Plot exact solution
plot(fine_t, fine_y, 'k-', 'LineWidth', 2, 'DisplayName', 'Exact Solution');
legend;
xlim([0, t_end]); ylim([-1, 1]);
title('Implicit Euler and Exact Solution');
xlabel('Time (t)'); ylabel('x(t)');

% Error Reduction Factor and Stability Analysis
fprintf('Errors and Reduction Factors:\n');
tableExplicitEuler = array2table([dt_values', errors_explicit', reduction_factors_explicit'], ...
    'VariableNames', {'Delta_t', 'Error', 'Reduction_Factor'});
tableImplicitEuler = array2table([dt_values', errors_implicit', reduction_factors_implicit'], ...
    'VariableNames', {'Delta_t', 'Error', 'Reduction_Factor'});

disp('Explicit Euler Results:');
disp(tableExplicitEuler);
disp('Implicit Euler Results:');
disp(tableImplicitEuler);

% Stability Table (Heuristic Stability Check)
stable_explicit = errors_explicit < 1;
stable_implicit = errors_implicit < 1;
stability_table = array2table([dt_values', stable_explicit', stable_implicit'], ...
    'VariableNames', {'Delta_t', 'Explicit_Stable', 'Implicit_Stable'});

disp('Stability Analysis Results:');
disp(stability_table);

% Exercise G: Van-der-Pol Oscillator (x'' − µ(1 - x^2)x' + x = 0)



% Properties for the differential equation
mu=4;
v_0 = [1;1];
t_end = 20;
f_v = @(v, t) [v(2); mu * (1 - v(1)^2) * v(2) - v(1)];
dt=0.1;
t = 0:dt:t_end;
[t_exp, y_exp] = explicit_euler(f_v, v_0, t);

% Plotting Van-der-Pol Oscillator over time.It does not converge with explicit Euler method.
figure("Name", "Van-der-Pol-Oscillator_Explicit", 'NumberTitle', 'off');
plot(0:0.1:5.9,y_exp(:,1:60), "LineWidth", 2);
title('Van-der-Pol-Oscillator-Explicit', 'FontSize', 14);


%Implicit Euler solution
dt=1;
f_v = @(v, t) [v(2); mu * (1 - v(1)^2) * v(2) - v(1)];
Df = @(v,t) [0, 1;  -2 * mu * v(1) * v(2) - 1, mu * (1 - v(1)^2)];
[t_imp, y_imp] = implicit_euler(f_v,Df,v_0, t, dt);
figure("Name", "Van-der-Pol-Oscillator_Impliciy", 'NumberTitle', 'off');
plot(0:0.1:t_end,y_imp, "LineWidth", 2);
title('Van-der-Pol-Oscillator', 'FontSize', 14);
figure("Name", "Van-der-Pol-Oscillator_Impliciy_2", 'NumberTitle', 'off');
plot(y_imp(1,:),y_imp(2,:));