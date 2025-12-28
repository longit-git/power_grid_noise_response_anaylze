function [t, y] = rungeKuttaSolver(odeFunc, tspan, y0, h)
    % odeFunc: Function handle for the system of ODEs, e.g., @(t, y) [dy1/dt; dy2/dt; ...]
    % tspan: Time span [t_start, t_end]
    % y0: Initial values of the variables in a column vector
    % h: Time step interval

    % Initialize variables
    t = tspan(1):h:tspan(2);
    n = length(t);
    y = zeros(length(y0), n);
    y(:, 1) = y0;

    % Runge-Kutta method
    for i = 1:n-1
        k1 = odeFunc(t(i), y(:, i));
        k2 = odeFunc(t(i) + h/2, y(:, i) + h*k1/2);
        k3 = odeFunc(t(i) + h/2, y(:, i) + h*k2/2);
        k4 = odeFunc(t(i) + h, y(:, i) + h*k3);

        y(:, i+1) = y(:, i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
