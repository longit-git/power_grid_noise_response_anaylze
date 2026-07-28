function [t, y] = rk4_solver(ode_func, tspan, y0, h)
% RK4_SOLVER  定步长四阶 Runge-Kutta 积分器
%
% 输入：
%   ode_func — @(t,y) 返回列向量 dy/dt
%   tspan    — [t_start, t_end]
%   y0       — 初始状态列向量
%   h        — 时间步长
%
% 输出：
%   t — 1×n 时间向量
%   y — length(y0)×n 状态矩阵

t = tspan(1) : h : tspan(2);
n = length(t);
y = zeros(length(y0), n);
y(:,1) = y0(:);

for i = 1:n-1
    k1 = ode_func(t(i),       y(:,i));
    k2 = ode_func(t(i) + h/2, y(:,i) + h/2 * k1);
    k3 = ode_func(t(i) + h/2, y(:,i) + h/2 * k2);
    k4 = ode_func(t(i) + h,   y(:,i) + h   * k3);
    y(:,i+1) = y(:,i) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end
end
