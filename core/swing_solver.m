function [t, y_der, evec, eval] = swing_solver(edgelist, noise, pos_noise, P, K, alpha, h)
% SWING_SOLVER  求解 swing equation（二阶 Kuramoto 模型），返回频率偏差序列
%
% 模型：
%   dω_i/dt = P_i(t) - α·ω_i + Σ_j K_ij·sin(θ_j - θ_i)
%   dθ_i/dt = ω_i
%
% 状态向量：y = [ω_1…ω_N, θ_1…θ_N]（长度 2N）
%
% 输入：
%   edgelist  — num_edges×2，每行 [i,j]
%   noise     — 长度 2T 的噪声列向量（半步长 h/2 分辨率：RK4 中点求值会用到全部 2T 点）
%   pos_noise — 噪声注入节点索引
%   P         — num_nodes×1 基础功率向量
%   K         — num_nodes×num_nodes 耦合强度矩阵
%   alpha     — 阻尼系数
%   h         — 时间步长
%
% 输出：
%   t      — 1×T 时间向量
%   y_der  — T×num_nodes 频率偏差矩阵
%   evec   — 平衡点线性化拉普拉斯矩阵特征向量
%   eval   — 对应特征值（对角矩阵）

num_nodes = size(K, 1);
T = length(noise) / 2;

%% 求平衡点（无噪声稳态）
y0   = zeros(2*num_nodes, 1);
opts = optimoptions('fsolve','Display','off','TolFun',1e-10,'TolX',1e-10);
y0   = fsolve(@(y) swing_ode(K, P, alpha, y), y0, opts);

%% 带噪声 ODE 积分
P_noisy    = @(t) inject_noise(P, noise, pos_noise, h, t);
ode_handle = @(t,y) swing_ode(K, P_noisy(t), alpha, y);
[t, y]     = rk4_solver(ode_handle, [0, T*h - h], y0, h);

%% 提取输出
y_der = y(1:num_nodes, :).';   % T×num_nodes

%% 平衡点处线性化拉普拉斯特征分解
theta_eq     = y0(num_nodes+1:end);
[evec, eval] = eig(lap_matrix(theta_eq, K));

end

%% ════════════════════════════════════════════════════════════
%  内部函数
%% ════════════════════════════════════════════════════════════

function dY = swing_ode(K, P, alpha, y)
% 向量化 swing equation
N   = length(P);
om  = y(1:N);
th  = y(N+1:2*N);
dth = th.' - th;                         % N×N：dth(i,j)=θ_j-θ_i
coupling = sum(K .* sin(dth), 2);        % N×1
dY = zeros(2*N, 1);
dY(1:N)     = P - alpha.*om + coupling;  % dω/dt
dY(N+1:2*N) = om;                        % dθ/dt
end

function P_out = inject_noise(P, noise, pos_noise, h, t)
idx              = round(t / (h/2)) + 1;
P_out            = P;
P_out(pos_noise) = P_out(pos_noise) + noise(idx);
end

function L = lap_matrix(theta, K)
dth = theta.' - theta;
C   = K .* cos(dth);
L   = diag(sum(C,2)) - C;
end
