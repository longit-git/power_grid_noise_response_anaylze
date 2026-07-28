%% compute_pd.m — Stage 2：计算概率密度缓存
%
% 依赖：Stage 1 输出的 y_der_<ib>_<ik>.mat 和 parameters.mat
%
% 输出：data/setting_<i>/pd_storage/pd_storage.mat
%   包含：
%     PD_storage   [num_b, num_k, num_nodes, 2, num_bin]
%                    (:,:,:,1,:) = bin centers
%                    (:,:,:,2,:) = counts（probability 归一化）
%     mu_storage   [num_b, num_k, num_nodes]  — 样本均值 mean(x)
%     sigma_storage[num_b, num_k, num_nodes]  — 样本标准差 std(x,1)（N 归一化）
%     edge_limit   [num_b, num_nodes, 2]      — 全局 bin 边界（供参考）
%     params       — 本 setting 参数结构体（从 parameters.mat 继承）
%
% 每个 setting 运行结束写 pd_storage_done.flag

function compute_pd(cfg, i_setting)

input_dir  = fullfile(cfg.path_data, sprintf('setting_%d', i_setting));
params     = load_params(input_dir);
num_nodes  = params.num_nodes;
num_b      = cfg.num_b;
num_k      = cfg.num_realizations;
num_bin    = cfg.num_bin;

flag_file = fullfile(input_dir, 'pd_storage', 'pd_storage_done.flag');
if check_done_flag(flag_file, 'pd_storage', i_setting)
    return
end

fprintf('\n[Setting %d] 开始计算 pd_storage（num_nodes=%d）\n', i_setting, num_nodes);

%% ── Stage 2-A：计算全局 bin 边界（跨所有实现取 min/max）────
% 目的：同一 b 值下所有实现用相同 bin，使直方图可直接比较（Hellinger 前提）
edge_limit = zeros(num_b, num_nodes, 2);
max1 = zeros(num_k, num_nodes);
min1 = zeros(num_k, num_nodes);

for i_b = 1:num_b
    for i_k = 1:num_k
        fn   = fullfile(input_dir, sprintf('y_der_%d_%d.mat', i_b, i_k));
        X    = load(fn, 'y_der_w').y_der_w;
        max1(i_k,:) = max(X, [], 1);
        min1(i_k,:) = min(X, [], 1);
    end
    edge_limit(i_b,:,1) = min(min1, [], 1);
    edge_limit(i_b,:,2) = max(max1, [], 1);
    fprintf('  Stage2-A: b index %d/%d done\n', i_b, num_b);
end

%% ── Stage 2-B：计算直方图 + fitdist ──────────────────────
PD_storage    = zeros(num_b, num_k, num_nodes, 2, num_bin);
mu_storage    = zeros(num_b, num_k, num_nodes);
sigma_storage = zeros(num_b, num_k, num_nodes);

parfor i_b = 1:num_b
    local_PD    = zeros(num_k, num_nodes, 2, num_bin);
    local_mu    = zeros(num_k, num_nodes);
    local_sigma = zeros(num_k, num_nodes);

    for i_k = 1:num_k
        fn = fullfile(input_dir, sprintf('y_der_%d_%d.mat', i_b, i_k));
        X  = load(fn, 'y_der_w').y_der_w;   % T×num_nodes

        for i_n = 1:num_nodes
            x_node = X(:, i_n);

            % 直方图（共享 bin 边界，probability 归一化）
            edges = linspace(edge_limit(i_b,i_n,1), edge_limit(i_b,i_n,2), num_bin+1);
            counts = histcounts(x_node, edges, 'Normalization', 'probability');
            centers = (edges(1:end-1) + edges(2:end)) / 2;

            local_PD(i_k, i_n, 1, :) = centers;
            local_PD(i_k, i_n, 2, :) = counts;

            % fitdist MLE 拟合（从原始数据点，精度最高）
            local_mu(i_k, i_n)    = mean(x_node);
            local_sigma(i_k, i_n) = std(x_node,1);
        end
    end

    PD_storage(i_b,:,:,:,:) = local_PD;
    mu_storage(i_b,:,:)     = local_mu;
    sigma_storage(i_b,:,:)  = local_sigma;

    fprintf('  Stage2-B: b index %d/%d done\n', i_b, num_b);
end

%% ── 保存 ──────────────────────────────────────────────────
out_dir = fullfile(input_dir, 'pd_storage');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

save(fullfile(out_dir, 'pd_storage.mat'), ...
    'PD_storage', 'mu_storage', 'sigma_storage', 'edge_limit', 'params', '-v7.3');

% 写完成标记
write_done_flag(flag_file);
fprintf('  [Setting %d] pd_storage 完成\n', i_setting);
end

%% ── 辅助：加载 parameters.mat ───────────────────────────
function params = load_params(input_dir)
s = load(fullfile(input_dir, 'parameters.mat'), 'params');
params = s.params;
end
