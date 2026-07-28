%% compute_ISE_vs_time.m — 计算 ISE 随观测时间的变化
%
% 用于复现 figure_4(b): ISE vs Time（不同 b 值下的曲线）
%
% 依赖：Stage 1 输出的 y_der_<ib>_<ik>.mat 和 parameters.mat
%
% 输出：
%   small_data_for_plotting/setting_<i>/ISE_vs_time.mat
%     ISE_vs_time_mean  [num_b, num_T]  — 对 realization 和节点取平均
%     ISE_vs_time_std   [num_b, num_T]  — 节点×realization 间标准差
%     T_list            [1, num_T]      — 观测时间列表 [s]
%     b                 [1, num_b]      — cfg.b
%
% 用法：
%   compute_ISE_vs_time(cfg, i_setting, T_list)
%   compute_ISE_vs_time(cfg, i_setting)   % 默认 T_list = [100, 200, 500, 1000]

function compute_ISE_vs_time(cfg, i_setting, T_list)

if nargin < 3 || isempty(T_list)
    T_list = [100, 200, 500, 1000];
end
T_list = sort(unique(T_list));

input_dir  = fullfile(cfg.path_data,       sprintf('setting_%d', i_setting));
output_dir = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting));

flag_file = fullfile(output_dir, 'ISE_vs_time_done.flag');
if check_done_flag(flag_file, 'ISE_vs_time', i_setting)
    return
end

params_file = fullfile(input_dir, 'parameters.mat');
if ~exist(params_file, 'file')
    error('[compute_ISE_vs_time] 找不到 %s，请先运行 Stage 1', params_file);
end

S_params = load(params_file, 'params');
num_nodes = S_params.params.num_nodes;
num_bin   = cfg.num_bin;
h         = cfg.h;

num_b = cfg.num_b;
num_k = cfg.num_realizations;
num_T = length(T_list);

ISE_vs_time_mean = zeros(num_b, num_T);
ISE_vs_time_std  = zeros(num_b, num_T);

fprintf('\n[Setting %d] 计算 ISE_vs_time ...\n', i_setting);

for i_b = 1:num_b
    for i_T = 1:num_T
        T = T_list(i_T);
        n_samples = round(T / h);
        if n_samples > cfg.N
            error('[compute_ISE_vs_time] T=%.0f s 超过可用数据长度 %.0f s', T, cfg.N*h);
        end

        ISE_all = zeros(num_k, num_nodes);

        for i_k = 1:num_k
            fn = fullfile(input_dir, sprintf('y_der_%d_%d.mat', i_b, i_k));
            if ~exist(fn, 'file')
                error('[compute_ISE_vs_time] 找不到 %s', fn);
            end

            X = load(fn, 'y_der_w').y_der_w;  % [N, num_nodes]
            X_T = X(1:n_samples, :);

            for i_n = 1:num_nodes
                x_node = X_T(:, i_n);

                mu    = mean(x_node);
                sigma = std(x_node, 1);  % 与 compute_pd 保持一致（N 标准差）

                x_min = min(x_node);
                x_max = max(x_node);
                edges = linspace(x_min, x_max, num_bin + 1);
                centers = (edges(1:end-1) + edges(2:end)) / 2;
                delta_x = centers(2) - centers(1);

                counts = histcounts(x_node, edges, 'Normalization', 'count');
                prob = counts / sum(counts);
                p_empirical = prob / delta_x;

                p_gaussian = gauss_pdf(centers, mu, sigma);

                % ISE = (δx)^2 / M * sum[ (p_empirical - p_gaussian)^2 ]
                ISE_all(i_k, i_n) = (delta_x^2 / num_bin) * ...
                                    sum((p_empirical - p_gaussian).^2);
            end
        end

        ISE_vs_time_mean(i_b, i_T) = mean(ISE_all(:));
        ISE_vs_time_std(i_b, i_T)  = std(ISE_all(:));
    end
    fprintf('  [Setting %d] ISE_vs_time b index %d/%d done\n', i_setting, i_b, num_b);
end

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
b = cfg.b;
save(fullfile(output_dir, 'ISE_vs_time.mat'), ...
    'ISE_vs_time_mean', 'ISE_vs_time_std', 'T_list', 'b');

%% 导出出图用小数据（tab 分隔 .txt）
% 格式：T_list.txt + ISE_vs_time_b<val>.txt（每个 b 一条曲线）
write_row_vector(fullfile(output_dir, 'T_list.txt'), T_list);
for i_b = 1:num_b
    fn = fullfile(output_dir, sprintf('ISE_vs_time_b%.2f.txt', b(i_b)));
    write_row_vector(fn, ISE_vs_time_mean(i_b, :));
end

write_done_flag(flag_file);
fprintf('  [Setting %d] ISE_vs_time 完成\n', i_setting);
end
