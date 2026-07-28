%% compute_HD_vs_M.m — 计算 Hellinger distance 随 realization 数量 M 的变化
%
% 用于复现 figure_3(b): H vs M（不同 b 值下的曲线）
%
% 依赖：Stage 2 输出的 pd_storage.mat（需要 PD_storage）
%
% 输出：
%   small_data_for_plotting/setting_<i>/HD_vs_M.mat
%     HD_vs_M_mean  [num_b, num_M]  — 对节点取平均后的 HD
%     HD_vs_M_std   [num_b, num_M]  — 节点间标准差
%     M_list        [1, num_M]      — realization 数量列表
%     b             [1, num_b]      — cfg.b
%
% 用法：
%   compute_HD_vs_M(cfg, i_setting, M_list)
%   compute_HD_vs_M(cfg, i_setting)   % 默认 M_list = 2:num_realizations

function compute_HD_vs_M(cfg, i_setting, M_list)

if nargin < 3 || isempty(M_list)
    M_list = 2:cfg.num_realizations;
end
M_list = unique(M_list);
M_list(M_list < 2) = [];

input_dir  = fullfile(cfg.path_data,       sprintf('setting_%d', i_setting));
output_dir = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting));

flag_file = fullfile(output_dir, 'HD_vs_M_done.flag');
if check_done_flag(flag_file, 'HD_vs_M', i_setting)
    return
end

pd_file = fullfile(input_dir, 'pd_storage', 'pd_storage.mat');
if ~exist(pd_file, 'file')
    error('[compute_HD_vs_M] 找不到 %s，请先运行 Stage 2 (compute_pd)', pd_file);
end

fprintf('\n[Setting %d] 计算 HD_vs_M ...\n', i_setting);
S = load(pd_file, 'PD_storage');

% PD_storage: [num_b, num_k, num_nodes, 2, num_bin]
[num_b, num_k, num_nodes, ~, ~] = size(S.PD_storage);
counts_all = squeeze(S.PD_storage(:,:,:,2,:));  % [num_b, num_k, num_nodes, num_bin]

M_list(M_list > num_k) = [];
num_M = length(M_list);

if num_M == 0
    error('[compute_HD_vs_M] M_list 为空或全部超出可用 realization 数量 %d', num_k);
end

HD_vs_M_mean = zeros(num_b, num_M);
HD_vs_M_std  = zeros(num_b, num_M);

for i_b = 1:num_b
    PD_b = squeeze(counts_all(i_b, :, :, :));  % [num_k, num_nodes, num_bin]

    for i_M = 1:num_M
        M = M_list(i_M);
        PD_M = PD_b(1:M, :, :);  % [M, num_nodes, num_bin]

        HD_node = zeros(num_nodes, 1);
        n_pairs = M * (M - 1) / 2;

        for i_n = 1:num_nodes
            PD_node = squeeze(PD_M(:, i_n, :));  % [M, num_bin]
            acc = 0;
            for i_k = 1:M-1
                for j_k = i_k+1:M
                    P = PD_node(i_k, :)';
                    Q = PD_node(j_k, :)';
                    acc = acc + hellinger(P, Q);
                end
            end
            HD_node(i_n) = acc / n_pairs;
        end

        HD_vs_M_mean(i_b, i_M) = mean(HD_node);
        HD_vs_M_std(i_b, i_M)  = std(HD_node);
    end
    fprintf('  [Setting %d] HD_vs_M b index %d/%d done\n', i_setting, i_b, num_b);
end

if ~exist(output_dir, 'dir'), mkdir(output_dir); end
b = cfg.b;
save(fullfile(output_dir, 'HD_vs_M.mat'), ...
    'HD_vs_M_mean', 'HD_vs_M_std', 'M_list', 'b');

%% 导出出图用小数据（tab 分隔 .txt）
% 格式：M_list.txt + HD_vs_M_b<val>.txt（每个 b 一条曲线）
write_row_vector(fullfile(output_dir, 'M_list.txt'), M_list);
for i_b = 1:num_b
    fn = fullfile(output_dir, sprintf('HD_vs_M_b%.2f.txt', b(i_b)));
    write_row_vector(fn, HD_vs_M_mean(i_b, :));
end

write_done_flag(flag_file);
fprintf('  [Setting %d] HD_vs_M 完成\n', i_setting);
end
