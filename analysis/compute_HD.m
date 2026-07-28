%% compute_HD.m — Stage 3a：从 pd_storage 计算 Hellinger distance
%
% 依赖：compute_pd 输出的 pd_storage.mat
%
% 输出：
%   small_data_for_plotting/setting_<i>/HD.mat   — HD [num_b, num_nodes]
%
% HD(b,n) = 平均 Hellinger distance，对所有实现对 (i_k, j_k) 取均值
%   HD = (1/C(K,2)) * Σ_{i<j} H(P_i, P_j)
%
%   H(P,Q) = ||sqrt(P)-sqrt(Q)||_2 / sqrt(2)
%
% 注意：这里的 P/Q 是 probability（直方图），不是 pdf，
%       Hellinger 公式对两者形式相同（归一化到1即可），与原版一致。

function compute_HD(cfg, i_setting)

input_dir  = fullfile(cfg.path_data,       sprintf('setting_%d', i_setting));
output_dir = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting));

flag_file = fullfile(output_dir, 'HD_done.flag');
if check_done_flag(flag_file, 'HD', i_setting)
    return
end

%% 加载 pd_storage（只需要 counts 那一层）
pd_file = fullfile(input_dir, 'pd_storage', 'pd_storage.mat');
fprintf('\n[Setting %d] 加载 pd_storage（HD 计算）...\n', i_setting);
S = load(pd_file, 'PD_storage');

% PD_storage: [num_b, num_k, num_nodes, 2, num_bin]
[num_b, num_k, num_nodes, ~, ~] = size(S.PD_storage);
% 取 counts 层（index=2）
counts_all = squeeze(S.PD_storage(:,:,:,2,:));
% 现在 counts_all: [num_b, num_k, num_nodes, num_bin]

n_pairs = num_k * (num_k - 1) / 2;   % C(K,2)

HD = zeros(num_b, num_nodes);

for i_b = 1:num_b
    % 取出当前 b 值的所有实现：[num_k, num_nodes, num_bin]
    PD_b = squeeze(counts_all(i_b, :, :, :));

    for i_n = 1:num_nodes
        % 取出当前节点所有实现的直方图：[num_k, num_bin]
        PD_node = squeeze(PD_b(:, i_n, :));

        acc = 0;
        for i_k = 1:num_k-1
            for j_k = i_k+1:num_k
                P = PD_node(i_k, :)';
                Q = PD_node(j_k, :)';
                acc = acc + hellinger(P, Q);
            end
        end
        HD(i_b, i_n) = acc / n_pairs;
    end
    fprintf('  [Setting %d] HD b index %d/%d done\n', i_setting, i_b, num_b);
end

%% 保存
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
save(fullfile(output_dir, 'HD.mat'), 'HD');

%% 导出出图用小数据（.mat + tab 分隔 .txt）
HD_mean = mean(HD, 2);            % [num_b x 1]
HD_std  = std(HD, 0, 2);          % [num_b x 1]
b_vec   = cfg.b(:);               % [num_b x 1]

sd_dir = fullfile(cfg.path_small_data, 'data_chi2_R2_H', sprintf('setting_%d', i_setting));
if ~exist(sd_dir, 'dir'), mkdir(sd_dir); end

% .mat
save(fullfile(sd_dir, 'HD_summary.mat'), 'b_vec', 'HD_mean', 'HD_std');

% .txt（tab 分隔）
write_col_vector(fullfile(sd_dir, 'b.txt'),        b_vec);
write_col_vector(fullfile(sd_dir, 'HD_mean.txt'),  HD_mean);
write_col_vector(fullfile(sd_dir, 'HD_std.txt'),   HD_std);

write_done_flag(flag_file);
fprintf('  [Setting %d] HD 完成\n', i_setting);
end
