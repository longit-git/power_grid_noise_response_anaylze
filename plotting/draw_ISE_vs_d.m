%% draw_ISE_vs_d.m — 绘制 ISE 随节点到扰动点图距离的变化（log 纵轴）
%
% 数据来源：data/setting_<i>/pd_storage/ISE_storage.mat（需先运行 compute_ise）
% 输出：figures/setting_<i>/setting_<i>_ISE_d.png
%
% 用法：
%   draw_ISE_vs_d(cfg, i_setting)

function draw_ISE_vs_d(cfg, i_setting)

root_dir = fileparts(fileparts(mfilename('fullpath')));

load(fullfile(cfg.path_key_data, cfg.settings(i_setting).graph + ".mat"), 'G');
load(fullfile(cfg.path_data, sprintf('setting_%d', i_setting), 'pd_storage', 'ISE_storage.mat'), 'ISE_storage');
load(fullfile(cfg.path_key_data, 'pos_noise.mat'), 'pos_noise');

num_nodes = numnodes(G);

%% 各节点到扰动点的图距离
dis = distances(G, pos_noise);
dis = dis(:);
max_dis = max(dis);

%% 对 realization 平均后，按距离分bin，对节点平均
ise = squeeze(mean(ISE_storage, 2));   % [num_b, num_nodes]
ise_vs_dd = zeros(cfg.num_b, max_dis+1, 2);
for i_n = 1:num_nodes
    ise_vs_dd(:, dis(i_n)+1, 1) = ise_vs_dd(:, dis(i_n)+1, 1) + ise(:, i_n);
    ise_vs_dd(:, dis(i_n)+1, 2) = ise_vs_dd(:, dis(i_n)+1, 2) + 1;
end
ise_vs_dm = ise_vs_dd(:, :, 1) ./ ise_vs_dd(:, :, 2);

%% 绘图（所有 b 值）
colortable = [0.90, 0.80, 0.40; 0.40, 0.70, 0.40; 0.19, 0.49, 0.57];
half = floor(cfg.num_b / 2);
colors = [color_gradient(colortable(1,:), colortable(2,:), half);
          color_gradient(colortable(2,:), colortable(3,:), cfg.num_b-half-1);
          colortable(end,:)];

fig = figure('Visible', 'off');
for i_b = 1:cfg.num_b
    plot(0:max_dis, ise_vs_dm(i_b, :), 'Color', colors(i_b, :), 'LineWidth', 2.5);
    hold on;
end
set(gca, 'YScale', 'log');
formatFigure(gca, 20, 1.2, [880 550], '');

%% 保存
save_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, sprintf('setting_%d_ISE_d.png', i_setting));
exportgraphics(fig, save_path, 'Resolution', 300);
close(fig);

fprintf('[draw_ISE_vs_d] Setting %d 已保存到 %s\n', i_setting, save_path);
end
