%% draw_HD_vs_d.m — 绘制 HD 随节点到扰动点图距离的变化
%
% 数据来源：data/setting_<i>/HD.mat（需先运行 compute_HD）
% 输出：figures/setting_<i>/setting_<i>_HD_d.png
%
% 用法：
%   draw_HD_vs_d(cfg, i_setting)

function draw_HD_vs_d(cfg, i_setting)

root_dir = fileparts(fileparts(mfilename('fullpath')));

load(fullfile(cfg.path_key_data, cfg.settings(i_setting).graph + ".mat"), 'G');

% HD.mat：当前版本 compute_HD 输出到 small_data_for_plotting，
% 旧版本输出到 data/ —— 优先读新位置，兼容旧数据目录
hd_file = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting), 'HD.mat');
if ~exist(hd_file, 'file')
    hd_file = fullfile(cfg.path_data, sprintf('setting_%d', i_setting), 'HD.mat');
    fprintf('[draw_HD_vs_d] 使用旧版 HD.mat 位置: %s\n', hd_file);
end
load(hd_file, 'HD');
load(fullfile(cfg.path_key_data, 'pos_noise.mat'), 'pos_noise');

num_nodes = numnodes(G);

%% 各节点到扰动点的图距离
dis = distances(G, pos_noise);
dis = dis(:);
max_dis = max(dis);

%% 按距离分bin，对节点平均
HD_vs_tt = zeros(cfg.num_b, max_dis+1, 2);
for i_n = 1:num_nodes
    HD_vs_tt(:, dis(i_n)+1, 1) = HD_vs_tt(:, dis(i_n)+1, 1) + HD(:, i_n);
    HD_vs_tt(:, dis(i_n)+1, 2) = HD_vs_tt(:, dis(i_n)+1, 2) + 1;
end
HD_vs_d = HD_vs_tt(:, :, 1) ./ HD_vs_tt(:, :, 2);

%% 绘图（b 取 5 个代表值）
colortable = [0.90, 0.80, 0.40; 0.40, 0.70, 0.40; 0.19, 0.49, 0.57];
half = floor(cfg.num_b / 2);
colors = [color_gradient(colortable(1,:), colortable(2,:), half);
          color_gradient(colortable(2,:), colortable(3,:), cfg.num_b-half-1);
          colortable(end,:)];

fig = figure('Visible', 'off');
for i_b = round(linspace(1, cfg.num_b, 5))
    plot(0:max_dis, HD_vs_d(i_b, :), 'Color', colors(i_b, :), 'LineWidth', 5);
    hold on;
end
fontsize(30, 'points');
ylim([0 0.3]);

%% 保存
save_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_path = fullfile(save_dir, sprintf('setting_%d_HD_d.png', i_setting));
exportgraphics(fig, save_path, 'Resolution', 300);
close(fig);

fprintf('[draw_HD_vs_d] Setting %d 已保存到 %s\n', i_setting, save_path);
end
