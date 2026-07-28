%% draw_ISE_vs_b.m — 绘制单 Y 轴 ISE(b) 图（带误差阴影）
%
% 数据来源：small_data_for_plotting/ISE/setting_<i>/pd_storage/
% 输出：figures/setting_<i>/setting_<i>_ISE_b.png
%
% 用法：
%   draw_ISE_vs_b(cfg, i_setting)

function draw_ISE_vs_b(cfg, i_setting)

root_dir = fileparts(fileparts(mfilename('fullpath')));

fig = figure('Visible', 'off');

% --- 加载数据 ---
ise_dir  = fullfile(cfg.path_small_data, 'ISE', sprintf('setting_%d', i_setting), 'pd_storage');
b        = load(fullfile(ise_dir, 'b_values.txt'));
ise_mean = load(fullfile(ise_dir, 'ise_mean.txt'));
ise_std  = load(fullfile(ise_dir, 'ise_std.txt'));

color1 = [0.37, 0.64, 0.93];  % 浅蓝色
plot_with_shaded_error(b, ise_mean, ise_std, color1);
ylim([0 5e-4]);
hold on;

% --- 图形美化 ---
xlabel('b');
ylabel('ISE');
title(sprintf('Setting %d', i_setting));
grid off;
fontsize(20, 'points');

% --- 保存图像为PNG (高DPI) ---
save_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

formatFigure(gca, 20, 1.2, [880 550], '');
save_path = fullfile(save_dir, sprintf('setting_%d_ISE_b.png', i_setting));
exportgraphics(fig, save_path, 'Resolution', 300);
close(fig);

fprintf('[draw_ISE_vs_b] Setting %d 已保存到 %s\n', i_setting, save_path);
end
