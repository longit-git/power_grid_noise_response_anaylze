%% draw_HD_vs_b.m — 绘制单 Y 轴 HD(b) 图（带误差阴影）
%
% 数据来源：small_data_for_plotting/data_chi2_R2_H/setting_<i>/
% 输出：figures/setting_<i>/setting_<i>_HD_b.png
%
% 用法：
%   draw_HD_vs_b(cfg, i_setting)

function draw_HD_vs_b(cfg, i_setting)

root_dir = fileparts(fileparts(mfilename('fullpath')));

fig = figure('Visible', 'off');

% --- 加载数据 ---
hd_dir = fullfile(cfg.path_small_data, 'data_chi2_R2_H', sprintf('setting_%d', i_setting));
b      = load(fullfile(hd_dir, 'b.txt'));
H_mean = load(fullfile(hd_dir, 'HD_mean.txt'));
H_std  = load(fullfile(hd_dir, 'HD_std.txt'));

color2 = [0.90, 0.48, 0.13];  % 橙色
plot_with_shaded_error(b, H_mean, H_std, color2);
ylim([0 0.4]);

% --- 图形美化 ---
xlabel('b');
ylabel('H');
title(sprintf('Setting %d', i_setting));
grid off;
fontsize(20, 'points');

% --- 保存图像为PNG (高DPI) ---
save_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

formatFigure(gca, 20, 1.2, [880 550], '');
save_path = fullfile(save_dir, sprintf('setting_%d_HD_b.png', i_setting));
exportgraphics(fig, save_path, 'Resolution', 300);
close(fig);

fprintf('[draw_HD_vs_b] Setting %d 已保存到 %s\n', i_setting, save_path);
end
