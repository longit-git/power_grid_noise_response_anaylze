%% draw_ise.m — 绘制 ISE 与 HD 的双 Y 轴图（带误差阴影）
%
% 每个 setting 生成一张图：
%   左 Y 轴：ISE（浅蓝色）
%   右 Y 轴：HD / Hellinger distance（橙色）
%
% 依赖（全部来自 small_data_for_plotting）：
%   ISE/setting_<i>/pd_storage/b_values.txt, ise_mean.txt, ise_std.txt
%   data_chi2_R2_H/setting_<i>/HD_mean.txt, HD_std.txt
%
% 用法：
%   draw_ise(cfg)
%   draw_ise(cfg, i_setting_list)   % 默认 1:length(cfg.settings)

function draw_ise(cfg, i_setting_list)

if nargin < 2
    i_setting_list = 1:length(cfg.settings);
end

root_dir = fileparts(fileparts(mfilename('fullpath')));

for i_setting = i_setting_list
    figure('Visible', 'off');

    % --- ISE 数据 ---
    ise_dir = fullfile(cfg.path_small_data, 'ISE', ...
                       sprintf('setting_%d', i_setting), 'pd_storage');
    b        = load(fullfile(ise_dir, 'b_values.txt'));
    ise_mean = load(fullfile(ise_dir, 'ise_mean.txt'));
    ise_std  = load(fullfile(ise_dir, 'ise_std.txt'));

    color1 = [0.37, 0.64, 0.93];  % 浅蓝色
    color2 = [0.90, 0.48, 0.13];  % 橙色

    % --- 左侧 Y 轴 (ISE) ---
    yyaxis left
    plot_with_shaded_error(b, ise_mean, ise_std, color1);
    ylim([0 4e-4]);
    hold on;
    ax = gca;
    ax.YColor = color1;
    ylabel('ISE');

    % --- 右侧 Y 轴 (HD) ---
    hd_dir = fullfile(cfg.path_small_data, 'data_chi2_R2_H', ...
                      sprintf('setting_%d', i_setting));
    HD_mean = load(fullfile(hd_dir, 'HD_mean.txt'));
    HD_std  = load(fullfile(hd_dir, 'HD_std.txt'));

    yyaxis right
    plot_with_shaded_error(b, HD_mean, HD_std, color2);
    ax.YColor = color2;
    ylim([0 0.3]);
    ylabel('HD');

    % --- 图形美化 ---
    xlabel('b');
    xlim([0 2]);
    title(sprintf('Setting %d', i_setting));
    grid on;
    fontsize(30, 'points');

    % --- 保存 ---
    save_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    save_path = fullfile(save_dir, sprintf('setting_%d_ise_H.png', i_setting));
    exportgraphics(gcf, save_path, 'Resolution', 300);
    close(gcf);
end

fprintf('[draw_ise] 已完成 %d 个 setting 的图片生成\n', length(i_setting_list));
end
