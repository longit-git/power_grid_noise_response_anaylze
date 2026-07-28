%% draw_ISE_vs_time.m — 复现 figure_4: response PDF 的 Gaussianity 分析
%
% (a) ISE vs b：Gaussian 拟合误差随噪声颜色 b 的变化，对节点和 realization 平均
% (b) ISE vs Time：不同 b 值下，ISE 随观测时间的变化（log-log）
%
% 依赖（全部来自 small_data_for_plotting）：
%   ISE/setting_<i>/pd_storage/b_values.txt, ise_mean.txt, ise_std.txt  → panel (a)
%   setting_<i>/T_list.txt, ISE_vs_time_b<val>.txt                       → panel (b)
%
% 用法：
%   draw_ISE_vs_time(cfg, i_setting)
%   draw_ISE_vs_time(cfg)   % 默认 i_setting = 10

function draw_ISE_vs_time(cfg, i_setting)

if nargin < 2, i_setting = 10; end

root_dir = fileparts(fileparts(mfilename('fullpath')));

b = cfg.b;

%% ── 颜色配置 ─────────────────────────────────────────────
n_color = 256;
cmap_points = [1.00, 0.90, 0.00;   % b=0  黄
               0.20, 0.60, 0.30;   % b=1  绿
               0.05, 0.15, 0.45];  % b=2  深蓝
customMap = interp1([0; 1; 2], cmap_points, linspace(0, 2, n_color)');

b_sel = [0, 0.5, 1, 1.5, 2];
[~, idx_b_sel] = ismember(b_sel, b);
idx_b_sel(idx_b_sel == 0) = [];

%% ═════════════════════════════════════════════════════════
%  Panel (a): ISE vs b
%% ═════════════════════════════════════════════════════════
sd_dir_a = fullfile(cfg.path_small_data, 'ISE', sprintf('setting_%d', i_setting), 'pd_storage');
b_a      = load(fullfile(sd_dir_a, 'b_values.txt'));
ise_mean = load(fullfile(sd_dir_a, 'ise_mean.txt'));
ise_std  = load(fullfile(sd_dir_a, 'ise_std.txt'));

fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);

subplot(1, 2, 1);
hold on;
fill([b_a', fliplr(b_a')], [ise_mean'+ise_std', fliplr(ise_mean'-ise_std')], ...
     [0.45, 0.63, 0.92], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(b_a, ise_mean, 'Color', [0.45, 0.63, 0.92], 'LineWidth', 3);
hold off;
xlabel('b');
ylabel('ISE');
xlim([0 2]);
ylim([0 max(ise_mean+ise_std) * 1.2]);
box on;
title('(a)');
fontsize(18, 'points');

% ── 左上角 semi-log inset ──
axes('Position', [0.18 0.62 0.12 0.25]);
semilogy(b_a, ise_mean, 'Color', [0.45, 0.63, 0.92], 'LineWidth', 1.5);
hold on;
fill([b_a', fliplr(b_a')], [ise_mean'+ise_std', fliplr(ise_mean'-ise_std')], ...
     [0.45, 0.63, 0.92], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;
xlim([0 2]);
ylim([1e-6 max(ise_mean) * 5]);
box on;
fontsize(10, 'points');

%% ═════════════════════════════════════════════════════════
%  Panel (b): ISE vs Time
%% ═════════════════════════════════════════════════════════
sd_dir_b = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting));
T_list   = load(fullfile(sd_dir_b, 'T_list.txt'));

subplot(1, 2, 2);
hold on;
for k = 1:length(idx_b_sel)
    i_b = idx_b_sel(k);
    fn = fullfile(sd_dir_b, sprintf('ISE_vs_time_b%.2f.txt', b(i_b)));
    ISE_curve = load(fn);
    color_idx = round((b(i_b) / 2) * (n_color - 1)) + 1;
    plot(T_list, ISE_curve, 'Color', customMap(color_idx, :), ...
         'LineWidth', 2.5, 'Marker', 'o', 'MarkerSize', 6);
end
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Time [s]');
ylabel('ISE');
xlim([min(T_list) max(T_list)]);
ylim([1e-6 1e-3]);
xticks(T_list);
box on;
title('(b)');
fontsize(18, 'points');

% 颜色条
colormap(customMap);
cb = colorbar('Location', 'eastoutside');
clim([0 2]);
cb.Ticks = b_sel;
cb.TickLabels = string(b_sel);
cb.Label.String = 'b';

%% ── 保存 ────────────────────────────────────────────────
fig_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
exportgraphics(fig, fullfile(fig_dir, 'figure_4_ISE_timespan.png'), 'Resolution', 300);
close(fig);

fprintf('[draw_ISE_vs_time] Setting %d 已保存到 %s\n', i_setting, fig_dir);
end
