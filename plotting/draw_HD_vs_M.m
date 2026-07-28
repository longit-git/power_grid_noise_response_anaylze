%% draw_HD_vs_M.m — 复现 figure_3: Hellinger distance 的 stationarity 分析
%
% (a) H vs b：30 次实现两两之间的平均 Hellinger distance，对节点平均
% (b) H vs M：不同 b 值下，Hellinger distance 随 realization 数量 M 的变化
%
% 依赖（全部来自 small_data_for_plotting）：
%   data_chi2_R2_H/setting_<i>/b.txt, HD_mean.txt, HD_std.txt  → panel (a)
%   setting_<i>/M_list.txt, HD_vs_M_b<val>.txt                  → panel (b)
%
% 用法：
%   draw_HD_vs_M(cfg, i_setting)
%   draw_HD_vs_M(cfg)   % 默认 i_setting = 10

function draw_HD_vs_M(cfg, i_setting)

if nargin < 2, i_setting = 10; end

root_dir = fileparts(fileparts(mfilename('fullpath')));

b = cfg.b;

%% ── 颜色配置 ─────────────────────────────────────────────
% 连续 colormap：b=0 黄色 → b=1 绿色 → b=2 深蓝
n_color = 256;
cmap_points = [1.00, 0.90, 0.00;   % b=0  黄
               0.20, 0.60, 0.30;   % b=1  绿
               0.05, 0.15, 0.45];  % b=2  深蓝
customMap = interp1([0; 1; 2], cmap_points, linspace(0, 2, n_color)');

b_sel = [0, 0.5, 1, 1.5, 2];
[~, idx_b_sel] = ismember(b_sel, b);
idx_b_sel(idx_b_sel == 0) = [];

%% ═════════════════════════════════════════════════════════
%  Panel (a): H vs b
%% ═════════════════════════════════════════════════════════
sd_dir_a = fullfile(cfg.path_small_data, 'data_chi2_R2_H', sprintf('setting_%d', i_setting));
b_a      = load(fullfile(sd_dir_a, 'b.txt'));
H_mean   = load(fullfile(sd_dir_a, 'HD_mean.txt'));
H_std    = load(fullfile(sd_dir_a, 'HD_std.txt'));

fig = figure('Visible', 'off', 'Position', [100 100 1200 500]);

subplot(1, 2, 1);
hold on;
fill([b_a', fliplr(b_a')], [H_mean'+H_std', fliplr(H_mean'-H_std')], ...
     [0.95 0.50 0.30], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(b_a, H_mean, 'Color', [0.95 0.50 0.30], 'LineWidth', 3);
hold off;
xlabel('b');
ylabel('H');
xlim([0 2]);
ylim([0 0.25]);
box on;
title('(a)');
fontsize(18, 'points');

% ── 左上角 semi-log inset ──
axes('Position', [0.18 0.62 0.12 0.25]);
semilogy(b_a, H_mean, 'Color', [0.95 0.50 0.30], 'LineWidth', 1.5);
hold on;
fill([b_a', fliplr(b_a')], [H_mean'+H_std', fliplr(H_mean'-H_std')], ...
     [0.95 0.50 0.30], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;
xlim([0 2]);
ylim([1e-2 0.25]);
box on;
fontsize(10, 'points');

%% ═════════════════════════════════════════════════════════
%  Panel (b): H vs M
%% ═════════════════════════════════════════════════════════
sd_dir_b = fullfile(cfg.path_small_data, sprintf('setting_%d', i_setting));
M_list   = load(fullfile(sd_dir_b, 'M_list.txt'));

subplot(1, 2, 2);
hold on;
for k = 1:length(idx_b_sel)
    i_b = idx_b_sel(k);
    fn = fullfile(sd_dir_b, sprintf('HD_vs_M_b%.2f.txt', b(i_b)));
    HD_curve = load(fn);
    color_idx = round((b(i_b) / 2) * (n_color - 1)) + 1;
    plot(M_list, HD_curve, 'Color', customMap(color_idx, :), 'LineWidth', 2.5);
end
hold off;
set(gca, 'XScale', 'log');
xlabel('M');
ylabel('H');
xlim([min(M_list) max(M_list)]);
ylim([0 0.5]);
xticks([2 10 20 30]);
box on;
title('(b)');
fontsize(18, 'points');

% 颜色条 / 图例
colormap(customMap);
cb = colorbar('Location', 'eastoutside');
clim([0 2]);
cb.Ticks = b_sel;
cb.TickLabels = string(b_sel);
cb.Label.String = 'b';

%% ── 保存 ────────────────────────────────────────────────
fig_dir = fullfile(root_dir, 'figures', sprintf('setting_%d', i_setting));
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
exportgraphics(fig, fullfile(fig_dir, 'figure_3_HD_stationarity.png'), 'Resolution', 300);
close(fig);

fprintf('[draw_HD_vs_M] Setting %d 已保存到 %s\n', i_setting, fig_dir);
end
