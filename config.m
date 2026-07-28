%% config.m — 全局参数配置
% 所有脚本通过 run(fullfile(root_dir,'config.m')) 加载
% 只改这一个文件，不改其他任何脚本

%% ── 时间序列 ──────────────────────────────────────────────
cfg.N    = 100000;        % 时间序列长度（采样点数）
cfg.h    = 0.01;          % 时间步长 [s]
cfg.Fs   = 1 / cfg.h;    % 采样率 [Hz]

%% ── 噪声 ──────────────────────────────────────────────────
cfg.b               = 0 : 0.1 : 2;   % 噪声颜色指数（0=白噪声, 2=布朗噪声）
cfg.noise_amplitude = 0.3;            % 噪声标准差目标值（相对 P_i 的比例）
cfg.noise_w         = 0.1;            % 均匀白噪声幅度参数

%% ── 统计实验 ──────────────────────────────────────────────
cfg.num_realizations = 30;    % 每个 b 值的重复实现次数
cfg.num_bin          = 35;    % 概率密度直方图 bin 数

%% ── 路径 ──────────────────────────────────────────────────
cfg.path_key_data       = './currently-using-key-data';
cfg.path_data           = './data';
cfg.path_small_data     = './small_data_for_plotting';

%% ── 并行 ──────────────────────────────────────────────────
cfg.num_workers = 6;   % parpool worker 数

%% ── 派生量（勿手动修改）──────────────────────────────────
cfg.num_b = length(cfg.b);

%% ── Settings 定义 ─────────────────────────────────────────
% 字段说明：
%   alpha    : 阻尼系数
%   graph    : 图文件名（currently-using-key-data/<graph>.mat）
%   P_mode   : 'default'  → +1/-1 交替各50%
%              '20_80'    → 20% P=+0.8 / 80% P=-0.2（总功率=0）
%              '40_60'    → 40% P=+0.6 / 60% P=-0.4（总功率=0）
%   K_mode   : 'fixed_30'   → 所有边 K=30
%              'rand_25_35' → 每条边 K~U[25,35]
%              'rand_20_40' → 每条边 K~U[20,40]

cfg.settings(1) = struct('alpha',0.5, 'graph','G20',              'P_mode','default', 'K_mode','fixed_30');
cfg.settings(2) = struct('alpha',1.5, 'graph','G20',              'P_mode','default', 'K_mode','fixed_30');
cfg.settings(3) = struct('alpha',1,   'graph','G_random_50nodes',  'P_mode','default', 'K_mode','fixed_30');
cfg.settings(4) = struct('alpha',1,   'graph','G_random_100nodes', 'P_mode','default', 'K_mode','fixed_30');
cfg.settings(5) = struct('alpha',1,   'graph','G20',              'P_mode','20_80',   'K_mode','fixed_30');
cfg.settings(6) = struct('alpha',1,   'graph','G20',              'P_mode','40_60',   'K_mode','fixed_30');
cfg.settings(7) = struct('alpha',1,   'graph','G20',              'P_mode','default', 'K_mode','rand_25_35');
cfg.settings(8) = struct('alpha',1,   'graph','G20',              'P_mode','default', 'K_mode','rand_20_40');
cfg.settings(9) = struct('alpha',1,   'graph','G_random_500nodes',              'P_mode','default', 'K_mode','fixed_30');
%备注：setting n的已存在映射不建议修改
cfg.settings(10) = struct('alpha',1,   'graph','G20',              'P_mode','default', 'K_mode','fixed_30');
cfg.settings(11) = struct('alpha',1,   'graph','G21',              'P_mode','default', 'K_mode','fixed_30');
cfg.settings(12) = struct('alpha',1,   'graph','G_BA',              'P_mode','default', 'K_mode','fixed_30');
cfg.settings(13) = struct('alpha',1,   'graph','G_chain_10_nodes',              'P_mode','default', 'K_mode','fixed_30');

