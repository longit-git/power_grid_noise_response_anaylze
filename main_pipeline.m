%% main_pipeline.m — 完整流程入口
%
% Stage 1 │ 数值模拟：生成有色噪声 + RK4 求解 swing equation → y_der / noise
% Stage 2 │ 概率密度缓存：直方图 + fitdist → pd_storage
% Stage 3 │ 指标计算 + 出图用小数据导出：HD、ISE → small_data_for_plotting
% Stage 4 │ HD_vs_M / ISE_vs_time 额外数据（主 setting）
% Stage 5 │ 生成非 VPM 图片
%
% 每个 stage 每个 setting 运行结束写 *_done.flag，重新运行自动跳过已完成项。
% 所有参数集中在 config.m，其他脚本不硬编码任何参数。
%
% 使用方式：
%   直接运行本脚本（跑全部 settings 的完整流程）
%   或注释掉不需要的 stage，单独重跑某一阶段
%
% 注意：
%   - VPM 相关计算完全由 main_vpm.m 负责，本脚本不涉及。
%   - create_random_power_grid.m / draw_network.m 为手动工具，不在本流程中。

clear; clc;

%% ── 初始化 ────────────────────────────────────────────────
root_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root_dir, 'core')));
addpath(genpath(fullfile(root_dir, 'analysis')));
addpath(genpath(fullfile(root_dir, 'plotting')));

run(fullfile(root_dir, 'config.m'));   % → cfg

load(fullfile(cfg.path_key_data, 'pos_noise.mat'));
cfg.pos_noise = pos_noise;

num_settings = length(cfg.settings);
PRIMARY_SETTING = 10;   % HD_vs_M / ISE_vs_time 图对应参数：alpha=1, G20, default P, K=30

fprintf('== powergrid-v2  |  %d settings  |  %d b值  |  %d 次实现 ==\n', ...
    num_settings, cfg.num_b, cfg.num_realizations);

%% 启动 parpool（若尚未运行）
if isempty(gcp('nocreate'))
    parpool(cfg.num_workers);
end

tic;

%% ════════════════════════════════════════════════════════════
%  Stage 1：数值模拟
%% ════════════════════════════════════════════════════════════
fprintf('\n======== Stage 1: 数值模拟 ========\n');

for i_setting = 1:num_settings
    s          = cfg.settings(i_setting);
    output_dir = fullfile(cfg.path_data, sprintf('setting_%d', i_setting));
    flag_file  = fullfile(output_dir, 'simulation_done.flag');

    if exist(flag_file, 'file')
        fprintf('[Setting %d] 模拟已完成，跳过\n', i_setting);
        continue
    end

    fprintf('\n[Setting %d/%d]  graph=%s  alpha=%.1f  P=%s  K=%s\n', ...
        i_setting, num_settings, s.graph, s.alpha, s.P_mode, s.K_mode);

    %% 加载图
    fn_G = fullfile(cfg.path_key_data, [s.graph, '.mat']);
    load(fn_G, 'G');
    edgelist  = G.Edges.EndNodes;
    num_nodes = numnodes(G);
    num_edges = size(edgelist, 1);

    %% 构建 P 和 K
    P = make_P(num_nodes, s.P_mode);
    K = make_K(edgelist, num_nodes, num_edges, s.K_mode);

    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    %% 遍历 b 值
    for i_b = 1:cfg.num_b
        b_val     = cfg.b(i_b);
        noise_all = zeros(2*cfg.N, cfg.num_realizations);
        y_der_all = zeros(cfg.N, num_nodes, cfg.num_realizations);

        parfor i_k = 1:cfg.num_realizations
            % 生成有色噪声并归一化
            nz = colored_noise(2*cfg.N, b_val, cfg.noise_w, cfg.Fs);
            nz = (cfg.noise_amplitude / std(nz)) * nz;

            % RK4 求解
            [~, y_der, ~, ~] = swing_solver( ...
                edgelist, nz, cfg.pos_noise, P, K, s.alpha, cfg.h);

            noise_all(:, i_k)    = nz;
            y_der_all(:, :, i_k) = y_der;
        end

        % 串行写盘（parfor 结束后）
        for i_k = 1:cfg.num_realizations
            y_der_w = y_der_all(:, :, i_k);  %#ok<NASGU>
            noise_w = noise_all(:, i_k);      %#ok<NASGU>
            save(fullfile(output_dir, sprintf('y_der_%d_%d.mat',  i_b, i_k)), 'y_der_w');
            save(fullfile(output_dir, sprintf('noise_%d_%d.mat',  i_b, i_k)), 'noise_w');
        end

        fprintf('  b=%.1f (%d/%d) done\n', b_val, i_b, cfg.num_b);
    end

    %% 保存参数记录
    params             = s;
    params.N           = cfg.N;
    params.h           = cfg.h;
    params.b           = cfg.b;
    params.num_b       = cfg.num_b;
    params.num_realizations = cfg.num_realizations;
    params.num_bin     = cfg.num_bin;
    params.pos_noise   = cfg.pos_noise;
    params.num_nodes   = num_nodes;
    params.num_edges   = num_edges;  %#ok<NASGU>
    save(fullfile(output_dir, 'parameters.mat'), 'params');

    fclose(fopen(flag_file, 'w'));
    fprintf('[Setting %d] 模拟完成\n', i_setting);
end

%% ════════════════════════════════════════════════════════════
%  Stage 2：概率密度缓存（pd_storage）
%% ════════════════════════════════════════════════════════════
fprintf('\n======== Stage 2: pd_storage ========\n');

for i_setting = 1:num_settings
    compute_pd(cfg, i_setting);
end

%% ════════════════════════════════════════════════════════════
%  Stage 3：HD / ISE 计算 + 出图用小数据导出
%% ════════════════════════════════════════════════════════════
fprintf('\n======== Stage 3: HD / ISE 计算与出图小数据导出 ========\n');

for i_setting = 1:num_settings
    compute_HD(cfg, i_setting);
    compute_ise(cfg, i_setting);
end

%% ════════════════════════════════════════════════════════════
%  Stage 4：HD_vs_M / ISE_vs_time 额外数据（主 setting）
%% ════════════════════════════════════════════════════════════
fprintf('\n======== Stage 4: HD_vs_M / ISE_vs_time 额外数据（Setting %d） ========\n', PRIMARY_SETTING);

compute_HD_vs_M(cfg, PRIMARY_SETTING);
compute_ISE_vs_time(cfg, PRIMARY_SETTING);

%% ════════════════════════════════════════════════════════════
%  Stage 5：生成非 VPM 图片
%% ════════════════════════════════════════════════════════════
fprintf('\n======== Stage 5: 图片生成 ========\n');

draw_ise(cfg);
draw_HD_vs_M(cfg, PRIMARY_SETTING);
draw_ISE_vs_time(cfg, PRIMARY_SETTING);
draw_HD_vs_b(cfg, PRIMARY_SETTING);
draw_ISE_vs_b(cfg, PRIMARY_SETTING);
draw_HD_vs_d(cfg, PRIMARY_SETTING);
draw_ISE_vs_d(cfg, PRIMARY_SETTING);

elapsed = toc;
fprintf('\n== 全部完成，总耗时 %.1f 秒（%.1f 分钟）==\n', elapsed, elapsed/60);
