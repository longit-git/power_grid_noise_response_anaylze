function compute_ise(cfg, i_setting)
%COMPUTE_ISE  计算 PD_storage 与 Gaussian 拟合之间的 Integrated Squared Error
%
%   用法：
%       compute_ise(cfg, i_setting)
%
%   输入：
%       cfg        — 配置结构体，需包含 cfg.path_data, cfg.path_small_data,
%                    cfg.num_bin, cfg.b
%       i_setting  — setting 编号
%
%   计算步骤：
%       1. 加载 pd_storage.mat (含 PD_storage, mu_storage, sigma_storage, edge_limit)
%       2. 将 probability → pdf (除以 bin width δx)
%       3. 计算 Gaussian pdf N(mu, sigma^2) 在 bin centers 处的值
%       4. ISE = (δx)^2 / M * sum_{j=1}^{M} [p'(x_j) - p(x_j)]^2
%       5. 保存 ISE_storage 为 .mat 和 .txt，并导出出图用小数据
%
%   输出：
%       data/setting_<i>/pd_storage/ISE_storage.mat  — ISE_storage 数组
%       data/setting_<i>/pd_storage/ISE_storage.txt  — 可读文本表格
%       data/setting_<i>/pd_storage/ISE_done.flag    — 完成标记（重跑自动跳过）
%       small_data_for_plotting/ISE/setting_<i>/pd_storage/ — 出图小数据

    num_bin = cfg.num_bin;
    M = num_bin;

    %% ── 加载 pd_storage.mat ──────────────────────────
    input_dir = fullfile(cfg.path_data, sprintf('setting_%d', i_setting), 'pd_storage');
    mat_file = fullfile(input_dir, 'pd_storage.mat');

    flag_file = fullfile(input_dir, 'ISE_done.flag');
    if check_done_flag(flag_file, 'ISE', i_setting)
        return
    end

    fprintf('[compute_ise] 加载 %s ...\n', mat_file);
    if ~exist(mat_file, 'file')
        error('文件不存在: %s', mat_file);
    end

    data = load(mat_file, 'PD_storage', 'mu_storage', 'sigma_storage', 'edge_limit');
    PD_storage    = data.PD_storage;      % [num_b, num_k, num_nodes, 2, num_bin]
    mu_storage    = data.mu_storage;      % [num_b, num_k, num_nodes]
    sigma_storage = data.sigma_storage;   % [num_b, num_k, num_nodes]
    edge_limit    = data.edge_limit;      % [num_b, num_nodes, 2]

    sz = size(PD_storage);
    num_b = sz(1);
    num_k = sz(2);
    num_nodes = sz(3);

    fprintf('[compute_ise] PD_storage 大小: %s\n', mat2str(sz));
    fprintf('[compute_ise] num_bin=%d, num_b=%d, num_k=%d, num_nodes=%d\n', ...
        num_bin, num_b, num_k, num_nodes);

    %% ── 预分配 ISE_storage ───────────────────────────
    ISE_storage = zeros(num_b, num_k, num_nodes);

    %% ── 逐点计算 ISE ─────────────────────────────────
    fprintf('[compute_ise] 开始计算 ISE ...\n');

    for i_b = 1:num_b
        for i_k = 1:num_k
            for i_n = 1:num_nodes
                % ── 提取 bin centers 和 probability ──
                centers = squeeze(PD_storage(i_b, i_k, i_n, 1, :));   % [num_bin x 1]
                prob    = squeeze(PD_storage(i_b, i_k, i_n, 2, :));   % [num_bin x 1]，和为1

                % ── 计算 bin width δx ──
                x_min = edge_limit(i_b, i_n, 1);
                x_max = edge_limit(i_b, i_n, 2);
                delta_x = (x_max - x_min) / num_bin;

                % ── probability → empirical pdf ──
                p_empirical = prob / delta_x;

                % ── Gaussian pdf (理论分布) ──
                mu    = mu_storage(i_b, i_k, i_n);
                sigma = sigma_storage(i_b, i_k, i_n);
                p_gaussian = gauss_pdf(centers, mu, sigma);

                % ── 计算 ISE ──
                % ISE = (δx)^2 / M * sum_{j=1}^{M} [p'(x_j) - p(x_j)]^2
                ISE_storage(i_b, i_k, i_n) = (delta_x^2 / M) * sum((p_empirical - p_gaussian).^2);
            end
        end

        fprintf('[compute_ise] b=%d/%d 完成\n', i_b, num_b);
    end

    %% ── 保存结果 ─────────────────────────────────────

    % ── .mat 文件 ──
    mat_out = fullfile(input_dir, 'ISE_storage.mat');
    save(mat_out, 'ISE_storage', '-v7.3');
    fprintf('[compute_ise] 已保存 .mat: %s\n', mat_out);

    % ── .txt 文件 ──
    % 格式：i_b  i_k  i_n  ISE
    txt_out = fullfile(input_dir, 'ISE_storage.txt');
    fid = fopen(txt_out, 'w');
    fprintf(fid, '# ISE_storage\n');
    fprintf(fid, '# 格式: i_b\ti_k\ti_n\tISE\n');
    fprintf(fid, '# num_b=%d, num_k=%d, num_nodes=%d, num_bin=%d\n', ...
            num_b, num_k, num_nodes, num_bin);
    fprintf(fid, '# ISE = (delta_x)^2 / M * sum[ (p_empirical - p_gaussian)^2 ]\n');
    fprintf(fid, '#\n');

    for i_b = 1:num_b
        for i_k = 1:num_k
            for i_n = 1:num_nodes
                fprintf(fid, '%d\t%d\t%d\t%.12e\n', ...
                        i_b, i_k, i_n, ISE_storage(i_b, i_k, i_n));
            end
        end
    end
    fclose(fid);
    fprintf('[compute_ise] 已保存 .txt: %s\n', txt_out);

    % ── 出图用小数据（.mat + tab 分隔 .txt）──
    ISE_2d    = reshape(ISE_storage, num_b, num_k * num_nodes);
    ise_mean  = mean(ISE_2d, 2);          % [num_b x 1]
    ise_std   = std(ISE_2d, 0, 2);        % [num_b x 1]
    b_values  = cfg.b(:);                 % [num_b x 1]

    sd_dir = fullfile(cfg.path_small_data, 'ISE', ...
                      sprintf('setting_%d', i_setting), 'pd_storage');
    if ~exist(sd_dir, 'dir'), mkdir(sd_dir); end

    % .mat
    save(fullfile(sd_dir, 'ISE_summary.mat'), 'b_values', 'ise_mean', 'ise_std');

    % .txt（tab 分隔）
    write_col_vector(fullfile(sd_dir, 'b_values.txt'), b_values);
    write_col_vector(fullfile(sd_dir, 'ise_mean.txt'), ise_mean);
    write_col_vector(fullfile(sd_dir, 'ise_std.txt'),  ise_std);
    fprintf('[compute_ise] 已保存出图小数据: %s\n', sd_dir);

    % ── 统计摘要 ──
    fprintf('[compute_ise] === ISE 统计摘要 ===\n');
    fprintf('  最小值: %.6e\n', min(ISE_storage(:)));
    fprintf('  最大值: %.6e\n', max(ISE_storage(:)));
    fprintf('  均值  : %.6e\n', mean(ISE_storage(:)));
    fprintf('  中位数: %.6e\n', median(ISE_storage(:)));
    fprintf('  标准差: %.6e\n', std(ISE_storage(:)));

    write_done_flag(flag_file);
    fprintf('[compute_ise] 全部完成！\n');
end
