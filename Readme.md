有色噪声驱动下电力网络 swing equation 数值模拟与分析（完整重构版）

## 目录结构

```
powergrid-v2/
├── config.m                      ← 所有参数 + 13 个 setting 定义（只改这里）
├── main_pipeline.m               ← 主流程入口：生成全部非 VPM 数据与图片
│
├── core/                         ← 纯函数库
│   ├── colored_noise.m           ← 1/f^b 有色噪声（频域滤波）
│   ├── make_P.m                  ← 节点功率向量（3种模式）
│   ├── make_K.m                  ← 耦合强度矩阵（3种模式）
│   ├── rk4_solver.m              ← 定步长 RK4 积分器
│   └── swing_solver.m            ← swing equation 求解
│
├── analysis/                     ← 分析脚本（被 main_pipeline 调用）
│   ├── compute_pd.m              ← Stage 2：直方图 + mu/sigma → pd_storage
│   ├── compute_HD.m              ← Stage 3a：Hellinger distance + 小数据导出
│   ├── compute_ise.m             ← Stage 3b：ISE + 小数据导出
│   ├── compute_HD_vs_M.m         ← Stage 4a：figure_3(b) H vs M
│   ├── compute_ISE_vs_time.m     ← Stage 4b：figure_4(b) ISE vs Time
│   ├── hellinger.m               ← 共享：Hellinger distance
│   ├── gauss_pdf.m               ← 共享：Gaussian pdf
│   ├── write_col_vector.m        ← 共享：txt 导出（tab 分隔）
│   ├── write_row_vector.m        ← 共享：txt 导出（tab 分隔）
│   ├── check_done_flag.m         ← 共享：done.flag 检查
│   ├── write_done_flag.m         ← 共享：done.flag 写入
│
├── plotting/                     ← 绘图函数（被 main_pipeline Stage 5 调用）
│   ├── draw_ise.m                ← ISE + HD 双 Y 轴图（每个 setting 一张）
│   ├── draw_HD_vs_M.m            ← figure_3 复现（HD vs b + HD vs M）
│   ├── draw_ISE_vs_time.m        ← figure_4 复现（ISE vs b + ISE vs Time）
│   ├── draw_HD_vs_b.m            ← 单轴 HD(b)（主 setting）
│   ├── draw_ISE_vs_b.m           ← 单轴 ISE(b)（主 setting）
│   ├── draw_HD_vs_d.m            ← HD vs 节点到扰动点的图距离（主 setting）
│   ├── draw_ISE_vs_d.m           ← ISE vs 图距离，log 纵轴（主 setting）
│   ├── plot_with_shaded_error.m  ← 共享：均值 + 误差阴影
│   ├── color_gradient.m          ← 共享：两色线性渐变
│   └── formatFigure.m            ← 共享：字体/边框/尺寸统一
│
├── figures/                      ← 只存放生成的图片（无任何代码）
│   ├── setting_<i>/*.png         ← main_pipeline Stage 5 输出
│   └── setting_<i>_thermal.png   ← main_vpm 输出（在 figures/ 根目录）
│
├── currently-using-key-data/     ← 图数据（从原项目复制）
│   ├── pos_noise.mat
│   ├── G20.mat
│   ├── G_random_50nodes.mat
│   ├── G_random_100nodes.mat
│   └── ...
│
├── data/                         ← 自动生成（原始大数据）
│   └── setting_<i>/
│       ├── y_der_<ib>_<ik>.mat
│       ├── noise_<ib>_<ik>.mat
│       ├── parameters.mat        ← 本 setting 完整参数
│       ├── simulation_done.flag  ← Stage 1 完成标记
│       └── pd_storage/
│           ├── pd_storage.mat    ← 直方图缓存 + mu/sigma
│           ├── ISE_storage.mat   ← ISE 原始数组
│           ├── ISE_storage.txt
│           └── pd_storage_done.flag
│
└── small_data_for_plotting/      ← 出图用小文件（.mat + tab 分隔 .txt）
    ├── setting_<i>/
    │   ├── HD.mat                ← HD [num_b x num_nodes]
    │   ├── HD_vs_M.mat           ← figure_3(b) 数据
    │   ├── ISE_vs_time.mat       ← figure_4(b) 数据
    │   ├── M_list.txt, HD_vs_M_b<val>.txt
    │   ├── T_list.txt, ISE_vs_time_b<val>.txt
    │   └── *_done.flag
    │
    ├── data_chi2_R2_H/setting_<i>/       ← HD 曲线小数据
    │   ├── HD_summary.mat
    │   ├── b.txt
    │   ├── HD_mean.txt
    │   └── HD_std.txt
    │
    └── ISE/setting_<i>/pd_storage/       ← ISE 曲线小数据
        ├── ISE_summary.mat
        ├── b_values.txt
        ├── ise_mean.txt
        └── ise_std.txt
```

## 使用方法

1. 确保`currently-using-key-data/` 下的 `.mat` 变量保存完整（包括graph对象，扰动节点序号）。
2. 按需修改 `config.m`（参数、setting 定义）。
3. MATLAB 中运行 `main_pipeline.m`：
   - 生成全部`data/` 和 `small_data_for_plotting/`。
   - Stage 5 生成 `figures/setting_<i>/` 下的所有图片（见下表）。
4. 如需偶尔生成新图或可视化网络：使用 `create_random_power_grid.m` / `draw_network.m`
   （后者需先手动把 `G` 载入工作区）。

重新运行时，已完成的 setting 自动跳过（靠 `*_done.flag`）。

### Stage 5 生成的图片

| 文件 | 函数 | setting 范围 |
|---|---|---|
| `setting_<i>_ise_H.png` | `draw_ise` | 全部 |
| `figure_3_HD_stationarity.png` | `draw_HD_vs_M` | 主 setting（默认 10） |
| `figure_4_ISE_timespan.png` | `draw_ISE_vs_time` | 主 setting |
| `setting_<i>_HD_b.png` | `draw_HD_vs_b` | 主 setting |
| `setting_<i>_ISE_b.png` | `draw_ISE_vs_b` | 主 setting |
| `setting_<i>_HD_d.png` | `draw_HD_vs_d` | 主 setting |
| `setting_<i>_ISE_d.png` | `draw_ISE_vs_d` | 主 setting |

## 数据流

```
[Stage 1] y_der / noise  （原始大数据，data/）
              ↓
[Stage 2] pd_storage.mat  （直方图 + mu/sigma，data/）
              ↓
[Stage 3] HD.mat（→ small_data）/ ISE_storage.mat（→ data/pd_storage）+ 出图小数据
              ↓
[Stage 4] HD_vs_M.mat / ISE_vs_time.mat（主 setting，默认 10）
              ↓
[Stage 5] figures/setting_<i>/*.png
```

## `small_data_for_plotting` 文件约定

所有出图用小文件均同时提供 `.mat` 与 tab 分隔 `.txt`：

| 用途 | 文件模式 | 示例 |
|---|---|---|
| 单条曲线 | `x.txt`, `y.txt` | `b.txt`, `HD_mean.txt` |
| 多条曲线 | `x.txt`, `y1.txt`, `y2.txt`, ... | `M_list.txt`, `HD_vs_M_b0.00.txt` |
| 带误差带曲线 | `x.txt`, `y_mean.txt`, `y_err.txt` | `b_values.txt`, `ise_mean.txt`, `ise_std.txt` |

## 13 个 Settings

| # | alpha | 网络 | P 分布 | K 分布 |
|---|---|---|---|---|
| 1 | 0.5 | G20 | 各50% | 固定30 |
| 2 | 1.5 | G20 | 各50% | 固定30 |
| 3 | 1 | G_random_50nodes | 各50% | 固定30 |
| 4 | 1 | G_random_100nodes | 各50% | 固定30 |
| 5 | 1 | G20 | 20%+ / 80%- | 固定30 |
| 6 | 1 | G20 | 40%+ / 60%- | 固定30 |
| 7 | 1 | G20 | 各50% | 均匀[25,35] |
| 8 | 1 | G20 | 各50% | 均匀[20,40] |
| 9 | 1 | G_random_500nodes | 各50% | 固定30 |
| 10 | 1 | G20 | 各50% | 固定30 |
| 11 | 1 | G21 | 各50% | 固定30 |
| 12 | 1 | G_BA | 各50% | 固定30 |
| 13 | 1 | G_chain_10_nodes | 各50% | 固定30 |

## 手动工具（不在主流程中）

| 文件 | 用途 | 状态 |
|---|---|---|
| `create_random_power_grid.m` | 生成随机电网图 | 可用；顶部 `clear`、参数硬编码在脚本内 |
| `draw_network.m` | 可视化网络 | 需先手动 `load` 图得到工作区变量 `G` |


