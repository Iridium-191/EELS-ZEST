# EELS-ZEST 中文说明

**Zero-loss Estimation and Sparse Thermometry**

EELS-ZEST 是面向 STEM-EELS 单谱和 datacube 的 MATLAB 工具箱，主要用于
零损峰估计、稀疏卷积信号重建、稀疏声子模态源强度分析和 attention
引导的 Bose-Einstein 温度拟合。

## v2.1 方法概览

1. **Stage 1 负责谱形重建。**  
   程序拟合完整 ZLP，并在能量轴对齐后使用峰心位于 `0 meV` 的总 fitted-ZLP
   作为卷积核，重建 Stokes 和 anti-Stokes 稀疏声子模态源强度。
2. **完整 datacube 的 Stage 1 使用分层 ZLP 模型。**  
   程序先以对齐、rank-limited SVD 后的 ZLP 谱形聚类，为每一类建立共享模板；
   soft membership 连续混合模板目标和允许的参数范围，避免不同类别边界发生
   硬切换。8 邻域 Huber 空间项只约束峰心、宽度、
   Pearson 形状和非对称性，不约束 ZLP 振幅或稀疏声子模态源强度。
3. **Stage 2 只用可信配对信号估温。**  
   attention 同时识别到的正负能量对参与 hard-BE 约束；其余能量点仍参与
   全谱拟合，但使用独立的软源强度系数，不主导温度。
4. **ZLP 差别识别 sharp interface，并控制单向温度补点。**  
   稳健的 ZLP-only edge score 标记 sharp interface，跨界边不再施加 Stage 1
   空间损失。Data-valid 温度完全固定，不做温度空间平滑；Data-invalid 像素
   只有在同侧至少存在两个有效邻居时才允许输出带标记的 ZLP 引导推断值，
   推断点不会反向影响有效锚点。
5. **诊断结果与温度一起输出。**  
   包括残差、hard/soft 声子模态源强度、BE 线性、temperature cost
   valley、不确定度、reliability 和 datacube 空间图。

默认配对分数为：

```text
pair_score(EM) = sqrt(att(+EM) * att(-EM))
```

`pair_score >= 0.35` 的能量对进入主要 hard-BE 集合；当配对数不足时，可从
`pair_score >= 0.10` 的候选中补充。温度诊断只使用 hard attention pairs。

## 运行环境

- MATLAB
- Optimization Toolbox
- 多 worker datacube 拟合需要 Parallel Computing Toolbox

## 快速开始

启动 GUI：

```matlab
addpath('EELS-ZEST');
app = run_eels_zest;
```

拟合单条谱线：

```matlab
settings = struct();
settings.run_stage2 = true;
settings.crop_range = [-250, 250];
settings.stage2_temperature_algorithm = 'grid_scan';
settings.stage1_hierarchical_enabled = true;
settings.stage2_spatial_map_enabled = true;

result = eels_fit_spectrum(energy_meV, spectrum, settings);
```

拟合随包提供的示例 datacube：

```matlab
settings = struct();
settings.run_stage2 = true;
settings.crop_range = [-250, 250];
settings.datacube_sample_strategy = 'all';
settings.datacube_parallel_enabled = true;
settings.datacube_parallel_workers = 16;
settings.test_output_dir = fullfile(pwd, 'eels_zest_output');

cuberes = eels_fit_datacube( ...
    fullfile('examples', 'data', 'example_data.mat'), settings);
```

也可以直接运行：

```matlab
run(fullfile('examples', 'zest_datacube_full_stage1.m'));
run(fullfile('examples', 'zest_datacube_full_stage2_parallel.m'));
```

## 推荐的 Datacube 流程

1. 载入 `.mat` 或 `.dm3` datacube。
2. 必要时先将 ZLP 对齐到 `0 meV`。
3. 进行低秩 SVD/PCA 降噪。
4. 裁剪到有效能区，通常为 `[-250, 250] meV`。
5. 对完整 datacube 运行分层 Stage 1，检查 cluster、membership、outlier、
   template deviation、残差、R2 和 NRMSE。
6. 检查 hard attention pairs 的 BE 线性和原始 temperature cost valley。
7. 以 Data-valid 为主要温度结果；对比 Raw 和可选的 ZLP-guided temperature，
   并检查 sharp-interface 与 Spatial-inferred mask 后再解释空间温度变化。

v2.1 GUI 会分别显示 template polish、每个空间颜色批次和 final
Stage 2 进度，并在有界并行批次之间检查取消请求。全量数值 MAT 和
summary 会先于任何可选绘图保存；大型诊断图默认不再自动导出。

## v2.1 Datacube 输出

- `T_map_raw`：逐像素原始温度，永不被空间结果覆盖；
- `T_map_data_valid`：排除边界 valley、单边 valley 和证据不足像素；
- `T_map_joint`：Data-valid 固定锚点加可选 ZLP 引导单向补点，不平滑有效温度；
- `T_uncertainty_joint_map`、`T_spatial_inferred_mask`、`T_joint_shift_map`；
- `cluster_label_map`、`cluster_membership_map`、`cluster_outlier_mask`；
- `template_deviation_map`、`zlp_spatial_loss_map`；
- `zlp_interface_score_map`、`zlp_interface_mask` 和方向 edge maps。

## 主要入口

| 文件 | 用途 |
| --- | --- |
| `run_eels_zest.m` | 启动 GUI |
| `eels_fit_spectrum.m` | 拟合单条谱线 |
| `eels_fit_datacube.m` | 拟合完整 datacube |
| `EELS_PearsonSparse_FullFit.m` | Stage 1/Stage 2 核心求解器与诊断 |
| `EELS_ZEST_App.m` | GUI 实现 |
| `EELS_ZEST_DataOps.m` | 稀疏声子源截断、重卷积、STA 与 EELS 导出 |
| `EELS_ZEST_PeakFit.m` | Gaussian、Lorentzian 和 Pseudo-Voigt 稳健多峰拟合 |
| `examples/` | 可直接运行的脚本 |
| `docs/example_data_figures/` | 示例图、CSV、MAT 数据与复现脚本 |
| `index.html` | 项目网页 |

## 输入格式

支持：

- 含 `.data` 和 `.ene` 字段的 MATLAB 对象或结构体；
- 含可识别 datacube 和能量轴变量的 `.mat` 文件；
- 通过 `DM3Import` 载入的 `.dm3` 文件。

推荐维度：

```text
Ny x Nx x Ne
```

其中 `Ne` 为能量通道数。

## v2.1 Data 页工具

- **Send signal to Data** 可按手动能量、Materials Project 声子 JSON
  或 Stage 2 最右 hard-attention 区间加拖尾截断 Stage 1/Stage 2 稀疏声子
  源强度，再使用同一 fitted-ZLP kernel 重卷积。最终谱保留自然卷积拖尾，
  不在 cutoff 处平切。导回的是扣除加性 ZLP 与背景后的声子重建谱；每像素
  fitted-ZLP kernel、Stage 1/Stage 2 温度图和有效掩膜会随 Data 状态写入
  `EELS.info`。因此 RL 解卷当前截断谱，STA 也能继续使用对应温度。
- **Multiply STA** 只保留正能量 loss 侧，并按
  `I(E)*E/2*[1-exp(-E/(kB*T))]` 乘 STA；温度可手动指定、使用 Stage 1
  soft 温度，或使用 Stage 2 温度并逐像素回退 Stage 1。
- **Stack Spectra** 使用非阻塞的 map 框选窗口；**Multi-peak fit** 支持
  ROI、x-line 和完整 datacube 的 Gaussian、Lorentzian、Pseudo-Voigt 稳健拟合。
- **Export data** 将 Data 页当前处理状态保存为 `EELS` 对象。其余导出按钮
  保存 MAT/CSV，图形内容还可保存 FIG/SVG/PNG/TIFF；写入 base workspace
  仅作为可选项。

## 示例图复现

示例数据位于：

```text
examples/data/example_data.mat
```

复现图和面板 CSV：

```matlab
cd(fullfile('docs', 'example_data_figures'));
export_example_figure_data;
plot_example_figures;
```

## 作者与引用

EELS-ZEST 由北京大学国际量子材料中心（ICQM）的 **Ruilin Mao** 开发。

如果 EELS-ZEST 对研究工作有所贡献，请按照 `CITATION.cff` 引用本软件。
GitHub 仓库会根据该文件提供 **Cite this repository** 引用入口。

## 许可证

Copyright (C) 2026 Ruilin Mao。

EELS-ZEST 采用 **GNU General Public License version 3 only**
（`GPL-3.0-only`）发布，完整条款见 `LICENSE`。欢迎研究使用者向作者通报，
并在相关论文、报告或数据成果中引用本软件。

## 文档

- `docs/GUI_TOOLBOX_GUIDE.md`：GUI 操作流程
- `docs/LOSS_AND_ATTENTION_THEORY_v1.3.md`：损失函数与 attention 方法
- `docs/theory.html`：公式和方法细节
- `index.html`：可浏览的项目说明

## 术语与 MATLAB 字段名

本文将拟合得到的非负系数称为**稀疏声子模态源强度**。它表示与 fitted-ZLP
卷积之前的源空间权重，不应直接表述为经过绝对标定的 phonon DOS。

为兼容现有 MATLAB 结果字段，代码中将 **phonon-mode intensity** 缩写为
**PM**。主要字段为 `PM_stokes`、`PM_antistokes` 和 `EM_final`。
