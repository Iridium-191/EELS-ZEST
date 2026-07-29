# EELS-ZEST 中文说明

**Zero-loss Estimation and Sparse Thermometry**

EELS-ZEST 是面向 STEM-EELS 单谱和 datacube 的 MATLAB 工具箱，主要用于
零损峰估计、稀疏卷积信号重建、稀疏声子模态源强度分析和 attention
引导的 Bose-Einstein 温度拟合。

## v2.0 方法概览

1. **Stage 1 负责谱形重建。**  
   程序拟合完整 ZLP，并在能量轴对齐后使用峰心位于 `0 meV` 的总 fitted-ZLP
   作为卷积核，重建 Stokes 和 anti-Stokes 稀疏声子模态源强度。
2. **Stage 2 只用可信配对信号估温。**  
   attention 同时识别到的正负能量对参与 hard-BE 约束；其余能量点仍参与
   全谱拟合，但使用独立的软源强度系数，不主导温度。
3. **诊断结果与温度一起输出。**  
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
5. 运行 Stage 1，检查 ZLP、重建信号、残差、R2 和 NRMSE。
6. 在 gain/loss 配对信息足够时运行 Stage 2。
7. 结合 BE 线性、temperature cost valley、reliability 和不确定度解释温度。

## 主要入口

| 文件 | 用途 |
| --- | --- |
| `run_eels_zest.m` | 启动 GUI |
| `eels_fit_spectrum.m` | 拟合单条谱线 |
| `eels_fit_datacube.m` | 拟合完整 datacube |
| `EELS_PearsonSparse_FullFit.m` | Stage 1/Stage 2 核心求解器与诊断 |
| `EELS_ZEST_App.m` | GUI 实现 |
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
