# EELS-ZEST 损失函数与注意力机制原理详解

本文档说明 EELS-ZEST v1.3 主线实现中 **ZLP + sparse convolution** 拟合、Stage 1 稀疏反卷、Stage 2 Bose--Einstein 温度约束、注意力窗口、误差评估与诊断图的数学结构。文档中的符号尽量与 `EELS_PearsonSparse_FullFit.m` 的变量命名一一对应，便于写文章、答辩或后续重构时引用。

> 版本范围：本文以 `pearson_sparse_fullfit/EELS_PearsonSparse_FullFit.m` 的主线 v1.3 实现为准。`version1.3_fullhardBE_beta` 只在 Stage 2 中改变 hard-BE EM 的选择策略，主体模型与损失函数结构相同。

---

## 1. 符号表

| 符号 | 代码变量 | 含义 |
|---|---|---|
| \(E_i\) | `ene`, `ene_corrected` | 第 \(i\) 个能量损失采样点，单位 meV |
| \(y_i\) | `raw_data_norm` | 归一化后的实验谱强度 |
| \(I_i\) | `raw_data` | 原始强度，未归一化 |
| \(\hat y_i\) | `fitted_spectrum_norm`, `stage2.model_norm` | 归一化拟合谱 |
| \(\hat I_i\) | `fitted_spectrum`, `stage2.model_norm * norm_val` | 原强度尺度下的拟合谱 |
| \(r_i\) | residual | 残差，\(r_i=y_i-\hat y_i\) 或 \(I_i-\hat I_i\) |
| \(w_i\) | `sqrtW.^2`, `weights_used` | 能量点 \(E_i\) 的数据拟合权重 |
| \(\sqrt{w_i}\) | `sqrtW` | 加权最小二乘中实际乘到残差上的权重 |
| \(Z(E)\) | `zlp_only`, `zlp_norm` | ZLP envelope |
| \(Z_k(E)\) | `zlp_components` | 第 \(k\) 个 ZLP 分量 |
| \(K(E)\) | `signal_kernel` | 用作稀疏卷积核的 ZLP kernel |
| \(\varepsilon_m\) | `EM_final`, `EM_hard`, `EM_soft` | 第 \(m\) 个 vibrational source / EM 中心，单位 meV，取正值 |
| \(p_m^{+}\) | loss/Stokes source intensity | loss/Stokes 侧 \(\varepsilon_m\) 的非负声子模态源强度系数 |
| \(p_m^{-}\) | gain/anti-Stokes source intensity | gain/anti-Stokes 侧 \(\varepsilon_m\) 的非负声子模态源强度系数 |
| \(T\) | `T_fit`, `T_soft` | 温度，单位 K |
| \(k_\mathrm{B}\) | hard-coded `0.08617333262` | Boltzmann 常数，单位 meV/K |
| \(q_m(T)\) | `q`, `q_hard`, `q_soft` | BE 详细平衡因子，\(q_m(T)=\exp[-\varepsilon_m/(k_\mathrm{B}T)]\) |
| \(\theta_Z\) | `zlp_params`, nonlinear parameters | ZLP 形状、宽度、偏移等非线性参数 |
| \(\mathbf{A}\) | `A` | 线性设计矩阵 |
| \(\mathbf{c}\) | `coef` | 线性非负系数向量 |
| \(\mathbf{R}\mathbf{c}\approx\mathbf{r}_0\) | `R`, `r0`, `Reg`, `reg_rhs` | 正则项写成的伪观测行 |
| \(\mathcal{A}(E)\) | `attention.feature_score` | Stage 2 feature attention score，范围 \([0,1]\) |
| \(W_\mathcal{A}(E)\) | `attention.weight` | 由注意力生成的 Stage 2 额外拟合权重 |
| \(\Omega_\mathcal{A}\) | `attention.active_mask`, `attention.bands_meV` | attention-selected feature window |
| \(\Omega_\mathrm{fit}\) | `fit_mask` | Stage 2 实际拟合能量窗口 |
| \(\Omega_\mathrm{tr}\) | `transition_mask` | ZLP core 与 temperature window 之间的 transition/soft 区域 |
| \(\Omega_\mathrm{tail}\) | tail mask | 高频拖尾约束区域 |
| \(\mathrm{CI}(E)\) | `confidence_interval` | 能量依赖的拟合置信区间 |

符号约定：

- 上标 \(+\) 表示 loss/Stokes 侧，能量位置 \(+\varepsilon_m\)。
- 上标 \(-\) 表示 gain/anti-Stokes 侧，能量位置 \(-\varepsilon_m\)。
- 所有源强度系数均被约束为非负：

\[
p_m^{+}\ge 0,\qquad p_m^{-}\ge 0.
\]

---

## 2. 核心物理模型

### 2.1 总谱模型

EELS-ZEST 的当前主线不再使用独立背景项。新结果中 background 列为空或为零，仅保留兼容字段。因此主要模型为

\[
\hat y(E)
= Z(E;\theta_Z) + V(E;\theta_Z,\mathbf{p}^{+},\mathbf{p}^{-},\boldsymbol{\varepsilon}),
\]

其中

\[
Z(E;\theta_Z)=\sum_{k=1}^{N_Z} a_k Z_k(E;\theta_k),
\]

\[
V(E)=
\sum_{m=1}^{M} p_m^{+}K(E-\varepsilon_m)
+\sum_{m=1}^{M} p_m^{-}K(E+\varepsilon_m).
\]

这里 \(Z(E)\) 是直接拟合的 zero-loss peak；\(V(E)\) 是用 ZLP kernel 卷积得到的振动信号。这个表达式的物理含义是：实验谱中每个振动损失/获得事件并不是一个 delta peak，而是被仪器零峰展宽卷积后的线形。因此 \(K(E)\) 应来自当前谱线的 ZLP 模型，而不是另造一个固定高斯核。

### 2.2 ZLP 分量

代码支持 asymmetric Pearson 与 asymmetric pseudo-Voigt。主线默认使用 asymmetric Pearson：

\[
Z_k(x)=
\left[
1+\frac{(x/s)^2}{m}
\right]^{-m},
\]

其中

\[
(s,m)=
\begin{cases}
(s_{L,k},m_{L,k}), & x<0,\\
(s_{R,k},m_{R,k}), & x\ge 0.
\end{cases}
\]

如果切换到 asymmetric pseudo-Voigt，则

\[
Z_k(x)=
\eta\,\frac{1}{1+(x/s)^2}
+(1-\eta)\exp\!\left[-\frac{1}{2}\left(\frac{x}{s}\right)^2\right],
\]

\[
(s,\eta)=
\begin{cases}
(s_{L,k},\eta_{L,k}), & x<0,\\
(s_{R,k},\eta_{R,k}), & x\ge 0.
\end{cases}
\]

ZLP 分量具有宽度上限、低频 guard 与 signal guard。物理目的不是让 ZLP 变得最灵活，而是让它足够描述零峰拖尾，同时避免超宽 ZLP 分量吞掉低频 acoustic/optical signal。

### 2.3 卷积核

代码中的 kernel 有两种策略：

\[
K(E)=
\begin{cases}
\sum_k Z_k(E), & \texttt{kernel\_mode = full\_zlp},\\
Z_{k^\star}(E), & \texttt{kernel\_mode = core\_component},
\end{cases}
\]

其中 \(k^\star\) 是平均宽度最小的 ZLP 分量：

\[
k^\star=\arg\min_k\frac{s_{L,k}+s_{R,k}}{2}.
\]

之后归一化：

\[
K(E)\leftarrow \frac{K(E)}{\sum_i K(E_i)}.
\]

v2.0 默认使用 `fitted_zlp`，使卷积核直接来自拟合后的总 ZLP/仪器响应。`core_component` 仍可作为调试选项，但不再作为默认值，因为单个窄分量的中心可能偏离总 ZLP 峰心，导致 Stokes/anti-Stokes 卷积峰产生系统性位移。

---

## 3. 线性设计矩阵

在固定非线性 ZLP 参数 \(\theta_Z\)、EM 网格 \(\boldsymbol{\varepsilon}\) 后，模型对线性系数是非负线性的。设计矩阵写为

\[
\mathbf{A}
=
\left[
\mathbf{Z}_1,\ldots,\mathbf{Z}_{N_Z},
\mathbf{S}_1,\ldots,\mathbf{S}_{M},
\mathbf{G}_1,\ldots,\mathbf{G}_{M}
\right],
\]

其中

\[
\mathbf{Z}_k(i)=Z_k(E_i),
\]

\[
\mathbf{S}_m(i)=K(E_i-\varepsilon_m),
\qquad
\mathbf{G}_m(i)=K(E_i+\varepsilon_m).
\]

对应线性系数为

\[
\mathbf{c}
=
\left[
a_1,\ldots,a_{N_Z},
p_1^{+},\ldots,p_M^{+},
p_1^{-},\ldots,p_M^{-}
\right]^\mathsf{T}.
\]

Stage 1 的拟合谱可以写为

\[
\hat{\mathbf{y}} = \mathbf{A}\mathbf{c}.
\]

非负约束为

\[
\mathbf{c}\ge 0.
\]

这个非负约束非常关键：它防止算法用正负振荡的源强度去“抵消”ZLP 或伪造谱线，从而使结果更接近非负 phonon-mode source intensity。

---

## 4. Stage 1 损失函数

Stage 1 的目标是 **shape-first sparse deconvolution**：先尽可能稳定地复现谱线并分离 ZLP 与 sparse vibrational response，不强迫所有源强度系数严格满足温度。它的基本优化问题为

\[
\min_{\theta_Z,\mathbf{c}\ge0,T_\mathrm{soft}}
\mathcal{L}_1
=
\mathcal{L}_\mathrm{data}
+\mathcal{L}_\mathrm{softBE}
+\mathcal{L}_\mathrm{smooth}
+\mathcal{L}_{2}
+\mathcal{L}_{1}
+\mathcal{L}_\mathrm{ZLPwide}
+\mathcal{L}_\mathrm{ZLPguard}
+\mathcal{L}_\mathrm{tail}.
\]

实际实现中，非线性变量 \(\theta_Z,T_\mathrm{soft}\) 外层迭代，固定它们后内层求一个非负线性问题。

### 4.1 数据残差项

令

\[
r_i=\hat y_i-y_i.
\]

加权数据残差为

\[
\mathcal{L}_\mathrm{data}
=
\sum_i w_i r_i^2.
\]

等价地，代码中用

\[
\mathbf{C}=
\begin{bmatrix}
\sqrt{\mathbf{W}}\mathbf{A}\\
\mathbf{R}
\end{bmatrix},
\qquad
\mathbf{d}=
\begin{bmatrix}
\sqrt{\mathbf{W}}\mathbf{y}\\
\mathbf{r}_0
\end{bmatrix},
\]

求解

\[
\min_{\mathbf{c}\ge0}
\|\mathbf{C}\mathbf{c}-\mathbf{d}\|_2^2+\boldsymbol{\lambda}_{1}^{\mathsf{T}}\mathbf{c},
\]

其中 \(\mathbf{R}\mathbf{c}\approx\mathbf{r}_0\) 表示所有正则项写成的伪观测行。

### 4.2 基础能量权重

主线权重不是简单均匀权重。能量权重由 ZLP core、shoulder、signal、高频和 tail 多个区域叠加：

\[
w_i =
w_\mathrm{core} C_i+
w_\mathrm{shoulder} S_i+
w_\mathrm{sig}(1-C_i)+
w_\mathrm{hf} H_i+
w_\mathrm{tail} T_i.
\]

其中用 sigmoid 平滑定义区域：

\[
C_i=\sigma\!\left(\frac{E_c-|E_i|}{\tau_c}\right),
\]

\[
S_i=
\sigma\!\left(\frac{|E_i|-E_c}{\tau_c}\right)
\sigma\!\left(\frac{E_h-|E_i|}{\tau_h}\right),
\]

\[
H_i=\sigma\!\left(\frac{|E_i|-E_h}{\tau_h}\right),
\]

\[
T_i=\sigma\!\left(\frac{|E_i|-E_t}{\tau_t}\right),
\]

\[
\sigma(x)=\frac{1}{1+\exp(-x)}.
\]

代码中大致对应：

| 数学符号 | 参数 |
|---|---|
| \(E_c\) | `core_downweight_width` |
| \(E_h\) | `highfreq_start_meV` |
| \(E_t\) | `tail_start_fraction * max(abs(E))` 与 `highfreq_start_meV` 的较大者 |
| \(w_\mathrm{core}\) | `core_weight` |
| \(w_\mathrm{shoulder}\) | `shoulder_weight` |
| \(w_\mathrm{sig}\) | `signal_weight` |
| \(w_\mathrm{hf}\) | `highfreq_weight` |
| \(w_\mathrm{tail}\) | `tail_weight` |

此外，为了避免强 ZLP core 完全支配目标函数，代码还会对局部低强度区域做动态补偿：

\[
B_i=
\left[
\frac{I_\mathrm{ref}}
\max\{I_\mathrm{local}(E_i),I_\mathrm{floor}\}}
\right]^{0.45},
\]

\[
B_i \leftarrow \mathrm{clip}(B_i,1,B_\mathrm{max}),
\]

\[
w_i\leftarrow
w_i\left[
1+\gamma(1-C_i)(B_i-1)
\right].
\]

物理意义：低频和高频振动峰远弱于 ZLP。如果直接用强度残差，ZLP core 的微小相对误差会淹没高频峰。因此权重函数主动降低 core 的相对统治力，提高 signal/tail 区域的小残差惩罚。

### 4.3 Stage 1 soft BE 正则

Stage 1 不硬性令 Stokes/anti-Stokes 满足温度，而是通过软惩罚项给出物理先验：

\[
q_m(T_\mathrm{soft})
=
\exp\left[-\frac{\varepsilon_m}{k_\mathrm{B}T_\mathrm{soft}}\right].
\]

Soft BE 残差为

\[
\Delta_m^\mathrm{BE}
=
p_m^{-}-q_m(T_\mathrm{soft})p_m^{+}.
\]

损失项为

\[
\mathcal{L}_\mathrm{softBE}
=
\lambda_\mathrm{BE}
\sum_m
\left(
p_m^{-}-q_m(T_\mathrm{soft})p_m^{+}
\right)^2.
\]

物理意义：如果 anti-Stokes 侧确实可靠，它会自然推动源强度接近 BE 比例；但如果 gain 侧弱、ZLP tail 未完全分离或存在噪声，Stage 1 仍允许源强度偏离 BE，以优先复现谱线形状。

### 4.4 源强度平滑正则

对 Stokes 和 anti-Stokes 的源强度分布分别使用二阶差分：

\[
\mathcal{L}_\mathrm{smooth}
=
\lambda_\mathrm{smooth}
\left(
\|\mathbf{D}_2\mathbf{p}^{+}\|_2^2
+
\|\mathbf{D}_2\mathbf{p}^{-}\|_2^2
\right).
\]

这里 \(\mathbf{D}_2\) 是沿 EM 轴的二阶差分矩阵。它抑制一个一个孤立尖峰，让稀疏声子模态源强度更接近连续 phonon source / DOS-like 分布。

### 4.5 源强度 L2 与 L1 稀疏正则

L2 正则：

\[
\mathcal{L}_2=
\lambda_2
\left(
\|\mathbf{p}^{+}\|_2^2+\|\mathbf{p}^{-}\|_2^2
\right).
\]

L1 稀疏项在内层线性问题中以线性代价加入：

\[
\mathcal{L}_1=
\boldsymbol{\lambda}_1^\mathsf{T}
\begin{bmatrix}
\mathbf{p}^{+}\\
\mathbf{p}^{-}
\end{bmatrix}.
\]

物理意义：EELS 振动谱通常可看作有限峰或有限谱带的卷积响应。L1 倾向减少不必要的源强度系数，L2 防止单个点过度爆炸；二者和源强度平滑一起控制“谱线毛刺”。

### 4.6 ZLP 宽分量惩罚

为了避免超宽 ZLP 分量吞掉真实低频信号，对较宽/靠后的 ZLP 分量幅度加 L2：

\[
\mathcal{L}_\mathrm{ZLPwide}
=
\lambda_\mathrm{wide}
\sum_{k\in\mathcal{K}_\mathrm{wide}}a_k^2.
\]

代码中由 `zlp_wide_component_l2` 和 `zlp_wide_component_start_index` 控制。

### 4.7 ZLP signal guard

ZLP signal guard 的作用是限制 ZLP 在振动信号区域继续抬高。设 guard 区域为

\[
\Omega_\mathrm{zg}=
\left\{
i:
E_\mathrm{zg}^{\min}\le |E_i|\le E_\mathrm{zg}^{\max}
\right\}.
\]

定义边界平滑 taper：

\[
g_\mathrm{lo}(E)=
\sigma\!\left(\frac{|E|-E_\mathrm{zg}^{\min}}{\tau_\mathrm{zg}}\right),
\]

\[
g_\mathrm{hi}(E)=
\sigma\!\left(\frac{E_\mathrm{zg}^{\max}-|E|}{\tau_\mathrm{zg}}\right).
\]

高能 tail release：

\[
g_\mathrm{rel}(E)
=
f_\mathrm{rel}
+
(1-f_\mathrm{rel})
\sigma\!\left(\frac{E_\mathrm{rel}-|E|}{\tau_\mathrm{rel}}\right).
\]

再乘上由谱线局部特征决定的 profile \(g_\mathrm{feat}(E)\)。最终

\[
g_\mathrm{ZLP}(E)
=
g_\mathrm{lo}(E)g_\mathrm{hi}(E)g_\mathrm{rel}(E)g_\mathrm{feat}(E).
\]

损失项为

\[
\mathcal{L}_\mathrm{ZLPguard}
=
\lambda_\mathrm{zg}
\sum_{i\in\Omega_\mathrm{zg}}
\left[
g_\mathrm{ZLP}(E_i) Z(E_i)
\right]^2.
\]

物理意义：如果某个能量区确实有 vibrational signal，ZLP 不应靠一个宽尾巴把它解释掉。guard 不是硬截断，而是软惩罚；它允许 ZLP tail 存在，但要求它在信号区保持克制。

### 4.8 低频 ZLP guard

低频 guard 是 signal guard 的更窄版本，用于 acoustic/low-frequency optical 区域：

\[
\Omega_\mathrm{lf}=
\{i:E_\mathrm{lf}^{\min}\le |E_i|\le E_\mathrm{lf}^{\max}\}.
\]

\[
\mathcal{L}_\mathrm{lowfreq}
=
\lambda_\mathrm{lf}
\sum_{i\in\Omega_\mathrm{lf}}
\left[
g_\mathrm{lf}(E_i)Z(E_i)
\right]^2.
\]

物理意义：低频区域最容易被 ZLP 宽分量吃掉，因此用较弱但专门的 guard 抑制 ZLP 低频吞峰。

### 4.9 高频拖尾 overshoot 约束

tail overshoot 不是普通正则，而是可作为不等式约束加入线性求解：

\[
\hat y_i \le y_i+\delta_\mathrm{tail},
\qquad i\in\Omega_\mathrm{tail}.
\]

其中

\[
\delta_\mathrm{tail}
=
\max
\left(
\delta_\mathrm{abs},
\eta_\mathrm{noise}\sigma_\mathrm{noise},
\eta_\mathrm{rel}A_\mathrm{tail}
\right).
\]

物理意义：高频尾部信号弱，模型稍微高于数据就会形成“抬尾”。这个约束防止 ZLP 或源强度在远端拖尾处系统性高估。

---

## 5. Stage 2 温度约束损失

Stage 2 在 Stage 1 的 ZLP 和源强度初值基础上估计温度。核心思想是：只在有足够 Stokes/anti-Stokes 证据的 EM 上硬施加 BE；其余区域保留 soft/free source intensity 以维护谱型。

### 5.1 hard / soft EM 分解

Stage 2 将 EM 分成两类：

\[
\mathcal{H}=\{\varepsilon_m:\text{hard-BE EM}\},
\]

\[
\mathcal{S}=\{\varepsilon_m:\text{soft/free EM}\}.
\]

hard-BE 部分强制

\[
p_m^{-}=q_m(T)p_m^{+},
\qquad
m\in\mathcal{H}.
\]

因此 hard-BE 的振动贡献为

\[
V_\mathrm{hard}(E;T)
=
\sum_{m\in\mathcal{H}}
b_m
\left[
K(E-\varepsilon_m)
+q_m(T)K(E+\varepsilon_m)
\right],
\]

其中 \(b_m=p_m^{+}\ge0\) 是 hard-BE 的唯一源强度系数。

soft 部分仍保持独立源强度：

\[
V_\mathrm{soft}(E)
=
\sum_{m\in\mathcal{S}}s_m^{+}K(E-\varepsilon_m)
+\sum_{m\in\mathcal{S}}s_m^{-}K(E+\varepsilon_m).
\]

Stage 2 总模型为

\[
\hat y(E;T)
=
Z(E)
+V_\mathrm{hard}(E;T)
+V_\mathrm{soft}(E).
\]

### 5.2 固定温度下的内层问题

对于某个给定温度 \(T\)，Stage 2 形成温度相关设计矩阵

\[
\mathbf{A}_2(T)
=
\left[
\mathbf{B}_\mathrm{hard}(T),
\mathbf{S}_\mathrm{soft},
\mathbf{G}_\mathrm{soft}
\right],
\]

其中

\[
\mathbf{B}_{m}(T)
=
\mathbf{S}_m+q_m(T)\mathbf{G}_m,
\qquad
m\in\mathcal{H}.
\]

拟合目标是 ZLP 扣除后的信号：

\[
\mathbf{y}_2=\mathbf{y}-\mathbf{Z}.
\]

内层求解：

\[
\min_{\mathbf{u}\ge0}
\left\|
\sqrt{\mathbf{W}_2}
\left[
\mathbf{A}_2(T)\mathbf{u}-\mathbf{y}_2
\right]
\right\|_2^2
+\left\|\mathbf{R}_2(T)\mathbf{u}-\mathbf{r}_2\right\|_2^2
+\boldsymbol{\lambda}_{1,2}^{\mathsf{T}}\mathbf{u}.
\]

### 5.3 Stage 2 温度外层 cost

对每个 \(T\)，内层求出最佳 \(\mathbf{u}(T)\)，然后计算

\[
\mathcal{C}(T)
=
\overline{
\left[
\sqrt{w_{2,i}}
(\hat y_i(T)-y_i)
\right]^2
}_{i\in\Omega_\mathrm{fit}}
\lambda_\mathrm{over}
\overline{
\left[
\sqrt{w_{2,i}}
\max(\hat y_i(T)-y_i,0)
\right]^2
}_{i\in\Omega_\mathrm{fit}}
\lambda_\mathrm{reg}
\overline{
\left[
\mathbf{R}_2(T)\mathbf{u}(T)-\mathbf{r}_2
\right]^2
}
\lambda_\mathrm{tr}\mathcal{L}_\mathrm{tr}
\lambda_\mathrm{hf}\mathcal{L}_\mathrm{hf}
\lambda_\mathrm{nbr}\mathcal{L}_\mathrm{nbr}.
\]

代码中对应：

| 数学项 | 代码项 |
|---|---|
| 数据残差 | `mean(rw.^2)` |
| 正向 overshoot | `stage2_overshoot_weight * mean(over.^2)` |
| 正则残差 | `stage2_regularization_cost_weight * mean(reg_res.^2)` |
| transition hard leakage | `stage2_transition_hard_leakage_weight * leak_score` |
| high-frequency overshoot | `stage2_highfreq_overshoot_weight * hf_over_score` |
| 邻域温度先验 | `stage2_neighbor_T_prior_weight * neighbor_prior_score` |

最终温度为

\[
T_\mathrm{fit}
=
\arg\min_T \mathcal{C}(T).
\]

代码支持 grid scan 和 coarse-to-continuous optimization。连续优化先粗扫，再用 `fminbnd` 在局部 bracket 中最小化 \(\mathcal{C}(T)\)。

### 5.4 Stage 2 正则项

Stage 2 的正则与 Stage 1 不完全相同，因为 hard-BE 源强度与 soft 源强度的角色不同。

#### hard 源强度先验

\[
\mathcal{L}_\mathrm{hard-prior}
=
\lambda_\mathrm{prior}
\sum_{m\in\mathcal{H}}
\left(
b_m-p_{m,1}^{+}
\right)^2.
\]

其中 \(p_{m,1}^{+}\) 是 Stage 1 插值得到的 Stokes 源强度。

#### hard anti-Stokes prior

代码还将 Stage 1 的 anti-Stokes 源强度作为弱先验：

\[
\mathcal{L}_\mathrm{anti-prior}
=
\lambda_\mathrm{anti}
\sum_{m\in\mathcal{H}}
c_m
\left(
q_m(T)b_m-p_{m,1}^{-}
\right)^2.
\]

其中 \(c_m\) 与 Stage 1 Stokes 源强度的相对强度有关，并有下限 `stage2_pair_prior_floor`。

#### soft 源强度先验

\[
\mathcal{L}_\mathrm{soft-prior}
=
0.35\lambda_\mathrm{prior}
\sum_{m\in\mathcal{S}}
\left[
(s_m^{+}-p_{m,1}^{+})^2
+(s_m^{-}-p_{m,1}^{-})^2
\right].
\]

#### soft BE

soft EM 不硬锁温度，但保留弱 BE 惩罚：

\[
\mathcal{L}_\mathrm{softBE,2}
=
\lambda_{\mathrm{softBE},2}
\sum_{m\in\mathcal{S}}
\left(
s_m^{-}-q_m(T)s_m^{+}
\right)^2.
\]

#### soft 源强度平滑

对 soft 源强度使用一阶和二阶局部差分：

\[
\mathcal{L}_\mathrm{first}
=
\lambda_\mathrm{D1}
\left(
\|\mathbf{D}_1\mathbf{s}^{+}\|_2^2
+\|\mathbf{D}_1\mathbf{s}^{-}\|_2^2
\right),
\]

\[
\mathcal{L}_\mathrm{smooth,2}
=
\lambda_\mathrm{D2}
\left(
\|\mathbf{D}_2\mathbf{s}^{+}\|_2^2
+\|\mathbf{D}_2\mathbf{s}^{-}\|_2^2
\right).
\]

局部差分只在 EM 间隔小于 `stage2_pm_smooth_max_gap_meV` 的邻近点之间建立，避免跨大能隙强行平滑。

---

## 6. 注意力机制

Stage 2 的 attention 不是机器学习网络，而是一个物理启发的谱线特征评分函数。它的作用包括：

1. 选择温度证据窗口；
2. 增强 Stokes/anti-Stokes 真实特征区域的拟合权重；
3. 选择 hard-BE EM；
4. 诊断谱线哪些区域支配损失。

### 6.1 broad attention window

设 Stage 2 温度硬窗口为

\[
\left[E_\mathrm{hard}^{\min},E_\mathrm{hard}^{\max}\right],
\]

attention 搜索窗口为

\[
\left[E_\mathcal{A}^{\min},E_\mathcal{A}^{\max}\right].
\]

代码保证 attention window 覆盖 hard window：

\[
E_\mathcal{A}^{\min}\le E_\mathrm{hard}^{\min},
\qquad
E_\mathcal{A}^{\max}\ge E_\mathrm{hard}^{\max}.
\]

定义 broad mask：

\[
B_i=
\mathbb{1}
\left(
E_\mathcal{A}^{\min}\le |E_i|\le E_\mathcal{A}^{\max}
\right).
\]

### 6.2 四个 feature channel

Stage 2 使用 Stage 1 结果构造四个特征通道。

#### signal channel

\[
S_i=
\mathrm{medianfilter}
\left(
|y_i-Z_i|
\right).
\]

它表示 ZLP 分离后剩余信号的强度。

#### gradient channel

\[
G_i=
\left|
\nabla
\mathrm{meanfilter}(S_i)
\right|.
\]

它强调峰边、shoulder、快速变化的光学支结构。

#### 源强度谱线通道

Stage 1 得到一组稀疏源强度 sticks。代码将这些 sticks 用高斯展宽成连续曲线：

\[
P_i
=
\sum_m
\frac{p_m^{+}}{p_\mathrm{max}}
\exp\left[
-\frac{1}{2}
\left(
\frac{E_i-\varepsilon_m}{\sigma_\mathcal{A}}
\right)^2
\right]
+
\sum_m
\frac{p_m^{-}}{p_\mathrm{max}}
\exp\left[
-\frac{1}{2}
\left(
\frac{E_i+\varepsilon_m}{\sigma_\mathcal{A}}
\right)^2
\right].
\]

这里 \(p_\mathrm{max}=\max(\mathbf{p}^{+},\mathbf{p}^{-})\)，\(\sigma_\mathcal{A}\) 对应 `stage2_attention_peak_sigma_meV`。

#### residual channel

\[
R_i=|y_i-\hat y_{1,i}|.
\]

它把 Stage 1 没有解释好的结构也纳入注意力。

### 6.3 robust unit scaling

每个 feature channel 会在 broad mask 内归一化到 \([0,1]\)。代码中的 `unit_scale` 可写成

\[
\mathcal{U}(x_i;\Omega)
=
\mathrm{clip}
\left[
\frac{x_i-\mathrm{median}_{j\in\Omega}(x_j)}
{Q_{0.98}(x_j:j\in\Omega)-\mathrm{median}_{j\in\Omega}(x_j)}
,0,1
\right].
\]

使用 98% 分位而不是最大值，是为了避免单个尖峰把全部 attention 压低。

### 6.4 feature score

主线 feature score 为

\[
\mathcal{A}_i
=
0.48\,\mathcal{U}(S_i)
+0.20\,\mathcal{U}(G_i)
+0.22\,\mathcal{U}(P_i)
+0.10\,\mathcal{U}(R_i),
\qquad i\in B.
\]

在 anti-Stokes 侧，算法额外计算

\[
\mathcal{A}_i^{-}
=
0.52\,\mathcal{U}(S_i)
+0.20\,\mathcal{U}(G_i)
+0.22\,\mathcal{U}(P_i)
+0.06\,\mathcal{U}(R_i).
\]

如果 anti 侧信号超过噪声阈值或源强度谱线足够强，则

\[
\mathcal{A}_i
\leftarrow
\max
\left[
\mathcal{A}_i,
\min
\left(
1,
\beta_-\mathcal{A}_i^{-}
\right)
\right],
\]

其中 \(\beta_-\) 对应 `stage2_antistokes_attention_boost`。

### 6.5 mirror symmetry

为了避免只看 Stokes 侧或只看 anti-Stokes 侧，attention 会镜像补偿：

\[
\mathcal{A}_i
\leftarrow
\max
\left[
\mathcal{A}_i,
0.65\,\mathcal{A}(-E_i)
\right].
\]

然后平滑并裁剪：

\[
\mathcal{A}_i
\leftarrow
\mathrm{clip}
\left(
\mathrm{movmean}(\mathcal{A}_i,5),0,1
\right).
\]

物理意义：温度来自 Stokes/anti-Stokes 比值。如果某一侧特征明显，另一侧对应能量也应被检查；镜像项让 attention 尽量成对。

### 6.6 active attention bands

active mask 定义为

\[
\Omega_\mathcal{A}
=
\left\{
i:
B_i=1,\,
\mathcal{A}_i\ge A_\mathrm{thr}
\right\},
\]

并进一步做镜像：

\[
\Omega_\mathcal{A}
\leftarrow
\Omega_\mathcal{A}\cup
\left\{
i:-E_i\in\Omega_\mathcal{A}
\right\}.
\]

其中 \(A_\mathrm{thr}\) 对应 `stage2_attention_threshold`。

### 6.7 attention weight

Stage 2 数据权重乘上 attention 权重：

\[
W_\mathcal{A}(E_i)
=
A_\mathrm{floor}
+g_\mathcal{A}\mathcal{A}_i.
\]

即

\[
\sqrt{w_{2,i}}
=
\sqrt{w_i}\sqrt{W_\mathcal{A}(E_i)}.
\]

其中 \(A_\mathrm{floor}\) 是 `stage2_attention_floor`，\(g_\mathcal{A}\) 是 `stage2_attention_gain`。

### 6.8 sigmoid fit transition

为了避免 attention window 边界处权重突然跳变，Stage 2 会对 `fit_mask` 做 sigmoid 过渡。若某个 band 为 \([E_a,E_b]\)，定义

\[
\alpha(E)
=
\sigma\left(\frac{E-E_a}{\tau_\alpha}\right)
\sigma\left(\frac{E_b-E}{\tau_\alpha}\right).
\]

多个 bands 取最大值：

\[
\alpha_\mathrm{fit}(E)
=
\max_b \alpha_b(E).
\]

然后数据权重再乘

\[
\sqrt{
w_\mathrm{floor}
+(1-w_\mathrm{floor})\alpha_\mathrm{fit}(E_i)
}.
\]

这一步就是为了避免“注意力边界处突然丢峰”或“刚出窗口拟合崩掉”的问题。

### 6.9 hard-BE EM 选择

主线 hard-BE EM 并不是所有 attention 点，而是从 Stage 1 源强度与 attention 中选择较可靠的配对 EM。对每个 EM：

\[
A_m^{+}=\mathcal{A}(+\varepsilon_m),
\qquad
A_m^{-}=\mathcal{A}(-\varepsilon_m).
\]

配对得分：

\[
P_m^\mathcal{A}
=
\sqrt{\max(A_m^{+},0)\max(A_m^{-},0)}.
\]

高能拖尾会 damp：

\[
D_m
=
d_\mathrm{floor}
+(1-d_\mathrm{floor})
\sigma\left(
\frac{E_\mathrm{damp}-\varepsilon_m}{\tau_\mathrm{damp}}
\right).
\]

hard score 为

\[
H_m=P_m^\mathcal{A}D_m.
\]

满足阈值、位于 Stage 2 window 内且间距足够的 EM 会成为 hard seed；再在 hard seed 附近加密：

\[
\varepsilon_m^\mathrm{dense}
=
\varepsilon_m^\mathrm{seed}+\Delta,
\qquad
\Delta\in[-h,h].
\]

主线这样做是为了避免把所有注意力点都硬锁 BE，防止谱型在 gain 侧弱信号、低频混叠或窗口边界处被温度约束带坏。

---

## 7. BE ratio 诊断图

诊断图中的 BE ratio panel 用源强度比值独立检查温度一致性。对每个可配对 EM：

\[
\ln\left(\frac{I_\mathrm{loss}}{I_\mathrm{gain}}\right)
\approx
\ln\left(\frac{p_m^{+}}{p_m^{-}}\right).
\]

BE 关系给出

\[
\frac{p_m^{-}}{p_m^{+}}
=
\exp\left[
-\frac{\varepsilon_m}{k_\mathrm{B}T}
\right],
\]

因此

\[
\ln\left(\frac{p_m^{+}}{p_m^{-}}\right)
=
\frac{\varepsilon_m}{k_\mathrm{B}T}.
\]

诊断图中过原点线性拟合：

\[
y_m=\beta\varepsilon_m,
\qquad
\beta=\frac{1}{k_\mathrm{B}T}.
\]

即

\[
T_\mathrm{ratio}
=
\frac{1}{k_\mathrm{B}\beta}.
\]

权重取

\[
\omega_m
\propto
\sqrt{p_m^{+}p_m^{-}},
\]

并归一化到 \([0,1]\)。现在诊断图会优先画 attention window 内所有可配对源强度，而不是只画 hard-BE 子集，这样更适合观察“整个注意力证据区域是否支持同一温度”。

---

## 8. 误差指标与置信区间

### 8.1 R2

\[
R^2
=
1-
\frac{
\sum_i (y_i-\hat y_i)^2
}{
\sum_i (y_i-\bar y)^2
}.
\]

它衡量总谱形复现程度，但强 ZLP 会使 \(R^2\) 很容易接近 1，因此不能单独作为温度可靠性的依据。

### 8.2 NRMSE

\[
\mathrm{NRMSE}
=
\frac{
\sqrt{\frac{1}{N}\sum_i(y_i-\hat y_i)^2}
}{
\max_i y_i-\min_i y_i
}.
\]

### 8.3 weighted NRMSE

\[
\mathrm{wNRMSE}
=
\frac{
\sqrt{
\frac{1}{N}\sum_i
\left[
\sqrt{w_i}(y_i-\hat y_i)
\right]^2
}
}{
\max_i(\sqrt{w_i}y_i)-\min_i(\sqrt{w_i}y_i)
}.
\]

wNRMSE 比普通 NRMSE 更关注被权重强调的 signal/high-frequency/tail 区。

### 8.4 high-frequency NRMSE

定义高频区域

\[
\Omega_\mathrm{hf}=\{i:|E_i|\ge E_\mathrm{hf}\}.
\]

\[
\mathrm{NRMSE}_\mathrm{hf}
=
\frac{
\sqrt{
\frac{1}{|\Omega_\mathrm{hf}|}
\sum_{i\in\Omega_\mathrm{hf}}
(y_i-\hat y_i)^2
}
}{
\max_{i\in\Omega_\mathrm{hf}}y_i
-
\min_{i\in\Omega_\mathrm{hf}}y_i
}.
\]

它用于发现“总 R2 很高但高频峰丢失或拖尾抬起”的问题。

### 8.5 能量依赖 CI

CI 是 energy-dependent confidence interval，不是温度误差。代码用局部 residual MAD 估计每个能量点的不确定性。

先取残差：

\[
\rho_i=I_i-\hat I_i.
\]

局部中心：

\[
m_i=
\mathrm{movmedian}_{j\in\mathcal{N}(i)}
(\rho_j).
\]

局部 MAD：

\[
\mathrm{MAD}_i=
\mathrm{movmedian}_{j\in\mathcal{N}(i)}
|\rho_j-m_j|.
\]

估计标准差：

\[
\sigma_i=
1.4826\,\mathrm{MAD}_i.
\]

加入噪声下限：

\[
\sigma_i
\leftarrow
\max
\left(
\sigma_i,
f_\mathrm{min}\max_j|I_j|
\right).
\]

置信区间半宽：

\[
\Delta I_i^\mathrm{CI}
=
z_\alpha\sigma_i.
\]

例如 95% 时 \(z_\alpha=1.96\)。最终

\[
\hat I_i^\mathrm{lower}
=
\hat I_i-\Delta I_i^\mathrm{CI},
\]

\[
\hat I_i^\mathrm{upper}
=
\hat I_i+\Delta I_i^\mathrm{CI}.
\]

物理意义：CI width 大的能量区域表示局部残差结构大、噪声大或模型不确定性强；它是谱线强度维度的不确定性，不是直接的温度不确定性。

### 8.6 Stage 2 温度不确定性

Stage 2 对温度扫描得到 \(\mathcal{C}(T)\)。温度不确定性来自 cost valley 的宽度。直观地：

- valley 窄：温度被 Stokes/anti-Stokes 比值强约束，uncertainty 小；
- valley 平：很多温度都能拟合谱形，uncertainty 大；
- 最优点在边界：uncertainty 被放大。

诊断中的 `T_uncertainty` 是这个 temperature cost profile 的估计宽度，单位 K。

### 8.7 reliability

Stage 2 reliability score 是经验组合：

\[
R_\mathrm{rel}
=
\mathrm{clip}
\left[
0.35R_\mathcal{A}
+0.30R_\mathrm{fit}
+0.35R_\mathrm{unc},
0,1
\right],
\]

其中

\[
R_\mathcal{A}
=
\min
\left(
1,\frac{f_\mathcal{A}}{0.18}
\right),
\]

\[
R_\mathrm{fit}
=
\min
\left[
1,
\max
\left(
0,
\frac{R^2-0.999}{0.001}
\right)
\right],
\]

\[
R_\mathrm{unc}
=
\exp
\left(
-\frac{\Delta T}{220}
\right).
\]

这里 \(f_\mathcal{A}\) 是 attention selected fraction，\(\Delta T\) 是 `T_uncertainty`。它不是严格统计置信度，而是 fit quality、attention evidence 和 temperature valley sharpness 的综合可靠性指标。

---

## 9. Loss budget

诊断图中的 loss budget 将总损失按能区与正则拆开。能区包括：

\[
\Omega_\mathrm{core}: |E|<E_c,
\]

\[
\Omega_\mathrm{shoulder}: E_c\le |E|<E_\mathrm{min},
\]

\[
\Omega_\mathrm{anti}: E\le -E_\mathrm{min},\ |E|<E_\mathrm{hf},
\]

\[
\Omega_\mathrm{stokes}: E\ge E_\mathrm{min},\ |E|<E_\mathrm{hf},
\]

\[
\Omega_\mathrm{hf}: E_\mathrm{hf}\le |E|<E_\mathrm{tail},
\]

\[
\Omega_\mathrm{tail}: |E|\ge E_\mathrm{tail}.
\]

每个区域的数据损失为

\[
\mathcal{L}_{\Omega}
=
\sum_{i\in\Omega}
w_i(y_i-\hat y_i)^2.
\]

额外项包括：

\[
\mathcal{L}_\mathrm{reg},
\quad
\mathcal{L}_\mathrm{ZLPguard},
\quad
\mathcal{L}_\mathrm{lowfreq},
\quad
\mathcal{L}_\mathrm{softBE},
\quad
\mathcal{L}_\mathrm{tail}.
\]

最终显示的是比例：

\[
f_j=
\frac{\mathcal{L}_j}{\sum_k\mathcal{L}_k}.
\]

这个图的用途是判断当前拟合的主要矛盾。例如：

- `zlp_core` 高：零峰 core 或 align 有问题；
- `stokes signal` / `anti signal` 高：温度窗口或稀疏模态网格还不够；
- `high frequency` 高：高频峰丢失或拖尾过高；
- `zlp_signal_guard` 高：ZLP 试图吃掉信号；
- `regularization` 高：谱形和正则先验冲突，可能过硬或 EM 太少。

---

## 10. 参数表

### 10.1 Stage 1 模型与 ZLP 参数

| 参数 | 默认值 | 作用 | 调参建议 |
|---|---:|---|---|
| `num_zlp_components` | `4` | ZLP 分量数 | 低频被 ZLP 吞掉时不要盲目增加；先检查宽度上限 |
| `zlp_model_type` | `asym_pearson` | ZLP 线形 | Pearson 更平滑；pVoigt 可用于特殊 ZLP 形状 |
| `zlp_component_width_max_meV` | `[3.3, 6.5, 12, 21]` | 各 ZLP 分量宽度上限 | 低频被吃时降低宽分量上限 |
| `zlp_wide_component_l2` | `6e-3` | 宽 ZLP 分量幅度惩罚 | 宽尾吃信号时增大 |
| `kernel_mode` | `fitted_zlp` | 卷积核来源 | 默认使用拟合总 ZLP；仅在明确需要窄核对照时改用 core_component |
| `em_min_meV` | `8` | 最小 EM | 太低会和 ZLP core 混叠 |
| `em_step_meV` | `6` | Stage 1 基础 EM 步长 | 谱线不光滑/丢峰时减小 |
| `em_lowfreq_dense_enabled` | `true` | 低频 EM 加密 | 低频 shoulder 重要时保留 |
| `stage1_fixed_EM` | `[]` | datacube 自动稀疏源强度后固定 EM | 用于全空间统一源强度维度 |

### 10.2 Stage 1 损失权重

| 参数 | 默认值 | 对应公式 | 作用 |
|---|---:|---|---|
| `core_weight` | `0.22` | \(w_\mathrm{core}\) | 降低 ZLP core 对总损失的统治 |
| `shoulder_weight` | `0.72` | \(w_\mathrm{shoulder}\) | 控制 core 外 shoulder |
| `signal_weight` | `1.45` | \(w_\mathrm{sig}\) | 强调主要振动信号区 |
| `highfreq_weight` | `2.2` | \(w_\mathrm{hf}\) | 强调高频弱峰 |
| `tail_weight` | `2.5` | \(w_\mathrm{tail}\) | 抑制高频拖尾抬升 |
| `lambda_be_soft` | `1e-3` | \(\lambda_\mathrm{BE}\) | Stage 1 soft BE |
| `lambda_pm_smooth` | `3e-3` | \(\lambda_\mathrm{smooth}\) | 源强度二阶平滑 |
| `lambda_pm_l2` | `2e-5` | \(\lambda_2\) | 源强度 L2 |
| `lambda_pm_l1` | `3e-5` | \(\lambda_1\) | 源强度稀疏 |
| `tail_overshoot_constraint` | `true` | tail inequality | 限制高频过拟合抬尾 |

### 10.3 ZLP guard 参数

| 参数 | 默认值 | 作用 |
|---|---:|---|
| `zlp_signal_guard_enabled` | `true` | 开启信号区 ZLP 抑制 |
| `zlp_signal_guard_start_meV` | `16` | signal guard 起点 |
| `zlp_signal_guard_end_meV` | `120` | signal guard 终点 |
| `zlp_signal_guard_weight` | `0.20` | guard 强度 |
| `zlp_lowfreq_guard_enabled` | `true` | 开启低频 guard |
| `zlp_lowfreq_guard_start_meV` | `14` | 低频 guard 起点 |
| `zlp_lowfreq_guard_end_meV` | `58` | 低频 guard 终点 |
| `zlp_lowfreq_guard_weight` | `0.05` | 低频 guard 强度 |

### 10.4 Stage 2 温度与 attention 参数

| 参数 | 默认值 | 作用 | 调参建议 |
|---|---:|---|---|
| `run_stage2` | `false` | 是否启用温度层 | 有 Stokes/anti-Stokes 证据时开启 |
| `stage2_window_meV` | `[40,150]` | hard-BE 候选温度窗口 | 太宽会纳入不可靠区域 |
| `stage2_attention_window_meV` | `[40,180]` | attention 搜索窗口 | 要覆盖高频峰，如 150 meV 峰 |
| `stage2_attention_threshold` | `0.18` | active attention 阈值 | 漏峰时降低，噪声多时提高 |
| `stage2_attention_gain` | `9` | attention 权重增益 | 特征区拟合不足时增大 |
| `stage2_antistokes_attention_boost` | `1.35` | anti 侧 attention boost | gain 侧弱但重要时增大 |
| `stage2_fit_sigmoid_taper_meV` | `10` | attention 边界过渡宽度 | 边界突然丢峰时增大 |
| `stage2_fit_sigmoid_min_weight` | `0.08` | 窗口外权重地板 | 窗口外谱形坏时增大 |
| `stage2_hard_em_threshold` | `0.30` | hard EM score 阈值 | 温度证据不足时降低会更激进 |
| `stage2_dense_em_step_meV` | `2` | hard seed 附近 EM 加密步长 | 不光滑时可降低 |
| `stage2_soft_continuum_enabled` | `true` | soft 区连续源强度补充 | 保护谱形 |

### 10.5 Stage 2 正则与温度可靠性

| 参数 | 默认值 | 作用 |
|---|---:|---|
| `stage2_pm_prior_weight` | `1e-3` | hard/soft 源强度靠近 Stage 1 初值 |
| `stage2_antistokes_prior_weight` | `2e-3` | hard 区 anti-Stokes 初值先验 |
| `stage2_soft_be_weight` | `1e-2` | soft 区 BE 软约束 |
| `stage2_pm_smooth_weight` | `4e-4` | soft 源强度二阶平滑 |
| `stage2_pm_firstdiff_weight` | `2e-4` | soft 源强度一阶平滑 |
| `stage2_pm_sparse_weight` | `8e-6` | Stage 2 源强度稀疏 |
| `stage2_overshoot_weight` | `0.35` | Stage 2 总体 overshoot 惩罚 |
| `stage2_highfreq_overshoot_weight` | `0.55` | 高频过拟合抬尾惩罚 |
| `stage2_ratio_min_pairs` | `5` | ratio 温度估计最少配对数 |
| `stage2_ratio_reliable_spread_K` | `140` | ratio 温度 IQR 可靠阈值 |

---

## 11. 物理解释：为什么要这样分层

### 11.1 Stage 1 的物理目标

Stage 1 的首要目标不是温度，而是谱形分解：

\[
\text{raw spectrum}
\rightarrow
\text{ZLP}
+
\text{sparse convolved vibrational response}.
\]

这一步的物理意义是通用 EELS 背底/零峰分离。即使数据没有可靠 anti-Stokes 信号，Stage 1 仍应能给出有用的 ZLP 扣除和稀疏声子模态源强度分布。它是 DOS-like 的源空间表示，但不是绝对标定的声子 DOS。

### 11.2 Stage 2 的物理目标

Stage 2 才引入温度。温度不是靠全谱强行拟合出来，而是靠 Stokes/anti-Stokes 中具有物理对应关系的能量窗口：

\[
\ln(I_\mathrm{loss}/I_\mathrm{gain})
\propto
\varepsilon.
\]

因此 Stage 2 既要硬约束一部分可信 EM，也要允许不可靠区域保持 soft/free，避免为了温度牺牲谱型。

### 11.3 attention 的物理意义

attention 的本质是一个“温度证据选择器”：

- \(S_i\) 找真实信号强度；
- \(G_i\) 找峰边和 shoulder；
- \(P_i\) 找 Stage 1 已经识别的稀疏模态位置；
- \(R_i\) 找 Stage 1 仍没解释好的结构；
- anti boost 防止 gain 侧弱信号被忽略；
- mirror symmetry 防止只看单侧峰。

这不是黑箱学习，而是物理启发的 feature engineering。

### 11.4 为什么主线不把所有 attention 都 hard BE

主线只把部分高置信 EM 设为 hard BE，因为实际数据中常见：

1. gain 侧信号接近噪声；
2. 低频区域与 ZLP tail 混叠；
3. 注意力窗口边界处 feature 不完整；
4. 某些峰只有 loss 侧明显；
5. ZLP 分解误差会伪造 anti 侧源强度。

如果整个 attention window 都硬 BE，温度线可能更漂亮，但谱型可能被带坏。`version1.3_fullhardBE_beta` 正是为了测试这个假设而保留的对照版本。

---

## 12. 可用于文章的方法段落草稿

本工作采用一种分层的 zero-loss peak estimation and sparse thermometry 方法。首先，将实验 EELS 谱表示为 zero-loss peak 与由 ZLP kernel 卷积的稀疏振动源谱之和：

\[
\hat I(E)=Z(E;\theta_Z)
+\sum_m p_m^{+}K(E-\varepsilon_m)
+\sum_m p_m^{-}K(E+\varepsilon_m).
\]

其中 \(Z(E)\) 由多分量不对称 Pearson 或 pseudo-Voigt 函数表示，卷积核 \(K(E)\) 由拟合得到的 ZLP 分量归一化得到。第一阶段通过非负加权最小二乘拟合 \(Z(E)\)、\(p_m^{+}\) 与 \(p_m^{-}\)，并引入源强度平滑、稀疏性、ZLP signal guard 和 soft Bose--Einstein penalty，以保证谱形复现和物理可解释性。第二阶段基于第一阶段的源强度与残差构造 attention score，用于选择 Stokes/anti-Stokes 均具有可靠特征的温度证据窗口。在这些高置信 EM 上施加 hard Bose--Einstein detailed-balance 关系 \(p_m^{-}=p_m^{+}\exp[-\varepsilon_m/(k_\mathrm{B}T)]\)，而窗口外或低置信区域保留独立的 soft 源强度，以避免温度约束破坏谱线形状。最终温度通过最小化温度相关的加权残差、正则残差、overshoot penalty 与邻域先验得到。拟合可靠性由 attention 覆盖率、谱线拟合优度和 temperature-cost valley 宽度共同评估。

---

## 13. 与代码字段的对应关系

本文前文统一使用“稀疏声子模态源强度”。为兼容 MATLAB 结果字段，软件在字段名中将 **phonon-mode intensity** 缩写为 **PM**；`PM_stokes` 和 `PM_antistokes` 分别保存 loss/Stokes 与 gain/anti-Stokes 侧的非负源强度系数。

| 文档概念 | 主要代码字段 |
|---|---|
| Stage 1 总拟合 | `result.fitted_spectrum`, `result.fitted_spectrum_norm` |
| Stage 1 ZLP | `result.zlp_only`, `result.stage1.zlp_norm` |
| Stage 1 sparse vib | `result.vib_only`, `result.stage1.vib_norm` |
| Stage 1 sparse source intensity | `result.EM_final`, `result.PM_stokes`, `result.PM_antistokes` |
| Stage 2 总拟合 | `result.stage2.model_norm * result.norm_val` |
| Stage 2 ZLP | `result.stage2.zlp_norm * result.norm_val` |
| Stage 2 vib | `result.stage2.vib_norm * result.norm_val` |
| Stage 2 hard vib | `result.stage2.vib_hard_norm * result.norm_val` |
| Stage 2 soft vib | `result.stage2.vib_soft_norm * result.norm_val` |
| Stage 2 hard EM | `result.stage2.EM_hard` |
| Stage 2 soft EM | `result.stage2.EM_soft` |
| attention score | `result.stage2.attention.feature_score` |
| attention bands | `result.stage2.attention.bands_meV` |
| CI | `result.confidence_interval` |
| loss budget | `result.loss_breakdown` |
| T uncertainty | `result.stage2.T_uncertainty` |
| reliability | `result.stage2.reliability_score` |
