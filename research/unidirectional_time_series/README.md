# Green–Laplace 单向波时间序列研究

本阶段完成项目隔离、实现审查和可重复基线，不开展新的完整物理推导。

更新：用户已选择输入为已知/已分离的 **eta1(t)**，并授权先做 OW3D 二阶和频试验，与 spectral MF12、VWA、Walker 比较。下文的初始化待定项保留为历史背景；当前执行和证据见 [二阶试验记录](SECOND_ORDER_PILOT.md)。本阶段不进行总波面反演。

## 输入问题的区别

1. **已知一阶谱**：给定母波幅值、波数和相位，在固定位置采样不同时间的束缚波。现有 `gl_spectral_surface(...,t)` 可以直接完成本基线。
2. **已分离 eta1(t)**：需要声明时间谱归一化、频率到波数映射、传播方向、水深、记录长度与采样，以及非空间网格对齐频率的处理；尚无本项目认证的直接时间序列入口。
3. **实测总波面 eta(t)**：还需要反演或分离自由一阶波和束缚分量；各时间 FFT 分量不能自动作为线性母波。是否研究该逆问题尚待决定。

## 可复用的实现及边界

`gl_spectral_coefficients` 准备一阶谱与线性色散频率；`gl_spectral_surface` 将一阶解析空间谱乘以 `exp(-i*omega*t)`，再调用现有 GL 图。统一接口输出 eta11/22/33 与真正的自由面势 psi11/22/33，后者保留自由面 Taylor 项。这里的求和不是完整三阶物理解：未包含差频或主谐波修正等未支持分量。

`gl_pure_sum_order4` 是独立、实验性、无 Stokes correction 的四阶接口，使用无量纲解析输入输出；不能通过统一接口的 order=4 调用。`paper/order4/` 包含既有 GL–WIT 复现材料；`symbolic/` 包含 Wolfram 源和冻结接口。保留历史符号身份不代表启用历史 Stokes 修正。

已审查根 README、CHANGELOG、依赖和无 Stokes 迁移说明、上述两个统一 API、最小示例、发布测试以及四阶论文 README。现有最小示例有 ky=[0,1]，属于既有发布回归；新的时间采样例专用 ky=0，不将原示例称为单向波。

## 本次合成例明确采用的约定

这些约定仅定义初始化夹具，不替用户决定正式实验输入：SI 单位，g=9.81 m/s²，h=1 m，无流，kx=[2,3] rad/m，ky=0，向 +x 传播；A=a+i*b，a=[0.010,0.006] m，b=[0,0.002] m。

`eta1(x,t)=Re sum A_j exp(i*k_j*x-i*omega_j*t)`，`omega_j=sqrt(g*k_j*tanh(k_j*h))` 为正角频率，f=omega/(2*pi)。因此正时间指数的 FFT 约定与此解析表示的符号关系必须显式处理，不能直接照搬空间谱。psi 单位为 m²/s。

空间域 Lx=Ly=2*pi m，Nx=32、Ny=4（保留现有二维执行器，沿 y 恒定），最高三阶和波数为 9，小于 x Nyquist=16；所有母波 kh>0.5。固定位置为网格点 x=0、y=0；t=[0,0.125,0.25,0.5] s，仅为稀疏调用检查，不作时间 FFT 分辨率或完整记录的声明。eta22_rank=6，三阶使用发布默认图。

## 运行和验证标准

从仓库根运行：

```matlab
restoredefaultpath; addpath(pwd); setup_green_laplace;
addpath('tests'); run_release_tests;
run('examples/run_minimal_example.m');
addpath('research/unidirectional_time_series');
run_known_spectrum_time_sampling;
```

小例逐分量检查：t=0 与同输入原 API 基线一致；非零 t 与先旋转一阶复幅值再在 t=0 调用一致（原始相对 L2、相对 Linf <1e-12）；线性 eta 和 psi 对显式线性表达式一致；沿 y 恒定，输出有限。原 API 不修改、不增补高阶参考场，不做任何对齐或拟合。结果写入 `results/unidirectional_time_series/`。这只是相位处理和 API 一致性验证，不是独立高阶物理精度认证，也不是 eta1(t) 重构或总波面反演算法。

正式输入还须明确：eta1 还是总 eta、实际水深/流条件、测点坐标和时间零点、传播方向、采样间隔/长度/缺测、相位和谱归一化、所需阶次及输出变量。均值和混合符号分量不在首阶段范围。不要为初始化自行去均值、滤波或拟合。

推荐下一项：确定输入类别后，先用同一已知母波合成 eta1(t)，在 MATLAB 中核查时间谱幅值/相位与线性色散映射，明确时间频率与空间 FFT 网格不一致的处理，再讨论直接时间序列执行器。若选择总 eta，先制定分离/逆问题的可辨识性与验证方案。任何新公式先用 Wolfram 精确认证，随后用 MATLAB 做同输入、同单位的独立场比较。
