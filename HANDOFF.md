# 初始化交接 — 2026-09-25

并发安排更新：暂不修改OW3D求解器。建议队列为4组/16个单相位运行，
短资源试跑后先4个、目标8个并发；12/16并发分别要求实测单例峰值内存
不超过40/30 GiB（600 GiB总预算，含25%余量），还需通过CPU/I/O检查。
当前尚无这套网格的峰值RSS实测，不能宣称已验证最大并发数。计划详见
`research/directional_wave_data/REMOTE_OW3D_PLAN.md` 的任务队列部分。

最新远程计算准备：Akp=0.12 结果已提交并推送至 `1038842`。2026-09-25
已检查远程资源和 MATLAB 启动；尚未创建远程研究目录或启动新 OW3D。
已确认原生 EP 支持只存 eta/自由面 phi。用户随后提出局部 kinematics 输出，
并接受三倍峰值频率每周期20点：现建议 0.2 s 输出、120--360 s 保存、
积分到360 s；先保留完整x和中心附近9条y网格线，以沿用既有正kx一阶
提取方式。四相位含所有原生 kinematics 字段约87.79 GiB，加稀疏EP约88.7 GiB。
若经验证改用局部时间分离，81x9小区域可减至约3.47 GiB；不能直接替换当前
空间投影。eta20需另查时间窗截断，裁剪记录还需修正初始谱的绝对时间相位。
无需改OW3D源码；所有原始数据及后处理留在远程。完整资源、输出依据和预算见
`research/directional_wave_data/REMOTE_OW3D_PLAN.md`。计划尚未执行；下文为
既有研究和初始化历史记录，旧“远程 MATLAB 未找到”条目已由此次核查更新。

## 当前研究进度：首个二阶时间序列试验已完成

最新用户范围：off-centerline测点的采样eta22峰值至少为同x中心线的1/3。
按此规则，旧半波长横向点仅约23%，更远横向点约0.04%，均不属于主要关注区。
Akp=0.02 的约±0.25λ、±0.35λ四点，峰值比例约70%/44%，主波组
eta22误差0.239%--0.335%、psi22误差0.189%--0.394%。现已按同规则完成
Akp=0.12：四点峰值比例约54%--71%，eta22误差2.767%--3.458%、psi22误差
1.735%--2.340%。强陡度数据的实际偏移为±52.734375/±70.3125 m，中心
(11250,4500) m；原始网格/域与弱陡度组不同。结果另存于
`results/directional_amplitude_gate_akp012/`，未覆盖Akp=0.02数据。
见 `research/directional_wave_data/AMPLITUDE_QUALIFIED_PROBES.md` 和
`results/directional_amplitude_gate/`。原远点压力测试保留，不作误差修正。

最新扩展：已执行13个方向性“算例—测点”组合，覆盖更远的 x±3λ、y±λ/±1.5λ，
kh=1/2/5、展宽5/15/25度和Akp=.02/.12，并加入方向性eta20。中心/沿x方向
eta22、psi22大多较好；远横向点发生明显失准，角度减半后仍未收敛。eta20
部分kh=1结果约2%--4%，kh=5某些结果约15%--57%，不能沿用eta22的精度结论。
独立频率带宽测试采用明确标注的合成谱与二阶MF12，不冒充OW3D证据。
完整范围、失败结果和解释见 `research/directional_wave_data/EXTENDED_SWEEP.md`；
结果在 `results/directional_sweep/` 和 `results/directional_bandwidth/`。未调参掩盖失败。

最新：方向性联合输入已做固定五点检查：中心、x方向约正负一个峰值波长、
y方向约正负半个峰值波长。7.5度方向离散下，主波组 eta22 误差为
0.1911%--0.6168%，psi22 为 0.1402%--0.8865%。使用同一初始谱和同一算法，
每点只输入其本地 eta1(t)，未按高阶误差选点或调参。每个主波组仍只有13个
原始保存时刻。见 `research/directional_wave_data/MULTIPROBE_TRIAL.md`，
图和汇总在 `results/directional_joint_input_multiprobe/`；中心旧结果未覆盖。

最新方向性试验：初始空间谱提供复方向权重，测点 eta1(t) 约束各频率的方向和。
test1、kh=1、展宽标签25度、Akp=0.02 的首例已执行；相位约定核查后，7.5度
方向离散的主波组 eta22/psi22 相对 L2 为 0.2668%/0.1402%。全窗误差为
19.74%/22.70%，不能把主波组结果表述为全记录认证；方向尾部存在高条件数。
详见 `research/directional_wave_data/JOINT_INPUT_PILOT.md`。只使用
`results/directional_joint_input/verified_convention/` 里的物理结果；上一级
保留的首轮未核查时间符号结果是无效调试记录。未进行高陡度方向性试验或新 OW3D 计算。

当前综合状态已整理到 `research/unidirectional_time_series/STATUS.md`。
新增自由面势 psi11/psi22 首轮试验，同一 eta1 输入、无 Stokes 修正；kh=1、
Alpha=1、Akp=0.02/0.12 的 GL2+2 psi22 主波组误差为 0.2597%/2.0860%，
二阶 MF12 为 0.05454%/1.9762%。psi20、psi33 尚未开展；不将自由面势
结果表述为任意垂向位置的体势验证。详见 SURFACE_POTENTIAL_PILOT.md。

最新推进：用户已授权 commit/push；研究分支已推送。完成非零差频 eta20 和
正和频 eta33 时间序列试验，仍为 Alpha=1、kh=1、Akp=0.02/0.12 的边界修正
原始数据。eta33 不计算、不绘制三阶 MF12；eta20 保留二阶 MF12 参考。
主波组误差：eta20 GL16 为 2.2463%/4.3139%；eta33 GL8 为 2.1583%/11.4886%。
eta20 是统一非零低频投影后的诊断，不包含严格零频均值；eta33 有明显求积阶数
依赖，尚不宣称 GL8 已收敛。详见 `research/unidirectional_time_series/ETA20_ETA33_TRIAL.md`。
图和字段位于 `results/unidirectional_time_series/eta20_eta33_boundary_alpha1_akp002/`
及 `..._akp012/`。公开 main 与原空间实现未修改。

最新：已比较 Alpha=1、kh=1 的 Akp=0.02 与 0.12，使用既有 ESC 项目的
原始测点 CSV，并明确区分旧边界污染批次和已有边界修正批次，未导入 ESC/
Stokes 修正。修正批次 Akp=0.12 主波组相对 OW3D 的 L2：GL12 3.8486%、
spectral MF12 3.8566%、VWA 3.3047%、Walker 8.2223%。旧批次的约 53% 差异
已复现，并有既有封闭水槽初始波包越界诊断作为来源解释。
详见 `research/unidirectional_time_series/ALPHA1_STEEPNESS.md`；主图在
`results/unidirectional_time_series/ow3d_boundary_kh1_alpha1_akp012/eta22_main_group.png`。
旧结果保留在 `ow3d_compact_kh1_alpha1_akp*`，不应作为首选有限水深验证数据。

后续按用户要求完成 kh=0.5、Alpha=1、Akp=0.02。结果在
`results/unidirectional_time_series/ow3d_kh0p5_alpha1_akp002/`。
GL12 与 spectral MF12 的相对 L2 差为 0.05252%，但相对 OW3D 全窗误差
分别为 23.4560% 和 23.4545%。OW3D 第二相位分量在主波包之后仍有明显振荡，
原因尚未确定；保留完整时间窗，不能把该差异简单归为 GL 误差。
原 kh=1 数据和图保持不变。详细对照和输入投影比例见 SECOND_ORDER_PILOT.md。

用户已明确正式输入为 eta1(t)，授权与 OW3D、spectral MF12、VWA、Walker 比较。
已实现独立的 GL 二阶逐对时间重构，保留公开空间接口不变；Wolfram 四项代数检查、
MATLAB 冻结核/空间接口/相位归一化检查通过。首例 kh=1、Alpha=1、Akp=0.02，
四相位共 1404 个原始快照已读取并哈希。全窗相对 OW3D 第二谐波记录的 L2：
GL6 0.2363%、GL8 0.1078%、GL12 0.09423%、spectral MF12 0.09420%、
VWA 0.8954%、Walker 11.50%。这是四相位谐波记录比较，不是严格扰动阶次分离认证。

执行说明、约定和限制见 `research/unidirectional_time_series/SECOND_ORDER_PILOT.md`；
图、原始时间序列、指标和哈希见 `results/unidirectional_time_series/ow3d_kh1_alpha1_akp002/`。
可运行 `plot_ow3d_eta22_pilot` 从保存的 MAT 结果重新绘图，无需再读取原始快照。
建议下一例保持规则不变，检查同水深同幅值的 Alpha=8。未推送、未部署远程、未修改历史数据。
下文初始化记录中的“唯一待回答的问题”已经解决。

## 当前目录入口（2026-09-25 清理后）

- 唯一的日常 GL 工作目录是 `C:\Users\spet5947\Documents\green-laplace-unidirectional-time-series`。在 Codex 中选择此目录；旧任务绑定的发布目录和 SWORD 路径不再作为工作入口。
- 本仓库的 `main` 保留公开版基线，`codex/unidirectional-time-series` 保留时间序列研究。GitHub `origin` 地址不变；目录清理没有合并或推送研究分支。
- 旧发布克隆、各论文仓库、SWORD worktree 和旧压缩备份统一存放在 `C:\Users\spet5947\Documents\Archive\GL-cleanup-20260925-025101\workspaces`，保留原目录名。
- 归档入口和移动清单在 `C:\Users\spet5947\Documents\Archive\GL-cleanup-20260925-025101`。历史目录作为只读资料；不要把论文、旧仓库或机器本地备份重新复制到公开代码目录。
- 下文记录的是初始化当时的路径和命令。需要查询旧资料时，先按归档中的 `move-plan.json` 查找新位置，不要依赖旧绝对路径。

## 原初始化记录

- 本地：`C:\Users\spet5947\Documents\green-laplace-unidirectional-time-series`，独立 Git clone（git-common-dir 为本仓库 `.git`），可作为本地 Codex 工作目录打开。
- 上游：<https://github.com/Xinyu01091002/green-laplace-wave-reconstruction.git>。
- 采用基线：`eaca1576854ad29ad047527903145271b5fddb4d`；初始化前及验证后均用 ls-remote 核对 main，未发现相对给定 SHA 的更新。
- 分支：`codex/unidirectional-time-series`。初始化提交用 `git log -1 --format=%H` 获取（交接文件随该提交保存，避免自引用 SHA）。没有 push。
- 原本地发布仓库 HEAD 同上，前后状态干净，未修改。历史来源和归档未访问或变更。
- 远程别名：`60.188.112.99:60093`。计划路径 `/home/lxy/green-laplace-unidirectional-time-series` 已检查不存在；本次未创建目录、未复制源码、未运行远程研究计算。后续需从明确本地提交部署并记录 SHA，不能复制未记录的工作目录。

## 实际执行命令

在原仓库检查 `git status --short`、`git rev-parse HEAD`、`git remote -v`，用 `Get-Item` 检查新路径（不存在），用 `Get-Command matlab,wolframscript,git,ssh,g++,clang++` 检查本机工具。初始化命令：

```powershell
git ls-remote https://github.com/Xinyu01091002/green-laplace-wave-reconstruction.git refs/heads/main
git clone https://github.com/Xinyu01091002/green-laplace-wave-reconstruction.git C:/Users/spet5947/Documents/green-laplace-unidirectional-time-series
# 以下在新仓库执行
git switch -c codex/unidirectional-time-series
New-Item -ItemType Directory -Force artifacts/initialization
matlab -batch "restoredefaultpath; addpath(pwd); setup_green_laplace; assert(~contains(path,'SWORD','IgnoreCase',true)); addpath('tests'); summary=run_release_tests; save('artifacts/initialization/release_summary.mat','summary'); set(groot,'defaultFigureVisible','off'); run('examples/run_minimal_example.m'); assert(startsWith(which('gl_spectral_surface'),pwd));" -logfile artifacts/initialization/release.log
matlab -batch "restoredefaultpath; addpath(pwd); setup_green_laplace; addpath('research/unidirectional_time_series'); report=run_known_spectrum_time_sampling; issues=checkcode('research/unidirectional_time_series/run_known_spectrum_time_sampling.m','-id'); disp(issues); assert(isempty(issues));" -logfile artifacts/initialization/time_sampling.log
git diff --check
git rev-parse --git-common-dir
git ls-remote origin refs/heads/main
wolframscript -code '$Version'
```

远程使用 `ssh -o BatchMode=yes -o ConnectTimeout=10 60.188.112.99:60093`，执行了：

```sh
pwd
if [ -e /home/lxy/green-laplace-unidirectional-time-series ]; then ls -la /home/lxy/green-laplace-unidirectional-time-series; else echo TARGET_ABSENT; fi
for t in matlab wolframscript gcc g++ cmake git; do command -v "$t"; done
uname -sr
g++ --version | head -1
cmake --version | head -1
wolframscript -code \$Version
ls -d /usr/local/MATLAB/*/bin/matlab /opt/MATLAB/*/bin/matlab /home/lxy/MATLAB/*/bin/matlab 2>/dev/null
test -f /usr/include/fftw3.h && echo FFTW_HEADER_PRESENT
command -v python3
command -v patch
timeout 15s wolframscript -code 1+1
```

本地初始化文件通过补丁创建。范围明确的提交命令：

```powershell
git add -- AGENTS.md HANDOFF.md research/unidirectional_time_series/README.md research/unidirectional_time_series/run_known_spectrum_time_sampling.m
git diff --cached --check
git commit -m "Initialize unidirectional time-series research and known-spectrum baseline"
git status --short
git log -1 --format=%H
```

## 实测结果与边界

- MATLAB R2022b：8/8 发布检查通过，涵盖线性前端、有序二阶比较、三阶自由面势、域检查、冻结接口、已有 eta20 诊断、无 Stokes 迁移和 Two-Scale 迁移。运行既有诊断不扩大本项目研究范围。
- 原最小示例成功：eta L2=5.35546377e-01，psi L2=1.13904433e+00；六个公开分量齐全。图与 MAT 文件在 `results/minimal_example/`。
- 新的两母波、四时刻单向采样通过：t=0 六分量逐元素完全一致；相位旋转等价检查最大相对 L2/Linf 均为 0；显式线性 eta/psi 检查最大相对 L2=9.772e-16；沿 y 恒定、有限值检查通过。Code Analyzer 无消息。数值/输入记录在 `results/unidirectional_time_series/known_spectrum_baseline.{mat,json}`。
- MATLAB 从默认路径启动并显式加载新仓库，无 SWORD 路径；运行 src/tests/examples 静态路径扫描仅发现测试禁止词及可选 MF12 占位路径。现有发布实现未改变。基线可独立脱离 SWORD 运行。
- 本地 Wolfram 14.0.0 可启动。远程 Linux 6.8.0-138-generic、g++ 13.3.0、CMake 3.28.3、Wolfram 15.0.1；Wolfram 简单求值返回 2，但报用户配置文件无法打开的警告。FFTW 头文件、Python3、patch、Git 可见；未编译或验证 FFTW 链接。
- 远程 MATLAB 未在 PATH 或上述常见目录找到，不能据此断言未安装。本机 g++/clang++ 未在 PATH 找到。未安装任何工具。
- MATLAB 重置默认路径时出现 Perl locale 回退警告，两个 batch 均 exit 0。
- 未重新运行 Wolfram 符号冻结/精确认证（只读取冻结接口并启动内核），未运行外部 MF12、完整四阶 GL–WIT 复现或 C++ 编译/计时；无新符号推导、无大规模计算。发布检查中的四阶 smoke 不等于完整四阶认证。

## 下一步

先确定正式输入类别；随后以同一母波合成 eta1(t)，核查时间谱归一化、相位、色散映射及非空间网格对齐问题，形成可验证的 MATLAB 输入契约。详见研究 README；当前采样例不宣称已完成直接时间序列算法或逆问题。

唯一待回答的问题：输入是已知/已分离的一阶时间序列 eta1(t)，还是测得的总波面 eta(t)？
