# 初始化交接 — 2026-09-25

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
