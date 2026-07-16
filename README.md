# JPsiPsi2s

J/ψ + ψ(2S) 关联产生四缪子分析代码。

## 常用命令

以下命令默认在 **`Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/`** 目录下执行，
除非另有注明。

### MC 模板分布对比

```bash
cd Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/

# RECO 级（默认）：读 raw NtupleMaker ntuple
python3 plot_mc_dists.py

# GEN 级
python3 plot_mc_dists.py --gen

# 模板拟合变量（--eff）：读 UnweightMC/ 下 rephrase MC 模式产出
python3 plot_mc_dists.py --eff

# 只跑指定 SPS 目录
python3 plot_mc_dists.py --sps-dir=ggpsi1psi2g_gpt0p8
```

输入：`/eos/home-l/leyao/26JP/Ntuple/SPS/<子目录>/*.root`（同目录多个文件自动 TChain 合并）。
DPS 输入：`/eos/home-l/leyao/26JP/Ntuple/DPS/filter3p5/v2_set*.root`。
`--eff` 模式输入：`UnweightMC/Unweight_*.root` + `UnweightMC/UnweightDPS.root`。
输出：`fig/mc_dists/<子目录>/`（`--eff` 模式 → `fig/mc_dists/eff/<channel>/`）。

### MC 模板产出（Unweight 文件）

编辑 `rephrase.cpp`，取消 `#define MC_MODE` 行的注释，设置 `MC_CHANNEL`：

```cpp
#define MC_MODE
string MC_CHANNEL = "DPS";                   // → UnweightMC/UnweightDPS.root
string MC_CHANNEL = "ggpsi1psi2g_gpt0p8";   // → UnweightMC/Unweight_ggpsi1psi2g_gpt0p8.root
```

然后运行（需 CMSSW 环境）：

```bash
cd /eos/home-l/leyao/26JP/crab_MC_Maker/CMSSW_10_6_20/src
cmssw-el7 && cmsenv
cd /eos/home-l/leyao/26JP/JPsiPsi2s/Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer
root -l -q -b rephrase.cpp
```

每次只跑一个 channel。恢复 Data 模式：重新注释 `#define MC_MODE`。

### 效率计算

```bash
# DPS 效率图（事件循环版本，兼容 v2 ntuple）
cd SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/
python3 efficiencyPlot.py              # 输出 plot/efficiency.txt + plot/*.png

# SPS 效率图
cd SKIM_tightfilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/
python3 efficiencyPlot.py

# 原版 GetEntries 备份
python3 efficiencyPlot_GetEntries.py    # 仅用于参考/对比
```

### 效率测试与性能对比

```bash
cd SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/
python3 test_v2_efficiency.py           # 事件循环（v2 兼容）
python3 test_v2_efficiency.py --old     # 原版 GetEntries（仅 Steps 1-3，用于计时对比）
```

### 上游同步检查

```bash
bash update/quick_sync.sh               # 快速检查 chensh 上游变更（~0.5s）
bash update/quick_sync.sh --diff        # 查看详细 diff
```

### 环境

```bash
# Python 脚本（EL9 直接跑）
python3 <script>.py

# ROOT 宏（需 CMSSW 环境，在 crab_MC_Maker 下 cmsenv）
cd /eos/home-l/leyao/26JP/crab_MC_Maker/CMSSW_10_6_20/src
cmssw-el7 && cmsenv
cd /eos/home-l/leyao/26JP/JPsiPsi2s/Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer
root -l -q -b <macro>.cpp
```

> **注意**：`Data/ULntuple18/CMSSW_10_6_20/` 不是 cmsrel 构造的 CMSSW 环境，不能在此做 `cmsenv`。

大多数脚本顶部为硬编码路径，运行前需确认输入 ROOT 文件路径是否存在。

---

## 项目简介

本项目是对 chensh 完整生产+分析工作区（`/eos/user/c/chensh/JPsiPsi2s`）做的**轻量级代码提取**：
仅保留 GEN acceptance、reco efficiency 分析宏、ntuple rootupler 及下游画图/拟合宏，不含 ROOT 文件、
CMSSW 编译产物及 MC 生产链（GEN/SIM/DIGI/RECO/HLT 配置与 CRAB 提交）。

分析末态为 J/ψ → μ⁺μ⁻ + ψ(2S) → μ⁺μ⁻，数据为 2018 UL（61.31 fb⁻¹）。

### 主分析链路

1. GEN acceptance: `GEN_nofilter/{SPS,DPS}/.../efficiencyPlot.py` → `plot/acceptance.txt`
2. 效率计算: `SKIM_tightfilter/{SPS,DPS}/.../efficiencyPlot.py` → `plot/efficiency.txt`
3. 混合效率: `merge.cpp` → `efficiency_SPS+DPS_6_4.txt`（SPS:DPS = 0.6:0.4）
4. 加权 (data): `rephrase.cpp` (Data 模式) → `WeightData.root`
5. MC 模板: `rephrase.cpp` (MC 模式) → `UnweightMC/Unweight_*.root`
6. 4D 拟合: `Fit_4D_tot.cpp` → `Model_4D_tot.root` + `fit_results_4D_tot.txt`
7. 微分截面: `Fit_4D_diff.cpp` → `fit_results_4D_diff.txt`
8. 模板拟合: `Tpl_Fit.cpp` → `fit_results_Tpl_Fit.txt`（f_DPS, σ_SPS, σ_DPS）

### 目录结构

- `GEN_nofilter/{SPS,DPS}/.../4mu_acc_eff/` — GEN acceptance（读 `GenAnalyzer/gen_tree`）
- `SKIM_tightfilter/{SPS,DPS}/.../4mu_acc_eff/` — 重建效率
- `SKIM_tightfilter/FourMuonAna/Onia/` — rootupler + 画图宏
- `Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/` — **下游拟合/截面提取主目录**
- `Data/NtupleMaker_repo/` — data ntuple maker 独立 repo

### 文档

- `AGENTS.md` — 项目长期规则
- `docs/analysis_project.md` — 分析蓝图
- `docs/selection_diff.md` — JP vs JJ selection 对齐
- `docs/sample_requirements.md` — MC 样本需求
- `docs/code_changelog.md` — 代码改动记录
- `docs/READ.md` — 历史 Q&A 记录
- `../docs/shared/` — 代码端 ↔ AN 端协作文档
