# JPsiPsi2s

J/ψ + ψ(2S) 关联产生四缪子分析代码。

## 常用命令

以下命令默认在 **`Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/`** 目录下执行。

### MC 模板分布对比

```bash
cd Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/

# 扫描 SPS/ 下所有子目录，每个子目录产生 5 张对比图（SPS vs DPS）
python3 plot_mc_dists.py

# 只跑指定 SPS 目录
python3 plot_mc_dists.py --sps-dir=ggpsi1psi2

# 使用 RECO 级变量（默认为 GEN 级）
python3 plot_mc_dists.py --reco
```

输入：`/eos/home-l/leyao/26JP/Ntuple/SPS/<子目录>/*.root`（同目录多个文件自动 TChain 合并）
DPS 输入固定为 `v2_set1.root`。
输出：`fig/mc_dists/<子目录名>/`，每子目录 5 张图（`delta_y / delta_phi / evt_mass / evt_y / evt_pt .png`）。

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
cmssw-el7                               # 切换到 EL7 + CMSSW 环境
cmsenv                                  # 在 CMSSW_*/src 下设置 ROOT/RooFit
```

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
3. Closure + 混合权重: `count.cpp` 读 acceptance/efficiency 文本表
4. Ntuple 产生: `MuMuGammaRootupler.cc` 或 `NtupleMaker.cc` → `rootuple/oniaTree`
5. 加权: `rephrase.cpp` → `WeightData*.root`
6. 4D 拟合: `Fit_4D_tot.cpp`（总截面）/ `Fit_4D_diff.cpp`（微分截面）
7. 模板拟合: `Tpl_Fit.cpp` / `Tplfit.cpp` → f_DPS → σ_SPS, σ_DPS

### 目录结构

- `GEN_nofilter/{SPS,DPS}/.../4mu_acc_eff/` — GEN acceptance（读 `GenAnalyzer/gen_tree`）
- `SKIM_tightfilter/{SPS,DPS}/.../4mu_acc_eff/` — 重建效率
- `SKIM_tightfilter/FourMuonAna/Onia/` — rootupler + 画图宏
- `Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/` — **下游拟合/截面提取主目录**
- `Data/NtupleMaker_repo/` — data ntuple maker 独立 repo

### 代码修改记录

v2 兼容性适配（2026-07-13）：`efficiencyPlot.py` ×2 已替换为单事件循环版本，兼容 v2 ntuple（RE 向量）
和旧 ntuple（GE 标量），自动检测格式。修复 PyROOT UChar_t 陷阱。详见 `docs/pyroot_uchar_trap.md`。

完整代码修改计划见 `docs/code_update_plan.md`，样本需求见 `docs/sample_requirements.md`。

### 文档

- `AGENTS.md` — 项目长期规则
- `docs/handoff.md` — 任务状态与下一步
- `docs/code_update_plan.md` — 代码修改计划
- `docs/sample_requirements.md` — MC 样本需求
- `docs/pyroot_uchar_trap.md` — PyROOT UChar_t 陷阱
- `docs/READ.md` — 历史 Q&A 记录
- `../docs/shared/` — 代码端 ↔ AN 端协作文档
