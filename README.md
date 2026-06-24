# JPsiPsi2s

J/psi + psi(2S) associated-production four-muon analysis code. This is a
lightweight extraction of the analysis macros from the full production+analysis
workspace: it keeps the GEN-level acceptance and reconstruction-efficiency
macros, the ntuple rootupler, and the downstream plotting macros. The MC
production chain (GEN/SIM/DIGI/RECO/HLT configs and CRAB submission) and all
ROOT/build artifacts are intentionally excluded.

The final state is four muons from one prompt J/psi -> mu mu and one
psi(2S) -> mu mu, measured in 2018 UL data.

## Main Workflow

1. Compute GEN acceptance maps in `GEN_nofilter/{SPS,DPS}/.../4mu_acc_eff/`
   with `efficiencyPlot.py`, producing 2D (pT, y) maps and `plot/acceptance.txt`.
2. Compute reconstruction/ID/Vtx/HLT efficiency maps in
   `SKIM_tightfilter/{SPS,DPS}/.../4mu_acc_eff/` with `efficiencyPlot.py`,
   producing `plot/efficiency.txt`.
3. Run the acceptance/efficiency closure check with `count.cpp`, which reads the
   text maps and computes per-event weights.
4. Produce the analysis ntuple from data/MC with the rootupler
   (`FourMuonAna/Onia/src/MuMuGammaRootupler.cc` or `NtupleMaker/.../NtupleMaker.cc`).
5. Read the analysis tree and draw distributions with
   `FourMuonAna/Onia/test/makeAllPlots/readTree2018.C`.

## Directory Layout

- `GEN_nofilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/`
  - SPS GEN acceptance maps (reads `GenAnalyzer/gen_tree`).
- `GEN_nofilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/`
  - DPS GEN acceptance maps.
- `SKIM_tightfilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/`
  - SPS reconstruction efficiency maps.
- `SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/`
  - DPS reconstruction efficiency maps.
- `SKIM_tightfilter/FourMuonAna/Onia/`
  - Ntuple rootupler (`MuMuGammaRootupler.cc`) and downstream plotting macros
    (`test/makeAllPlots/`).
- `SKIM_tightfilter/SPS/ULPythia2018/CMSSW_10_6_20/src/NtupleMaker/NtupleMaker/`
  - Standalone ntuple maker (`NtupleMaker.cc`).

## Key Macros

`efficiencyPlot.py` bins events in (pT, y) for the J/psi and the psi(2S),
applies the muon acceptance/quality cuts, and fills 2D efficiency histograms and
their statistical errors; results are written as a flat text map under `plot/`.

`count.cpp` loads a text map (`plot/acceptance.txt` or `efficiency.txt`) and
computes per-candidate weights as the product of the J/psi and psi(2S) cell
weights, used for the acceptance/efficiency closure test.

`efficiencyFactors_JJ.py`, `plot2dsigeffs_JJ.py`, `LoadData_JJ.py`,
`sample_shortnames_JJ.py` and `tdrStyle.py` make up the shared plotting
framework for the 2D acceptance/efficiency maps.

`MuMuGammaRootupler.cc` is the EDAnalyzer that turns RECO/MiniAOD events into the
flat analysis ntuple (the four-muon `oniaTree`).

`readTree2018.C` / `plot_compare.C` read the analysis tree and draw comparison
distributions.

## Environment

The ROOT/RooFit macros and CMSSW analyzers are old-style and are intended for
the lxplus EL7/CMSSW environment:

```bash
cmssw-el7
cmsenv   # inside the relevant CMSSW_*/src
```

Most macros use hardcoded input paths; check the settings at the top of each
script before running.
