#!/usr/bin/env python3
"""
SKIM efficiency map: DPS (single event loop, v2 ntuple compatible).

Replaces the original GetEntries-per-bin approach (2603 full-tree scans) with a
single event loop.  Compatible with both v2 NtupleMaker (RE vectors for Steps 4-5)
and old NtupleMaker (GE scalars).

Performance: ~1300x faster than GetEntries at production scale.
Reference: docs/READ.md (2026-07-13), test_v2_efficiency.py in same directory.
"""

import ROOT
from array import array
import os
from math import sqrt
import numpy as np

# ── matplotlib / PyROOT compatibility ───────────────────────
# ROOT 6.40+ disallows 'from ROOT import *'; import ROOT and prefix.
# The local tdrStyle.py may use the old pattern — import it conditionally.
try:
    from tdrStyle import *       # old-style local copy (may still work)
except (ImportError, SyntaxError):
    pass                         # setTDRStyle() defined locally below

# ── Config ──────────────────────────────────────────────────
ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin  = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])
NPT = len(ptBin) - 1   # 19
NY  = len(yBin) - 1    # 10

outPlotDir = "plot/"
if not os.path.exists(outPlotDir):
    os.system("mkdir -p " + outPlotDir)

# Input files (hardcoded paths — adjust for your environment)
inFileDir = "pTHat4/"
inFileList = ["Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 41)]
# + ["2/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]
# + ["3/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]
# + ["4/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]

# ── PyROOT UChar_t helper ───────────────────────────────────
# ROOT ≥ 6.30 on EL9 returns UChar_t values as Python str('\x00')/'\\x01'.
# bool('\x00') is True (non-empty string), so we must convert explicitly.
def uchar_bool(val):
    """Convert PyROOT UChar_t to bool safely across ROOT versions."""
    if val is None:
        return False
    if isinstance(val, str):
        return ord(val[0]) != 0
    if isinstance(val, bytes):
        return val[0] != 0
    return bool(val)

# ── Bin finding ─────────────────────────────────────────────
def find_pt_bin(pt_val):
    if pt_val <= ptBin[0] or pt_val >= ptBin[-1]:
        return -1
    for i in range(NPT):
        if ptBin[i] <= pt_val < ptBin[i + 1]:
            return i
    return -1

def find_y_bin(y_val):
    if y_val <= yBin[0] or y_val >= yBin[-1]:
        return -1
    for i in range(NY):
        if yBin[i] <= y_val < yBin[i + 1]:
            return i
    return -1

# ── Open files & detect format ──────────────────────────────
inFiles = {}
inTrees = []
for inFileName in inFileList:
    path = inFileDir + inFileName
    if not os.path.exists(path):
        print(f"Warning: {path} not found, skipping")
        continue
    f = ROOT.TFile(path, "READ")
    inFiles[inFileName] = f
    t = f.Get("rootuple/oniaTree")
    if not t:
        print(f"Warning: no rootuple/oniaTree in {path}, skipping")
        continue
    inTrees.append(t)
    print(f"  {path}  ({t.GetEntries()} entries)")

if not inTrees:
    raise RuntimeError("No input files found!")

# Detect ntuple format: old (has GEevt_matchTrg) vs v2 (RE vectors only)
has_GE_matchTrg = bool(inTrees[0].GetListOfBranches().FindObject("GEevt_matchTrg"))
has_GE_muPtMax  = bool(inTrees[0].GetListOfBranches().FindObject("GEevt_muPtMax"))
is_v2 = not (has_GE_matchTrg and has_GE_muPtMax)
print(f"Ntuple format: {'v2 (RE vectors)' if is_v2 else 'old (GE scalars)'}")

# Enable branches (GE: Steps 1-4, RE: Steps 4-5 for v2)
for t in inTrees:
    t.SetBranchStatus("*", 0)
    for b in ["GEevt_valid", "GEevt_passAcc", "GEevt_fourMuMass",
              "GEJpsi_pt", "GEJpsi_y", "GEJpsi_matchGEN", "GEJpsi_passID", "GEJpsi_passVtx",
              "GEpsi2S_pt", "GEpsi2S_y", "GEpsi2S_matchGEN", "GEpsi2S_passID", "GEpsi2S_passVtx",
              "GEevt_passHLT"]:
        t.SetBranchStatus(b, 1)
    if is_v2:
        for b in ["REevt_fourMuMass", "REevt_matchTrg", "REmu_pt"]:
            t.SetBranchStatus(b, 1)
    else:
        for b in ["GEevt_matchTrg", "GEevt_muPtMax"]:
            t.SetBranchStatus(b, 1)

# ── Accumulator matrices (numpy, filled in event loop) ──────
mBin_Jpsi   = np.zeros((NY, NPT), dtype=int)
mReco_Jpsi  = np.zeros((NY, NPT), dtype=int)
mId_Jpsi    = np.zeros((NY, NPT), dtype=int)
mVtx_Jpsi   = np.zeros((NY, NPT), dtype=int)
mBin_psi2S  = np.zeros((NY, NPT), dtype=int)
mReco_psi2S = np.zeros((NY, NPT), dtype=int)
mId_psi2S   = np.zeros((NY, NPT), dtype=int)
mVtx_psi2S  = np.zeros((NY, NPT), dtype=int)

mVtx_evt = np.zeros((NPT, NPT), dtype=int)
mHlt_evt = np.zeros((NPT, NPT), dtype=int)
mTrg_evt = np.zeros((NPT, NPT), dtype=int)

# ── Single event loop ───────────────────────────────────────
print("Running event loop ...")
n_total = 0
for inTree in inTrees:
    nEvt = inTree.GetEntries()
    for i in range(nEvt):
        inTree.GetEntry(i)
        n_total += 1

        # Entry cut (all GE flags are UChar_t)
        if not (uchar_bool(inTree.GEevt_valid) and
                uchar_bool(inTree.GEevt_passAcc) and
                inTree.GEevt_fourMuMass > 7.5):
            continue

        jpsi_pt  = inTree.GEJpsi_pt
        jpsi_y   = inTree.GEJpsi_y
        psi2S_pt = inTree.GEpsi2S_pt
        psi2S_y  = inTree.GEpsi2S_y

        iPt_j = find_pt_bin(jpsi_pt)
        iY_j  = find_y_bin(jpsi_y)
        iPt_p = find_pt_bin(psi2S_pt)
        iY_p  = find_y_bin(psi2S_y)

        psi2S_in_fid = (iPt_p >= 0 and iY_p >= 0)
        jpsi_in_fid  = (iPt_j >= 0 and iY_j >= 0)

        # ── Steps 1-3: single dimuon efficiency ──
        if iPt_j >= 0 and iY_j >= 0 and psi2S_in_fid:
            mBin_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(inTree.GEJpsi_matchGEN):
                mReco_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(inTree.GEJpsi_passID):
                mId_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(inTree.GEJpsi_passVtx):
                mVtx_Jpsi[iY_j][iPt_j] += 1

        if iPt_p >= 0 and iY_p >= 0 and jpsi_in_fid:
            mBin_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(inTree.GEpsi2S_matchGEN):
                mReco_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(inTree.GEpsi2S_passID):
                mId_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(inTree.GEpsi2S_passVtx):
                mVtx_psi2S[iY_p][iPt_p] += 1

        # ── Steps 4-5: event-level efficiency ──
        if iPt_j < 0 or iY_j < 0 or iPt_p < 0 or iY_p < 0:
            continue
        if not (uchar_bool(inTree.GEJpsi_passVtx) and
                uchar_bool(inTree.GEpsi2S_passVtx)):
            continue

        mVtx_evt[iPt_p][iPt_j] += 1

        hlt_ok = uchar_bool(inTree.GEevt_passHLT)
        if hlt_ok:
            mHlt_evt[iPt_p][iPt_j] += 1

        # Trigger match: RE vectors (v2) or GE scalars (old)
        trg_ok = False
        if is_v2:
            n_combo = inTree.REevt_fourMuMass.size()
            for j in range(n_combo):
                if inTree.REevt_matchTrg[j]:
                    mu_pt_max = max(inTree.REmu_pt) if inTree.REmu_pt.size() > 0 else 0
                    if mu_pt_max > 5.5:
                        trg_ok = True
                    break
        else:
            trg_ok = (uchar_bool(inTree.GEevt_matchTrg) and
                      inTree.GEevt_muPtMax > 5.5)
        if trg_ok:
            mTrg_evt[iPt_p][iPt_j] += 1

print(f"  Processed {n_total} events total")

# ── Fill ROOT histograms from matrices ──────────────────────
# Same histogram names as original for backward compatibility.
hReco_Jpsi  = ROOT.TH2D("hReco_Jpsi",  "hReco_Jpsi",  NPT, ptBin, NY, yBin)
hReco_psi2S = ROOT.TH2D("hReco_psi2S", "hReco_psi2S", NPT, ptBin, NY, yBin)
hId_Jpsi    = ROOT.TH2D("hId_Jpsi",    "hId_Jpsi",    NPT, ptBin, NY, yBin)
hId_psi2S   = ROOT.TH2D("hId_psi2S",   "hId_psi2S",   NPT, ptBin, NY, yBin)
hVtx_Jpsi   = ROOT.TH2D("hVtx_Jpsi",   "hVtx_Jpsi",   NPT, ptBin, NY, yBin)
hVtx_psi2S  = ROOT.TH2D("hVtx_psi2S",  "hVtx_psi2S",  NPT, ptBin, NY, yBin)
dReco_Jpsi  = ROOT.TH2D("dReco_Jpsi",  "dReco_Jpsi",  NPT, ptBin, NY, yBin)
dReco_psi2S = ROOT.TH2D("dReco_psi2S", "dReco_psi2S", NPT, ptBin, NY, yBin)
dId_Jpsi    = ROOT.TH2D("dId_Jpsi",    "dId_Jpsi",    NPT, ptBin, NY, yBin)
dId_psi2S   = ROOT.TH2D("dId_psi2S",   "dId_psi2S",   NPT, ptBin, NY, yBin)
dVtx_Jpsi   = ROOT.TH2D("dVtx_Jpsi",   "dVtx_Jpsi",   NPT, ptBin, NY, yBin)
dVtx_psi2S  = ROOT.TH2D("dVtx_psi2S",  "dVtx_psi2S",  NPT, ptBin, NY, yBin)
hBin_Jpsi   = ROOT.TH2D("hBin_Jpsi",   "hBin_Jpsi",   NPT, ptBin, NY, yBin)
hBin_psi2S  = ROOT.TH2D("hBin_psi2S",  "hBin_psi2S",  NPT, ptBin, NY, yBin)
hHlt        = ROOT.TH2D("hHlt",        "hHlt",        NPT, ptBin, NPT, ptBin)
hTrg        = ROOT.TH2D("hTrg",        "hTrg",        NPT, ptBin, NPT, ptBin)
dHlt        = ROOT.TH2D("dHlt",        "dHlt",        NPT, ptBin, NPT, ptBin)
dTrg        = ROOT.TH2D("dTrg",        "dTrg",        NPT, ptBin, NPT, ptBin)
hBin_evt    = ROOT.TH2D("hBin_evt",    "hBin_evt",    NPT, ptBin, NPT, ptBin)

nBin_Jpsi_max, nBin_psi2S_max, nBin_evt_max = 0, 0, 0

# ── Fill per-meson histograms (Steps 1-3) ───────────────────
for ptIdx in range(NPT):
    for yIdy in range(NY):
        nBin_Jpsi  = int(mBin_Jpsi[yIdy][ptIdx])
        nReco_Jpsi = int(mReco_Jpsi[yIdy][ptIdx])
        nId_Jpsi   = int(mId_Jpsi[yIdy][ptIdx])
        nVtx_Jpsi  = int(mVtx_Jpsi[yIdy][ptIdx])
        nBin_psi2S  = int(mBin_psi2S[yIdy][ptIdx])
        nReco_psi2S = int(mReco_psi2S[yIdy][ptIdx])
        nId_psi2S   = int(mId_psi2S[yIdy][ptIdx])
        nVtx_psi2S  = int(mVtx_psi2S[yIdy][ptIdx])

        print(f"\033[7mProcessing pt-y bin: {ptBin[ptIdx]:.0f}, {yBin[yIdy]:.1f}\033[0m")
        print(f"[\033[32mJ/psi\033[0m] Vtx/Id/Reco/Total = {nVtx_Jpsi}/{nId_Jpsi}/{nReco_Jpsi}/{nBin_Jpsi}")
        print(f"[\033[31mpsi(2S)\033[0m] Vtx/Id/Reco/Total = {nVtx_psi2S}/{nId_psi2S}/{nReco_psi2S}/{nBin_psi2S}")

        ibin = hReco_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        hBin_Jpsi.SetBinContent(ibin, nBin_Jpsi)
        hBin_psi2S.SetBinContent(ibin, nBin_psi2S)
        if nBin_Jpsi > nBin_Jpsi_max:
            nBin_Jpsi_max = nBin_Jpsi
        if nBin_psi2S > nBin_psi2S_max:
            nBin_psi2S_max = nBin_psi2S

        if nBin_Jpsi:
            hReco_Jpsi.SetBinContent(ibin, nReco_Jpsi / nBin_Jpsi)
            dReco_Jpsi.SetBinContent(ibin, sqrt((nBin_Jpsi - nReco_Jpsi) * nReco_Jpsi / nBin_Jpsi) / nBin_Jpsi)
        if nReco_Jpsi:
            hId_Jpsi.SetBinContent(ibin, nId_Jpsi / nReco_Jpsi)
            dId_Jpsi.SetBinContent(ibin, sqrt((nReco_Jpsi - nId_Jpsi) * nId_Jpsi / nReco_Jpsi) / nReco_Jpsi)
        if nId_Jpsi:
            hVtx_Jpsi.SetBinContent(ibin, nVtx_Jpsi / nId_Jpsi)
            dVtx_Jpsi.SetBinContent(ibin, sqrt((nId_Jpsi - nVtx_Jpsi) * nVtx_Jpsi / nId_Jpsi) / nId_Jpsi)

        if nBin_psi2S:
            hReco_psi2S.SetBinContent(ibin, nReco_psi2S / nBin_psi2S)
            dReco_psi2S.SetBinContent(ibin, sqrt((nBin_psi2S - nReco_psi2S) * nReco_psi2S / nBin_psi2S) / nBin_psi2S)
        if nReco_psi2S:
            hId_psi2S.SetBinContent(ibin, nId_psi2S / nReco_psi2S)
            dId_psi2S.SetBinContent(ibin, sqrt((nReco_psi2S - nId_psi2S) * nId_psi2S / nReco_psi2S) / nReco_psi2S)
        if nId_psi2S:
            hVtx_psi2S.SetBinContent(ibin, nVtx_psi2S / nId_psi2S)
            dVtx_psi2S.SetBinContent(ibin, sqrt((nId_psi2S - nVtx_psi2S) * nVtx_psi2S / nId_psi2S) / nId_psi2S)

# ── Fill event-level histograms (Steps 4-5) ─────────────────
for ptIdy in range(NPT):
    for ptIdx in range(NPT):
        nVtx = int(mVtx_evt[ptIdy][ptIdx])
        nHlt = int(mHlt_evt[ptIdy][ptIdx])
        nTrg = int(mTrg_evt[ptIdy][ptIdx])

        print(f"\033[7mProcessing pt-pt bin: {ptBin[ptIdx]:.0f}, {ptBin[ptIdy]:.0f}\033[0m")

        ibin = hHlt.GetBin(ptIdx + 1, ptIdy + 1)
        hBin_evt.SetBinContent(ibin, nVtx)
        if nVtx > nBin_evt_max:
            nBin_evt_max = nVtx
        if nVtx:
            hHlt.SetBinContent(ibin, nHlt / nVtx)
            dHlt.SetBinContent(ibin, sqrt((nVtx - nHlt) * nHlt / nVtx) / nVtx)
        if nHlt:
            hTrg.SetBinContent(ibin, nTrg / nHlt)
            dTrg.SetBinContent(ibin, sqrt((nHlt - nTrg) * nTrg / nHlt) / nHlt)

# ── Plotting (same as original) ─────────────────────────────
def setTDRStyle():
    """Minimal TDR style if local tdrStyle.py import fails."""
    try:
        from tdrStyle import setTDRStyle as _set
        _set()
    except (ImportError, NameError):
        ROOT.gROOT.SetStyle("Plain")
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetPadTopMargin(0.10)
        ROOT.gStyle.SetPadRightMargin(0.20)

width, height = 3000, 1200
squarish = False
zRange = 1.

def plotSave(histo, titles, name):
    c = ROOT.TCanvas(name, name, width, height)
    c.cd()
    c.SetTopMargin(0.10)
    c.SetRightMargin(0.20)
    histo.GetXaxis().SetTitle(titles[0])
    histo.GetYaxis().SetTitle(titles[1])
    histo.GetZaxis().SetTitle(titles[2])
    histo.GetXaxis().SetTitleOffset(1.0)
    histo.GetYaxis().SetTitleOffset(0.8)
    histo.GetZaxis().SetTitleOffset(0.6)
    histo.GetZaxis().SetRangeUser(0., zRange)
    histo.Draw("colzTEXT" + ("45" if squarish else ""))
    histo.GetXaxis().SetNdivisions(505)
    if width == height:
        histo.GetYaxis().SetNdivisions(505)
    c.SaveAs(outPlotDir + name + ".png")
    del c

setTDRStyle()
ROOT.gStyle.SetPaintTextFormat("1.0f")

# Per-meson efficiency plots
plotSave(hReco_Jpsi,  ["pT_{J/#psi}", "rapidity_{J/#psi}", "#varepsilon_{RECO}^{J/#psi}"], "Eff_RECO_Jpsi")
plotSave(hId_Jpsi,    ["pT_{J/#psi}", "rapidity_{J/#psi}", "#varepsilon_{ID}^{J/#psi}"],   "Eff_ID_Jpsi")
plotSave(hVtx_Jpsi,   ["pT_{J/#psi}", "rapidity_{J/#psi}", "#varepsilon_{Dimuon}^{J/#psi}"], "Eff_Dimuon_Jpsi")
plotSave(hReco_psi2S, ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#varepsilon_{RECO}^{#psi(2S)}"], "Eff_RECO_psi2S")
plotSave(hId_psi2S,   ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#varepsilon_{ID}^{#psi(2S)}"],   "Eff_ID_psi2S")
plotSave(hVtx_psi2S,  ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#varepsilon_{Dimuon}^{#psi(2S)}"], "Eff_Dimuon_psi2S")
plotSave(dReco_Jpsi,  ["pT_{J/#psi}", "rapidity_{J/#psi}", "#delta#varepsilon_{RECO}^{J/#psi}"], "EffErr_RECO_Jpsi")
plotSave(dId_Jpsi,    ["pT_{J/#psi}", "rapidity_{J/#psi}", "#delta#varepsilon_{ID}^{J/#psi}"],   "EffErr_ID_Jpsi")
plotSave(dVtx_Jpsi,   ["pT_{J/#psi}", "rapidity_{J/#psi}", "#delta#varepsilon_{Dimuon}^{J/#psi}"], "EffErr_Dimuon_Jpsi")
plotSave(dReco_psi2S, ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#delta#varepsilon_{RECO}^{#psi(2S)}"], "EffErr_RECO_psi2S")
plotSave(dId_psi2S,   ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#delta#varepsilon_{ID}^{#psi(2S)}"],   "EffErr_ID_psi2S")
plotSave(dVtx_psi2S,  ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "#delta#varepsilon_{Dimuon}^{#psi(2S)}"], "EffErr_Dimuon_psi2S")

# Distribution plots
zRange = nBin_Jpsi_max
ROOT.gStyle.SetPaintTextFormat("1.0f")
plotSave(hBin_Jpsi,  ["pT_{J/#psi}", "rapidity_{J/#psi}", "N_{J/#psi}"],  "N_distri_Jpsi")
zRange = nBin_psi2S_max
plotSave(hBin_psi2S, ["pT_{#psi(2S)}", "rapidity_{#psi(2S)}", "N_{#psi(2S)}"], "N_distri_psi2S")

# Event-level plots (square)
squarish, height = True, 2500
zRange = nBin_evt_max
plotSave(hBin_evt, ["pT_{J/#psi}", "pT_{#psi(2S)}", "N_{event}"], "N_distri_dimuon_pT")
ROOT.gStyle.SetPaintTextFormat("1.2f")
zRange = 1.
plotSave(hHlt, ["pT_{J/#psi}", "pT_{#psi(2S)}", "#varepsilon_{HLT}"],     "Eff_HLT")
plotSave(dHlt, ["pT_{J/#psi}", "pT_{#psi(2S)}", "#delta#varepsilon_{HLT}"], "EffErr_HLT")
plotSave(hTrg, ["pT_{J/#psi}", "pT_{#psi(2S)}", "#varepsilon_{FourMu}"],  "Eff_FourMu")
plotSave(dTrg, ["pT_{J/#psi}", "pT_{#psi(2S)}", "#delta#varepsilon_{FourMu}"], "EffErr_FourMu")

# ── Write efficiency.txt (same format as original) ──────────
with open(outPlotDir + "efficiency.txt", 'w') as f:
    n_Jpsi  = int(np.sum(mBin_Jpsi))
    n_psi2S = int(np.sum(mBin_psi2S))
    f.write(f"{NPT} {NY} {n_Jpsi} {n_psi2S}\n")
    f.write(' '.join(str(x) for x in ptBin) + '\n')
    f.write(' '.join(str(x) for x in yBin) + '\n')
    # J/psi: Vtx, Total
    for i in range(NY):
        for j in range(NPT):
            f.write(f"{int(mVtx_Jpsi[i][j])} {int(mBin_Jpsi[i][j])} ")
        f.write('\n')
    # psi(2S): Vtx, Total
    for i in range(NY):
        for j in range(NPT):
            f.write(f"{int(mVtx_psi2S[i][j])} {int(mBin_psi2S[i][j])} ")
        f.write('\n')
    # Event: Trg, Vtx
    for i in range(NPT):
        for j in range(NPT):
            f.write(f"{int(mTrg_evt[i][j])} {int(mVtx_evt[i][j])} ")
        f.write('\n')

print(f"\nDone. Output: {outPlotDir}efficiency.txt + plots")
