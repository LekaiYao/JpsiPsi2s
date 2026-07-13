#!/usr/bin/env python3
"""
Test script: run efficiency calculation on v2 ntuple (single event loop).
Compares performance: single event loop vs original GetEntries-per-bin approach.

Usage:
    python3 test_v2_efficiency.py             # event loop mode (v2 compatible)
    python3 test_v2_efficiency.py --old       # GetEntries mode (Steps 1-3 only, timing ref)
"""

import ROOT
from array import array
import os, sys, time
from math import sqrt
import numpy as np

# ── Config ──────────────────────────────────────────────────
INFILE  = "/eos/home-l/leyao/26JP/Ntuple/DPS/filter3p5/v2_test1_nt2.root"
TREE_NAME = "rootuple/oniaTree"
OUTDIR = "plot_v2_test/"
USE_OLD_METHOD = "--old" in sys.argv

ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin  = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])
NPT = len(ptBin) - 1   # 19
NY  = len(yBin) - 1    # 10

if not os.path.exists(OUTDIR):
    os.system("mkdir -p " + OUTDIR)

# ── Open file ───────────────────────────────────────────────
f = ROOT.TFile(INFILE, "READ")
t = f.Get(TREE_NAME)
nEntries = t.GetEntries()
print(f"Tree: {TREE_NAME}, entries: {nEntries}")

# Enable branches we need
t.SetBranchStatus("*", 0)
for b in ["GEevt_valid", "GEevt_passAcc", "GEevt_fourMuMass",
          "GEJpsi_pt", "GEJpsi_y", "GEJpsi_matchGEN", "GEJpsi_passID", "GEJpsi_passVtx",
          "GEpsi2S_pt", "GEpsi2S_y", "GEpsi2S_matchGEN", "GEpsi2S_passID", "GEpsi2S_passVtx",
          "GEevt_passHLT"]:
    t.SetBranchStatus(b, 1)
for b in ["REevt_fourMuMass", "REevt_matchTrg", "REmu_pt"]:
    t.SetBranchStatus(b, 1)

# ── Histograms & matrices ───────────────────────────────────
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

# ── Helper: safe UChar_t conversion (PyROOT returns str '\x00'/'\\x01') ──
def uchar_bool(val):
    """Convert PyROOT UChar_t value to bool. Handles str/bytes/int."""
    if isinstance(val, str):
        return ord(val) != 0 if len(val) == 1 else bool(val)
    return bool(val)

# ── Helper: find bin index ──────────────────────────────────
def find_pt_bin(pt_val):
    if pt_val <= ptBin[0] or pt_val >= ptBin[-1]:
        return -1
    for i in range(NPT):
        if ptBin[i] <= pt_val < ptBin[i+1]:
            return i
    return -1

def find_y_bin(y_val):
    if y_val <= yBin[0] or y_val >= yBin[-1]:
        return -1
    for i in range(NY):
        if yBin[i] <= y_val < yBin[i+1]:
            return i
    return -1

# ── Timing ──────────────────────────────────────────────────
t_start = time.time()

if USE_OLD_METHOD:
    print("\n=== Using GetEntries method (Steps 1-3 only, for timing) ===\n")
    entryCut = "GEevt_valid && GEevt_passAcc && GEevt_fourMuMass > 7.5 && "
    fidJpsiCut = f"GEJpsi_pt > {ptBin[0]} && GEJpsi_pt < {ptBin[-1]} && GEJpsi_y > {yBin[0]} && GEJpsi_y < {yBin[-1]} && "
    fidpsi2SCut = f"GEpsi2S_pt > {ptBin[0]} && GEpsi2S_pt < {ptBin[-1]} && GEpsi2S_y > {yBin[0]} && GEpsi2S_y < {yBin[-1]} && "

    t_old_start = time.time()
    n_calls = 0
    for ptIdx in range(NPT):
        for yIdy in range(NY):
            bcJ = entryCut + fidpsi2SCut + f"GEJpsi_pt > {ptBin[ptIdx]} && GEJpsi_pt < {ptBin[ptIdx+1]} && GEJpsi_y > {yBin[yIdy]} && GEJpsi_y < {yBin[yIdy+1]}"
            bcP = entryCut + fidJpsiCut + f"GEpsi2S_pt > {ptBin[ptIdx]} && GEpsi2S_pt < {ptBin[ptIdx+1]} && GEpsi2S_y > {yBin[yIdy]} && GEpsi2S_y < {yBin[yIdy+1]}"
            mBin_Jpsi[yIdy][ptIdx]  = int(t.GetEntries(bcJ))
            mReco_Jpsi[yIdy][ptIdx] = int(t.GetEntries(bcJ + " && GEJpsi_matchGEN"))
            mId_Jpsi[yIdy][ptIdx]   = int(t.GetEntries(bcJ + " && GEJpsi_passID"))
            mVtx_Jpsi[yIdy][ptIdx]  = int(t.GetEntries(bcJ + " && GEJpsi_passVtx"))
            mBin_psi2S[yIdy][ptIdx]  = int(t.GetEntries(bcP))
            mReco_psi2S[yIdy][ptIdx] = int(t.GetEntries(bcP + " && GEpsi2S_matchGEN"))
            mId_psi2S[yIdy][ptIdx]   = int(t.GetEntries(bcP + " && GEpsi2S_passID"))
            mVtx_psi2S[yIdy][ptIdx]  = int(t.GetEntries(bcP + " && GEpsi2S_passVtx"))
            n_calls += 8
    t_old = time.time() - t_old_start
    print(f"GetEntries method: {n_calls} calls in {t_old:.2f}s "
          f"(~{n_calls/t_old:.0f} calls/s, ~{t_old/n_calls*1000:.1f} ms/call)")

else:
    print("\n=== Using single event loop (v2 compatible) ===\n")

    t_loop_start = time.time()
    n_entry_pass = 0
    n_vtx_both = 0
    n_hlt = 0
    n_trg = 0

    for i in range(nEntries):
        t.GetEntry(i)

        # Entry cut (all GEevt_* flags are UChar_t → use uchar_bool)
        if not (uchar_bool(t.GEevt_valid) and uchar_bool(t.GEevt_passAcc) and t.GEevt_fourMuMass > 7.5):
            continue
        n_entry_pass += 1

        jpsi_pt  = t.GEJpsi_pt
        jpsi_y   = t.GEJpsi_y
        psi2S_pt = t.GEpsi2S_pt
        psi2S_y  = t.GEpsi2S_y

        # ── Steps 1-3: single dimuon efficiency (all GE flags are UChar_t) ──
        # Per original code: J/psi denominator requires psi(2S) in global fiducial, and vice versa
        iPt_j = find_pt_bin(jpsi_pt)
        iY_j  = find_y_bin(jpsi_y)
        iPt_p = find_pt_bin(psi2S_pt)
        iY_p  = find_y_bin(psi2S_y)

        psi2S_in_fid = (iPt_p >= 0 and iY_p >= 0)
        jpsi_in_fid  = (iPt_j >= 0 and iY_j >= 0)

        if iPt_j >= 0 and iY_j >= 0 and psi2S_in_fid:
            mBin_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(t.GEJpsi_matchGEN):
                mReco_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(t.GEJpsi_passID):
                mId_Jpsi[iY_j][iPt_j] += 1
            if uchar_bool(t.GEJpsi_passVtx):
                mVtx_Jpsi[iY_j][iPt_j] += 1

        if iPt_p >= 0 and iY_p >= 0 and jpsi_in_fid:
            mBin_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(t.GEpsi2S_matchGEN):
                mReco_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(t.GEpsi2S_passID):
                mId_psi2S[iY_p][iPt_p] += 1
            if uchar_bool(t.GEpsi2S_passVtx):
                mVtx_psi2S[iY_p][iPt_p] += 1

        # ── Steps 4-5: event-level efficiency ──
        if iPt_j < 0 or iY_j < 0 or iPt_p < 0 or iY_p < 0:
            continue
        if not (uchar_bool(t.GEJpsi_passVtx) and uchar_bool(t.GEpsi2S_passVtx)):
            continue

        mVtx_evt[iPt_p][iPt_j] += 1
        n_vtx_both += 1

        # HLT: GEevt_passHLT is UChar_t
        hlt_ok = uchar_bool(t.GEevt_passHLT)
        if hlt_ok:
            mHlt_evt[iPt_p][iPt_j] += 1
            n_hlt += 1

        # Trigger match: RE vectors (v2 only)
        trg_ok = False
        n_combo = t.REevt_fourMuMass.size()
        for j in range(n_combo):
            if t.REevt_matchTrg[j]:
                mu_pt_max = max(t.REmu_pt) if t.REmu_pt.size() > 0 else 0
                if mu_pt_max > 5.5:
                    trg_ok = True
                break
        if trg_ok:
            mTrg_evt[iPt_p][iPt_j] += 1
            n_trg += 1

    t_loop = time.time() - t_loop_start
    print(f"Event loop: {nEntries} events in {t_loop:.2f}s "
          f"({nEntries/t_loop:.0f} evt/s)")
    print(f"  Entry pass: {n_entry_pass}")
    print(f"  Both Vtx:   {n_vtx_both}")
    print(f"  HLT pass:   {n_hlt}")
    print(f"  Trg match:  {n_trg}")

# ── Print summary ───────────────────────────────────────────
print(f"\n{'='*70}")
print(f"Efficiency Summary ({'GetEntries' if USE_OLD_METHOD else 'event loop'})")
print(f"{'='*70}")

total_j = np.sum(mBin_Jpsi)
total_p = np.sum(mBin_psi2S)
print(f"\n--- Steps 1-3 (single dimuon) ---")
print(f"  J/psi:   total={total_j:5d}")
print(f"  psi(2S): total={total_p:5d}")

if not USE_OLD_METHOD:
    total_evt = np.sum(mVtx_evt)
    total_hlt = np.sum(mHlt_evt)
    total_trg = np.sum(mTrg_evt)
    print(f"\n--- Steps 4-5 (event-level) ---")
    print(f"  nVtx (both dimuons): {total_evt}")
    print(f"  nHlt:                {total_hlt}  (eff = {total_hlt/max(1,total_evt)*100:.1f}%)")
    print(f"  nTrg:                {total_trg}  (eff = {total_trg/max(1,total_hlt)*100:.1f}%)")

# ── Per-bin spot check ──────────────────────────────────────
print(f"\n--- Per-bin spot check (first 3x3 J/psi bins) ---")
for iy in range(min(3, NY)):
    for ip in range(min(3, NPT)):
        b = mBin_Jpsi[iy][ip]
        r = mReco_Jpsi[iy][ip]
        d = mId_Jpsi[iy][ip]
        v = mVtx_Jpsi[iy][ip]
        print(f"  (pT={ptBin[ip]:.0f}-{ptBin[ip+1]:.0f}, y={yBin[iy]:.1f}-{yBin[iy+1]:.1f}): "
              f"tot={b:4d} rec={r:4d} id={d:4d} vtx={v:4d}")

# ── Performance analysis ────────────────────────────────────
if not USE_OLD_METHOD:
    n_bins_steps123 = NPT * NY   # 190
    n_bins_steps45  = NPT * NPT  # 361
    n_calls_old = n_bins_steps123 * 8 + n_bins_steps45 * 3  # 2603

    t_g = time.time()
    t.GetEntries("GEevt_valid")
    t_one_getentries = time.time() - t_g

    print(f"\n{'='*70}")
    print(f"Performance Analysis")
    print(f"{'='*70}")
    print(f"  Single GetEntries (2494 evt): {t_one_getentries*1000:.1f} ms")
    print(f"  Old method total calls:       {n_calls_old}")
    print(f"  Old method estimated:         {n_calls_old * t_one_getentries:.1f}s")
    print(f"  Event loop actual:            {t_loop:.2f}s")
    if t_one_getentries > 0 and t_loop > 0:
        speedup = n_calls_old * t_one_getentries / t_loop
        print(f"  Speedup (event loop):         ~{speedup:.0f}x")

    # Project to production scale
    print(f"\n  --- Projected to production scale ---")
    for scale, label in [(10, "v2_set1 (227k evt)"),
                          (45, "1 old file (1M evt)"),
                          (1800, "40 old files (40M evt)")]:
        est_old = n_calls_old * t_one_getentries * scale
        est_new = t_loop * scale
        print(f"  {label:30s}: old={est_old/60:6.1f}min  new={est_new:6.1f}s  "
              f"ratio={est_old/max(est_new,0.001):.0f}x")

# ── Write efficiency.txt ────────────────────────────────────
if not USE_OLD_METHOD:
    with open(OUTDIR + "efficiency.txt", 'w') as fout:
        n, m = NPT, NY
        n_Jpsi = np.sum(mBin_Jpsi)
        n_psi2S = np.sum(mBin_psi2S)
        fout.write(f"{n} {m} {n_Jpsi} {n_psi2S}\n")
        fout.write(' '.join(str(x) for x in ptBin) + '\n')
        fout.write(' '.join(str(x) for x in yBin) + '\n')
        for i in range(m):
            for j in range(n):
                fout.write(f"{mVtx_Jpsi[i][j]} {mBin_Jpsi[i][j]} ")
            fout.write('\n')
        for i in range(m):
            for j in range(n):
                fout.write(f"{mVtx_psi2S[i][j]} {mBin_psi2S[i][j]} ")
            fout.write('\n')
        for i in range(n):
            for j in range(n):
                fout.write(f"{mTrg_evt[i][j]} {mVtx_evt[i][j]} ")
            fout.write('\n')
    print(f"\nWrote {OUTDIR}efficiency.txt")

f.Close()
print(f"\nTotal runtime: {time.time() - t_start:.1f}s")
