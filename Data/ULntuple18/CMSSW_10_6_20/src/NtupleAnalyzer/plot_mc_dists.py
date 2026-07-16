#!/usr/bin/env python3
"""
Scan SPS samples vs DPS: plot GEN-level differential distributions with ratio panel.

For each SPS subdirectory under SPS_BASE containing .root files, produce
5 overlaid + ratio plots (delta_y, delta_phi, evt_mass, evt_y, evt_pt).

Usage:
    python3 plot_mc_dists.py                     # scan all SPS subdirs
    python3 plot_mc_dists.py --sps-dir=ggpsi1psi2  # single SPS folder only
    python3 plot_mc_dists.py --gen               # GEN-level (opt-in)
    python3 plot_mc_dists.py --eff               # template-fit variables from UnweightMC/
Output: fig/mc_dists/<sps_folder_name>/
"""

import ROOT, os, sys, glob, numpy as np
from array import array

# ── Config ──────────────────────────────────────────────────
SPS_BASE  = "/eos/home-l/leyao/26JP/Ntuple/SPS"
DPS_FILE  = "/eos/home-l/leyao/26JP/Ntuple/DPS/filter3p5/v2_set1.root"
TREE_NAME = "rootuple/oniaTree"
BASE_OUT  = "fig/mc_dists/"
USE_RECO  = "--gen" not in sys.argv and "--eff" not in sys.argv   # default: RECO level
USE_GEN   = "--gen" in sys.argv
USE_EFF   = "--eff" in sys.argv

# Which SPS folder(s) to process
SPS_DIR = None
for a in sys.argv[1:]:
    if a.startswith("--sps-dir="):
        SPS_DIR = a.split("=", 1)[1]

if USE_EFF:
    SPS_BASE  = "UnweightMC"
    DPS_FILE  = "UnweightMC/UnweightDPS.root"
    TREE_NAME = "data"

os.makedirs(BASE_OUT, exist_ok=True)

# ── PyROOT UChar_t helper ───────────────────────────────────
def uchar_bool(val):
    if isinstance(val, str): return ord(val[0]) != 0 if len(val) == 1 else bool(val)
    if isinstance(val, bytes): return val[0] != 0
    return bool(val)

# ── Bin definitions: fine binning for trend visualisation ───
BINS = {
    "delta_y": {
        "title": "|#Deltay(J/#psi, #psi(2S))|",
        "edges": np.linspace(0, 4.0, 41),
        "xlabel": "|#Deltay|",
    },
    "delta_phi": {
        "title": "|#Delta#phi(J/#psi, #psi(2S))|/#pi",
        "edges": np.linspace(0, 1.0, 51),
        "xlabel": "|#Delta#phi| / #pi",
    },
    "evt_mass": {
        "title": "M_{inv}(J/#psi, #psi(2S))",
        "edges": np.linspace(7.5, 107.5, 51),
        "xlabel": "M_{inv} [GeV]",
    },
    "evt_y": {
        "title": "|y(J/#psi, #psi(2S))|",
        "edges": np.linspace(0, 2.0, 41),
        "xlabel": "|y|",
    },
    "evt_pt": {
        "title": "p_{T}(J/#psi, #psi(2S))",
        "edges": np.linspace(0, 80, 41),
        "xlabel": "p_{T} [GeV]",
    },
}

# ── Branch list ─────────────────────────────────────────────
GE_BRANCHES = [
    "GEevt_valid", "GEevt_passAcc", "GEevt_fourMuMass",
    "GEJpsi_pt", "GEJpsi_eta", "GEJpsi_phi", "GEJpsi_mass",
    "GEpsi2S_pt", "GEpsi2S_eta", "GEpsi2S_phi", "GEpsi2S_mass",
]
RECO_BRANCHES = [
    "REJpsi_pt", "REJpsi_y", "REJpsi_phi", "REJpsi_mass",
    "REpsi2S_pt", "REpsi2S_y", "REpsi2S_phi", "REpsi2S_mass",
    "REevt_fourMuMass", "REevt_matchTrg", "REmu_pt",
    "REevt_passHLT",
] if USE_RECO else []

# ── Extract distributions ───────────────────────────────────
def extract(files, label):
    """files: single path or list of paths.  Returns {vals, n_pass} or None."""
    if isinstance(files, str):
        files = [files]
    chain = ROOT.TChain(TREE_NAME)
    for p in files:
        chain.Add(p)
    n_total = chain.GetEntries()
    if n_total == 0:
        print(f"  ERROR: 0 entries for {label}")
        return None

    vals = {k: [] for k in BINS}

    if USE_EFF:
        # ── Eff mode: read pre-computed template-fit variables from data tree ──
        chain.SetBranchStatus("*", 0)
        chain.SetBranchStatus("delta_y", 1);   chain.SetBranchStatus("delta_phi", 1)
        chain.SetBranchStatus("evt_mass", 1);  chain.SetBranchStatus("evt_y", 1)
        chain.SetBranchStatus("evt_pt", 1)
        for i in range(n_total):
            chain.GetEntry(i)
            vals["delta_y"].append(chain.delta_y)
            vals["delta_phi"].append(chain.delta_phi / np.pi)  # tree stores in rad
            vals["evt_mass"].append(chain.evt_mass)
            vals["evt_y"].append(chain.evt_y)
            vals["evt_pt"].append(chain.evt_pt)
        n_pass = n_total
    else:
        # ── GEN/RECO mode: read raw ntuple, compute four-vectors ──
        chain.SetBranchStatus("*", 0)
        for b in GE_BRANCHES + RECO_BRANCHES:
            chain.SetBranchStatus(b, 1)

        n_pass = 0
        for i in range(n_total):
            chain.GetEntry(i)
            if not (uchar_bool(chain.GEevt_valid) and uchar_bool(chain.GEevt_passAcc)
                    and chain.GEevt_fourMuMass > 7.5):
                continue

            if USE_RECO and chain.REevt_fourMuMass.size() > 0:
                j_pt = chain.REJpsi_pt[0];   j_y = chain.REJpsi_y[0]
                j_phi = chain.REJpsi_phi[0]; j_m = chain.REJpsi_mass[0]
                p_pt = chain.REpsi2S_pt[0];  p_y = chain.REpsi2S_y[0]
                p_phi = chain.REpsi2S_phi[0]; p_m = chain.REpsi2S_mass[0]
                j_eta, p_eta = j_y, p_y
            else:
                j_pt = chain.GEJpsi_pt;   j_eta = chain.GEJpsi_eta
                j_phi = chain.GEJpsi_phi; j_m   = chain.GEJpsi_mass
                p_pt = chain.GEpsi2S_pt;  p_eta = chain.GEpsi2S_eta
                p_phi = chain.GEpsi2S_phi; p_m  = chain.GEpsi2S_mass

            if not (10 < j_pt < 40 and -2 < j_eta < 2): continue
            if not (10 < p_pt < 40 and -2 < p_eta < 2): continue

            if USE_RECO:
                n_combo = chain.REevt_fourMuMass.size()
                if n_combo == 0: continue
                has_hlt = any(chain.REevt_passHLT[j] for j in range(n_combo))
                has_trg = any(chain.REevt_matchTrg[j] for j in range(n_combo))
                mu_pt_max = max(chain.REmu_pt) if chain.REmu_pt.size() > 0 else 0
                if not (has_hlt and has_trg and mu_pt_max > 5.5): continue

            jv = ROOT.TLorentzVector()
            jv.SetPtEtaPhiM(j_pt, j_eta, j_phi, j_m)
            pv = ROOT.TLorentzVector()
            pv.SetPtEtaPhiM(p_pt, p_eta, p_phi, p_m)
            sv = jv + pv

            vals["delta_y"].append(abs(jv.Rapidity() - pv.Rapidity()))
            vals["delta_phi"].append(
                abs(ROOT.TVector2.Phi_mpi_pi(jv.Phi() - pv.Phi())) / np.pi)
            vals["evt_mass"].append(sv.M())
            vals["evt_y"].append(abs(sv.Rapidity()))
            vals["evt_pt"].append(sv.Pt())
            n_pass += 1

    print(f"  {label}: {n_pass} / {n_total} events")
    return {"label": label, "n": n_pass, "vals": vals}

# ── Plot ─────────────────────────────────────────────────────
COLORS = {"DPS": ROOT.kRed + 1}
FILLS  = {"DPS": ROOT.kRed - 9}
SPS_COLORS = [ROOT.kBlue + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1,
              ROOT.kCyan + 1, ROOT.kOrange + 1]

def plot_with_ratio(dps_data, sps_list, varname, cfg, outdir):
    """Top: normalised overlay of DPS + all SPS samples.
       Bottom: each SPS/DPS ratio centred at 1."""
    edges = cfg["edges"]
    nbins = len(edges) - 1
    xlabel = cfg["xlabel"]
    edge_arr = array('d', edges)

    # ── DPS histogram ─────────────────────────────────────
    dps_arr = np.array(dps_data["vals"][varname])
    h_dps = ROOT.TH1D(f"hd_{varname}", f";{xlabel};Events", nbins, edge_arr)
    h_dps_norm = ROOT.TH1D(f"hdn_{varname}", f";{xlabel};dN/d{xlabel} (norm.)",
                           nbins, edge_arr)
    for v in dps_arr:
        h_dps.Fill(v)
        h_dps_norm.Fill(v)
    if h_dps_norm.Integral() > 0:
        h_dps_norm.Scale(1.0 / h_dps_norm.Integral(), "width")

    # ── SPS histograms ────────────────────────────────────
    hs_norm = {}
    hs_raw  = {}
    for idx, sps in enumerate(sps_list):
        lab = sps["label"]
        arr = np.array(sps["vals"][varname])
        hn = ROOT.TH1D(f"hn_{varname}_{idx}",
                       f";{xlabel};dN/d{xlabel} (norm.)",
                       nbins, edge_arr)
        hr = ROOT.TH1D(f"hr_{varname}_{idx}",
                       f";{xlabel};Events", nbins, edge_arr)
        for v in arr:
            hn.Fill(v)
            hr.Fill(v)
        if hn.Integral() > 0:
            hn.Scale(1.0 / hn.Integral(), "width")
        hs_norm[lab] = hn
        hs_raw[lab]  = hr

    # ── Ratio histograms ──────────────────────────────────
    h_ratios = {}
    max_hi = 1.0
    min_lo = 1.0
    for lab, hr in hs_raw.items():
        hr_ratio = hr.Clone(f"ratio_{varname}_{lab}")
        hr_ratio.Divide(h_dps)
        h_ratios[lab] = hr_ratio
        for j in range(1, nbins + 1):
            v = hr_ratio.GetBinContent(j)
            e = hr_ratio.GetBinError(j)
            if v > 0:
                hi = v + e
                lo = v - e
                if hi > max_hi: max_hi = hi
                if lo < min_lo: min_lo = lo
    ratio_hi = 1.2 * max_hi
    ratio_lo = 0.8 * min_lo

    # ── Canvas: two stacked pads ──────────────────────────
    cname = f"c_{os.path.basename(outdir)}_{varname}"
    c = ROOT.TCanvas(cname, cfg["title"], 900, 700)
    pad_top = ROOT.TPad(f"pt_{cname}", "top", 0, 0.30, 1, 1)
    pad_bot = ROOT.TPad(f"pb_{cname}", "bot", 0, 0, 1, 0.30)
    pad_top.SetBottomMargin(0.02); pad_top.SetTopMargin(0.08)
    pad_top.SetRightMargin(0.05)
    pad_bot.SetTopMargin(0.02); pad_bot.SetBottomMargin(0.30)
    pad_bot.SetRightMargin(0.05)
    pad_top.Draw(); pad_bot.Draw()

    # ── Top pad: normalised overlay ───────────────────────
    pad_top.cd()
    # Find global maximum across all histograms
    peak = h_dps_norm.GetMaximum()
    for hn in hs_norm.values():
        m = hn.GetMaximum()
        if m > peak: peak = m
    h_dps_norm.SetMinimum(0.0)
    h_dps_norm.SetMaximum(1.2 * peak)

    leg = ROOT.TLegend(0.68, 0.75, 0.88, 0.88)
    leg.SetTextSize(0.035)
    # DPS
    h_dps_norm.SetLineColor(COLORS["DPS"])
    h_dps_norm.SetLineWidth(2)
    h_dps_norm.Draw("HIST")
    leg.AddEntry(h_dps_norm, "DPS", "l")
    # SPS samples
    for idx, (lab, hn) in enumerate(hs_norm.items()):
        color = SPS_COLORS[idx % len(SPS_COLORS)]
        hn.SetLineColor(color)
        hn.SetLineWidth(2)
        hn.SetLineStyle(idx + 1)
        hn.Draw("HIST SAME")
        leg.AddEntry(hn, lab, "l")
    leg.Draw()

    # ── Bottom pad: ratios ────────────────────────────────
    pad_bot.cd()
    first = True
    for idx, (lab, hr) in enumerate(h_ratios.items()):
        color = SPS_COLORS[idx % len(SPS_COLORS)]
        hr.SetTitle(f";{xlabel};SPS/DPS")
        hr.SetLineColor(color)
        hr.SetLineWidth(2)
        hr.SetLineStyle(idx + 1)
        hr.SetMarkerStyle(20 + idx)
        hr.SetMarkerSize(0.4)
        hr.SetMinimum(ratio_lo)
        hr.SetMaximum(ratio_hi)
        hr.GetYaxis().SetNdivisions(505)
        hr.GetYaxis().SetTitleSize(0.12)
        hr.GetYaxis().SetTitleOffset(0.35)
        hr.GetYaxis().SetLabelSize(0.10)
        hr.GetXaxis().SetTitleSize(0.14)
        hr.GetXaxis().SetTitleOffset(0.80)
        hr.GetXaxis().SetLabelSize(0.10)
        if first:
            hr.Draw("HIST E")
            first = False
        else:
            hr.Draw("HIST E SAME")

    # Reference line at y=1
    line = ROOT.TLine(edges[0], 1.0, edges[-1], 1.0)
    line.SetLineColor(ROOT.kGray + 2)
    line.SetLineStyle(2)
    line.Draw("SAME")

    c.SaveAs(outdir + f"/{varname}.png")
    print(f"    {varname}.png")


# ── Main ─────────────────────────────────────────────────────
if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    if USE_EFF:
        # ── Eff mode: read from UnweightMC/ ──────────────────
        BASE_OUT = "fig/mc_dists/eff/"
        os.makedirs(BASE_OUT, exist_ok=True)

        # Find SPS files: UnweightMC/Unweight_*.root (not DPS)
        all_unweight = sorted(glob.glob(os.path.join(SPS_BASE, "Unweight_*.root")))
        sps_to_run = []
        for p in all_unweight:
            name = os.path.basename(p)
            if name == "UnweightDPS.root": continue
            # Extract channel name: Unweight_<channel>.root -> <channel>
            ch = name[len("Unweight_"):-len(".root")]
            sps_to_run.append((ch, [p]))

        if not sps_to_run:
            print("No Unweight_*.root files found in %s" % SPS_BASE)
            sys.exit(1)

        # Load DPS
        print("Loading DPS (eff mode)...")
        dps = extract([DPS_FILE], "DPS")
        if dps is None:
            print("ERROR: cannot load DPS from %s" % DPS_FILE)
            sys.exit(1)
        print()

        for ch_name, sps_files in sps_to_run:
            tag = "Unweight_" + ch_name
            outdir = os.path.join(BASE_OUT, ch_name)
            os.makedirs(outdir, exist_ok=True)
            print("=== %s ===" % tag)
            sps = extract(sps_files, tag)
            if sps is None:
                print("  SKIP %s: extraction failed\n" % tag)
                continue
            print("  Plotting template-fit variables vs DPS...")
            for vname, cfg in BINS.items():
                plot_with_ratio(dps, [sps], vname, cfg, outdir)
            print("  -> %s/\n" % outdir)

    else:
        # ── GEN/RECO mode (original) ─────────────────────────
        if SPS_DIR:
            sps_dirs = [SPS_DIR]
        else:
            sps_dirs = sorted(
                d for d in os.listdir(SPS_BASE)
                if os.path.isdir(os.path.join(SPS_BASE, d)))

        sps_to_run = []
        for d in sps_dirs:
            full = os.path.join(SPS_BASE, d)
            roots = sorted(glob.glob(os.path.join(full, "*.root")))
            if roots:
                sps_to_run.append((d, roots))
            else:
                print(f"Skip {d}: no .root files")

        if not sps_to_run:
            print("No SPS samples found.")
            sys.exit(1)

        print("Loading DPS...")
        dps = extract(DPS_FILE, "DPS")
        if dps is None:
            print("ERROR: cannot load DPS")
            sys.exit(1)
        print()

        for folder_name, sps_files in sps_to_run:
            tag = f"SPS_{folder_name}"
            level = "RECO" if USE_RECO else "GEN"
            outdir = os.path.join(BASE_OUT, level, folder_name)
            os.makedirs(outdir, exist_ok=True)
            print(f"=== {tag} ({len(sps_files)} file(s)) ===")

            sps = extract(sps_files, tag)
            if sps is None:
                print(f"  SKIP {tag}: extraction failed\n")
                continue

            print(f"  Plotting vs DPS at {level} level...")

            for vname, cfg in BINS.items():
                plot_with_ratio(dps, [sps], vname, cfg, outdir)

            print(f"  -> {outdir}/\n")

    print(f"Done. Plots in {BASE_OUT}*/")
