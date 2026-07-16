#!/usr/bin/env python3
"""Compare our WeightData.root vs chensh WeightData_prime.root on fit variables."""
import ROOT, numpy as np, os
from array import array

OURS   = "WeightData.root"
CHENSH = "/eos/user/c/chensh/JPsiPsi2s/Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/WeightData_prime.root"
OUTDIR = "fig/"
os.makedirs(OUTDIR, exist_ok=True)

VARS = [
    ("Jpsi_mass",   2.95, 3.25, 30, "M(J/#psi) [GeV]"),
    ("psi2S_mass",  3.35, 4.05, 30, "M(#psi(2S)) [GeV]"),
    ("Jpsi_ctau",  -0.03, 0.16, 30, "c#tau(J/#psi) [cm]"),
    ("psi2S_ctau", -0.03, 0.16, 30, "c#tau(#psi(2S)) [cm]"),
]

f_ours   = ROOT.TFile(OURS)
f_chensh = ROOT.TFile(CHENSH)
t_ours   = f_ours.Get("data")
t_chensh = f_chensh.Get("data")
print("Our entries: %d   Chensh entries: %d" % (t_ours.GetEntries(), t_chensh.GetEntries()))

for vname, xmin, xmax, nbins, title in VARS:
    edge_arr = array('d', np.linspace(xmin, xmax, nbins+1))

    h_ours   = ROOT.TH1D("ho_"+vname, ";%s;Weighted events" % title, nbins, edge_arr)
    h_chensh = ROOT.TH1D("hc_"+vname, ";%s;Weighted events" % title, nbins, edge_arr)
    h_ours_n   = ROOT.TH1D("hon_"+vname, ";%s;dN/d%s (norm.)" % (title, title), nbins, edge_arr)
    h_chensh_n = ROOT.TH1D("hcn_"+vname, ";%s;dN/d%s (norm.)" % (title, title), nbins, edge_arr)

    # Need to read branches
    t_ours.SetBranchStatus("*", 0)
    t_ours.SetBranchStatus(vname, 1); t_ours.SetBranchStatus("evt_weight", 1)
    t_chensh.SetBranchStatus("*", 0)
    t_chensh.SetBranchStatus(vname, 1); t_chensh.SetBranchStatus("evt_weight", 1)

    for i in range(t_ours.GetEntries()):
        t_ours.GetEntry(i)
        h_ours.Fill(getattr(t_ours, vname), t_ours.evt_weight)
        h_ours_n.Fill(getattr(t_ours, vname), t_ours.evt_weight)
    for i in range(t_chensh.GetEntries()):
        t_chensh.GetEntry(i)
        h_chensh.Fill(getattr(t_chensh, vname), t_chensh.evt_weight)
        h_chensh_n.Fill(getattr(t_chensh, vname), t_chensh.evt_weight)

    if h_ours_n.Integral() > 0:   h_ours_n.Scale(1.0/h_ours_n.Integral(), "width")
    if h_chensh_n.Integral() > 0: h_chensh_n.Scale(1.0/h_chensh_n.Integral(), "width")

    # Ratio
    h_ratio = h_ours_n.Clone("ratio_"+vname)
    h_ratio.Divide(h_chensh_n)
    max_dev = max(abs(h_ratio.GetBinContent(j)-1.0) for j in range(1, nbins+1) if h_ratio.GetBinContent(j)>0)
    if max_dev < 0.05: max_dev = 0.5
    ratio_lo, ratio_hi = 1.0-1.1*max_dev, 1.0+1.1*max_dev

    # Canvas
    c = ROOT.TCanvas("c_"+vname, title, 900, 700)
    ptop = ROOT.TPad("pt_"+vname, "top", 0, 0.30, 1, 1)
    pbot = ROOT.TPad("pb_"+vname, "bot", 0, 0, 1, 0.30)
    ptop.SetBottomMargin(0.02); ptop.SetTopMargin(0.08); ptop.SetRightMargin(0.05)
    pbot.SetTopMargin(0.02); pbot.SetBottomMargin(0.30); pbot.SetRightMargin(0.05)
    ptop.Draw(); pbot.Draw()

    ptop.cd()
    h_ours.SetLineColor(ROOT.kBlue+1);   h_ours.SetLineWidth(2)
    h_chensh.SetLineColor(ROOT.kRed+1);  h_chensh.SetLineWidth(2)
    peak = max(h_ours.GetMaximum(), h_chensh.GetMaximum())
    h_ours.SetMinimum(0); h_ours.SetMaximum(1.2*peak)
    h_ours.Draw("HIST"); h_chensh.Draw("HIST SAME")
    leg = ROOT.TLegend(0.68, 0.78, 0.88, 0.88)
    leg.AddEntry(h_ours, "Ours (%d)" % t_ours.GetEntries(), "l")
    leg.AddEntry(h_chensh, "Chensh (%d)" % t_chensh.GetEntries(), "l")
    leg.Draw()

    pbot.cd()
    h_ratio.SetTitle(";%s;Ours / Chensh" % title)
    h_ratio.SetLineColor(ROOT.kViolet+2); h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerStyle(20); h_ratio.SetMarkerSize(0.5)
    h_ratio.SetMinimum(ratio_lo); h_ratio.SetMaximum(ratio_hi)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(0.12); h_ratio.GetYaxis().SetTitleOffset(0.35)
    h_ratio.GetYaxis().SetLabelSize(0.10)
    h_ratio.GetXaxis().SetTitleSize(0.14); h_ratio.GetXaxis().SetTitleOffset(0.80)
    h_ratio.GetXaxis().SetLabelSize(0.10)
    h_ratio.Draw("HIST E")
    ROOT.TLine(xmin, 1.0, xmax, 1.0).Draw("SAME")

    c.SaveAs(OUTDIR+"compare_"+vname+".png")
    print("  Saved %scompare_%s.png" % (OUTDIR, vname))

f_ours.Close(); f_chensh.Close()
print("Done.")
