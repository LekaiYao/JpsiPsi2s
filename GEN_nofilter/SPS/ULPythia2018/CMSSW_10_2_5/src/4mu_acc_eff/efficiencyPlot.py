from ROOT import *
from array import array
import os
from tdrStyle import *
from math import sqrt
import numpy as np

# Specify bin settings for efficiency maps.
ptBin = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])
yBin = array('d', [-2, -1.75, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.75, 2])

# Specify output directory.
outPlotDir = "plot/"
if not os.path.exists(outPlotDir):
    os.system("mkdir " + outPlotDir)

# Handle input files and extract trees.
inFileDir = "./"
inFileList = ["SPS_2018_Psi2SJ_{}.root".format(i) for i in range(1, 11)]
inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("GenAnalyzer/gen_tree"))
for inTree in inTrees:
    print(inTree)

# Define histograms, including efficiency(h) and error(d).
# x-axis: pT, y-axis: y
hEta_Jpsi = TH2D("hEta_Jpsi", "hEta_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hEta_psi2S = TH2D("hEta_psi2S", "hEta_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hPt_Jpsi = TH2D("hPt_Jpsi", "hPt_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hPt_psi2S = TH2D("hPt_psi2S", "hPt_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dEta_Jpsi = TH2D("dEta_Jpsi", "dEta_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dEta_psi2S = TH2D("dEta_psi2S", "dEta_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dPt_Jpsi = TH2D("dPt_Jpsi", "dPt_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dPt_psi2S = TH2D("dPt_psi2S", "dPt_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hBin_Jpsi = TH2D("hBin_Jpsi", "hBin_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(J/psi)
hBin_psi2S = TH2D("hBin_psi2S", "hBin_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(psi(2S))

# Define matrices
mBin_Jpsi, mPt_Jpsi, mBin_psi2S, mPt_psi2S = tuple([np.zeros((len(yBin) - 1, len(ptBin) - 1), dtype=int) for i in range(4)])

# Loop on all rectangular bins, calculate efficiency and fill histograms.

rJpsi_pt, rJpsi_y = "GENjpsi_pt[GEN_pair_id[0].first]", "GENjpsi_y[GEN_pair_id[0].first]"
rJpsi_mupt, rJpsi_mueta = "GENjpsi_mu_pt[GEN_pair_id[0].first]", "GENjpsi_mu_eta[GEN_pair_id[0].first]"
rpsi2S_pt, rpsi2S_y = "GENpsi2s_pt[GEN_pair_id[0].second]", "GENpsi2s_y[GEN_pair_id[0].second]"
rpsi2S_mupt, rpsi2S_mueta = "GENpsi2s_mu_pt[GEN_pair_id[0].second]", "GENpsi2s_mu_eta[GEN_pair_id[0].second]"
etaCut_Jpsi, etaCut_psi2S = " && abs({}[0]) < 2.4 && abs({}[1]) < 2.4".format(rJpsi_mueta, rJpsi_mueta), " && abs({}[0]) < 2.4 && abs({}[1]) < 2.4".format(rpsi2S_mueta, rpsi2S_mueta)
# evtMassCut = "GENevt_mass > 7.5 && "
nBin_Jpsi_max, nBin_psi2S_max = 0, 0
for ptIdx in range(len(ptBin) - 1):
    for yIdy in range(len(yBin) - 1):
        binCut_Jpsi = "{} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi_pt, ptBin[ptIdx], rJpsi_pt, ptBin[ptIdx + 1], rJpsi_y, yBin[yIdy], rJpsi_y, yBin[yIdy + 1])
        binCut_psi2S = "{} > {} && {} < {} && {} > {} && {} < {}".format(rpsi2S_pt, ptBin[ptIdx], rpsi2S_pt, ptBin[ptIdx + 1], rpsi2S_y, yBin[yIdy], rpsi2S_y, yBin[yIdy + 1])
        fidCut_Jpsi = " && {} > {} && {} < {} && {} > {} && {} < {}".format(rJpsi_pt, ptBin[0], rJpsi_pt, ptBin[-1], rJpsi_y, yBin[0], rJpsi_y, yBin[-1])
        fidCut_psi2S = " && {} > {} && {} < {} && {} > {} && {} < {}".format(rpsi2S_pt, ptBin[0], rpsi2S_pt, ptBin[-1], rpsi2S_y, yBin[0], rpsi2S_y, yBin[-1])
        nBin_Jpsi, nBin_psi2S, nEta_Jpsi, nEta_psi2S, nPt_Jpsi, nPt_psi2S = (0.,) * 6
        for inTree in inTrees:
            nBin_Jpsi += inTree.GetEntries(binCut_Jpsi + fidCut_psi2S)# evtMassCut + 
            nEta_Jpsi += inTree.GetEntries(binCut_Jpsi + fidCut_psi2S + etaCut_Jpsi)
            nPt_Jpsi += inTree.GetEntries(binCut_Jpsi + fidCut_psi2S + etaCut_Jpsi + " && {}[0] > 3.5 && {}[1] > 3.5".format(rJpsi_mupt, rJpsi_mupt))
            nBin_psi2S += inTree.GetEntries(binCut_psi2S + fidCut_Jpsi)
            nEta_psi2S += inTree.GetEntries(binCut_psi2S + fidCut_Jpsi + etaCut_psi2S)
            nPt_psi2S += inTree.GetEntries(binCut_psi2S + fidCut_Jpsi + etaCut_psi2S + " && {}[0] > 3.5 && {}[1] > 3.5".format(rpsi2S_mupt, rpsi2S_mupt))
        print("\033[7mProcessing pt-y bin: {}, {}\033[0m".format(ptBin[ptIdx], yBin[yIdy]))
        print("[\033[32mJ/psi\033[0m] PtEta/Eta/Total = %d/%d/%d" % (nPt_Jpsi, nEta_Jpsi, nBin_Jpsi))
        print("[\033[31mpsi(2S)\033[0m] PtEta/Eta/Total = %d/%d/%d" % (nPt_psi2S, nEta_psi2S, nBin_psi2S))
        bin = hEta_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        hBin_Jpsi.SetBinContent(bin, nBin_Jpsi)
        hBin_psi2S.SetBinContent(bin, nBin_psi2S)
        if nBin_Jpsi > nBin_Jpsi_max:
            nBin_Jpsi_max = nBin_Jpsi
        if nBin_psi2S > nBin_psi2S_max:
            nBin_psi2S_max = nBin_psi2S
        if nBin_Jpsi:
            hEta_Jpsi.SetBinContent(bin, nEta_Jpsi / nBin_Jpsi)
            dEta_Jpsi.SetBinContent(bin, sqrt((nBin_Jpsi - nEta_Jpsi) * nEta_Jpsi / nBin_Jpsi) / nBin_Jpsi)
        if nEta_Jpsi:
            hPt_Jpsi.SetBinContent(bin, nPt_Jpsi / nEta_Jpsi)
            dPt_Jpsi.SetBinContent(bin, sqrt((nEta_Jpsi - nPt_Jpsi) * nPt_Jpsi / nEta_Jpsi) / nEta_Jpsi)
        if nBin_psi2S:
            hEta_psi2S.SetBinContent(bin, nEta_psi2S / nBin_psi2S)
            dEta_psi2S.SetBinContent(bin, sqrt((nBin_psi2S - nEta_psi2S) * nEta_psi2S / nBin_psi2S) / nBin_psi2S)
        if nEta_psi2S:
            hPt_psi2S.SetBinContent(bin, nPt_psi2S / nEta_psi2S)
            dPt_psi2S.SetBinContent(bin, sqrt((nEta_psi2S - nPt_psi2S) * nPt_psi2S / nEta_psi2S) / nEta_psi2S)
        mBin_Jpsi[yIdy][ptIdx] = int(nBin_Jpsi)
        mPt_Jpsi[yIdy][ptIdx] = int(nPt_Jpsi)
        mBin_psi2S[yIdy][ptIdx] = int(nBin_psi2S)
        mPt_psi2S[yIdy][ptIdx] = int(nPt_psi2S)

# Plot [histo] with [titles(x, y, z)] and save as [name].
width, height = 3000, 1200
zRange = 1
def plotSave(histo, titles, name):
    c = TCanvas(name, name, width, height)
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
    histo.Draw("colzTEXT")
    histo.GetXaxis().SetNdivisions(505)
    # latex2 = TLatex()
    # latex2.SetNDC()
    # latex2.SetTextSize(0.5 * c.GetTopMargin())
    # latex2.SetTextFont(62)
    # latex2.SetTextAlign(11) 
    # latex2.DrawLatex(0.18, 0.92, "CMS")
    # latex2.SetTextSize(0.4 * c.GetTopMargin())
    # latex2.SetTextFont(52)
    # latex2.SetTextAlign(11)
    # latex2.DrawLatex(0.25, 0.92, "Preliminary")
    # latex2.SetTextFont(42)
    # latex2.SetTextSize(0.25 * c.GetTopMargin())
    # latex2.DrawLatex(0.40, 0.92, "SPS_Jpsi+psi(2S)")
    # latex2.DrawLatex(0.72, 0.92, " #sqrt{s} = 13 TeV")
    c.SaveAs(outPlotDir + name + ".png")
    del c

# Draw and save every plot.
setTDRStyle()
plotSave(hEta_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "acc_Eta_Jpsi"], "Acc_Eta_Jpsi")
plotSave(hEta_psi2S, ["pT_psi2S", "rapidity_psi2S", "acc_Eta_psi2S"], "Acc_Eta_psi2S")
plotSave(hPt_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "acc_Pt_Jpsi"], "Acc_Pt_Jpsi")
plotSave(hPt_psi2S, ["pT_psi2S", "rapidity_psi2S", "acc_Pt_psi2S"], "Acc_Pt_psi2S")
plotSave(dEta_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "accErr_Eta_Jpsi"], "AccErr_Eta_Jpsi")
plotSave(dEta_psi2S, ["pT_psi2S", "rapidity_psi2S", "accErr_Eta_psi2S"], "AccErr_Eta_psi2S")
plotSave(dPt_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "accErr_Pt_Jpsi"], "AccErr_Pt_Jpsi")
plotSave(dPt_psi2S, ["pT_psi2S", "rapidity_psi2S", "accErr_Pt_psi2S"], "AccErr_Pt_psi2S")
ROOT.gStyle.SetPaintTextFormat("1.0f");
zRange = nBin_Jpsi_max
plotSave(hBin_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "N_Jpsi"], "N_distri_Jpsi")
zRange = nBin_psi2S_max
plotSave(hBin_psi2S, ["pT_psi2S", "rapidity_psi2S", "N_psi2S"], "N_distri_psi2S")

# Output efficiency file format:
# n(number of pT bins), m(number of y bins)
# (n+1) double matrix: pT bin boundary
# (m+1) double matrix: y bin boundary
# m*n integer matrix: N1(number of accepted J/psi), N0(number of initial J/psi)
# m*n integer matrix: N1(number of accepted psi(2S)), N0(number of initial psi(2S))
with open(outPlotDir + "acceptance.txt", 'w') as f:
    n, m = len(ptBin) - 1, len(yBin) - 1
    f.write("{} {}\n".format(n, m))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mPt_Jpsi[i][j]) + ' ' + str(mBin_Jpsi[i][j]) + ' '
        l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mPt_psi2S[i][j]) + ' ' + str(mBin_psi2S[i][j]) + ' '
        l += '\n'
    f.write(l)
