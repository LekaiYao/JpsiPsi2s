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
inFileDir = "pTHat4/"
inFileList = ["Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 41)]\
# + ["2/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]\
# + ["3/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]\
# + ["4/Ntuple_2018_DPS_{}.root".format(i) for i in range(1, 11)]
inFiles = {}
inTrees = []
for inFileName in inFileList:
    inFiles[inFileName] = TFile(inFileDir + inFileName, "READ")
    inTrees.append(inFiles[inFileName].Get("rootuple/oniaTree"))
for inTree in inTrees:
    print(inTree)

# Define histograms, including efficiency(h) and error(d).
# x-axis: pT, y-axis: y
hReco_Jpsi = TH2D("hReco_Jpsi", "hReco_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hReco_psi2S = TH2D("hReco_psi2S", "hReco_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hId_Jpsi = TH2D("hId_Jpsi", "hId_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hId_psi2S = TH2D("hId_psi2S", "hId_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hVtx_Jpsi = TH2D("hVtx_Jpsi", "hVtx_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hVtx_psi2S = TH2D("hVtx_psi2S", "hVtx_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dReco_Jpsi = TH2D("dReco_Jpsi", "dReco_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dReco_psi2S = TH2D("dReco_psi2S", "dReco_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dId_Jpsi = TH2D("dId_Jpsi", "dId_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dId_psi2S = TH2D("dId_psi2S", "dId_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dVtx_Jpsi = TH2D("dVtx_Jpsi", "dVtx_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
dVtx_psi2S = TH2D("dVtx_psi2S", "dVtx_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)
hBin_Jpsi = TH2D("hBin_Jpsi", "hBin_Jpsi", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(J/psi)
hBin_psi2S = TH2D("hBin_psi2S", "hBin_psi2S", len(ptBin) - 1, ptBin, len(yBin) - 1, yBin)# data point distribution(psi(2S))
# x-axis: J/psi pT, y-axis: psi(2S) pT
hHlt = TH2D("hHlt", "hHlt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
hTrg = TH2D("hTrg", "hTrg", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
dHlt = TH2D("dHlt", "dHlt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
dTrg = TH2D("dTrg", "dTrg", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)
hBin_evt = TH2D("hBin_evt", "hBin_evt", len(ptBin) - 1, ptBin, len(ptBin) - 1, ptBin)# data point distribution

# Define matrices
mBin_Jpsi, mVtx_Jpsi, mBin_psi2S, mVtx_psi2S = tuple([np.zeros((len(yBin) - 1, len(ptBin) - 1), dtype=int) for i in range(4)])
mVtx_evt, mTrg_evt = tuple([np.zeros((len(ptBin) - 1, len(ptBin) - 1), dtype=int) for i in range(2)])

# Loop on all rectangular bins, calculate efficiency and fill histograms.
entryCut = "GEevt_valid && GEevt_passAcc && GEevt_fourMuMass > 7.5 && "# fabs(GEJpsi_phi - GEpsi2S_phi) > 0.6 && fabs(GEJpsi_phi - GEpsi2S_phi) < 1.2 && 
fidJpsiCut = "GEJpsi_pt > {} && GEJpsi_pt < {} && GEJpsi_y > {} && GEJpsi_y < {} && ".format(ptBin[0], ptBin[-1], yBin[0], yBin[-1])
fidpsi2SCut = "GEpsi2S_pt > {} && GEpsi2S_pt < {} && GEpsi2S_y > {} && GEpsi2S_y < {} && ".format(ptBin[0], ptBin[-1], yBin[0], yBin[-1])
fidYCut = "GEJpsi_y > {} && GEJpsi_y < {} && GEpsi2S_y > {} && GEpsi2S_y < {} && ".format(yBin[0], yBin[-1], yBin[0], yBin[-1])
nBin_Jpsi_max, nBin_psi2S_max, nBin_evt_max = 0, 0, 0
for ptIdx in range(len(ptBin) - 1):
    for yIdy in range(len(yBin) - 1):
        binCut_Jpsi = entryCut + fidpsi2SCut + "GEJpsi_pt > {} && GEJpsi_pt < {} && GEJpsi_y > {} && GEJpsi_y < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], yBin[yIdy], yBin[yIdy + 1])
        binCut_psi2S = entryCut + fidJpsiCut + "GEpsi2S_pt > {} && GEpsi2S_pt < {} && GEpsi2S_y > {} && GEpsi2S_y < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], yBin[yIdy], yBin[yIdy + 1])
        nBin_Jpsi, nBin_psi2S, nReco_Jpsi, nReco_psi2S, nId_Jpsi, nId_psi2S, nVtx_Jpsi, nVtx_psi2S = (0.,) * 8
        for inTree in inTrees:
            nBin_Jpsi += inTree.GetEntries(binCut_Jpsi)
            nReco_Jpsi += inTree.GetEntries(binCut_Jpsi + " && GEJpsi_matchGEN")
            nId_Jpsi += inTree.GetEntries(binCut_Jpsi + " && GEJpsi_passID")
            nVtx_Jpsi += inTree.GetEntries(binCut_Jpsi + " && GEJpsi_passVtx")
            nBin_psi2S += inTree.GetEntries(binCut_psi2S)
            nReco_psi2S += inTree.GetEntries(binCut_psi2S + " && GEpsi2S_matchGEN")
            nId_psi2S += inTree.GetEntries(binCut_psi2S + " && GEpsi2S_passID")
            nVtx_psi2S += inTree.GetEntries(binCut_psi2S + " && GEpsi2S_passVtx")
        print("\033[7mProcessing pt-y bin: {}, {}\033[0m".format(ptBin[ptIdx], yBin[yIdy]))
        print("[\033[32mJ/psi\033[0m] Vtx/Id/Reco/Total = %d/%d/%d/%d" % (nVtx_Jpsi, nId_Jpsi, nReco_Jpsi, nBin_Jpsi))
        print("[\033[31mpsi(2S)\033[0m] Vtx/Id/Reco/Total = %d/%d/%d/%d" % (nVtx_psi2S, nId_psi2S, nReco_psi2S, nBin_psi2S))
        bin = hReco_Jpsi.GetBin(ptIdx + 1, yIdy + 1)
        hBin_Jpsi.SetBinContent(bin, nBin_Jpsi)
        hBin_psi2S.SetBinContent(bin, nBin_psi2S)
        if nBin_Jpsi > nBin_Jpsi_max:
            nBin_Jpsi_max = nBin_Jpsi
        if nBin_psi2S > nBin_psi2S_max:
            nBin_psi2S_max = nBin_psi2S
        if nBin_Jpsi:
            hReco_Jpsi.SetBinContent(bin, nReco_Jpsi / nBin_Jpsi)
            dReco_Jpsi.SetBinContent(bin, sqrt((nBin_Jpsi - nReco_Jpsi) * nReco_Jpsi / nBin_Jpsi) / nBin_Jpsi)
        if nReco_Jpsi:
            hId_Jpsi.SetBinContent(bin, nId_Jpsi / nReco_Jpsi)
            dId_Jpsi.SetBinContent(bin, sqrt((nReco_Jpsi - nId_Jpsi) * nId_Jpsi / nReco_Jpsi) / nReco_Jpsi)
        if nId_Jpsi:
            hVtx_Jpsi.SetBinContent(bin, nVtx_Jpsi / nId_Jpsi)
            dVtx_Jpsi.SetBinContent(bin, sqrt((nId_Jpsi - nVtx_Jpsi) * nVtx_Jpsi / nId_Jpsi) / nId_Jpsi)
        if nBin_psi2S:
            hReco_psi2S.SetBinContent(bin, nReco_psi2S / nBin_psi2S)
            dReco_psi2S.SetBinContent(bin, sqrt((nBin_psi2S - nReco_psi2S) * nReco_psi2S / nBin_psi2S) / nBin_psi2S)
        if nReco_psi2S:
            hId_psi2S.SetBinContent(bin, nId_psi2S / nReco_psi2S)
            dId_psi2S.SetBinContent(bin, sqrt((nReco_psi2S - nId_psi2S) * nId_psi2S / nReco_psi2S) / nReco_psi2S)
        if nId_psi2S:
            hVtx_psi2S.SetBinContent(bin, nVtx_psi2S / nId_psi2S)
            dVtx_psi2S.SetBinContent(bin, sqrt((nId_psi2S - nVtx_psi2S) * nVtx_psi2S / nId_psi2S) / nId_psi2S)
        mBin_Jpsi[yIdy][ptIdx] = int(nBin_Jpsi)
        mVtx_Jpsi[yIdy][ptIdx] = int(nVtx_Jpsi)
        mBin_psi2S[yIdy][ptIdx] = int(nBin_psi2S)
        mVtx_psi2S[yIdy][ptIdx] = int(nVtx_psi2S)
    for ptIdy in range(len(ptBin) - 1):
        binCut = entryCut + fidYCut + "GEJpsi_pt > {} && GEJpsi_pt < {} && GEpsi2S_pt > {} && GEpsi2S_pt < {}".format(ptBin[ptIdx], ptBin[ptIdx + 1], ptBin[ptIdy], ptBin[ptIdy + 1])
        nVtx, nHlt, nTrg = (0.,) * 3
        for inTree in inTrees:
            nVtx += inTree.GetEntries(binCut + " && GEJpsi_passVtx && GEpsi2S_passVtx")
            nHlt += inTree.GetEntries(binCut + " && GEevt_passHLT")
            nTrg += inTree.GetEntries(binCut + " && GEevt_matchTrg && GEevt_muPtMax > 5.5")# or: GEevt_L1muPtMax
        print("\033[7mProcessing pt-pt bin: {}, {}\033[0m".format(ptBin[ptIdx], ptBin[ptIdy]))
        # print("[\033[36mevent\033[0m] Trg/Hlt/Vtx = %d/%d/%d" % (nTrg, nHlt, nVtx))
        bin = hHlt.GetBin(ptIdx + 1, ptIdy + 1)
        hBin_evt.SetBinContent(bin, nVtx)
        if nVtx > nBin_evt_max:
            nBin_evt_max = nVtx
        if nVtx:
            hHlt.SetBinContent(bin, nHlt / nVtx)
            dHlt.SetBinContent(bin, sqrt((nVtx - nHlt) * nHlt / nVtx) / nVtx)
        if nHlt:
            hTrg.SetBinContent(bin, nTrg / nHlt)
            dTrg.SetBinContent(bin, sqrt((nHlt - nTrg) * nTrg / nHlt) / nHlt)
        mVtx_evt[ptIdy][ptIdx] = int(nVtx)
        mTrg_evt[ptIdy][ptIdx] = int(nHlt)#nTrg

# Plot [histo] with [titles(x, y, z)] and save as [name].
width, height = 3000, 1200
squarish = False
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
    histo.Draw("colzTEXT" + ("45" if squarish else ""))
    histo.GetXaxis().SetNdivisions(505)
    if squarish:
        histo.GetYaxis().SetNdivisions(505)
    if width == height:
        histo.GetYaxis().SetNdivisions(505)
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
plotSave(hReco_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_RECO_Jpsi"], "Eff_RECO_Jpsi")
plotSave(hId_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_ID_Jpsi"], "Eff_ID_Jpsi")
plotSave(hVtx_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "eff_Dimuon_Jpsi"], "Eff_Dimuon_Jpsi")
plotSave(hReco_psi2S, ["pT_psi2S", "rapidity_psi2S", "eff_RECO_psi2S"], "Eff_RECO_psi2S")
plotSave(hId_psi2S, ["pT_psi2S", "rapidity_psi2S", "eff_ID_psi2S"], "Eff_ID_psi2S")
plotSave(hVtx_psi2S, ["pT_psi2S", "rapidity_psi2S", "eff_Dimuon_psi2S"], "Eff_Dimuon_psi2S")
plotSave(dReco_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_RECO_Jpsi"], "EffErr_RECO_Jpsi")
plotSave(dId_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_ID_Jpsi"], "EffErr_ID_Jpsi")
plotSave(dVtx_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "effErr_Dimuon_Jpsi"], "EffErr_Dimuon_Jpsi")
plotSave(dReco_psi2S, ["pT_psi2S", "rapidity_psi2S", "effErr_RECO_psi2S"], "EffErr_RECO_psi2S")
plotSave(dId_psi2S, ["pT_psi2S", "rapidity_psi2S", "effErr_ID_psi2S"], "EffErr_ID_psi2S")
plotSave(dVtx_psi2S, ["pT_psi2S", "rapidity_psi2S", "effErr_Dimuon_psi2S"], "EffErr_Dimuon_psi2S")
ROOT.gStyle.SetPaintTextFormat("1.0f");
zRange = nBin_Jpsi_max
plotSave(hBin_Jpsi, ["pT_Jpsi", "rapidity_Jpsi", "N_Jpsi"], "N_distri_Jpsi")
zRange = nBin_psi2S_max
plotSave(hBin_psi2S, ["pT_psi2S", "rapidity_psi2S", "N_psi2S"], "N_distri_psi2S")
squarish, height = True, 2500
zRange = nBin_evt_max
plotSave(hBin_evt, ["pT_Jpsi", "pT_psi2S", "N_event"], "N_distri_dimuon_pT")
ROOT.gStyle.SetPaintTextFormat("1.2f");
zRange = 1.
plotSave(hHlt, ["pT_Jpsi", "pT_psi2S", "eff_HLT"], "Eff_HLT")
plotSave(dHlt, ["pT_Jpsi", "pT_psi2S", "effErr_HLT"], "EffErr_HLT")
plotSave(hTrg, ["pT_Jpsi", "pT_psi2S", "eff_FourMu"], "Eff_FourMu")
plotSave(dTrg, ["pT_Jpsi", "pT_psi2S", "effErr_FourMu"], "EffErr_FourMu")

# Output efficiency file format:
# n(number of pT bins), m(number of y bins)
# (n+1) double matrix: pT bin boundary
# (m+1) double matrix: y bin boundary
# m*(2n) integer matrix: N1(number of J/psi with vertex), N0(number of accepted J/psi)
# m*(2n) integer matrix: N1(number of psi(2S) with vertex), N0(number of accepted psi(2S))
# n*(2n) integer matrix: N1(number of trigger-matched events), N0(number of events with vertices)
with open(outPlotDir + "efficiency.txt", 'w') as f:
    n_Jpsi, n_psi2S = 0, 0
    n, m = len(ptBin) - 1, len(yBin) - 1
    for i in range(m):
        for j in range(n):
            n_Jpsi += mBin_Jpsi[i][j]
            n_psi2S += mBin_psi2S[i][j]
    f.write("{} {} {} {}\n".format(n, m, n_Jpsi, n_psi2S))
    l = ""
    for pt in ptBin:
        l += str(pt) + ' '
    l += '\n'
    for y in yBin:
        l += str(y) + ' '
    l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mVtx_Jpsi[i][j]) + ' ' + str(mBin_Jpsi[i][j]) + ' '
        l += '\n'
    for i in range(m):
        for j in range(n):
            l += str(mVtx_psi2S[i][j]) + ' ' + str(mBin_psi2S[i][j]) + ' '
        l += '\n'
    for i in range(n):
        for j in range(n):
            l += str(mTrg_evt[i][j]) + ' ' + str(mVtx_evt[i][j]) + ' '
        l += '\n'
    f.write(l)
