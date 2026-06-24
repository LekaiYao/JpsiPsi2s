import optparse
import numpy as np
from ROOT import *
from LoadData_JJ import *
from tdrStyle import setTDRStyle
from tqdm import tqdm

setTDRStyle()

# def parseOptions():
#     usage = ('usage: %prog [options]\n' + '%prog -h for help')
#     parser = optparse.OptionParser(usage)
#     # input options
#     parser.add_option('', '--obsBins', dest='OBSBINS', type='string', default='',
#         help='Bin boundaries for the diff. measurement separated by "|", e.g. as "|0|50|100|", use the defalut if empty string')
#     # store options and arguments as global variables
#     global opt, args
#     (opt, args) = parser.parse_args()
# # parse the arguments and options
# global opt, args
# parseOptions()
# obs_bins_pt = opt.OBSBINS.strip('|').split('|')
# obs_bins_pt = array('d', [float(obs_bins_pt[i]) for i in range(len(obs_bins_pt))])
obs_bins_pt = array('d', [i for i in range(10, 24)] + [24, 26, 28, 30, 35, 40])

gen_mu_eta_cut = 2.4
gen_mu_pt_cut = 3.5

eta_acc_Jpsi = [TH1D("DPS_eta_acc_Jpsi", "DPS_eta_acc_Jpsi", len(obs_bins_pt) - 1, obs_bins_pt), TH1D("SPS_eta_acc_Jpsi", "SPS_eta_acc_Jpsi", len(obs_bins_pt) - 1, obs_bins_pt)]
pteta_acc_Jpsi = [TH1D("DPS_pteta_acc_Jpsi", "DPS_pteta_acc_Jpsi", len(obs_bins_pt) - 1, obs_bins_pt), TH1D("SPS_pteta_acc_Jpsi", "SPS_pteta_acc_Jpsi", len(obs_bins_pt) - 1, obs_bins_pt)]
eta_acc_psi2S = [TH1D("DPS_eta_acc_psi2S", "DPS_eta_acc_psi2S", len(obs_bins_pt) - 1, obs_bins_pt), TH1D("SPS_eta_acc_psi2S", "SPS_eta_acc_psi2S", len(obs_bins_pt) - 1, obs_bins_pt)]
pteta_acc_psi2S = [TH1D("DPS_pteta_acc_psi2S", "DPS_pteta_acc_psi2S", len(obs_bins_pt) - 1, obs_bins_pt), TH1D("SPS_pteta_acc_psi2S", "SPS_pteta_acc_psi2S", len(obs_bins_pt) - 1, obs_bins_pt)]

for bin in tqdm(range(len(obs_bins_pt) - 1)):
    ptmin, ptmax = obs_bins_pt[bin], obs_bins_pt[bin + 1]
    fid_Jpsi = " && GENjpsi_pt[0] > {} && GENjpsi_pt[0] < {} && abs(GENjpsi_y[0]) < 2".format(10, 40)
    fid_psi2S = " && GENpsi2s_pt[0] > {} && GENpsi2s_pt[0] < {} && abs(GENpsi2s_y[0]) < 2".format(10, 40)
    binsel_Jpsi = "GENjpsi_mass && GENpsi2s_mass && GENjpsi_pt[0] > {} && GENjpsi_pt[0] < {} && abs(GENjpsi_y[0]) < 2".format(ptmin, ptmax) + fid_psi2S
    binsel_psi2S = "GENjpsi_mass && GENpsi2s_mass && GENpsi2s_pt[0] > {} && GENpsi2s_pt[0] < {} && abs(GENpsi2s_y[0]) < 2".format(ptmin, ptmax) + fid_Jpsi
    ##
    ntot, neta, npteta = [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
    for sample in SamplesMC:
        samplename = sample[:-5]
        i = 1 if sample[0] == 'S' else 0
        ntot[i] += Tree[samplename].GetEntries(binsel_Jpsi)
        neta[i] += Tree[samplename].GetEntries(binsel_Jpsi + " && abs(GENjpsi_mu_eta[0]) < {}".format(gen_mu_eta_cut))
        npteta[i] += Tree[samplename].GetEntries(binsel_Jpsi +
            " && abs(GENjpsi_mu_eta[0]) < {} && GENjpsi_mu_pt[0] > {}".format(gen_mu_eta_cut, gen_mu_pt_cut))
    for i in range(2):
        eta_acc_Jpsi[i].SetBinContent(bin + 1, neta[i] / ntot[i])
        eta_acc_Jpsi[i].SetBinError(bin + 1, np.sqrt(neta[i] * (ntot[i] - neta[i]) / ntot[i]) / ntot[i])
        pteta_acc_Jpsi[i].SetBinContent(bin + 1, npteta[i] / neta[i])
        pteta_acc_Jpsi[i].SetBinError(bin + 1, np.sqrt(npteta[i] * (neta[i] - npteta[i]) / neta[i]) / neta[i])
    ##
    ntot, neta, npteta = [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
    for sample in SamplesMC:
        samplename = sample[:-5]
        i = 1 if sample[0] == 'S' else 0
        ntot[i] += Tree[samplename].GetEntries(binsel_psi2S)
        neta[i] += Tree[samplename].GetEntries(binsel_psi2S + " && abs(GENpsi2s_mu_eta[0]) < {}".format(gen_mu_eta_cut))
        npteta[i] += Tree[samplename].GetEntries(binsel_psi2S +
            " && abs(GENpsi2s_mu_eta[0]) < {} && GENpsi2s_mu_pt[0] > {}".format(gen_mu_eta_cut, gen_mu_pt_cut))
    for i in range(2):
        eta_acc_psi2S[i].SetBinContent(bin + 1, neta[i] / ntot[i])
        eta_acc_psi2S[i].SetBinError(bin + 1, np.sqrt(neta[i] * (ntot[i] - neta[i]) / ntot[i]) / ntot[i])
        pteta_acc_psi2S[i].SetBinContent(bin + 1, npteta[i] / neta[i])
        pteta_acc_psi2S[i].SetBinError(bin + 1, np.sqrt(npteta[i] * (neta[i] - npteta[i]) / neta[i]) / neta[i])


def plotnsave(hists, name):
    c = TCanvas("c", "c", 2000, 1500)
    c.cd()
    c.SetTopMargin(0.10)
    c.SetRightMargin(0.20)
    for hist in hists:
        hist.GetXaxis().SetTitle('pT(gen.)')
        hist.GetYaxis().SetTitle('Acc.(gen.)')
        hist.GetXaxis().SetTitleOffset(1.0)
        hist.GetYaxis().SetTitleOffset(0.8)
    hists[0].SetLineColor(kRed)
    hists[0].Draw("")
    hists[1].SetLineColor(kBlue)
    hists[1].Draw("same")
    for hist in hists:
        hist.GetXaxis().SetNdivisions(505)
    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.5 * c.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right  
    latex2.DrawLatex(0.18, 0.92, "CMS")
    latex2.SetTextSize(0.4 * c.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(15)
    latex2.DrawLatex(0.23, 0.92, "Preliminary")
    latex2.SetTextSize(0.25*c.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.DrawLatex(0.4, 0.92, name)
    latex2.DrawLatex(0.72, 0.92, " #sqrt{s} = 13 TeV")
    c.SaveAs(name + ".png")
    del c

plotnsave(eta_acc_Jpsi, "Acc_eta_Jpsi")
plotnsave(pteta_acc_Jpsi, "Acc_pt_eta_Jpsi")
plotnsave(eta_acc_psi2S, "Acc_eta_psi2S")
plotnsave(pteta_acc_psi2S, "Acc_pt_eta_psi2S")
