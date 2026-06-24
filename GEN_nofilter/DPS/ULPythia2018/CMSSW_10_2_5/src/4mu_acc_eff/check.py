import sys, os, string, re, pwd, commands, ast, optparse, shlex, time
from array import array
from math import *
from decimal import *
from sample_shortnames_JJ import *
from ROOT import *
from LoadData_JJ import *

List = []
for long, short in sample_shortnames.iteritems():
    List.append(long)

Histos = {}
Histos["check"] = TH1D("check", "check", 25, 0, 10000)
Histos["check"].Sumw2()
Histos["check_2"] = TH1D("check_2", "check_2", 25, 0, 10000)
Histos["check_2"].Sumw2()

genweight = "1.0"

obs_gen_low_pt = 10
obs_gen_high_pt = 15
obs_gen_low_y = 0.5
obs_gen_high_y = 1.0
gen_mu_eta_cut = 2.4
gen_mu_pt_cut = 3.5

nEvents = {}
nEvents["check"] = 0.0
nEvents["check_2"] = 0.0

for Sample in List:
    cut_charm_gen_ya = "(GENjpsi_eta[0] > "+str(obs_gen_low_y)+ " && GENjpsi_eta[0] < " + str(obs_gen_high_y) + ")"
    # cut_charm_gen_yb = "(GENpsi2s_eta[0] > "+str(obs_gen_low_y)+ " && GENpsi2s_eta[0] < " + str(obs_gen_high_y) + ")"
    cut_charm_gen_pta = "(GENjpsi_pt[0] > "+str(obs_gen_low_pt)+ " && GENjpsi_pt[0] < " + str(obs_gen_high_pt) + ")"  
    # cut_charm_gen_ptb = "(GENpsi2s_pt[0] > "+str(obs_gen_low_pt)+ " && GENpsi2s_pt[0] < " + str(obs_gen_high_pt) + ")"  
    cut_mu_gen_eta_a = "(abs(GENjpsi_mu_eta[0]) < "+str(gen_mu_eta_cut)+")"
    cut_mu_gen_pt_a = "(GENjpsi_mu_pt[0] > "+str(gen_mu_pt_cut)+")"
    # cut_mu_gen_eta_b = "(abs(GENpsi2s_mu_eta[0]) < "+str(gen_mu_eta_cut)+ " && abs(GENpsi2s_mu_eta[1]) < "+str(gen_mu_eta_cut)+")"
    # cut_mu_gen_pt_b = "(GENpsi2s_mu_pt[0] > "+str(gen_mu_pt_cut)+ " && GENpsi2s_mu_pt[1] > "+str(gen_mu_pt_cut)+")"
    cut_mu_status_a = "(GEN_pair_id[0].first)"
    # cut_mu_status_b = "(GENpsi2s_mu_pt[0])"

    #fs_a
    Tree[Sample].Draw("GENjpsi_mass[0]>> check","("+genweight+")*("+cut_charm_gen_ya+" && "+cut_charm_gen_pta+" && "+cut_mu_status_a+")","goff")
    #fid_eta_a
    Tree[Sample].Draw("GENjpsi_mass[0]>> check_2","("+genweight+")*("+cut_charm_gen_ya+" && "+cut_charm_gen_pta+" && "+cut_mu_gen_eta_a+")","goff")

    nEvent = Histos["check"].Integral()
    print Sample, nEvent
    nEvents["check"] += nEvent

    nEvent = Histos["check_2"].Integral()
    print Sample, nEvent
    nEvents["check_2"] += nEvent

print nEvents["check"]
