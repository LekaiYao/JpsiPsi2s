from ROOT import *
from array import array
import os
dirDPS = '/eos/home-c/chensh/JPsiPsi2s/GEN_nofilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff'
dirSPS = '/eos/home-c/chensh/JPsiPsi2s/GEN_nofilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff'
SamplesMC = ['6/DPS_2018_Psi2SJ_{}.root'.format(i) for i in range(1, 9)] \
+ ['8/DPS_2018_Psi2SJ_{}.root'.format(i) for i in range(1, 6)] \
+ ['8/DPS_2018_Psi2SJ_{}.root'.format(i) for i in range(7, 11)] \
+ ['8/{}/DPS_2018_Psi2SJ_{}.root'.format(i, j) for i in range(1, 6) for j in range(1, 21)] \
+ ['SPS_2018_Psi2SJ_%d.root' % i for i in range(1, 11)]
RootFile = {}
Tree = {}
for i in range(len(SamplesMC)):
    sample = SamplesMC[i][:-5]
    if i < 117:
        RootFile[sample] = TFile(dirDPS + '/' + SamplesMC[i], "READ")
    else:
	    RootFile[sample] = TFile(dirSPS + '/' + SamplesMC[i], "READ")
    Tree[sample] = RootFile[sample].Get("GenAnalyzer/gen_tree")
    if not Tree[sample]:
        print sample+' has no gen. tree'
    else:
        print sample
