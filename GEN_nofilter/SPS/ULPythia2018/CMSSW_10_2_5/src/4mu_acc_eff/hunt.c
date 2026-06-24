#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
using namespace std;
// pair id control
// #define usePairId
//
#ifdef usePairId
    #define jpsiId _GEN_pair_id->at(0).first
    #define psi2sId _GEN_pair_id->at(0).second
#else
    #define jpsiId 0
    #define psi2sId 0
#endif

void hunt() {
    // read in files and trees
    int nFile = 1;
    TFile **files = new TFile *[nFile + 1];
    TTree **trees = new TTree *[nFile + 1];
    for(int i = 0; i < nFile; i++) {
        files[i] = new TFile(("SPS_JJ_withGenFilter_" + std::to_string(i + 1) + ".root").c_str(), "READ");
        trees[i] = (TTree *)files[i]->Get("GenAnalyzer/gen_tree");
    }
    // plot [J/psi, psi(2S)][mass, pT, eta], muon[pT, eta]
    int nHist = 20;
    TH1D **histos = new TH1D *[nHist + 1];
    histos[0] = new TH1D("J/psi_mass", "J/#psi_mass", 100, 3.0962, 3.0976);
    histos[2] = new TH1D("J/psi_mass_extra", "J/#psi_mass_extra", 100, 3.0962, 3.0976);
    histos[1] = new TH1D("psi(2S)_mass", "#psi(2S)_mass", 100, 3.6848, 3.6874);
    histos[3] = new TH1D("psi(2S)_mass_extra", "#psi(2S)_mass_extra", 100, 3.6848, 3.6874);//mass
    histos[4] = new TH1D("J/psi_pt", "J/#psi_pt", 100, 0, 50.0);
    histos[6] = new TH1D("J/psi_pt_extra", "J/#psi_pt_extra", 100, 0, 50.0);
    histos[5] = new TH1D("psi(2S)_pt", "#psi(2S)_pt", 100, 0, 50.0);
    histos[7] = new TH1D("psi(2S)_pt_extra", "#psi(2S)_pt_extra", 100, 0, 50.0);//pT
    histos[8] = new TH1D("J/psi_eta", "J/#psi_eta", 100, -3.0, 3.0);
    histos[10] = new TH1D("J/psi_eta_extra", "J/#psi_eta_extra", 100, -3.0, 3.0);
    histos[9] = new TH1D("psi(2S)_eta", "#psi(2S)_eta", 100, -3.0, 3.0);
    histos[11] = new TH1D("psi(2S)_eta_extra", "#psi(2S)_eta_extra", 100, -3.0, 3.0);//eta
    histos[12] = new TH1D("J/psi_mu_pt", "J/#psi_mu_pt", 100, 0, 25.0);
    histos[14] = new TH1D("J/psi_mu_pt_extra", "J/#psi_mu_pt_extra", 100, 0, 25.0);
    histos[13] = new TH1D("psi(2S)_mu_pt", "#psi(2S)_mu_pt", 100, 0, 25.0);
    histos[15] = new TH1D("psi(2S)_mu_pt_extra", "#psi(2S)_mu_pt_extra", 100, 0, 25.0);//pT
    histos[16] = new TH1D("J/psi_mu_eta", "J/#psi_mu_eta", 100, -3.0, 3.0);
    histos[18] = new TH1D("J/psi_mu_eta_extra", "J/#psi_mu_eta_extra", 100, -3.0, 3.0);
    histos[17] = new TH1D("psi(2S)_mu_eta", "#psi(2S)_mu_eta", 100, -3.0, 3.0);
    histos[19] = new TH1D("psi(2S)_mu_eta_extra", "#psi(2S)_mu_eta_extra", 100, -3.0, 3.0);//eta
    for(int i = 0; i < nFile; i++) {
        Int_t nEntry = trees[0]->GetEntries();
        std::vector<double> *_GENjpsi_mass = 0, *_GENpsi2s_mass = 0;
        std::vector<double> *_GENjpsi_pt = 0, *_GENpsi2s_pt = 0;
        std::vector<double> *_GENjpsi_eta = 0, *_GENpsi2s_eta = 0;
        std::vector<std::vector<double>> *_GENjpsi_mu_pt = 0, *_GENpsi2s_mu_pt = 0;
        std::vector<std::vector<double>> *_GENjpsi_mu_eta = 0, *_GENpsi2s_mu_eta = 0;
        TBranch *jpsi_mass_branch = trees[0]->GetBranch("GENjpsi_mass");
        TBranch *psi2s_mass_branch = trees[0]->GetBranch("GENpsi2s_mass");
        TBranch *jpsi_pt_branch = trees[0]->GetBranch("GENjpsi_pt");
        TBranch *psi2s_pt_branch = trees[0]->GetBranch("GENpsi2s_pt");
        TBranch *jpsi_eta_branch = trees[0]->GetBranch("GENjpsi_eta");
        TBranch *psi2s_eta_branch = trees[0]->GetBranch("GENpsi2s_eta");
        TBranch *jpsi_mu_pt_branch = trees[0]->GetBranch("GENjpsi_mu_pt");
        TBranch *psi2s_mu_pt_branch = trees[0]->GetBranch("GENpsi2s_mu_pt");
        TBranch *jpsi_mu_eta_branch = trees[0]->GetBranch("GENjpsi_mu_eta");
        TBranch *psi2s_mu_eta_branch = trees[0]->GetBranch("GENpsi2s_mu_eta");
        jpsi_mass_branch->SetAddress(&_GENjpsi_mass);
        psi2s_mass_branch->SetAddress(&_GENpsi2s_mass);
        jpsi_pt_branch->SetAddress(&_GENjpsi_pt);
        psi2s_pt_branch->SetAddress(&_GENpsi2s_pt);
        jpsi_eta_branch->SetAddress(&_GENjpsi_eta);
        psi2s_eta_branch->SetAddress(&_GENpsi2s_eta);
        jpsi_mu_pt_branch->SetAddress(&_GENjpsi_mu_pt);
        psi2s_mu_pt_branch->SetAddress(&_GENpsi2s_mu_pt);
        jpsi_mu_eta_branch->SetAddress(&_GENjpsi_mu_eta);
        psi2s_mu_eta_branch->SetAddress(&_GENpsi2s_mu_eta);
        #ifdef usePairId
        std::vector<std::pair<int, int>> *_GEN_pair_id = 0;
        TBranch *pair_id_branch = trees[0]->GetBranch("GEN_pair_id");
        pair_id_branch->SetAddress(&_GEN_pair_id);
        #endif
        for(Int_t j = 0; j < nEntry; j++) {
            Long64_t jEntry = trees[0]->LoadTree(j);
            jpsi_mass_branch->GetEntry(jEntry);
            psi2s_mass_branch->GetEntry(jEntry);
            jpsi_pt_branch->GetEntry(jEntry);
            psi2s_pt_branch->GetEntry(jEntry);
            jpsi_eta_branch->GetEntry(jEntry);
            psi2s_eta_branch->GetEntry(jEntry);
            jpsi_mu_pt_branch->GetEntry(jEntry);
            psi2s_mu_pt_branch->GetEntry(jEntry);
            jpsi_mu_eta_branch->GetEntry(jEntry);
            psi2s_mu_eta_branch->GetEntry(jEntry);
            #ifdef usePairId
            pair_id_branch->GetEntry(jEntry);
            #endif
            histos[0]->Fill(_GENjpsi_mass->at(0));
            histos[1]->Fill(_GENpsi2s_mass->at(0));
            histos[4]->Fill(_GENjpsi_pt->at(0));
            histos[5]->Fill(_GENpsi2s_pt->at(0));
            histos[8]->Fill(_GENjpsi_eta->at(0));
            histos[9]->Fill(_GENpsi2s_eta->at(0));
            histos[12]->Fill(_GENjpsi_mu_pt->at(jpsiId).at(0));
            histos[12]->Fill(_GENjpsi_mu_pt->at(jpsiId).at(1));
            histos[13]->Fill(_GENpsi2s_mu_pt->at(psi2sId).at(0));
            histos[13]->Fill(_GENpsi2s_mu_pt->at(psi2sId).at(1));
            histos[16]->Fill(_GENjpsi_mu_eta->at(jpsiId).at(0));
            histos[16]->Fill(_GENjpsi_mu_eta->at(jpsiId).at(1));
            histos[17]->Fill(_GENpsi2s_mu_eta->at(psi2sId).at(0));
            histos[17]->Fill(_GENpsi2s_mu_eta->at(psi2sId).at(1));
            for(int k = _GENjpsi_mass->size() - 1; k > 0; k--) {
                histos[2]->Fill(_GENjpsi_mass->at(k));
                histos[6]->Fill(_GENjpsi_pt->at(k));
                histos[10]->Fill(_GENjpsi_eta->at(k));
                histos[14]->Fill(_GENjpsi_mu_pt->at(k).at(0));
                histos[14]->Fill(_GENjpsi_mu_pt->at(k).at(1));
                histos[18]->Fill(_GENjpsi_mu_eta->at(k).at(0));
                histos[18]->Fill(_GENjpsi_mu_eta->at(k).at(1));
            }
            for(int k = _GENpsi2s_mass->size() - 1; k > 0; k--) {
                histos[3]->Fill(_GENpsi2s_mass->at(k));
                histos[7]->Fill(_GENpsi2s_pt->at(k));
                histos[11]->Fill(_GENpsi2s_eta->at(k));
                histos[15]->Fill(_GENpsi2s_mu_pt->at(k).at(0));
                histos[15]->Fill(_GENpsi2s_mu_pt->at(k).at(1));
                histos[19]->Fill(_GENpsi2s_mu_eta->at(k).at(0));
                histos[19]->Fill(_GENpsi2s_mu_eta->at(k).at(1));
            }
        }
    }
    // save histograms
    for(int i = 0; i < 5; i++) {
        TCanvas *c = new TCanvas();
        c->Divide(2, 2);
        for(int j = 0; j < 4; j++) {
            c->cd(j + 1);
            histos[4 * i + j]->Draw();
        }
        c->SaveAs(("HuntOutput_" + std::to_string(i) + ".pdf").c_str());
    }

    // TCanvas *c = new TCanvas();
    // c->Divide(2, 2);
    // c->cd(1);
    // gPad->SetLogy();
    // t->Draw("GENjpsi_mass[0]");
    // t->Draw("GENjpsi_mass[1]>>h1", "GENjpsi_ins > 1", "same");
    // TH1D *h1 = (TH1D *)gROOT->FindObject("h1");
    // h1->SetLineColor(kRed);
    // c->SaveAs("HuntOutput.pdf");
    return;
}