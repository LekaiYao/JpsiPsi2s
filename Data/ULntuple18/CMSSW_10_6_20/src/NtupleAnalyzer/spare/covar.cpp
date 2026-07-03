// Author: Shiyang CHEN
// Description: Covariance check for DPS Monte-Carlo
// Implementation: This program should be run with ROOT no lower than v6.30
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TScatter.h"
using namespace std;
// I/O settings area
// DPS
#define N_DIR 1
string prefix[N_DIR] = {""};
string infix = "/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/Ntuple/Ntuple_2018_DPS";
int suffix[N_DIR] = {10};
string type = "DPS";
// SPS
// #define N_DIR 1
// string prefix[N_DIR] = {""};
// string infix = "/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/Ntuple/Ntuple_2018_SPS";
// int suffix[N_DIR] = {1};
// string type = "SPS";

class Process {
    private:
    vector<Double_t> Jpsi_pt, Jpsi_eta, Jpsi_phi;
    vector<Double_t> psi2S_pt, psi2S_eta, psi2S_phi;
    void readTree(string& fileName) {
        // Handle root file
        TFile file(fileName.c_str(), "READ");
        TTree *inTree = (TTree *)file.Get("rootuple/oniaTree");
        int nEvent = inTree->GetEntries();
        // Define input tree variables
        Double_t GEJpsi_pt = 0, GEJpsi_eta = 0, GEJpsi_phi = 0;
        Double_t GEpsi2S_pt = 0, GEpsi2S_eta = 0, GEpsi2S_phi = 0;
        Double_t GEevt_fourMuMass = 0;
        UChar_t GEevt_valid = 0, GEevt_passAcc = 0, GEJpsi_passID = 0, GEpsi2S_passID = 0, GEJpsi_passVtx = 0, GEpsi2S_passVtx = 0,
        GEevt_passHLT = 0, GEevt_matchTrg = 0;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("GEJpsi_pt", &GEJpsi_pt);
        inTree->SetBranchAddress("GEJpsi_eta", &GEJpsi_eta);
        inTree->SetBranchAddress("GEJpsi_phi", &GEJpsi_phi);
        inTree->SetBranchAddress("GEpsi2S_pt", &GEpsi2S_pt);
        inTree->SetBranchAddress("GEpsi2S_eta", &GEpsi2S_eta);
        inTree->SetBranchAddress("GEpsi2S_phi", &GEpsi2S_phi);
        inTree->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
        inTree->SetBranchAddress("GEevt_valid", &GEevt_valid);
        inTree->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
        inTree->SetBranchAddress("GEJpsi_passID", &GEJpsi_passID);
        inTree->SetBranchAddress("GEpsi2S_passID", &GEpsi2S_passID);
        inTree->SetBranchAddress("GEJpsi_passVtx", &GEJpsi_passVtx);
        inTree->SetBranchAddress("GEpsi2S_passVtx", &GEpsi2S_passVtx);
        inTree->SetBranchAddress("GEevt_passHLT", &GEevt_passHLT);
        inTree->SetBranchAddress("GEevt_matchTrg", &GEevt_matchTrg);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            if(!GEevt_valid) continue;
            if(!GEevt_passAcc) continue;
            // if(!GEJpsi_passID || !GEpsi2S_passID) continue;
            // if(!GEJpsi_passVtx || !GEpsi2S_passVtx) continue;
            if(fabs(GEJpsi_eta) > 2 || fabs(GEpsi2S_eta) > 2) continue;
            if(GEJpsi_pt < 10 || GEpsi2S_pt < 10 || GEJpsi_pt > 40 || GEpsi2S_pt > 40) continue;
            if(GEevt_fourMuMass < 7.5) continue;
            // if(!GEevt_passHLT) continue;
            // bool matchTrg = false;
            // for(int j = 0; j < (int)REevt_matchTrg->size(); j++) {
            //     if(REevt_matchTrg->at(j)) {
            //         matchTrg = true;
            //         break;
            //     }
            // }
            // if(!matchTrg) continue;
            // if(!GEevt_matchTrg) continue;
            Jpsi_pt.push_back(GEJpsi_pt);
            Jpsi_eta.push_back(GEJpsi_eta);
            Jpsi_phi.push_back(GEJpsi_phi);
            psi2S_pt.push_back(GEpsi2S_pt);
            psi2S_eta.push_back(GEpsi2S_eta);
            psi2S_phi.push_back(GEpsi2S_phi);
            hPt.Fill(GEJpsi_pt, GEpsi2S_pt, 1);
            hEta.Fill(GEJpsi_eta, GEpsi2S_eta, 1);
            hPhi.Fill(GEJpsi_phi, GEpsi2S_phi, 1);
            hDEta.Fill(fabs(GEJpsi_eta - GEpsi2S_eta));
            hDPhi.Fill(fabs(GEJpsi_phi - GEpsi2S_phi));
            totEvent++;
        }
        sPt = new TScatter(totEvent, &Jpsi_pt[0], &psi2S_pt[0]);
        sEta = new TScatter(totEvent, &Jpsi_eta[0], &psi2S_eta[0]);
        sPhi = new TScatter(totEvent, &Jpsi_phi[0], &psi2S_phi[0]);
        corr_pt = hPt.GetCorrelationFactor();
        corr_eta = hEta.GetCorrelationFactor();
        corr_phi = hPhi.GetCorrelationFactor();
        cout<<"Finished.  >> "<<endl;
        return;
    }
    public:
    int totEvent = 0;
    TScatter *sPt, *sEta, *sPhi;
    TH2D hPt = TH2D("hPt", "hPt", 10, 0., 40., 10, 0., 40.);
    TH2D hEta = TH2D("hEta", "hEta", 10, -2.4, 2.4, 10, -2.4, 2.4);
    TH2D hPhi = TH2D("hPhi", "hPhi", 10, -3.15, 3.15, 10, -3.15, 3.15);
    TH1D hDEta = TH1D("hDEta", "hDEta", 6, 0, 2.4);
    TH1D hDPhi = TH1D("hDPhi", "hDPhi", 5, 0, 3);
    Double_t corr_pt, corr_eta, corr_phi;
    void loopOn() {
        for(int i = 0; i < N_DIR; i++) {
            for(int j = 1; j <= suffix[i]; j++) {
                string fileName = prefix[i] + infix + prefix[i] + "_" + to_string(j) + ".root";
                readTree(fileName);
            }
        }
    }
    void showCovar() {
        cout<<"Covariance in pT: "<<corr_pt<<endl;
        cout<<"Covariance in eta: "<<corr_eta<<endl;
        cout<<"Covariance in phi: "<<corr_phi<<endl;
        return;
    }
    void saveFig() {
        TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
        sPt->SetTitle((type + " Pt Scatter corr=" + to_string(corr_pt) + ";p_T(J/#psi);p_T(#psi(2S))").c_str());
        sPt->GetXaxis()->SetRangeUser(0., 40.);
        sPt->GetYaxis()->SetRangeUser(0., 40.);
        sPt->SetMarkerStyle(8);
        sPt->SetMarkerSize(0.3);
        sPt->Draw("scat");
        canvas->SaveAs(("plots/Covar_pt_" + type + ".png").c_str());
        delete canvas;
        canvas = new TCanvas("canvas", "canvas", 1000, 1000);
        sEta->SetTitle((type + " Eta Scatter corr=" + to_string(corr_eta) + ";#eta(J/#psi);#eta(#psi(2S))").c_str());
        sEta->GetXaxis()->SetRangeUser(-2.4, 2.4);
        sEta->GetYaxis()->SetRangeUser(-2.4, 2.4);
        sEta->SetMarkerStyle(8);
        sEta->SetMarkerSize(0.3);
        sEta->Draw("scat");
        canvas->SaveAs(("plots/Covar_eta_" + type + ".png").c_str());
        delete canvas;
        canvas = new TCanvas("canvas", "canvas", 1000, 1000);
        sPhi->SetTitle((type + " Phi Scatter corr=" + to_string(corr_phi) + ";#phi(J/#psi);#phi(#psi(2S))").c_str());
        sPhi->GetXaxis()->SetRangeUser(-3.15, 3.15);
        sPhi->GetYaxis()->SetRangeUser(-3.15, 3.15);
        sPhi->SetMarkerStyle(8);
        sPhi->SetMarkerSize(0.3);
        sPhi->Draw("scat");
        canvas->SaveAs(("plots/Covar_phi_" + type + ".png").c_str());
        delete canvas;
        canvas = new TCanvas("canvas", "canvas", 1000, 1000);
        hDEta.SetTitle((type + " DEta GEN level;#Delta#eta").c_str());
        hDEta.GetYaxis()->SetRangeUser(0, 1.1 * hDEta.GetMaximum());
        hDEta.Draw();
        canvas->SaveAs(("plots/GEN_DEta_" + type + ".png").c_str());
        delete canvas;
        canvas = new TCanvas("canvas", "canvas", 1000, 1000);
        hDPhi.SetTitle((type + " DPhi GEN level;#Delta#phi").c_str());
        hDPhi.GetYaxis()->SetRangeUser(0, 1.1 * hDPhi.GetMaximum());
        hDPhi.Draw();
        canvas->SaveAs(("plots/GEN_DPhi_" + type + ".png").c_str());
        delete canvas;
    }
};

void covar() {
    Process process;
    process.loopOn();
    process.showCovar();
    process.saveFig();
    cout<<"Valid event: "<<process.totEvent<<endl;
    return;
}
