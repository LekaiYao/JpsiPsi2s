// Author: Shiyang CHEN
// Description: Kinematics investigation of GEN/RECO level DPS sample
// Implementation: 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
using namespace std;
// #define RECO

void loadFile(vector<string>& filenames) {
    // // No filter DPS datasets - GEN level
    // string infix = "/DPS_2018_Psi2SJ_nofilter_";
    // for(int i = 1; i <= 10; i++) filenames.push_back("1" + infix + to_string(i) + ".root");
    // for(int i = 11; i <= 20; i++) filenames.push_back("2" + infix + to_string(i) + ".root");
    // for(int i = 21; i <= 30; i++) filenames.push_back("3" + infix + to_string(i) + ".root");
    // for(int i = 31; i <= 45; i++) filenames.push_back("4" + infix + to_string(i) + ".root");
    // for(int i = 46; i <= 60; i++) filenames.push_back("5" + infix + to_string(i) + ".root");

    // // Filter DPS datasets without pTHat - RECO level
    // string infix = "/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/Ntuple/Ntuple_2018_DPS_";
    // for(int i = 1; i <= 25; i++) filenames.push_back(infix + to_string(i) + ".root");

    // // Filter DPS datasets with pTHat=2 - GEN level
    // string infix = "pTHat2/DPS_2018_Psi2SJ_HO_";
    // for(int i = 1; i <= 10; i++) filenames.push_back(infix + to_string(i) + ".root");

    // // Filter DPS datasets with pTHat=4 - GEN level
    // string infix = "pTHat4/DPS_2018_Psi2SJ_HO_";
    // for(int i = 1; i <= 10; i++) filenames.push_back(infix + to_string(i) + ".root");

    // GEN only DPS - GEN level
    string prefix = "8/DPS_2018_Psi2SJ_";
    for(int i = 1; i <= 10; i++) {
        if(i == 6) continue;
        filenames.push_back(prefix + to_string(i) + ".root");
    }
    for(int i = 1; i <= 5; i++) {
        prefix = "8/" + to_string(i) + "/DPS_2018_Psi2SJ_";
        for(int j = 1; j <= 20; j++) filenames.push_back(prefix + to_string(j) + ".root");
    }

    // // Event mixing DPS - GEN level
    // filenames.push_back("MixedDPS.root");
}
void vs() {
    vector<string> filenames;
    loadFile(filenames);
    #ifdef RECO
    TChain *chain = new TChain("rootuple/oniaTree");
    #else
    TChain *chain = new TChain("GenAnalyzer/gen_tree");
    // TChain *chain = new TChain("data");
    #endif
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        chain->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    // Define histograms
    TH1D Jpsi_mu_pt("Jpsi_mu_pt", "Jpsi_mu_pt", 50, 3.3, 9.3);
    TH1D Jpsi_mu_eta("Jpsi_mu_eta", "Jpsi_mu_eta", 50, -2.4, 2.4);
    TH1D psi2S_mu_pt("psi2S_mu_pt", "psi2S_mu_pt", 50, 3.3, 9.3);
    TH1D psi2S_mu_eta("psi2S_mu_eta", "psi2S_mu_eta", 50, -2.4, 2.4);// 1~4
    TH1D Jpsi_pt("Jpsi_pt", "Jpsi_pt", 50, 10, 20);
    TH1D Jpsi_y("Jpsi_y", "Jpsi_y", 50, -2.4, 2.4);
    // TH1D Jpsi_eta("Jpsi_eta", "Jpsi_eta", 50, -2.4, 2.4);
    TH1D psi2S_pt("psi2S_pt", "psi2S_pt", 50, 10, 20);
    TH1D psi2S_y("psi2S_y", "psi2S_y", 50, -2.4, 2.4);// 5~8
    // TH1D psi2S_eta("psi2S_eta", "psi2S_eta", 50, -2.4, 2.4);// 5~8
    TH2D Jpsi_mu_pt_2D("Jpsi_mu_pt_2D", "Jpsi_mu_pt_2D", 50, 3.3, 9.3, 50, 3.3, 9.3);
    TH2D psi2S_mu_pt_2D("psi2S_mu_pt_2D", "psi2S_mu_pt_2D", 50, 3.3, 9.3, 50, 3.3, 9.3);
    TH2D Jpsi_psi2S_pt_2D("Jpsi_psi2S_pt_2D", "Jpsi_psi2S_pt_2D", 50, 10, 20, 50, 10, 20);
    TH1D Delta_y("Delta_y", "Delta_y", 50, 0, 2.5);// 9~12
    // TH1D Delta_eta("Delta_eta", "Delta_eta", 50, 0, 2.5);// 9~12
    // Fill histograms
    #ifdef RECO
    string cut = "GEevt_valid && GEmu_pt[0] > 3.3 && GEmu_pt[1] > 3.3 && GEmu_pt[2] > 3.3 && GEmu_pt[3] > 3.3 && "
        "fabs(GEmu_eta[0]) < 2.5 && fabs(GEmu_eta[1]) < 2.5 && fabs(GEmu_eta[2]) < 2.5 && fabs(GEmu_eta[3]) < 2.5 && "
        "GEJpsi_pt > 10 && GEJpsi_pt < 40 && GEpsi2S_pt > 10 && GEpsi2S_pt < 40 && "
        "fabs(GEJpsi_eta) < 2 && fabs(GEpsi2S_eta) < 2";
    chain->Draw("GEmu_pt[0] >> Jpsi_mu_pt", cut.c_str());chain->Draw("GEmu_pt[1] >> Jpsi_mu_pt", cut.c_str());
    chain->Draw("GEmu_eta[0] >> Jpsi_mu_eta", cut.c_str());chain->Draw("GEmu_eta[1] >> Jpsi_mu_eta", cut.c_str());
    chain->Draw("GEmu_pt[2] >> psi2S_mu_pt", cut.c_str());chain->Draw("GEmu_pt[3] >> psi2S_mu_pt", cut.c_str());
    chain->Draw("GEmu_eta[2] >> psi2S_mu_eta", cut.c_str());chain->Draw("GEmu_eta[3] >> psi2S_mu_eta", cut.c_str());// 1~4
    chain->Draw("GEJpsi_pt >> Jpsi_pt", cut.c_str());
    chain->Draw("GEJpsi_eta >> Jpsi_eta", cut.c_str());
    chain->Draw("GEpsi2S_pt >> psi2S_pt", cut.c_str());
    chain->Draw("GEpsi2S_eta >> psi2S_eta", cut.c_str());// 5~8
    chain->Draw("GEmu_pt[0]:GEmu_pt[1] >> Jpsi_mu_pt_2D", cut.c_str());
    chain->Draw("GEmu_pt[2]:GEmu_pt[3] >> psi2S_mu_pt_2D", cut.c_str());
    chain->Draw("GEJpsi_pt:GEpsi2S_pt >> Jpsi_psi2S_pt_2D", cut.c_str());
    chain->Draw("fabs(GEJpsi_eta-GEpsi2S_eta) >> Delta_eta", cut.c_str());// 9~12
    #else
    // string cut = "";
    string cut = 
        // "GENjpsi_mu_pt[0][0] > 3.3 && GENjpsi_mu_pt[0][1] > 3.3 && GENpsi2s_mu_pt[0][0] > 3.3 && GENpsi2s_mu_pt[0][1] > 3.3 && "
        // "fabs(GENjpsi_mu_eta[0][0]) < 2.5 && fabs(GENjpsi_mu_eta[0][1]) < 2.5 && fabs(GENpsi2s_mu_eta[0][0]) < 2.5 && fabs(GENpsi2s_mu_eta[0][1]) < 2.5 && "
        "GENjpsi_pt[0] > 10 && GENjpsi_pt[0] < 40 && GENpsi2s_pt[0] > 10 && GENpsi2s_pt[0] < 40 && "
        "fabs(GENjpsi_y[0]) < 2 && fabs(GENpsi2s_y[0]) < 2";
        // "fabs(GENjpsi_eta[0]) < 2 && fabs(GENpsi2s_eta[0]) < 2";
    chain->Draw("GENjpsi_mu_pt[0] >> Jpsi_mu_pt", cut.c_str());
    chain->Draw("GENjpsi_mu_eta[0] >> Jpsi_mu_eta", cut.c_str());
    chain->Draw("GENpsi2s_mu_pt[0] >> psi2S_mu_pt", cut.c_str());
    chain->Draw("GENpsi2s_mu_eta[0] >> psi2S_mu_eta", cut.c_str());// 1~4
    chain->Draw("GENjpsi_pt[0] >> Jpsi_pt", cut.c_str());
    chain->Draw("GENjpsi_y[0] >> Jpsi_y", cut.c_str());
    // chain->Draw("GENjpsi_eta[0] >> Jpsi_eta", cut.c_str());
    chain->Draw("GENpsi2s_pt[0] >> psi2S_pt", cut.c_str());
    chain->Draw("GENpsi2s_y[0] >> psi2S_y", cut.c_str());// 5~8
    // chain->Draw("GENpsi2s_eta[0] >> psi2S_eta", cut.c_str());// 5~8
    chain->Draw("GENjpsi_mu_pt[0][0]:GENjpsi_mu_pt[0][1] >> Jpsi_mu_pt_2D", cut.c_str());
    chain->Draw("GENpsi2s_mu_pt[0][0]:GENpsi2s_mu_pt[0][1] >> psi2S_mu_pt_2D", cut.c_str());
    chain->Draw("GENjpsi_pt[0]:GENpsi2s_pt[0] >> Jpsi_psi2S_pt_2D", cut.c_str());
    chain->Draw("fabs(GENjpsi_y[0]-GENpsi2s_y[0]) >> Delta_y", cut.c_str());// 9~12
    // chain->Draw("fabs(GENjpsi_eta[0]-GENpsi2s_eta[0]) >> Delta_eta", cut.c_str());// 9~12
    #endif
    // Draw histograms
    TCanvas *outCanvas = new TCanvas("Kinematics", "Kinematics", 3200, 2400);
    outCanvas->Divide(4, 3);
    outCanvas->cd(1);
    Jpsi_mu_pt.Draw();
    outCanvas->cd(2);
    Jpsi_mu_eta.Draw();
    outCanvas->cd(3);
    psi2S_mu_pt.Draw();
    outCanvas->cd(4);
    psi2S_mu_eta.Draw();// 1~4
    outCanvas->cd(5);
    Jpsi_pt.Draw();
    outCanvas->cd(6);
    Jpsi_y.Draw();
    // Jpsi_eta.Draw();
    outCanvas->cd(7);
    psi2S_pt.Draw();
    outCanvas->cd(8);
    psi2S_y.Draw();// 5~8
    // psi2S_eta.Draw();// 5~8
    outCanvas->cd(9);
    Jpsi_mu_pt_2D.Draw();
    outCanvas->cd(10);
    psi2S_mu_pt_2D.Draw();
    outCanvas->cd(11);
    Jpsi_psi2S_pt_2D.Draw();
    outCanvas->cd(12);
    Delta_y.Draw();// 9~12
    // Delta_eta.Draw();// 9~12
    outCanvas->SaveAs("Kine_nofilter_cut.png");
    delete outCanvas;
    return;
}
