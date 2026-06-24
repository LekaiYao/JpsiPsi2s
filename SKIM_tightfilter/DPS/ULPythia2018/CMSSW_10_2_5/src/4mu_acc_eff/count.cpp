#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    string prefix = "pTHat4/Ntuple_2018_DPS_";
    for(int j = 1; j <= 40; j++) filenames.push_back(prefix + to_string(j) + ".root");
}

vector<Double_t> eff_pt, eff_y;
Double_t **nBin_Jpsi, **nVtx_Jpsi, **nBin_psi2S, **nVtx_psi2S, **nVtx_evt, **nTrg_evt;
Double_t calWeight(Double_t Jpsi_pt, Double_t Jpsi_y, Double_t psi2S_pt, Double_t psi2S_y) {
    int i = upper_bound(eff_y.begin(), eff_y.end(), Jpsi_y) - eff_y.begin() - 1;
    int j = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt) - eff_pt.begin() - 1;
    int k = upper_bound(eff_y.begin(), eff_y.end(), psi2S_y) - eff_y.begin() - 1;
    int l = upper_bound(eff_pt.begin(), eff_pt.end(), psi2S_pt) - eff_pt.begin() - 1;
    Double_t w = 1;
    if(nVtx_Jpsi[i][j] != 0) w *= nBin_Jpsi[i][j] / nVtx_Jpsi[i][j];
    if(nVtx_psi2S[k][l] != 0) w *= nBin_psi2S[k][l] / nVtx_psi2S[k][l];
    if(nTrg_evt[l][j] != 0) w *= nVtx_evt[l][j] / nTrg_evt[l][j];
    return w;
}
void loadEff() {
    string line;
    // Save acc in arrays
    ifstream effFile("plot/efficiency.txt");
    // ifstream effFile("/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/efficiency_SPS+2DPS.txt");
    if(!effFile.is_open()) return;
    int eff_ptBin = 0, eff_yBin = 0, lineCnt = 0;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            // Double_t x, y;
            // iss>>eff_ptBin>>eff_yBin>>x>>y;
            iss>>eff_ptBin>>eff_yBin;
            nBin_Jpsi = new Double_t*[eff_yBin];
            nVtx_Jpsi = new Double_t*[eff_yBin];
            nBin_psi2S = new Double_t*[eff_yBin];
            nVtx_psi2S = new Double_t*[eff_yBin];
            nVtx_evt = new Double_t*[eff_ptBin];
            nTrg_evt = new Double_t*[eff_ptBin];
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= eff_ptBin; i++) {
                iss>>pt;
                eff_pt.push_back(pt);
            }
        }else if(lineCnt == 2) {
            Double_t y;
            for(int i = 0; i <= eff_yBin; i++) {
                iss>>y;
                eff_y.push_back(y);
            }
        }else if(lineCnt <= 2 + eff_yBin) {
            int i = lineCnt - 3;
            nVtx_Jpsi[i] = new Double_t[eff_ptBin];
            nBin_Jpsi[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j];
        }else if(lineCnt <= 2 + 2 * eff_yBin) {
            int i = lineCnt - 3 - eff_yBin;
            nVtx_psi2S[i] = new Double_t[eff_ptBin];
            nBin_psi2S[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_psi2S[i][j]>>nBin_psi2S[i][j];
        }else if(lineCnt <= 2 + 2 * eff_yBin + eff_ptBin) {
            int i = lineCnt - 3 - 2 * eff_yBin;
            nVtx_evt[i] = new Double_t[eff_ptBin];
            nTrg_evt[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nTrg_evt[i][j]>>nVtx_evt[i][j];
        }else break;
        lineCnt++;
    }
    effFile.close();
    return;
}

void count() {
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("rootuple/oniaTree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    loadEff();
    // loop on tree entries
    UChar_t GEevt_valid = 0, GEevt_passAcc = 0;
    Double_t GEJpsi_pt = 0, GEJpsi_y = 0, GEpsi2S_pt = 0, GEpsi2S_y = 0, GEevt_fourMuMass = 0;
    vector<Double_t> *REmu_pt = 0, *REmu_eta = 0;
    vector<Double_t> *REJpsi_pt = 0, *REJpsi_y = 0, *REpsi2S_pt = 0, *REpsi2S_y = 0;
    vector<int> *REJpsi_muId1 = 0, *REpsi2S_muId1 = 0;
    vector<int> *REevt_JpsiId = 0, *REevt_psi2SId = 0;
    vector<bool> *REevt_matchTrg = 0;
    ch->SetBranchAddress("GEevt_valid", &GEevt_valid);
    ch->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
    ch->SetBranchAddress("GEJpsi_pt", &GEJpsi_pt);
    ch->SetBranchAddress("GEJpsi_y", &GEJpsi_y);
    ch->SetBranchAddress("GEpsi2S_pt", &GEpsi2S_pt);
    ch->SetBranchAddress("GEpsi2S_y", &GEpsi2S_y);
    ch->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
    ch->SetBranchAddress("REmu_pt", &REmu_pt);
    ch->SetBranchAddress("REmu_eta", &REmu_eta);
    ch->SetBranchAddress("REJpsi_pt", &REJpsi_pt);
    ch->SetBranchAddress("REJpsi_y", &REJpsi_y);
    ch->SetBranchAddress("REJpsi_muId1", &REJpsi_muId1);
    ch->SetBranchAddress("REpsi2S_pt", &REpsi2S_pt);
    ch->SetBranchAddress("REpsi2S_y", &REpsi2S_y);
    ch->SetBranchAddress("REpsi2S_muId1", &REpsi2S_muId1);
    ch->SetBranchAddress("REevt_JpsiId", &REevt_JpsiId);
    ch->SetBranchAddress("REevt_psi2SId", &REevt_psi2SId);
    ch->SetBranchAddress("REevt_matchTrg", &REevt_matchTrg);
    Int_t nEntry = ch->GetEntries(), nPassAcc = 0, nMatchTrg = 0;
    Double_t totWeight = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(!GEevt_valid || !GEevt_passAcc) continue;
        if(GEJpsi_pt < 10 || GEJpsi_pt > 40 || GEpsi2S_pt < 10 || GEpsi2S_pt > 40) continue;
        if(fabs(GEJpsi_y) > 2 || fabs(GEpsi2S_y) > 2) continue;
        if(GEevt_fourMuMass < 7.5) continue;
        nPassAcc++;
        for(int j = (int)REevt_matchTrg->size() - 1; j >= 0; j--) {
            if(!REevt_matchTrg->at(j)) continue;
            int JpsiId = REevt_JpsiId->at(j), psi2SId = REevt_psi2SId->at(j);
            if(REJpsi_pt->at(JpsiId) > 40 || REJpsi_pt->at(JpsiId) < 10) continue;
            if(REpsi2S_pt->at(psi2SId) > 40 || REpsi2S_pt->at(psi2SId) < 10) continue;
            if(REmu_pt->at(REJpsi_muId1->at(JpsiId)) < 5.5 && REmu_pt->at(REpsi2S_muId1->at(psi2SId)) < 5.5) continue;
            nMatchTrg++;
            totWeight += calWeight(REJpsi_pt->at(JpsiId), REJpsi_y->at(JpsiId), REpsi2S_pt->at(psi2SId), REpsi2S_y->at(psi2SId));
            break;
        }
    }
    cout<<"\nnEvent="<<nPassAcc<<"\nnPassAcc="<<nMatchTrg<<"\ntotWeight="<<totWeight<<endl;
    return;
}