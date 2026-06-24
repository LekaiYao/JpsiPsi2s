#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    string prefix = "8/DPS_2018_Psi2SJ_";
    for(int i = 1; i <= 10; i++) {
        if(i == 6) continue;
        filenames.push_back(prefix + to_string(i) + ".root");
    }
    for(int i = 1; i <= 5; i++) {
        prefix = "8/" + to_string(i) + "/DPS_2018_Psi2SJ_";
        for(int j = 1; j <= 20; j++) filenames.push_back(prefix + to_string(j) + ".root");
    }
}

void covar() {
    // Configuration
    int nVarJpsi = 2, nVarpsi2S = 2;
    string varJpsi[] = {"GENjpsi_pt", "GENjpsi_y"}, varpsi2S[] = {"GENpsi2s_pt", "GENpsi2s_y"};
    vector<Double_t> **brJpsi = new vector<Double_t>*[nVarJpsi], **brpsi2S = new vector<Double_t>*[nVarpsi2S];
    vector<Double_t> *evt_mass = 0;
    for(int i = 0; i < nVarJpsi; i++) brJpsi[i] = 0;
    for(int i = 0; i < nVarpsi2S; i++) brpsi2S[i] = 0;
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("GenAnalyzer/gen_tree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    // loop on tree entries
    for(int i = 0; i < nVarJpsi; i++) ch->SetBranchAddress(varJpsi[i].c_str(), &(brJpsi[i]));
    for(int i = 0; i < nVarpsi2S; i++) ch->SetBranchAddress(varpsi2S[i].c_str(), &(brpsi2S[i]));
    ch->SetBranchAddress("GENevt_mass", &evt_mass);
    // Unified version
    double *sumJpsi = new double[nVarJpsi], *sumpsi2S = new double[nVarpsi2S];
    double *sum2Jpsi = new double[nVarJpsi], *sum2psi2S = new double[nVarpsi2S];
    double *sumCovar = new double[nVarJpsi * nVarpsi2S];
    for(int j = 0; j < nVarJpsi; j++) {
        sumJpsi[j] = 0;
        sum2Jpsi[j] = 0;
    }
    for(int j = 0; j < nVarpsi2S; j++) {
        sumpsi2S[j] = 0;
        sum2psi2S[j] = 0;
    }
    for(int j = 0; j < nVarJpsi; j++) {
        for(int k = 0; k < nVarpsi2S; k++) sumCovar[j * nVarpsi2S + k] = 0;
    }
    Int_t nEntry = ch->GetEntries(), nEvent = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(brJpsi[0]->empty() || brpsi2S[0]->empty()) continue;
        if(brJpsi[0]->at(0) < 10 || brJpsi[0]->at(0) > 40 || brpsi2S[0]->at(0) < 10 || brpsi2S[0]->at(0) > 40) continue;
        if(fabs(brJpsi[1]->at(0)) > 2 || fabs(brpsi2S[1]->at(0)) > 2) continue;
        if(evt_mass->at(0) < 7.5) continue;
        for(int j = 0; j < nVarJpsi; j++) {
            sumJpsi[j] += brJpsi[j]->at(0);
            sum2Jpsi[j] += brJpsi[j]->at(0) * brJpsi[j]->at(0);
        }
        for(int j = 0; j < nVarpsi2S; j++) {
            sumpsi2S[j] += brpsi2S[j]->at(0);
            sum2psi2S[j] += brpsi2S[j]->at(0) * brpsi2S[j]->at(0);
        }
        for(int j = 0; j < nVarJpsi; j++) {
            for(int k = 0; k < nVarpsi2S; k++) sumCovar[j * nVarpsi2S + k] += brJpsi[j]->at(0) * brpsi2S[k]->at(0);
        }
        nEvent++;
    }
    for(int j = 0; j < nVarJpsi; j++) {
        sumJpsi[j] /= nEvent;
        sum2Jpsi[j] /= nEvent;
    }
    for(int j = 0; j < nVarpsi2S; j++) {
        sumpsi2S[j] /= nEvent;
        sum2psi2S[j] /= nEvent;
    }
    for(int j = 0; j < nVarJpsi; j++) {
        for(int k = 0; k < nVarpsi2S; k++) sumCovar[j * nVarpsi2S + k] /= nEvent;
    }
    cout<<endl<<"Total event number: "<<nEvent<<endl;
    for(int i = 0; i < nVarJpsi; i++) {
        for(int j = 0; j < nVarpsi2S; j++) cout<<(sumCovar[i * nVarpsi2S + j] - sumJpsi[i] * sumpsi2S[j]) / sqrt((sum2Jpsi[i] - sumJpsi[i] * sumJpsi[i]) * (sum2psi2S[j] - sumpsi2S[j] * sumpsi2S[j]))<<' ';
        cout<<endl;
    }
    return;
}
