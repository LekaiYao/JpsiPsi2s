#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    // SPS
    // filenames.push_back("WeightSPS.root");
    // DPS
    // filenames.push_back("WeightDPS.root");
    // Bdecay
    filenames.push_back("WeightBdecay.root");
}

void corr() {
    // Configuration
    int nVar = 4;
    string varName[] = {"Jpsi_mass", "psi2S_mass", "Jpsi_ctau", "psi2S_ctau"};
    Double_t *var = new Double_t[nVar], w;
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("data");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    // loop on tree entries
    for(int i = 0; i < nVar; i++) ch->SetBranchAddress(varName[i].c_str(), &(var[i]));
    ch->SetBranchAddress("evt_weight", &w);
    // Unified version
    double *sum = new double[nVar], **sumt = new double*[nVar];
    for(int i = 0; i < nVar; i++) {
        sum[i] = 0;
        sumt[i] = new double[nVar];
        for(int j = 0; j < nVar; j++) sumt[i][j] = 0;
    }
    Int_t nEntry = ch->GetEntries();
    double sumw = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        for(int i = 0; i < nVar; i++) {
            sum[i] += var[i] * w;
            for(int j = 0; j <= i; j++) sumt[i][j] += var[i] * var[j] * w;
        }
        sumw += w;
    }
    for(int i = 0; i < nVar; i++) {
        sum[i] /= sumw;
        for(int j = 0; j <= i; j++) sumt[i][j] /= sumw;
    }
    cout<<endl<<"Total event number: "<<nEntry<<endl<<"Total event weight: "<<sumw<<endl;
    for(int i = 0; i < nVar; i++) {
        for(int j = 0; j <= i; j++) cout<<(sumt[i][j] - sum[i] * sum[j]) / sqrt((sumt[i][i] - sum[i] * sum[i]) * (sumt[j][j] - sum[j] * sum[j]))<<' ';
        cout<<endl;
    }
    return;
}
