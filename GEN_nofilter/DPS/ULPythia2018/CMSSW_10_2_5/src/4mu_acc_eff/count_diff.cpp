#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;
#define PI 3.14159265359

void loadFile(vector<string>& filenames) {
    filenames.push_back("MixedDPS.root");
}
vector<Double_t> acc_pt, acc_y;
Double_t **nGen_Jpsi, **nAcc_Jpsi, **nGen_psi2S, **nAcc_psi2S;
Double_t calWeight(Double_t Jpsi_pt, Double_t Jpsi_y, Double_t psi2S_pt, Double_t psi2S_y) {
    int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
    int k = upper_bound(acc_y.begin(), acc_y.end(), psi2S_y) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), psi2S_pt) - acc_pt.begin() - 1;// k-psi(2S) y, l-psi(2S) pT
    Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_psi2S[k][l] / nAcc_psi2S[k][l];
    return w;
}
void loadAcc(bool statis) {
    string line;
    Double_t nJpsi, npsi2S;
    // Save acc in arrays
    // ifstream accFile("/eos/home-c/chensh/JPsiPsi2s/GEN_nofilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/acceptance_SPS+2DPS.txt");
    // ifstream accFile("plot/acceptance.txt");
    ifstream accFile("/eos/home-c/chensh/JPsiPsi2s/GEN_nofilter/SPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/plot/acceptance.txt");
    if(!accFile.is_open()) return;
    int acc_ptBin = 0, acc_yBin = 0, lineCnt = 0;
    while(getline(accFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>acc_ptBin>>acc_yBin;
            nGen_Jpsi = new Double_t*[acc_yBin];
            nAcc_Jpsi = new Double_t*[acc_yBin];
            nGen_psi2S = new Double_t*[acc_yBin];
            nAcc_psi2S = new Double_t*[acc_yBin];
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= acc_ptBin; i++) {
                iss>>pt;
                acc_pt.push_back(pt);
            }
        }else if(lineCnt == 2) {
            Double_t y;
            for(int i = 0; i <= acc_yBin; i++) {
                iss>>y;
                acc_y.push_back(y);
            }
        }else if(lineCnt <= 2 + acc_yBin) {
            int i = lineCnt - 3;
            nAcc_Jpsi[i] = new Double_t[acc_ptBin];
            nGen_Jpsi[i] = new Double_t[acc_ptBin];
            for(int j = 0; j < acc_ptBin; j++) {
                iss>>nAcc_Jpsi[i][j]>>nGen_Jpsi[i][j];
                nJpsi += nGen_Jpsi[i][j];
            }
        }else if(lineCnt <= 2 + 2 * acc_yBin) {
            int i = lineCnt - 3 - acc_yBin;
            nAcc_psi2S[i] = new Double_t[acc_ptBin];
            nGen_psi2S[i] = new Double_t[acc_ptBin];
            for(int j = 0; j < acc_ptBin; j++) {
                iss>>nAcc_psi2S[i][j]>>nGen_psi2S[i][j];
                npsi2S += nGen_psi2S[i][j];
            }
        }else break;
        lineCnt++;
    }
    accFile.close();
    if(statis) cout<<"Input acc map: "<<nJpsi<<' '<<npsi2S<<' '<<(nJpsi == npsi2S)<<endl;
    return;
}

void count_diff() {
    // Initialization
    int nVar = 5, nBin[] = {5, 5, 6, 6, 5}, *xBin = new int[nVar];
    string varName[] = {"evt_mass", "evt_y", "evt_pt", "delta_y", "delta_phi"};
    vector<vector<double>> binning = {
        {7.5, 17.5, 27.5, 37.5, 47.5, 107.5},
        {0, 0.4, 0.8, 1.2, 1.6, 2.0},
        {0, 8, 16, 24, 32, 40, 80},
        {0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0},
        {0, 0.2 * PI, 0.4 * PI, 0.6 * PI, 0.8 * PI, PI}
    };
    vector<vector<int>> nEvent, nPassAcc;
    vector<vector<double>> totWeight;
    for(int i = 0; i < nVar; i++) {
        nEvent.push_back(vector<int>());
        nPassAcc.push_back(vector<int>());
        totWeight.push_back(vector<double>());
        for(int j = 0; j < nBin[i]; j++) {
            nEvent[i].push_back(0);
            nPassAcc[i].push_back(0);
            totWeight[i].push_back(0);
        }
    }
    int nEvent0 = 0, nPassAcc0 = 0;
    double totWeight0 = 0;
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("data");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    loadAcc(false);
    // loop on tree entries
    vector<Double_t> *Jpsi_pt = 0, *Jpsi_y = 0, *psi2S_pt = 0, *psi2S_y = 0;
    Double_t *var = new Double_t[nVar];
    bool evt_acc;
    ch->SetBranchAddress("GENjpsi_pt", &Jpsi_pt);
    ch->SetBranchAddress("GENjpsi_y", &Jpsi_y);
    ch->SetBranchAddress("GENpsi2s_pt", &psi2S_pt);
    ch->SetBranchAddress("GENpsi2s_y", &psi2S_y);
    for(int i = 0; i < nVar; i++) ch->SetBranchAddress(varName[i].c_str(), &var[i]);
    ch->SetBranchAddress("evt_acc", &evt_acc);
    Int_t nEntry = ch->GetEntries();
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        if(!(i & 511)) cout<<"Processing: "<<(i+1)<<"th entry\r";
        // if(Jpsi_pt->empty() || psi2S_pt->empty()) continue;
        // if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || psi2S_pt->at(0) < 10 || psi2S_pt->at(0) > 40) continue;
        // if(fabs(Jpsi_y->at(0)) > 2 || fabs(psi2S_y->at(0)) > 2) continue;
        // if(evt_mass->at(0) < 7.5) continue;
        nEvent0++;
        if(evt_acc) {
            nPassAcc0++;
            totWeight0 += calWeight(Jpsi_pt->at(0), Jpsi_y->at(0), psi2S_pt->at(0), psi2S_y->at(0));
        }
        for(int j = 0; j < nVar; j++) {
            xBin[j] = upper_bound(binning[j].begin(), binning[j].end(), var[j]) - binning[j].begin() - 1;
            nEvent[j][xBin[j]]++;
            if(!evt_acc) continue;
            nPassAcc[j][xBin[j]]++;
            totWeight[j][xBin[j]] += calWeight(Jpsi_pt->at(0), Jpsi_y->at(0), psi2S_pt->at(0), psi2S_y->at(0));
        }
    }
    cout<<"\nnEvent="<<nEvent0<<"\nnPassAcc="<<nPassAcc0<<"\ntotWeight="<<totWeight0<<endl;
    for(int i = 0; i < nVar; i++) {
        cout<<varName[i]<<' ';
        for(int j = 0; j < nBin[i]; j++) cout<<(nEvent[i][j] - totWeight[i][j]) / totWeight[i][j]<<" ";
        cout<<endl;
    }
    return;
}