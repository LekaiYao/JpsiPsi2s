#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
using namespace std;

void loadFile(vector<string>& filenames) {
    string prefix = "SPS_2018_Psi2SJ_";
    for(int i = 1; i <= 1; i++) filenames.push_back(prefix + to_string(i) + ".root");
}

vector<Double_t> acc_pt, acc_y;
Double_t **nGen_Jpsi, **nAcc_Jpsi, **nGen_psi2S, **nAcc_psi2S;
Double_t calWeight(Double_t Jpsi_pt, Double_t Jpsi_y, Double_t psi2S_pt, Double_t psi2S_y) {
    int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
    int k = upper_bound(acc_y.begin(), acc_y.end(), psi2S_y) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), psi2S_pt) - acc_pt.begin() - 1;// k-psi(2S) y, l-psi(2S) pT
    Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_psi2S[k][l] / nAcc_psi2S[k][l];
    return w;
}
Double_t calWeight_Jpsi(Double_t Jpsi_pt, Double_t Jpsi_y) {
    int i = upper_bound(acc_y.begin(), acc_y.end(), Jpsi_y) - acc_y.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt) - acc_pt.begin() - 1;// i-J/psi y, j-J/psi pT
    Double_t w = nGen_Jpsi[i][j] / nAcc_Jpsi[i][j];
    return w;
}
Double_t calWeight_psi2S(Double_t psi2S_pt, Double_t psi2S_y) {
    int k = upper_bound(acc_y.begin(), acc_y.end(), psi2S_y) - acc_y.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), psi2S_pt) - acc_pt.begin() - 1;// k-psi(2S) y, l-psi(2S) pT
    Double_t w = nGen_psi2S[k][l] / nAcc_psi2S[k][l];
    return w;
}
void loadAcc(bool statis) {
    string line;
    Double_t nJpsi, npsi2S;
    // Save acc in arrays
    // ifstream accFile("acceptance_SPS+2DPS.txt");
    ifstream accFile("plot/acceptance.txt");
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

void count() {
    // Handle input files and construct TChain
    vector<string> filenames;
    loadFile(filenames);
    TChain *ch = new TChain("GenAnalyzer/gen_tree");
    for(int i = 0; i < (int)filenames.size(); i++) {
        cout<<"Processing "<<(i + 1)<<"th file: "<<filenames[i];
        ch->Add(filenames[i].c_str());
        cout<<". Done!"<<'\n';
    }
    loadAcc(false);
    // loop on tree entries
    vector<Double_t> *Jpsi_pt = 0, *Jpsi_y = 0, *psi2S_pt = 0, *psi2S_y = 0;
    vector<vector<Double_t>> *Jpsi_mu_pt = 0, *Jpsi_mu_eta = 0, *psi2S_mu_pt = 0, *psi2S_mu_eta = 0;
    ch->SetBranchAddress("GENjpsi_pt", &Jpsi_pt);
    ch->SetBranchAddress("GENjpsi_y", &Jpsi_y);
    ch->SetBranchAddress("GENpsi2s_pt", &psi2S_pt);
    ch->SetBranchAddress("GENpsi2s_y", &psi2S_y);
    ch->SetBranchAddress("GENjpsi_mu_pt", &Jpsi_mu_pt);
    ch->SetBranchAddress("GENjpsi_mu_eta", &Jpsi_mu_eta);
    ch->SetBranchAddress("GENpsi2s_mu_pt", &psi2S_mu_pt);
    ch->SetBranchAddress("GENpsi2s_mu_eta", &psi2S_mu_eta);
    // Unified version
    Int_t nEntry = ch->GetEntries(), nEvent = 0, nPassAcc = 0;
    Double_t totWeight = 0;
    for(int i = 0; i < nEntry; i++) {
        ch->GetEntry(i);
        cout<<"Processing: "<<(i+1)<<"th entry\r";
        if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || psi2S_pt->at(0) < 10 || psi2S_pt->at(0) > 40) continue;
        if(fabs(Jpsi_y->at(0)) > 2 || fabs(psi2S_y->at(0)) > 2) continue;
        nEvent++;
        if(Jpsi_mu_pt->at(0)[0] < 3.5 || Jpsi_mu_pt->at(0)[1] < 3.5 || psi2S_mu_pt->at(0)[0] < 3.5 || psi2S_mu_pt->at(0)[1] < 3.5) continue;
        if(fabs(Jpsi_mu_eta->at(0)[0]) > 2.4 || fabs(Jpsi_mu_eta->at(0)[1]) > 2.4 || fabs(psi2S_mu_eta->at(0)[0]) > 2.4 || fabs(psi2S_mu_eta->at(0)[1]) > 2.4) continue;
        nPassAcc++;
        totWeight += calWeight(Jpsi_pt->at(0), Jpsi_y->at(0), psi2S_pt->at(0), psi2S_y->at(0));
    }
    cout<<"\nnEvent="<<nEvent<<"\nnPassAcc="<<nPassAcc<<"\ntotWeight="<<totWeight<<endl;
    // Separated version
    // Int_t nEntry = ch->GetEntries(), nEvent[] = {0, 0}, nPassAcc[] = {0, 0};
    // Double_t totWeight[] = {0, 0};
    // for(int i = 0; i < nEntry && i < 60000; i++) {
    //     ch->GetEntry(i);
    //     cout<<"Processing: "<<(i+1)<<"th entry\r";
    //     if(Jpsi_pt->empty() || psi2S_pt->empty()) continue;
    //     if(Jpsi_pt->at(0) < 10 || Jpsi_pt->at(0) > 40 || psi2S_pt->at(0) < 10 || psi2S_pt->at(0) > 40) continue;
    //     if(fabs(Jpsi_y->at(0)) > 2 || fabs(psi2S_y->at(0)) > 2) continue;
    //     nEvent[0]++;
    //     if(Jpsi_mu_pt->at(0)[0] > 3.5 && Jpsi_mu_pt->at(0)[1] > 3.5 && fabs(Jpsi_mu_eta->at(0)[0]) < 2.4 && fabs(Jpsi_mu_eta->at(0)[1]) < 2.4) {
    //         nPassAcc[0]++;
    //         totWeight[0] += calWeight_Jpsi(Jpsi_pt->at(0), Jpsi_y->at(0));
    //     }
    //     nEvent[1]++;
    //     if(psi2S_mu_pt->at(0)[0] > 3.5 && psi2S_mu_pt->at(0)[1] > 3.5 && fabs(psi2S_mu_eta->at(0)[0]) < 2.4 && fabs(psi2S_mu_eta->at(0)[1]) < 2.4) {
    //         nPassAcc[1]++;
    //         totWeight[1] += calWeight_psi2S(psi2S_pt->at(0), psi2S_y->at(0));
    //     }
    // }
    // cout<<"\n[J/psi]\nnEvent="<<nEvent[0]<<"\nnPassAcc="<<nPassAcc[0]<<"\ntotWeight="<<totWeight[0];
    // cout<<"\n[psi(2S)]\nnEvent="<<nEvent[1]<<"\nnPassAcc="<<nPassAcc[1]<<"\ntotWeight="<<totWeight[1]<<endl;
    return;
}
