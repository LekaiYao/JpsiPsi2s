#include <iostream>
#include <fstream>
using namespace std;

vector<Double_t> eff_pt, eff_eta;
Double_t **nBin_Jpsi, **nVtx_Jpsi, **nBin_psi2S, **nVtx_psi2S, **nVtx_evt, **nTrg_evt;

void merge() {
    string line;
    // Read SPS efficiency map
    ifstream effFile("plot/efficiency.txt");
    if(!effFile.is_open()) return;
    int eff_ptBin = 0, eff_etaBin = 0, lineCnt = 0, n_Jpsi_SPS, n_psi2S_SPS;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>eff_ptBin>>eff_etaBin>>n_Jpsi_SPS>>n_psi2S_SPS;
            if(n_Jpsi_SPS != n_psi2S_SPS) {
                cout<<"Invalid Jpsi/psi(2S) number!\n";
                return;
            }
            nBin_Jpsi = new Double_t*[eff_etaBin];
            nVtx_Jpsi = new Double_t*[eff_etaBin];
            nBin_psi2S = new Double_t*[eff_etaBin];
            nVtx_psi2S = new Double_t*[eff_etaBin];
            nVtx_evt = new Double_t*[eff_ptBin];
            nTrg_evt = new Double_t*[eff_ptBin];
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= eff_ptBin; i++) {
                iss>>pt;
                eff_pt.push_back(pt);
            }
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= eff_etaBin; i++) {
                iss>>eta;
                eff_eta.push_back(eta);
            }
        }else if(lineCnt <= 2 + eff_etaBin) {
            int i = lineCnt - 3;
            nVtx_Jpsi[i] = new Double_t[eff_ptBin];
            nBin_Jpsi[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j];
        }else if(lineCnt <= 2 + 2 * eff_etaBin) {
            int i = lineCnt - 3 - eff_etaBin;
            nVtx_psi2S[i] = new Double_t[eff_ptBin];
            nBin_psi2S[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_psi2S[i][j]>>nBin_psi2S[i][j];
        }else if(lineCnt <= 2 + 2 * eff_etaBin + eff_ptBin) {
            int i = lineCnt - 3 - 2 * eff_etaBin;
            nVtx_evt[i] = new Double_t[eff_ptBin];
            nTrg_evt[i] = new Double_t[eff_ptBin];
            for(int j = 0; j < eff_ptBin; j++) iss>>nTrg_evt[i][j]>>nVtx_evt[i][j];
        }else break;
        lineCnt++;
    }
    effFile.close();
    // Add DPS efficiency map with a factor
    effFile = ifstream("/eos/home-c/chensh/JPsiPsi2s/SKIM_tightfilter/DPS/ULPythia2018/CMSSW_10_2_5/src/4mu_acc_eff/plot/efficiency.txt");
    if(!effFile.is_open()) return;
    int n_Jpsi_DPS, n_psi2S_DPS;
    eff_ptBin = 0;
    eff_etaBin = 0;
    lineCnt = 0;
    Double_t fac_DPS = 2;
    while(getline(effFile, line)) {
        istringstream iss(line);
        if(lineCnt == 0) {
            iss>>eff_ptBin>>eff_etaBin>>n_Jpsi_DPS>>n_psi2S_DPS;
            if(n_Jpsi_DPS != n_psi2S_DPS) {
                cout<<"Invalid Jpsi/psi(2S) number!\n";
                return;
            }
            fac_DPS *= n_Jpsi_SPS / n_Jpsi_DPS * 0.223810 / 0.339468;
        }else if(lineCnt == 1) {
            Double_t pt;
            for(int i = 0; i <= eff_ptBin; i++) iss>>pt;
        }else if(lineCnt == 2) {
            Double_t eta;
            for(int i = 0; i <= eff_etaBin; i++) iss>>eta;
        }else if(lineCnt <= 2 + eff_etaBin) {
            int i = lineCnt - 3;
            Double_t nVtx_Jpsi_ij, nBin_Jpsi_ij;
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nVtx_Jpsi_ij>>nBin_Jpsi_ij;
                nVtx_Jpsi[i][j] += fac_DPS * nVtx_Jpsi_ij;
                nBin_Jpsi[i][j] += fac_DPS * nBin_Jpsi_ij;
            }
        }else if(lineCnt <= 2 + 2 * eff_etaBin) {
            int i = lineCnt - 3 - eff_etaBin;
            Double_t nVtx_psi2S_ij, nBin_psi2S_ij;
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nVtx_psi2S_ij>>nBin_psi2S_ij;
                nVtx_psi2S[i][j] += fac_DPS * nVtx_psi2S_ij;
                nBin_psi2S[i][j] += fac_DPS * nBin_psi2S_ij;
        }
        }else if(lineCnt <= 2 + 2 * eff_etaBin + eff_ptBin) {
            int i = lineCnt - 3 - 2 * eff_etaBin;
            Double_t nVtx_evt_ij, nTrg_evt_ij;
            for(int j = 0; j < eff_ptBin; j++) {
                iss>>nTrg_evt_ij>>nVtx_evt_ij;
                nVtx_evt[i][j] += fac_DPS * nVtx_evt_ij;
                nTrg_evt[i][j] += fac_DPS * nTrg_evt_ij;
        }
        }else break;
        lineCnt++;
    }
    effFile.close();
    // Output merged efficiency map
    ofstream outFile("efficiency_SPS+2DPS.txt");
    outFile<<eff_ptBin<<' '<<eff_etaBin<<endl;
    for(int i = 0; i <= eff_ptBin; i++) outFile<<eff_pt[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i <= eff_etaBin; i++) outFile<<eff_eta[i]<<' ';
    outFile<<'\n';
    for(int i = 0; i < eff_etaBin; i++) {
        for(int j = 0; j < eff_ptBin; j++) outFile<<nVtx_Jpsi[i][j]<<' '<<nBin_Jpsi[i][j]<<' ';
        outFile<<'\n';
    }
    for(int i = 0; i < eff_etaBin; i++) {
        for(int j = 0; j < eff_ptBin; j++) outFile<<nVtx_psi2S[i][j]<<' '<<nBin_psi2S[i][j]<<' ';
        outFile<<'\n';
    }
    for(int i = 0; i < eff_ptBin; i++) {
        for(int j = 0; j < eff_ptBin; j++) outFile<<nTrg_evt[i][j]<<' '<<nVtx_evt[i][j]<<' ';
        outFile<<'\n';
    }
    outFile.close();
    return;
}
