// Author: Shiyang CHEN
// Description: Provides event-mixing based on processed data or SPS sample.
// Implementation: 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
using namespace std;

// class Process {
//     private:
//     vector<Double_t> acc_pt, acc_eta, eff_pt, eff_eta;
//     int **nGen_Jpsi, **nAcc_Jpsi, **nGen_psi2S, **nAcc_psi2S;
//     int **nBin_Jpsi, **nVtx_Jpsi, **nBin_psi2S, **nVtx_psi2S, **nVtx_evt, **nTrg_evt;
//     public:
//     Double_t calWeight(Double_t Jpsi_pt, Double_t Jpsi_eta, Double_t psi2S_pt, Double_t psi2S_eta) {
//         int i = upper_bound(acc_eta.begin(), acc_eta.end(), Jpsi_eta) - acc_eta.begin() - 1, j = upper_bound(acc_pt.begin(), acc_pt.end(), Jpsi_pt) - acc_pt.begin() - 1;// i-J/psi eta, j-J/psi pT
//         int k = upper_bound(acc_eta.begin(), acc_eta.end(), psi2S_eta) - acc_eta.begin() - 1, l = upper_bound(acc_pt.begin(), acc_pt.end(), psi2S_pt) - acc_pt.begin() - 1;// k-psi(2S) eta, l-psi(2S) pT
//         Double_t w = (Double_t)nGen_Jpsi[i][j] / nAcc_Jpsi[i][j] * nGen_psi2S[k][l] / nAcc_psi2S[k][l];
//         i = upper_bound(eff_eta.begin(), eff_eta.end(), Jpsi_eta) - eff_eta.begin() - 1;
//         j = upper_bound(eff_pt.begin(), eff_pt.end(), Jpsi_pt) - eff_pt.begin() - 1;
//         k = upper_bound(eff_eta.begin(), eff_eta.end(), psi2S_eta) - eff_eta.begin() - 1;
//         l = upper_bound(eff_pt.begin(), eff_pt.end(), psi2S_pt) - eff_pt.begin() - 1;
//         w *= (Double_t)nBin_Jpsi[i][j] / nVtx_Jpsi[i][j] * nBin_psi2S[k][l] / nVtx_psi2S[k][l] * nVtx_evt[l][j] / nTrg_evt[l][j];
//         return w;
//     }
//     Process() {
//         // Save acc&eff in arrays
//         ifstream accFile("acceptance.txt");
//         if(!accFile.is_open()) return;
//         int acc_ptBin = 0, acc_etaBin = 0, eff_ptBin = 0, eff_etaBin = 0, lineCnt = 0;
//         string line;
//         while(getline(accFile, line)) {
//             istringstream iss(line);
//             if(lineCnt == 0) {
//                 iss>>acc_ptBin>>acc_etaBin;
//                 nGen_Jpsi = new int*[acc_etaBin];
//                 nAcc_Jpsi = new int*[acc_etaBin];
//                 nGen_psi2S = new int*[acc_etaBin];
//                 nAcc_psi2S = new int*[acc_etaBin];
//             }else if(lineCnt == 1) {
//                 Double_t pt;
//                 for(int i = 0; i <= acc_ptBin; i++) {
//                     iss>>pt;
//                     acc_pt.push_back(pt);
//                 }
//             }else if(lineCnt == 2) {
//                 Double_t eta;
//                 for(int i = 0; i <= acc_etaBin; i++) {
//                     iss>>eta;
//                     acc_eta.push_back(eta);
//                 }
//             }else if(lineCnt <= 2 + acc_etaBin) {
//                 int i = lineCnt - 3;
//                 nAcc_Jpsi[i] = new int[acc_ptBin];
//                 nGen_Jpsi[i] = new int[acc_ptBin];
//                 for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_Jpsi[i][j]>>nGen_Jpsi[i][j];
//             }else if(lineCnt <= 2 + 2 * acc_etaBin) {
//                 int i = lineCnt - 3 - acc_etaBin;
//                 nAcc_psi2S[i] = new int[acc_ptBin];
//                 nGen_psi2S[i] = new int[acc_ptBin];
//                 for(int j = 0; j < acc_ptBin; j++) iss>>nAcc_psi2S[i][j]>>nGen_psi2S[i][j];
//             }else break;
//             lineCnt++;
//         }
//         accFile.close();
//         ifstream effFile("efficiency.txt");
//         if(!effFile.is_open()) return;
//         lineCnt = 0;
//         while(getline(effFile, line)) {
//             istringstream iss(line);
//             if(lineCnt == 0) {
//                 iss>>eff_ptBin>>eff_etaBin;
//                 nBin_Jpsi = new int*[eff_etaBin];
//                 nVtx_Jpsi = new int*[eff_etaBin];
//                 nBin_psi2S = new int*[eff_etaBin];
//                 nVtx_psi2S = new int*[eff_etaBin];
//                 nVtx_evt = new int*[eff_ptBin];
//                 nTrg_evt = new int*[eff_ptBin];
//             }else if(lineCnt == 1) {
//                 Double_t pt;
//                 for(int i = 0; i <= eff_ptBin; i++) {
//                     iss>>pt;
//                     eff_pt.push_back(pt);
//                 }
//             }else if(lineCnt == 2) {
//                 Double_t eta;
//                 for(int i = 0; i <= eff_etaBin; i++) {
//                     iss>>eta;
//                     eff_eta.push_back(eta);
//                 }
//             }else if(lineCnt <= 2 + eff_etaBin) {
//                 int i = lineCnt - 3;
//                 nVtx_Jpsi[i] = new int[eff_ptBin];
//                 nBin_Jpsi[i] = new int[eff_ptBin];
//                 for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_Jpsi[i][j]>>nBin_Jpsi[i][j];
//             }else if(lineCnt <= 2 + 2 * eff_etaBin) {
//                 int i = lineCnt - 3 - eff_etaBin;
//                 nVtx_psi2S[i] = new int[eff_ptBin];
//                 nBin_psi2S[i] = new int[eff_ptBin];
//                 for(int j = 0; j < eff_ptBin; j++) iss>>nVtx_psi2S[i][j]>>nBin_psi2S[i][j];
//             }else if(lineCnt <= 2 + 2 * eff_etaBin + eff_ptBin) {
//                 int i = lineCnt - 3 - 2 * eff_etaBin;
//                 nVtx_evt[i] = new int[eff_ptBin];
//                 nTrg_evt[i] = new int[eff_ptBin];
//                 for(int j = 0; j < eff_ptBin; j++) iss>>nTrg_evt[i][j]>>nVtx_evt[i][j];
//             }else break;
//             lineCnt++;
//         }
//         effFile.close();
//         return;
//     }
// };
// vector<Double_t> JpsiEta, JpsiPhi, JpsiCtau, JpsiPt, JpsiMass;
// vector<Double_t> psi2SEta, psi2SPhi, psi2SCtau, psi2SPt, psi2SMass;
// vector<Double_t> evtWeight;
// int nEvent;
// void read() {
//     // TFile inFile("WeightData.root", "READ");
//     // TFile inFile("WeightSPS.root", "READ");
//     TFile inFile("WeightDPS.root", "READ");
//     TTree *inTree = (TTree *)inFile.Get("data");
//     // Define variables
//     Double_t Jpsi_eta, Jpsi_phi, Jpsi_ctau, Jpsi_pt, Jpsi_mass;
//     Double_t psi2S_eta, psi2S_phi, psi2S_ctau, psi2S_pt, psi2S_mass;
//     Double_t evt_weight;//, evt_mass;
//     inTree->SetBranchAddress("Jpsi_eta", &Jpsi_eta);
//     inTree->SetBranchAddress("Jpsi_phi", &Jpsi_phi);
//     inTree->SetBranchAddress("Jpsi_ctau", &Jpsi_ctau);
//     inTree->SetBranchAddress("Jpsi_pt", &Jpsi_pt);
//     inTree->SetBranchAddress("Jpsi_mass", &Jpsi_mass);
//     inTree->SetBranchAddress("psi2S_eta", &psi2S_eta);
//     inTree->SetBranchAddress("psi2S_phi", &psi2S_phi);
//     inTree->SetBranchAddress("psi2S_ctau", &psi2S_ctau);
//     inTree->SetBranchAddress("psi2S_pt", &psi2S_pt);
//     inTree->SetBranchAddress("psi2S_mass", &psi2S_mass);
//     inTree->SetBranchAddress("evt_weight", &evt_weight);
//     nEvent = inTree->GetEntries();
//     for(int i = 0; i < nEvent; i++) {
//         inTree->GetEntry(i);
//         JpsiEta.push_back(Jpsi_eta);
//         JpsiPhi.push_back(Jpsi_phi);
//         JpsiCtau.push_back(Jpsi_ctau);
//         JpsiPt.push_back(Jpsi_pt);
//         JpsiMass.push_back(Jpsi_mass);
//         psi2SEta.push_back(psi2S_eta);
//         psi2SPhi.push_back(psi2S_phi);
//         psi2SCtau.push_back(psi2S_ctau);
//         psi2SPt.push_back(psi2S_pt);
//         psi2SMass.push_back(psi2S_mass);
//         evtWeight.push_back(evt_weight);
//     }
//     return;
// }
// void mix(Int_t num = 1000) {
//     read();
//     Process process;
//     cout<<nEvent<<" events in pool."<<endl;
//     TFile file("MixingDPS.root", "RECREATE");
//     TTree *outTree = new TTree("data", "data");
//     // Define output tree variables
//     Double_t Jpsi_eta, Jpsi_phi;
//     Double_t psi2S_eta, psi2S_phi;
//     Double_t evt_weight;
//     // Set output tree branch address
//     outTree->Branch("Jpsi_eta", &Jpsi_eta, "Jpsi_eta/D");
//     outTree->Branch("Jpsi_phi", &Jpsi_phi, "Jpsi_phi/D");
//     outTree->Branch("psi2S_eta", &psi2S_eta, "psi2S_eta/D");
//     outTree->Branch("psi2S_phi", &psi2S_phi, "psi2S_phi/D");
//     outTree->Branch("evt_weight", &evt_weight, "evt_weight/D");
//     // Randomly composite J/psi and psi(2S)
//     for(int i = 0; i < num; i++) {
//         int j = gRandom->Integer(nEvent), k = gRandom->Integer(nEvent);
//         // while(fabs(JpsiCtau[j]) > 0.02) j = gRandom->Integer(nEvent);
//         // while(k == j) k = gRandom->Integer(nEvent);// || fabs(psi2SCtau[k]) > 0.02
//         TLorentzVector JpsiLV, psi2SLV;
//         JpsiLV.SetPtEtaPhiM(JpsiPt[j], JpsiEta[j], JpsiPhi[j], JpsiMass[j]);
//         psi2SLV.SetPtEtaPhiM(psi2SPt[k], psi2SEta[k], psi2SPhi[k], psi2SMass[k]);
//         if((JpsiLV + psi2SLV).M() < 7.5) continue;
//         Jpsi_eta = JpsiEta[j];
//         Jpsi_phi = JpsiPhi[j];
//         psi2S_eta = psi2SEta[k];
//         psi2S_phi = psi2SPhi[k];
//         evt_weight = process.calWeight(JpsiPt[j], JpsiEta[j], psi2SPt[k], psi2SEta[k]);
//         outTree->Fill();
//         cout<<"Picking "<<j<<", "<<k<<"   \r";
//     }
//     outTree->Write();
//     cout<<endl;
//     return;
// }

vector<Double_t> JpsiMass, JpsiCtau, psi2SMass, psi2SCtau, evtWeight;
int nSrc = 2;
string files[] = {"WeightSPS.root", "WeightDPS.root", "WeightBdecay.root"};
// int ns[] = {5000, 10000, 15000};
double wlist[] = {490.72, 981.43, 3066.57}, scale = 100;
void read() {
    srand(time(0));
    delete gRandom;
    gRandom = new TRandom3(rand());
    for(int i = 0; i < nSrc; i++) {
        string file = files[i];
        double w = wlist[i] * scale;
        cout<<"Demanding "<<w<<" initial events from "<<file<<endl;
        TFile f(file.c_str(), "READ");
        TTree *t = (TTree *)f.Get("data");
        Double_t Jpsi_mass, Jpsi_ctau, psi2S_mass, psi2S_ctau, evt_weight;
        t->SetBranchAddress("Jpsi_mass", &Jpsi_mass);
        t->SetBranchAddress("Jpsi_ctau", &Jpsi_ctau);
        t->SetBranchAddress("psi2S_mass", &psi2S_mass);
        t->SetBranchAddress("psi2S_ctau", &psi2S_ctau);
        t->SetBranchAddress("evt_weight", &evt_weight);
        int nEvent = t->GetEntries();
        double sumw = 0, actw = 0;
        for(int j = 0; j < nEvent; j++) {
            t->GetEntry(j);
            sumw += evt_weight;
        }
        cout<<"[event number: "<<nEvent<<"] [sum of weights: "<<sumw<<"]"<<endl;
        // double r = (double)nEvent / sumw, p = w / (double)nEvent;
        double r = 1, p = w / sumw;
        if(p > 1) cout<<"No enough events from "<<file<<endl;
        for(int j = 0; j < nEvent; j++) {
            // int x = gRandom->Integer(nEvent);
            // if(x > n) continue;
            if(gRandom->Rndm() > p) continue;
            t->GetEntry(j);
            JpsiMass.push_back(Jpsi_mass);
            JpsiCtau.push_back(Jpsi_ctau);
            psi2SMass.push_back(psi2S_mass);
            psi2SCtau.push_back(psi2S_ctau);
            evtWeight.push_back(evt_weight * r);
            actw += evt_weight;
        }
        cout<<"[actual weight: "<<actw * r<<"]"<<endl;
    }
}
void mix() {
    read();
    TFile file("WeightMix.root", "RECREATE");
    TTree *outTree = new TTree("data", "data");
    Double_t Jpsi_mass, Jpsi_ctau, psi2S_mass, psi2S_ctau, evt_weight;
    outTree->Branch("Jpsi_mass", &Jpsi_mass, "Jpsi_mass/D");
    outTree->Branch("Jpsi_ctau", &Jpsi_ctau, "Jpsi_ctau/D");
    outTree->Branch("psi2S_mass", &psi2S_mass, "psi2S_mass/D");
    outTree->Branch("psi2S_ctau", &psi2S_ctau, "psi2S_ctau/D");
    outTree->Branch("evt_weight", &evt_weight, "evt_weight/D");
    int nEvent = JpsiMass.size();
    for(int i = 0; i < nEvent; i++) {
        Jpsi_mass = JpsiMass[i];
        Jpsi_ctau = JpsiCtau[i];
        psi2S_mass = psi2SMass[i];
        psi2S_ctau = psi2SCtau[i];
        evt_weight = evtWeight[i];
        outTree->Fill();
    }
    outTree->Write();
    return;
}