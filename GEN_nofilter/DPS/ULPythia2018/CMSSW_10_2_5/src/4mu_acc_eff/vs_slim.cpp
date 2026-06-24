// Author: Shiyang CHEN
// Description: NtupleAnalyzer for 2018 MINIAOD/Monte-Carlo
// Implementation: Take Ntuple as input, apply physical seletion and acceptance&efficiency correction.
// Output root file contains following information: event(mass, weight), dimuon(mass, pt, eta, phi, ctau), muon(pt)
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
using namespace std;
#define PI 3.14159265359

class Process {
    private:
    void readTree(string& fileName) {
        // Handle root file
        TFile file(fileName.c_str(), "READ");
        TTree *inTree = (TTree *)file.Get("GenAnalyzer/gen_tree");
        int nEvent = inTree->GetEntries();
        // Define input tree variables
        vector<Double_t> *GENjpsi_pt = 0, *GENjpsi_eta = 0, *GENjpsi_y = 0, *GENjpsi_phi = 0, *GENjpsi_mass = 0, *GENjpsi_ctau = 0;
        vector<Double_t> *GENpsi2s_pt = 0, *GENpsi2s_eta = 0, *GENpsi2s_y = 0, *GENpsi2s_phi = 0, *GENpsi2s_mass = 0, *GENpsi2s_ctau = 0;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("GENjpsi_pt", &GENjpsi_pt);
        inTree->SetBranchAddress("GENjpsi_eta", &GENjpsi_eta);
        inTree->SetBranchAddress("GENjpsi_y", &GENjpsi_y);
        inTree->SetBranchAddress("GENjpsi_phi", &GENjpsi_phi);
        inTree->SetBranchAddress("GENjpsi_mass", &GENjpsi_mass);
        inTree->SetBranchAddress("GENjpsi_ctau", &GENjpsi_ctau);

        inTree->SetBranchAddress("GENpsi2s_pt", &GENpsi2s_pt);
        inTree->SetBranchAddress("GENpsi2s_eta", &GENpsi2s_eta);
        inTree->SetBranchAddress("GENpsi2s_y", &GENpsi2s_y);
        inTree->SetBranchAddress("GENpsi2s_phi", &GENpsi2s_phi);
        inTree->SetBranchAddress("GENpsi2s_mass", &GENpsi2s_mass);
        inTree->SetBranchAddress("GENpsi2s_ctau", &GENpsi2s_ctau);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            // Apply fiducial cut
            if(GENjpsi_mass->empty() || GENpsi2s_mass->empty()) continue;
            if(GENjpsi_pt->at(0) > 40 || GENjpsi_pt->at(0) < 10 || GENpsi2s_pt->at(0) > 40 || GENpsi2s_pt->at(0) < 10) continue;
            if(fabs(GENjpsi_y->at(0)) > 2 || fabs(GENpsi2s_y->at(0)) > 2) continue;

            TLorentzVector JpsiLV, psi2SLV;
            JpsiLV.SetPtEtaPhiM(GENjpsi_pt->at(0), GENjpsi_eta->at(0), GENjpsi_phi->at(0), GENjpsi_mass->at(0));
            psi2SLV.SetPtEtaPhiM(GENpsi2s_pt->at(0), GENpsi2s_eta->at(0), GENpsi2s_phi->at(0), GENpsi2s_mass->at(0));
            evt_mass.push_back((JpsiLV + psi2SLV).M());
            delta_phi.push_back(PI - fabs(fabs(GENjpsi_phi->at(0) - GENpsi2s_phi->at(0)) - PI));
            delta_y.push_back(fabs(GENjpsi_y->at(0) - GENpsi2s_y->at(0)));
            evt_pt.push_back((JpsiLV + psi2SLV).Pt());
            evt_y.push_back(fabs((JpsiLV + psi2SLV).Rapidity()));

            Jpsi_mass.push_back(GENjpsi_mass->at(0));
            Jpsi_ctau.push_back(GENjpsi_ctau->at(0));
            Jpsi_pt.push_back(GENjpsi_pt->at(0));
            Jpsi_y.push_back(GENjpsi_y->at(0));

            psi2S_mass.push_back(GENpsi2s_mass->at(0));
            psi2S_ctau.push_back(GENpsi2s_ctau->at(0));
            psi2S_pt.push_back(GENpsi2s_pt->at(0));
            psi2S_y.push_back(GENpsi2s_y->at(0));
            totEvent++;
        }
        cout<<"Finished.  >> "<<endl;
        return;
    }

    public:
    int totEvent = 0;
    vector<Double_t> Jpsi_mass, Jpsi_ctau, Jpsi_pt, Jpsi_y;
    vector<Double_t> psi2S_mass, psi2S_ctau, psi2S_pt, psi2S_y;
    vector<Double_t> evt_mass, evt_y, evt_pt, delta_y, delta_phi;
    
    void loopOn() {
        string prefix = "8/DPS_2018_Psi2SJ_";
        for(int i = 1; i <= 10; i++) {
            if(i == 6) continue;
            string filename = prefix + to_string(i) + ".root";
            readTree(filename);
        }
        for(int i = 1; i <= 5; i++) {
            prefix = "8/" + to_string(i) + "/DPS_2018_Psi2SJ_";
            for(int j = 1; j <= 20; j++) {
                string filename = prefix + to_string(j) + ".root";
                readTree(filename);
            }
        }
    }
};

string outFile = "UnweightDPS.root";
void slim() {
    Process process;
    process.loopOn();
    TFile file(outFile.c_str(), "RECREATE");
    TTree *outTree = new TTree("data", "data");
    // Define output tree variables
    Double_t Jpsi_mass, Jpsi_ctau, Jpsi_pt, Jpsi_y;
    Double_t psi2S_mass, psi2S_ctau, psi2S_pt, psi2S_y;
    Double_t evt_mass, evt_y, evt_pt, delta_y, delta_phi;
    Double_t evt_weight;
    // Set output tree SetBranchAddress address
    outTree->Branch("Jpsi_mass", &Jpsi_mass, "Jpsi_mass/D");
    outTree->Branch("Jpsi_ctau", &Jpsi_ctau, "Jpsi_ctau/D");
    outTree->Branch("Jpsi_pt", &Jpsi_pt, "Jpsi_pt/D");
    outTree->Branch("Jpsi_y", &Jpsi_y, "Jpsi_y/D");

    outTree->Branch("psi2S_mass", &psi2S_mass, "psi2S_mass/D");
    outTree->Branch("psi2S_ctau", &psi2S_ctau, "psi2S_ctau/D");
    outTree->Branch("psi2S_pt", &psi2S_pt, "psi2S_pt/D");
    outTree->Branch("psi2S_y", &psi2S_y, "psi2S_y/D");

    outTree->Branch("evt_mass", &evt_mass, "evt_mass/D");
    outTree->Branch("evt_y", &evt_y, "evt_y/D");
    outTree->Branch("evt_pt", &evt_pt, "evt_pt/D");
    outTree->Branch("delta_y", &delta_y, "delta_y/D");
    outTree->Branch("delta_phi", &delta_phi, "delta_phi/D");

    outTree->Branch("evt_weight", &evt_weight, "evt_weight/D");
    for(int i = 0; i < process.totEvent; i++) {
        Jpsi_mass = process.Jpsi_mass[i];
        Jpsi_ctau = process.Jpsi_ctau[i];
        Jpsi_pt = process.Jpsi_pt[i];
        Jpsi_y = process.Jpsi_y[i];

        psi2S_mass = process.psi2S_mass[i];
        psi2S_ctau = process.psi2S_ctau[i];
        psi2S_pt = process.psi2S_pt[i];
        psi2S_y = process.psi2S_y[i];
        
        evt_mass = process.evt_mass[i];
        evt_y = process.evt_y[i];
        evt_pt = process.evt_pt[i];
        delta_y = process.delta_y[i];
        delta_phi = process.delta_phi[i];

        evt_weight = 1;
        outTree->Fill();
    }
    outTree->Write();
    cout<<"Done! Total event number: "<<process.totEvent<<endl;
    return;
}
