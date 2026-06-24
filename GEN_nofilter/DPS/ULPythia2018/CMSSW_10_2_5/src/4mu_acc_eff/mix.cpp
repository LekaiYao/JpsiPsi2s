// Author: Shiyang CHEN
// Description: Event mixing for 2018 DPS Monte-Carlo
// Implementation: Randomly composite J/psi and psi(2S) to make a new event.
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
        vector<Double_t> *GENjpsi_pt = 0, *GENjpsi_eta = 0, *GENjpsi_phi = 0, *GENjpsi_mass = 0, *GENjpsi_y = 0;
        vector<Double_t> *GENpsi2s_pt = 0, *GENpsi2s_eta = 0, *GENpsi2s_phi = 0, *GENpsi2s_mass = 0, *GENpsi2s_y = 0;
        vector<vector<Double_t>> *GENjpsi_mu_pt = 0, *GENjpsi_mu_eta = 0, *GENpsi2s_mu_pt = 0, *GENpsi2s_mu_eta = 0;
        // Set input tree SetBranchAddress address
        inTree->SetBranchAddress("GENjpsi_pt", &GENjpsi_pt);
        inTree->SetBranchAddress("GENjpsi_eta", &GENjpsi_eta);
        inTree->SetBranchAddress("GENjpsi_phi", &GENjpsi_phi);
        inTree->SetBranchAddress("GENjpsi_mass", &GENjpsi_mass);
        inTree->SetBranchAddress("GENjpsi_y", &GENjpsi_y);
        inTree->SetBranchAddress("GENjpsi_mu_pt", &GENjpsi_mu_pt);
        inTree->SetBranchAddress("GENjpsi_mu_eta", &GENjpsi_mu_eta);

        inTree->SetBranchAddress("GENpsi2s_pt", &GENpsi2s_pt);
        inTree->SetBranchAddress("GENpsi2s_eta", &GENpsi2s_eta);
        inTree->SetBranchAddress("GENpsi2s_phi", &GENpsi2s_phi);
        inTree->SetBranchAddress("GENpsi2s_mass", &GENpsi2s_mass);
        inTree->SetBranchAddress("GENpsi2s_y", &GENpsi2s_y);
        inTree->SetBranchAddress("GENpsi2s_mu_pt", &GENpsi2s_mu_pt);
        inTree->SetBranchAddress("GENpsi2s_mu_eta", &GENpsi2s_mu_eta);
        // Loop on input tree entries
        cout<<fileName<<" has events: "<<nEvent<<endl;
        for(int i = 0; i < nEvent; i++) {
            cout<<"Processing No."<<i<<'\r';
            inTree->GetEntry(i);
            // Apply fiducial cut
            if(!GENjpsi_mass->empty() && GENjpsi_pt->at(0) < 40 && GENjpsi_pt->at(0) > 10 && fabs(GENjpsi_y->at(0)) < 2) {
                Jpsi_pt.push_back(GENjpsi_pt->at(0));
                Jpsi_eta.push_back(GENjpsi_eta->at(0));
                Jpsi_phi.push_back(GENjpsi_phi->at(0));
                Jpsi_mass.push_back(GENjpsi_mass->at(0));
                Jpsi_y.push_back(GENjpsi_y->at(0));
                if(GENjpsi_mu_pt->at(0)[0] > 3.5 && GENjpsi_mu_pt->at(0)[1] > 3.5 && fabs(GENjpsi_mu_eta->at(0)[0]) < 2.4 && fabs(GENjpsi_mu_eta->at(0)[1]) < 2.4) Jpsi_acc.push_back(true);
                else Jpsi_acc.push_back(false);
                Jpsi_mu_pt.push_back(GENjpsi_mu_pt->at(0));
                Jpsi_mu_eta.push_back(GENjpsi_mu_eta->at(0));
            }
            if(!GENpsi2s_mass->empty() && GENpsi2s_pt->at(0) < 40 && GENpsi2s_pt->at(0) > 10 && fabs(GENpsi2s_y->at(0)) < 2) {
                psi2S_pt.push_back(GENpsi2s_pt->at(0));
                psi2S_eta.push_back(GENpsi2s_eta->at(0));
                psi2S_phi.push_back(GENpsi2s_phi->at(0));
                psi2S_mass.push_back(GENpsi2s_mass->at(0));
                psi2S_y.push_back(GENpsi2s_y->at(0));
                if(GENpsi2s_mu_pt->at(0)[0] > 3.5 && GENpsi2s_mu_pt->at(0)[1] > 3.5 && fabs(GENpsi2s_mu_eta->at(0)[0]) < 2.4 && fabs(GENpsi2s_mu_eta->at(0)[1]) < 2.4) psi2S_acc.push_back(true);
                else psi2S_acc.push_back(false);
                psi2S_mu_pt.push_back(GENpsi2s_mu_pt->at(0));
                psi2S_mu_eta.push_back(GENpsi2s_mu_eta->at(0));
            }
        }
        cout<<"Finished.  >> "<<endl;
        return;
    }

    public:
    vector<Double_t> Jpsi_pt, Jpsi_eta, Jpsi_phi, Jpsi_mass, Jpsi_y, Jpsi_acc;
    vector<Double_t> psi2S_pt, psi2S_eta, psi2S_phi, psi2S_mass, psi2S_y, psi2S_acc;
    vector<vector<Double_t>> Jpsi_mu_pt, Jpsi_mu_eta, psi2S_mu_pt, psi2S_mu_eta;
    
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
        prefix = "9/1/DPS_2018_Psi2SJ_";
        for(int i = 1; i <= 20; i++) {
            string filename = prefix + to_string(i) + ".root";
            readTree(filename);
        }
    }
};

string outFile = "MixedDPS.root";
void mix(int goal = 200000) {
    Process process;
    process.loopOn();
    TFile file(outFile.c_str(), "RECREATE");
    // TTree *outTree = new TTree("GenAnalyzer/gen_tree", "gen_tree");
    TTree *outTree = new TTree("data", "data");
    // Define output tree variables
    vector<Double_t> Jpsi_pt, Jpsi_y, psi2S_pt, psi2S_y;
    vector<vector<Double_t>> Jpsi_mu_pt, Jpsi_mu_eta, psi2S_mu_pt, psi2S_mu_eta;
    Double_t evt_mass, evt_y, evt_pt, delta_y, delta_phi;
    bool evt_acc;
    // Set output tree SetBranchAddress address
    outTree->Branch("GENjpsi_pt", &Jpsi_pt);
    outTree->Branch("GENjpsi_y", &Jpsi_y);
    outTree->Branch("GENpsi2s_pt", &psi2S_pt);
    outTree->Branch("GENpsi2s_y", &psi2S_y);

    outTree->Branch("GENjpsi_mu_pt", &Jpsi_mu_pt);
    outTree->Branch("GENjpsi_mu_eta", &Jpsi_mu_eta);
    outTree->Branch("GENpsi2s_mu_pt", &psi2S_mu_pt);
    outTree->Branch("GENpsi2s_mu_eta", &psi2S_mu_eta);
    
    outTree->Branch("evt_mass", &evt_mass, "evt_mass/D");
    outTree->Branch("evt_y", &evt_y, "evt_y/D");
    outTree->Branch("evt_pt", &evt_pt, "evt_pt/D");
    outTree->Branch("delta_y", &delta_y, "delta_y/D");
    outTree->Branch("delta_phi", &delta_phi, "delta_phi/D");
    outTree->Branch("evt_acc", &evt_acc, "evt_acc/B");
    int nEvt = 1, nJpsi = process.Jpsi_pt.size(), npsi2S = process.psi2S_pt.size();
    srand(time(0));
    delete gRandom;
    gRandom = new TRandom3(rand());
    while(nEvt <= goal) {
        int xJpsi = gRandom->Rndm() * nJpsi, xpsi2S = gRandom->Rndm() * npsi2S;
        TLorentzVector JpsiLV, psi2SLV;
        JpsiLV.SetPtEtaPhiM(process.Jpsi_pt[xJpsi], process.Jpsi_eta[xJpsi], process.Jpsi_phi[xJpsi], process.Jpsi_mass[xJpsi]);
        psi2SLV.SetPtEtaPhiM(process.psi2S_pt[xpsi2S], process.psi2S_eta[xpsi2S], process.psi2S_phi[xpsi2S], process.psi2S_mass[xpsi2S]);
        evt_mass = (JpsiLV + psi2SLV).M();
        if(evt_mass < 7.5) continue;
        Jpsi_pt.clear();
        Jpsi_y.clear();
        psi2S_pt.clear();
        psi2S_y.clear();
        Jpsi_mu_pt.clear();
        Jpsi_mu_eta.clear();
        psi2S_mu_pt.clear();
        psi2S_mu_eta.clear();
        Jpsi_pt.push_back(process.Jpsi_pt[xJpsi]);
        Jpsi_y.push_back(process.Jpsi_y[xJpsi]);
        psi2S_pt.push_back(process.psi2S_pt[xpsi2S]);
        psi2S_y.push_back(process.psi2S_y[xpsi2S]);
        Jpsi_mu_pt.push_back(process.Jpsi_mu_pt[xJpsi]);
        Jpsi_mu_eta.push_back(process.Jpsi_mu_eta[xJpsi]);
        psi2S_mu_pt.push_back(process.psi2S_mu_pt[xpsi2S]);
        psi2S_mu_eta.push_back(process.psi2S_mu_eta[xpsi2S]);
        evt_y = fabs((JpsiLV + psi2SLV).Rapidity());
        evt_pt = (JpsiLV + psi2SLV).Pt();
        delta_y = fabs(process.Jpsi_y[xJpsi] - process.psi2S_y[xpsi2S]);
        delta_phi = PI - fabs(fabs(process.Jpsi_phi[xJpsi] - process.psi2S_phi[xpsi2S]) - PI);
        evt_acc = process.Jpsi_acc[xJpsi] && process.psi2S_acc[xpsi2S];
        outTree->Fill();
        nEvt++;
    }
    outTree->Write();
    cout<<"Done! Generate "<<(nEvt - 1)<<" events from ("<<nJpsi<<", "<<npsi2S<<")."<<endl;
    return;
}
