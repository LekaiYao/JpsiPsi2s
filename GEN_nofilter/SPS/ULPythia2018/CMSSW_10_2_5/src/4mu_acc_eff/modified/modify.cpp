#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

void process(string fn) {
    TFile f(fn.c_str(), "update");
    TTree *t = (TTree *)f.Get("GenAnalyzer/gen_tree");
    int nEvt = t->GetEntries();
    Double_t evt_mass;
    vector<Double_t> *Jpsi_pt = 0, *Jpsi_eta = 0, *Jpsi_phi = 0, *Jpsi_mass = 0;
    vector<Double_t> *psi2S_pt = 0, *psi2S_eta = 0, *psi2S_phi = 0, *psi2S_mass = 0;
    t->Branch("evt_mass", &evt_mass, "evt_mass/D");
    t->SetBranchAddress("GENjpsi_pt", &Jpsi_pt);
    t->SetBranchAddress("GENpsi2s_pt", &psi2S_pt);
    t->SetBranchAddress("GENjpsi_eta", &Jpsi_eta);
    t->SetBranchAddress("GENpsi2s_eta", &psi2S_eta);
    t->SetBranchAddress("GENjpsi_phi", &Jpsi_phi);
    t->SetBranchAddress("GENpsi2s_phi", &psi2S_phi);
    t->SetBranchAddress("GENjpsi_mass", &Jpsi_mass);
    t->SetBranchAddress("GENpsi2s_mass", &psi2S_mass);
    for(int i = 0; i < nEvt; i++) {
        cout<<"Processing "<<i<<"th...\r";
        t->GetEntry(i);
        TLorentzVector l1, l2;
        l1.SetPtEtaPhiM(Jpsi_pt->at(0), Jpsi_eta->at(0), Jpsi_phi->at(0), Jpsi_mass->at(0));
        l2.SetPtEtaPhiM(psi2S_pt->at(0), psi2S_eta->at(0), psi2S_phi->at(0), psi2S_mass->at(0));
        evt_mass = (l1 + l2).M();
        t->Fill();
    }
    cout<<endl;
    t->Write("", TObject::kOverwrite);
    return;
}

void modify() {
    int n = 1;
    for(int i = 1; i <= n; i++) process("SPS_2018_Psi2SJ_" + to_string(i) + ".root");
    return;
}