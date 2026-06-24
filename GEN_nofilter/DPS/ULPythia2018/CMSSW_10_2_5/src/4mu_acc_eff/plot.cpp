// Author: Shiyang CHEN
// Description: 
// Implementation: 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
using namespace std;

void plot() {
    // Configuration
    const string prefix = "7/DPS_2018_Psi2SJ_";
    const int infix = 16, nVars = 2;
    const string vars[] = {"GENjpsi_ctau", "GENpsi2S_ctau"};
    const double uppers[] = {0.16, 0.16}, lowers[] = {-0.03, -0.03};
    const bool setYLog[] = {true, true};
    // Primary loop
    for(int i = 0; i < nVars; i++) {
        string var = vars[i];
        const double upX = uppers[i], lowX = lowers[i];
        // Secondary loop
        TH1D *hs = new TH1D(var.c_str(), var.c_str(), 100, lowX, upX);
        int sumn = 0, sume = 0;
        for(int j = 1; j <= infix; j++) {
            string file = prefix + to_string(j) + ".root";
            TFile inFile(file.c_str(), "READ");
            TTree *inTree = (TTree *)inFile.Get("GenAnalyzer/gen_tree");
            vector<Double_t> *val = 0, *jpsi_pt = 0, *jpsi_y = 0, *psi2S_pt = 0, *psi2S_y = 0;
            inTree->SetBranchAddress(var.c_str(), &val);
            inTree->SetBranchAddress("GENjpsi_pt", &jpsi_pt);
            inTree->SetBranchAddress("GENjpsi_y", &jpsi_y);
            inTree->SetBranchAddress("GENpsi2s_pt", &psi2S_pt);
            inTree->SetBranchAddress("GENpsi2s_y", &psi2S_y);
            int nEvent = inTree->GetEntries();
            sume += nEvent;
            for(int k = 0; k < nEvent; k++) {
                inTree->GetEntry(k);
                if(val->empty()) continue;
                int id[] = {-1, -1};
                for(int l = 0; l < jpsi_pt->size(); l++) {
                    if(jpsi_pt->at(l) < 10 || jpsi_pt->at(l) > 40) continue;
                    if(jpsi_y->at(l) < -2 || jpsi_y->at(l) > 2) continue;
                    id[0] = l;
                    break;
                }
                for(int l = 0; l < psi2S_pt->size(); l++) {
                    if(psi2S_pt->at(l) < 10 || psi2S_pt->at(l) > 40) continue;
                    if(psi2S_y->at(l) < -2 || psi2S_y->at(l) > 2) continue;
                    id[1] = l;
                    break;
                }
                if(id[0] == -1 || id[1] == -1) continue;
                hs->Fill(val->at(id[i]));
                sumn += 1;
            } 
        }
        // hs->Scale(1.0 / sumn);
        double maxHs = hs->GetMaximum();
        TCanvas *c = new TCanvas("c", var.c_str(), 2000, 1000);
        gStyle->SetOptStat(0);
        if(setYLog[i]) {
            gPad->SetLogy();
            hs->GetYaxis()->SetRangeUser(0.5, 1.1 * maxHs);
        }else hs->GetYaxis()->SetRangeUser(0, 1.1 * maxHs);
        hs->Draw();
        c->SaveAs(("plots/" + var + ".png").c_str());
        cout<<"Found "<<sume<<" events in "<<var<<" dimension, choosing "<<sumn<<" events\n";
        delete c;
        delete hs;
    }
}