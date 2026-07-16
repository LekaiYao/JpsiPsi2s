// Author: Shiyang CHEN
// Description: For each variables, compare its distribution among datasets.
// Implementation: ...
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
using namespace std;

void para_plot() {
    // Configuration
    const int nFiles = 4, nVars = 4;
    const string files[] = {"WeightData.root", "WeightSPS.root", "WeightDPS.root", "WeightBdecay.root"};
    // Setting-1 nVars=4
    // const string vars[] = {"Jpsi_pt", "psi2S_pt", "Jpsi_y", "psi2S_y"};
    // const double uppers[] = {40, 40, 2, 2}, lowers[] = {10, 10, -2, -2};
    // const bool setYLog[] = {false, false, false, false};
    // Setting-2 nVars=4
    const string vars[] = {"Jpsi_mass", "Jpsi_ctau", "psi2S_mass", "psi2S_ctau"};
    const double uppers[] = {3.25, 0.16, 4.05, 0.16}, lowers[] = {2.95, -0.04,  3.35, -0.04};
    const bool setYLog[] = {false, true, false, true};
    // Setting-3 nVars=2
    // const string vars[] = {"evt_vtxProb", "evt_mass"};
    // const double uppers[] = {0.1, 100}, lowers[] = {0, 0};
    // const bool setYLog[] = {true, false};
    const unsigned int colors[] = {kBlue, kViolet, kRed, kBlack};
    // Primary loop
    TH1D **hs = new TH1D*[nFiles], **whs = new TH1D*[nFiles];
    for(int i = 0; i < nVars; i++) {
        string var = vars[i];
        const double upper = uppers[i], lower = lowers[i];
        // Secondary loop
        double maxHs = 0, maxWhs = 0;
        for(int j = 0; j < nFiles; j++) {
            string file = files[j];
            hs[j] = new TH1D(var.c_str(), var.c_str(), 100, lower, upper);
            whs[j] = new TH1D(var.c_str(), (var + "-weighted").c_str(), 100, lower, upper);
            TFile inFile(file.c_str(), "READ");
            TTree *inTree = (TTree *)inFile.Get("data");
            Double_t val, w, sumn = 0, sumw = 0, vtx = 0;
            inTree->SetBranchAddress(var.c_str(), &val);
            inTree->SetBranchAddress("evt_weight", &w);
            // if(j == 2) inTree->SetBranchAddress("evt_vtxProb", &vtx);
            int nEvent = inTree->GetEntries();
            for(int k = 0; k < nEvent; k++) {
                inTree->GetEntry(k);
                // if(j == 2 && vtx == 0) continue;
                hs[j]->Fill(val);
                whs[j]->Fill(val, w);
                sumn += 1;
                sumw += w;
            }
            hs[j]->SetLineColor(colors[j]);
            hs[j]->Scale(1.0 / sumn);
            whs[j]->SetLineColor(colors[j]);
            whs[j]->Scale(1.0 / sumw);
            if(hs[j]->GetMaximum() > maxHs) maxHs = hs[j]->GetMaximum();
            if(whs[j]->GetMaximum() > maxWhs) maxWhs = whs[j]->GetMaximum();
        }
        TCanvas *c = new TCanvas("c", var.c_str(), 2000, 1000);
        gStyle->SetOptStat(0);
        c->Divide(2, 1);
        c->cd(1);
        if(setYLog[i]) {
            gPad->SetLogy();
            hs[0]->GetYaxis()->SetRangeUser(1e-4, 1.1 * maxHs);
        }else hs[0]->GetYaxis()->SetRangeUser(0, 1.1 * maxHs);
        hs[0]->Draw();
        for(int j = 1; j < nFiles; j++) hs[j]->Draw("HIST SAME");
        c->cd(2);
        TLegend l(.6, .6, .9, .9);
        if(setYLog[i]) {
            gPad->SetLogy();
            whs[0]->GetYaxis()->SetRangeUser(1e-4, 1.1 * maxWhs);
        }else whs[0]->GetYaxis()->SetRangeUser(0, 1.1 * maxWhs);
        whs[0]->Draw();
        l.AddEntry(whs[0], files[0].c_str(), "l");
        for(int j = 1; j < nFiles; j++) {
            whs[j]->Draw("HIST SAME");
            l.AddEntry(whs[j], files[j].c_str(), "l");
        }
        l.DrawClone();
        c->SaveAs(("plots/" + var + ".png").c_str());
        delete c;
        for(int j = 0; j < nFiles; j++) {
            delete hs[j];
            delete whs[j];
        }
    }
    delete [] hs;
    delete [] whs;
}
