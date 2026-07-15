// plot_mc_templates.C
// Compare SPS vs DPS GEN-level distributions for the 5 template-fit variables.
// Purpose: assess whether SPS/DPS shapes are distinguishable enough for f_DPS extraction.
//
// Usage: root -l -q 'plot_mc_templates.C'

#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include <iostream>
#include <cmath>

void plot_mc_templates() {
    gStyle->SetOptStat(0);

    // ── Input files ──────────────────────────────────────────
    const char* spsPath = "/eos/home-l/leyao/26JP/Ntuple/SPS/ggpsi1psi2/ggpsi1psi2.root";
    const char* dpsPath = "/eos/home-l/leyao/26JP/Ntuple/DPS/filter3p5/v2_set1.root";

    TFile *fsps = TFile::Open(spsPath);
    TFile *fdps = TFile::Open(dpsPath);
    if (!fsps || !fdps) { cerr << "File not found" << endl; return; }

    TTree *tsps = (TTree*)fsps->Get("rootuple/oniaTree");
    TTree *tdps = (TTree*)fdps->Get("rootuple/oniaTree");
    if (!tsps || !tdps) { cerr << "Tree not found" << endl; return; }

    // ── Branch addresses ─────────────────────────────────────
    UChar_t GEevt_valid, GEevt_passAcc;
    Double_t GEevt_fourMuMass;
    Double_t GEJpsi_pt, GEJpsi_eta, GEJpsi_phi, GEJpsi_mass;
    Double_t GEpsi2S_pt, GEpsi2S_eta, GEpsi2S_phi, GEpsi2S_mass;

    for (auto* t : {tsps, tdps}) {
        t->SetBranchAddress("GEevt_valid", &GEevt_valid);
        t->SetBranchAddress("GEevt_passAcc", &GEevt_passAcc);
        t->SetBranchAddress("GEevt_fourMuMass", &GEevt_fourMuMass);
        t->SetBranchAddress("GEJpsi_pt", &GEJpsi_pt);
        t->SetBranchAddress("GEJpsi_eta", &GEJpsi_eta);
        t->SetBranchAddress("GEJpsi_phi", &GEJpsi_phi);
        t->SetBranchAddress("GEJpsi_mass", &GEJpsi_mass);
        t->SetBranchAddress("GEpsi2S_pt", &GEpsi2S_pt);
        t->SetBranchAddress("GEpsi2S_eta", &GEpsi2S_eta);
        t->SetBranchAddress("GEpsi2S_phi", &GEpsi2S_phi);
        t->SetBranchAddress("GEpsi2S_mass", &GEpsi2S_mass);
    }

    // ── Binning (matching AN §7.1 / Tpl_Fit.cpp) ────────────
    const int nDeltaY = 6;
    double binsDeltaY[nDeltaY+1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0};
    const int nDeltaPhi = 5;
    double binsDeltaPhi[nDeltaPhi+1] = {0.0, 0.2*TMath::Pi(), 0.4*TMath::Pi(), 0.6*TMath::Pi(), 0.8*TMath::Pi(), 1.0*TMath::Pi()};
    const int nMass = 5;
    double binsMass[nMass+1] = {7.5, 17.5, 27.5, 37.5, 47.5, 107.5};
    const int nY = 5;
    double binsY[nY+1] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0};
    const int nPt = 6;
    double binsPt[nPt+1] = {0, 8, 16, 24, 32, 40, 80};

    // ── Histograms ───────────────────────────────────────────
    TH1D *hSPS_dY   = new TH1D("hSPS_dY",   ";|#Deltay|;Events (norm.)", nDeltaY, binsDeltaY);
    TH1D *hDPS_dY   = new TH1D("hDPS_dY",   ";|#Deltay|;Events (norm.)", nDeltaY, binsDeltaY);
    TH1D *hSPS_dPhi = new TH1D("hSPS_dPhi", ";|#Delta#phi|;Events (norm.)", nDeltaPhi, binsDeltaPhi);
    TH1D *hDPS_dPhi = new TH1D("hDPS_dPhi", ";|#Delta#phi|;Events (norm.)", nDeltaPhi, binsDeltaPhi);
    TH1D *hSPS_mass = new TH1D("hSPS_mass", ";M_{inv}^{J/#psi-#psi(2S)} [GeV];Events (norm.)", nMass, binsMass);
    TH1D *hDPS_mass = new TH1D("hDPS_mass", ";M_{inv}^{J/#psi-#psi(2S)} [GeV];Events (norm.)", nMass, binsMass);
    TH1D *hSPS_y    = new TH1D("hSPS_y",    ";|y^{J/#psi-#psi(2S)}|;Events (norm.)", nY, binsY);
    TH1D *hDPS_y    = new TH1D("hDPS_y",    ";|y^{J/#psi-#psi(2S)}|;Events (norm.)", nY, binsY);
    TH1D *hSPS_pt   = new TH1D("hSPS_pt",   ";p_{T}^{J/#psi-#psi(2S)} [GeV];Events (norm.)", nPt, binsPt);
    TH1D *hDPS_pt   = new TH1D("hDPS_pt",   ";p_{T}^{J/#psi-#psi(2S)} [GeV];Events (norm.)", nPt, binsPt);

    // ── Event loop ───────────────────────────────────────────
    struct { TTree* t; TH1D *h[5]; const char* label; } samples[] = {
        {tsps, {hSPS_dY, hSPS_dPhi, hSPS_mass, hSPS_y, hSPS_pt}, "SPS (HELAC-Onia LO)"},
        {tdps, {hDPS_dY, hDPS_dPhi, hDPS_mass, hDPS_y, hDPS_pt}, "DPS (Pythia8)"}
    };

    for (auto& s : samples) {
        int nEvt = s.t->GetEntries();
        int nPass = 0;
        for (int i = 0; i < nEvt; i++) {
            s.t->GetEntry(i);
            // Entry cut
            if (!GEevt_valid || !GEevt_passAcc || GEevt_fourMuMass <= 7.5) continue;
            // Fiducial
            if (GEJpsi_pt <= 10 || GEJpsi_pt >= 40 || fabs(GEJpsi_eta) >= 2.0) continue;  // using eta≈y at high pT
            if (GEpsi2S_pt <= 10 || GEpsi2S_pt >= 40 || fabs(GEpsi2S_eta) >= 2.0) continue;

            // Build four-vectors
            TLorentzVector jpsi, psi2S;
            jpsi.SetPtEtaPhiM(GEJpsi_pt, GEJpsi_eta, GEJpsi_phi, GEJpsi_mass);
            psi2S.SetPtEtaPhiM(GEpsi2S_pt, GEpsi2S_eta, GEpsi2S_phi, GEpsi2S_mass);
            TLorentzVector sys = jpsi + psi2S;

            // Differential variables
            double dY   = fabs(jpsi.Rapidity() - psi2S.Rapidity());
            double dPhi = fabs(TVector2::Phi_mpi_pi(jpsi.Phi() - psi2S.Phi()));
            double mass = sys.M();
            double Y    = fabs(sys.Rapidity());
            double pt   = sys.Pt();

            s.h[0]->Fill(dY);
            s.h[1]->Fill(dPhi);
            s.h[2]->Fill(mass);
            s.h[3]->Fill(Y);
            s.h[4]->Fill(pt);
            nPass++;
        }
        cout << s.label << ": " << nPass << " events pass fiducial cuts" << endl;
    }

    // ── Normalise to unit area ───────────────────────────────
    for (auto* h : {hSPS_dY, hDPS_dY, hSPS_dPhi, hDPS_dPhi, hSPS_mass, hDPS_mass,
                    hSPS_y, hDPS_y, hSPS_pt, hDPS_pt}) {
        if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
    }

    // ── Style ────────────────────────────────────────────────
    hSPS_dY->SetLineColor(kBlue);   hSPS_dY->SetLineWidth(2);
    hDPS_dY->SetLineColor(kRed);    hDPS_dY->SetLineWidth(2);
    hSPS_dPhi->SetLineColor(kBlue); hSPS_dPhi->SetLineWidth(2);
    hDPS_dPhi->SetLineColor(kRed);  hDPS_dPhi->SetLineWidth(2);
    hSPS_mass->SetLineColor(kBlue); hSPS_mass->SetLineWidth(2);
    hDPS_mass->SetLineColor(kRed);  hDPS_mass->SetLineWidth(2);
    hSPS_y->SetLineColor(kBlue);    hSPS_y->SetLineWidth(2);
    hDPS_y->SetLineColor(kRed);     hDPS_y->SetLineWidth(2);
    hSPS_pt->SetLineColor(kBlue);   hSPS_pt->SetLineWidth(2);
    hDPS_pt->SetLineColor(kRed);    hDPS_pt->SetLineWidth(2);

    // ── Draw ─────────────────────────────────────────────────
    const char* titles[] = {"|#Deltay|", "|#Delta#phi|", "M_{inv}", "|y|", "p_{T}"};
    TH1D* hSPS_arr[] = {hSPS_dY, hSPS_dPhi, hSPS_mass, hSPS_y, hSPS_pt};
    TH1D* hDPS_arr[] = {hDPS_dY, hDPS_dPhi, hDPS_mass, hDPS_y, hDPS_pt};
    const char* fnames[] = {"delta_y", "delta_phi", "evt_mass", "evt_y", "evt_pt"};

    TCanvas *c = new TCanvas("c", "c", 1000, 700);
    for (int i = 0; i < 5; i++) {
        hDPS_arr[i]->Draw("HIST");
        hSPS_arr[i]->Draw("HIST SAME");
        // Adjust ranges
        if (hSPS_arr[i]->GetMaximum() > hDPS_arr[i]->GetMaximum())
            hSPS_arr[i]->SetMaximum(hSPS_arr[i]->GetMaximum() * 1.3);
        else
            hDPS_arr[i]->SetMaximum(hDPS_arr[i]->GetMaximum() * 1.3);

        TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(hSPS_arr[i], "SPS (HELAC-Onia)", "l");
        leg->AddEntry(hDPS_arr[i], "DPS (Pythia8)", "l");
        leg->Draw();

        c->SaveAs(Form("fig/mc_template_%s.png", fnames[i]));
        delete leg;
    }

    cout << "\nPlots saved to fig/mc_template_*.png" << endl;

    fsps->Close();
    fdps->Close();
}
