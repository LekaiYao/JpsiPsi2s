#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

using namespace std;
using namespace RooFit;

constexpr long double operator"" _PI(long double f) {
    return 3.14159265359 * f;
}

void Tplfit() {
    // Differential cross-section results
    TH1D *hDiff = new TH1D("Diff Xsection", "Diff Xsection", 5, 0., 1.0_PI);
    for(double i = 0.2; i < 1.0_PI; i += 0.2_PI) hDiff->Fill(i);
    RooRealVar evt_dPhi("evt_dPhi", "evt_dPhi", 0., 1.0_PI);
    RooDataHist dh_Diff("dh_Diff", "dh_Diff", evt_dPhi, hDiff);
    double xSec[] = {34.87, 9.74, 15.98, 8.06, 15.34}, xSecErr[] = {1.81, 1.04, 1.17, 1.04, 1.37};
    // double xSec[] = {22.26, 12.43, 10.73}, xSecErr[] = {3.98, 2.13, 2.26};
    for(int i = 0; i < 5; i++) {
        dh_Diff.get(i);
        dh_Diff.set(xSec[i], xSecErr[i]);
    }
    // SPS event yield
    TH1D *hSPS = new TH1D("SPS histo", "SPS histo", 5, 0., 1.0_PI);
    TFile spsFile("../UnweightSPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data");
    Double_t sps_Jpsi_phi, sps_psi2S_phi, sps_weight;
    spsTree->SetBranchAddress("delta_phi", &sps_Jpsi_phi);
    // spsTree->SetBranchAddress("psi2S_phi", &sps_psi2S_phi);
    spsTree->SetBranchAddress("evt_weight", &sps_weight);
    int nSPSEvt = spsTree->GetEntries();
    Double_t sps_sum = 0;
    // for(double i = 0; i < 2.6; i += 0.6) {
    //     string sel = "fabs(Jpsi_phi - psi2S_phi) < " + to_string(i + 0.6) + " && fabs(Jpsi_phi - psi2S_phi) > " + to_string(i);
    //     int nBin = spsTree->GetEntries(sel.c_str());
    //     for(int j = 0; j < nBin; j++) hSPS->Fill(i + 0.2);
    // }
    for(int i = 0; i < nSPSEvt; i++) {
        spsTree->GetEntry(i);
        hSPS->Fill(sps_Jpsi_phi, sps_weight);
        sps_sum += sps_weight;
    }
    hSPS->Scale(1.0 / sps_sum);
    RooDataHist dh_SPS("dh_SPS", "dh_SPS", evt_dPhi, hSPS);
    RooHistPdf pdf_SPS("pdf_SPS", "pdf_SPS", evt_dPhi, dh_SPS);
    // DPS event yield
    TH1D *hDPS = new TH1D("DPS histo", "DPS histo", 5, 0., 1.0_PI);
    TFile dpsFile("../UnweightDPS.root", "READ");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    Double_t dps_Jpsi_phi, dps_psi2S_phi, dps_weight;
    dpsTree->SetBranchAddress("delta_phi", &dps_Jpsi_phi);
    // dpsTree->SetBranchAddress("psi2S_phi", &dps_psi2S_phi);
    dpsTree->SetBranchAddress("evt_weight", &dps_weight);
    int nDPSEvt = dpsTree->GetEntries();
    Double_t dps_sum = 0;
    // for(double i = 0; i < 2.6; i += 0.6) {
    //     string sel = "fabs(Jpsi_phi - psi2S_phi) < " + to_string(i + 0.6) + " && fabs(Jpsi_phi - psi2S_phi) > " + to_string(i);
    //     int nBin = dpsTree->GetEntries(sel.c_str());
    //     for(int j = 0; j < nBin; j++) hDPS->Fill(i + 0.2);
    // }
    for(int i = 0; i < nDPSEvt; i++) {
        dpsTree->GetEntry(i);
        hDPS->Fill(dps_Jpsi_phi, dps_weight);
        dps_sum += dps_weight;
    }
    hDPS->Scale(1.0 / dps_sum);
    RooDataHist dh_DPS("dh_DPS", "dh_DPS", evt_dPhi, hDPS);
    RooHistPdf pdf_DPS("pdf_DPS", "pdf_DPS", evt_dPhi, dh_DPS);
    // Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
    RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(pdf_SPS, pdf_DPS), frac);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
    // Template fit
    pdf_all.fitTo(dh_Diff);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1300);
    RooPlot* frame = evt_dPhi.frame();
    dh_Diff.plotOn(frame, Name("Data"));
    pdf_all.plotOn(frame, LineColor(kBlack), Name("All"));
    pdf_all.plotOn(frame, Components(pdf_SPS), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    pdf_all.plotOn(frame, Components(pdf_DPS), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle("#Delta#phi");
    frame->SetYTitle("#frac{#partial#sigma}{#partial#Delta#phi}");
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs("fits/Diff_Xsec_DPhi.png");
    pdf_all.getVariables()->Print("v");
    cout<<frac.getValV()<<"+/-"<<frac.getError();
    return;
}
