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

void Tplfit() {
    // Differential cross-section results
    TH1D *hDiff = new TH1D("Diff Xsection", "Diff Xsection", 6, 0., 2.4);
    for(double i = 0.2; i < 2.4; i += 0.4) hDiff->Fill(i);
    RooRealVar evt_dEta("evt_dEta", "evt_dEta", 0., 2.4);
    RooDataHist dh_Diff("dh_Diff", "dh_Diff", evt_dEta, hDiff);
    double xSec[] = {57.56, 24.00, 13.7, 11.62, 10.06, 6.82}, xSecErr[] = {8.89, 4.85, 3.05, 2.68, 2.37, 2.15};
    for(int i = 0; i < 6; i++) {
        dh_Diff.get(i);
        dh_Diff.set(xSec[i], xSecErr[i]);
    }
    // SPS event yield
    TH1D *hSPS = new TH1D("SPS histo", "SPS histo", 6, 0., 2.4);
    TFile spsFile("../WeightSPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data");
    Double_t sps_Jpsi_eta, sps_psi2S_eta, sps_weight;
    spsTree->SetBranchAddress("Jpsi_eta", &sps_Jpsi_eta);
    spsTree->SetBranchAddress("psi2S_eta", &sps_psi2S_eta);
    spsTree->SetBranchAddress("evt_weight", &sps_weight);
    int nSPSEvt = spsTree->GetEntries();
    Double_t sps_sum = 0;
    // for(double i = 0; i < 2.2; i += 0.4) {
    //     string sel = "fabs(Jpsi_eta - psi2S_eta) < " + to_string(i + 0.4) + " && fabs(Jpsi_eta - psi2S_eta) > " + to_string(i);
    //     int nBin = spsTree->GetEntries(sel.c_str());
    //     for(int j = 0; j < nBin; j++) hSPS->Fill(i + 0.2);
    // }
    for(int i = 0; i < nSPSEvt; i++) {
        spsTree->GetEntry(i);
        hSPS->Fill(fabs(sps_Jpsi_eta - sps_psi2S_eta), sps_weight);
        sps_sum += sps_weight;
    }
    hSPS->Scale(1.0 / sps_sum);
    RooDataHist dh_SPS("dh_SPS", "dh_SPS", evt_dEta, hSPS);
    RooHistPdf pdf_SPS("pdf_SPS", "pdf_SPS", evt_dEta, dh_SPS);
    // DPS event yield
    TH1D *hDPS = new TH1D("DPS histo", "DPS histo", 6, 0., 2.4);
    TFile dpsFile("../WeightDPS.root", "READ");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    Double_t dps_Jpsi_eta, dps_psi2S_eta, dps_weight;
    dpsTree->SetBranchAddress("Jpsi_eta", &dps_Jpsi_eta);
    dpsTree->SetBranchAddress("psi2S_eta", &dps_psi2S_eta);
    dpsTree->SetBranchAddress("evt_weight", &dps_weight);
    int nDPSEvt = dpsTree->GetEntries();
    Double_t dps_sum = 0;
    // for(double i = 0; i < 2.2; i += 0.4) {
    //     string sel = "fabs(Jpsi_eta - psi2S_eta) < " + to_string(i + 0.4) + " && fabs(Jpsi_eta - psi2S_eta) > " + to_string(i);
    //     int nBin = dpsTree->GetEntries(sel.c_str());
    //     for(int j = 0; j < nBin; j++) hDPS->Fill(i + 0.2);
    // }
    for(int i = 0; i < nDPSEvt; i++) {
        dpsTree->GetEntry(i);
        hDPS->Fill(fabs(dps_Jpsi_eta - dps_psi2S_eta), dps_weight);
        dps_sum += dps_weight;
    }
    hDPS->Scale(1.0 / dps_sum);
    RooDataHist dh_DPS("dh_DPS", "dh_DPS", evt_dEta, hDPS);
    RooHistPdf pdf_DPS("pdf_DPS", "pdf_DPS", evt_dEta, dh_DPS);
    // Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
    RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(pdf_SPS, pdf_DPS), frac);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
    // Template fit
    pdf_all.fitTo(dh_Diff);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1000);
    RooPlot* frame = evt_dEta.frame();
    dh_Diff.plotOn(frame, Name("Data"));
    pdf_all.plotOn(frame, LineColor(kBlack), Name("All"));
    pdf_all.plotOn(frame, Components(pdf_SPS), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    pdf_all.plotOn(frame, Components(pdf_DPS), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle("#Delta#eta");
    frame->SetYTitle("#frac{#partial#sigma}{#partial#Delta#eta}");
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs("fits/Diff_Xsec_DEta.pdf");
    pdf_all.getVariables()->Print("v");
    return;
}