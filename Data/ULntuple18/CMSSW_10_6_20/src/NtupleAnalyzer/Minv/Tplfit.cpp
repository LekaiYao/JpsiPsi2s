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
    TH1D *hDiff = new TH1D("Diff Xsection", "Diff Xsection", 6, 7.5, 67.5);
    for(int i = 10; i < 67.5; i+=10) hDiff->Fill(i);
    RooRealVar evt_mass("evt_mass", "evt_mass", 7.5, 67.5);
    RooDataHist dh_Diff("dh_Diff", "dh_Diff", evt_mass, hDiff);
    double xSec[] = {529.46, 116.305, 197.713, 61.2989, 42.2767, 0.881255}, xSecErr[] = {28.9165, 16.0057, 17.6876, 12.5414, 9.09304, 5.66849};
    for(int i = 0; i < 6; i++) {
        dh_Diff.get(i);
        dh_Diff.set(xSec[i] / 29.23 / 10, xSecErr[i] / 29.23 / 10);
    }
    // SPS event yield
    TH1D *hSPS = new TH1D("SPS histo", "SPS histo", 6, 7.5, 67.5);
    TFile spsFile("../WeightSPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data");
    Double_t sps_mass, sps_weight;
    spsTree->SetBranchAddress("evt_mass", &sps_mass);
    spsTree->SetBranchAddress("evt_weight", &sps_weight);
    int nSPSEvt = spsTree->GetEntries();
    Double_t sps_sum = 0;
    for(int i = 0; i < nSPSEvt; i++) {
        spsTree->GetEntry(i);
        hSPS->Fill(sps_mass, sps_weight);
        sps_sum += sps_weight;
    }
    hSPS->Scale(1.0 / sps_sum);
    RooDataHist dh_SPS("dh_SPS", "dh_SPS", evt_mass, hSPS);
    RooHistPdf pdf_SPS("pdf_SPS", "pdf_SPS", evt_mass, dh_SPS);
    // DPS event yield
    TH1D *hDPS = new TH1D("DPS histo", "DPS histo", 6, 7.5, 67.5);
    TFile dpsFile("../WeightDPS.root", "READ");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    Double_t dps_mass, dps_weight;
    dpsTree->SetBranchAddress("evt_mass", &dps_mass);
    dpsTree->SetBranchAddress("evt_weight", &dps_weight);
    int nDPSEvt = dpsTree->GetEntries();
    Double_t dps_sum = 0;
    for(int i = 0; i < nDPSEvt; i++) {
        dpsTree->GetEntry(i);
        hDPS->Fill(dps_mass, dps_weight);
        dps_sum += dps_weight;
    }
    hDPS->Scale(1.0 / dps_sum);
    RooDataHist dh_DPS("dh_DPS", "dh_DPS", evt_mass, hDPS);
    RooHistPdf pdf_DPS("pdf_DPS", "pdf_DPS", evt_mass, dh_DPS);
    // Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
    RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(pdf_SPS, pdf_DPS), frac);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
    // Template fit
    pdf_all.fitTo(dh_Diff);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    RooPlot* frame = evt_mass.frame();
    dh_Diff.plotOn(frame);
    pdf_all.plotOn(frame, LineColor(kBlack));
    pdf_all.plotOn(frame, Components(pdf_SPS), LineColor(kBlue), LineStyle(kDashed));
    pdf_all.plotOn(frame, Components(pdf_DPS), LineColor(kRed), LineStyle(kDotted));
    frame->Draw();
    canvas->SaveAs("fits/Diff_Xsec_Minv.pdf");
    pdf_all.getVariables()->Print("v");
    return;
}