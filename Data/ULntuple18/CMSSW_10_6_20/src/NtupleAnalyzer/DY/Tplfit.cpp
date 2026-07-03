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
    int nDY = 6;
    double xDY[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0};
    TH1D *hDiff = new TH1D("Diff Xsection", "Diff Xsection", nDY, xDY);
    // for(int i = 0; i < nDY; i++) hDiff->Fill(xDY[i] + 0.1);
    RooRealVar evt_dY("evt_dY", "evt_dY", 0., 4.);
    RooDataHist dh_Diff("dh_Diff", "dh_Diff", evt_dY, hDiff);
    // double xSec[] = {55.05, 26.73, 13.33, 12.72, 7.18, 2.23}, xSecErr[] = {8.89, 5.05, 2.94, 2.82, 1.88, 0.6};
    double xSec[] = {16.12, 7.58, 1, 5.50, 1, 0}, xSecErr[] = {4.97, 3.12, 0.5, 2.8, 0.5, 0};
    for(int i = 0; i < nDY; i++) {
        dh_Diff.get(i);
        dh_Diff.set(xSec[i], xSecErr[i]);
    }
    // SPS event yield
    TH1D *hSPS = new TH1D("SPS histo", "SPS histo", nDY, xDY);
    TFile spsFile("../UnweightSPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data");
    Double_t sps_delta_y, sps_weight, sps_Jpsi_pt, sps_Jpsi_y, sps_psi2S_pt, sps_psi2S_y;
    spsTree->SetBranchAddress("delta_y", &sps_delta_y);
    spsTree->SetBranchAddress("evt_weight", &sps_weight);
    spsTree->SetBranchAddress("Jpsi_pt", &sps_Jpsi_pt);
    spsTree->SetBranchAddress("Jpsi_y", &sps_Jpsi_y);
    spsTree->SetBranchAddress("psi2S_pt", &sps_psi2S_pt);
    spsTree->SetBranchAddress("psi2S_y", &sps_psi2S_y);
    int nSPSEvt = spsTree->GetEntries();
    Double_t sps_sum = 0;
    for(int i = 0; i < nSPSEvt; i++) {
        spsTree->GetEntry(i);
        if(sps_Jpsi_pt < 20 || sps_psi2S_pt < 20 || fabs(sps_Jpsi_y) > 1.2 || fabs(sps_psi2S_y) > 1.2) continue;
        if(sps_delta_y < 2) {
            hSPS->Fill(sps_delta_y, sps_weight);
            sps_sum += sps_weight;
        }else {
            hSPS->Fill(sps_delta_y, sps_weight / 5);
            sps_sum += sps_weight;// / 5;
        }
        
    }
    hSPS->Scale(1.0 / sps_sum);
    RooDataHist dh_SPS("dh_SPS", "dh_SPS", evt_dY, hSPS);
    RooHistPdf pdf_SPS("pdf_SPS", "pdf_SPS", evt_dY, dh_SPS);
    // DPS event yield
    TH1D *hDPS = new TH1D("DPS histo", "DPS histo", nDY, xDY);
    TFile dpsFile("../UnweightDPS.root", "READ");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    Double_t dps_delta_y, dps_weight, dps_Jpsi_pt, dps_Jpsi_y, dps_psi2S_pt, dps_psi2S_y;
    dpsTree->SetBranchAddress("delta_y", &dps_delta_y);
    dpsTree->SetBranchAddress("evt_weight", &dps_weight);
    dpsTree->SetBranchAddress("Jpsi_pt", &dps_Jpsi_pt);
    dpsTree->SetBranchAddress("Jpsi_y", &dps_Jpsi_y);
    dpsTree->SetBranchAddress("psi2S_pt", &dps_psi2S_pt);
    dpsTree->SetBranchAddress("psi2S_y", &dps_psi2S_y);
    int nDPSEvt = dpsTree->GetEntries();
    Double_t dps_sum = 0;
    for(int i = 0; i < nDPSEvt; i++) {
        dpsTree->GetEntry(i);
        if(dps_Jpsi_pt < 20 || dps_psi2S_pt < 20 || fabs(dps_Jpsi_y) > 1.2 || fabs(dps_psi2S_y) > 1.2) continue;
        if(dps_delta_y < 2) {
            hDPS->Fill(dps_delta_y, dps_weight);
            dps_sum += dps_weight;
        }else {
            hDPS->Fill(dps_delta_y, dps_weight / 5);
            dps_sum += dps_weight;// / 5;
        }
    }
    hDPS->Scale(1.0 / dps_sum);
    RooDataHist dh_DPS("dh_DPS", "dh_DPS", evt_dY, hDPS);
    RooHistPdf pdf_DPS("pdf_DPS", "pdf_DPS", evt_dY, dh_DPS);
    // Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar n_sps_dps("n_sps_dps", "n_sps_dps", 10, 1, 1e4);
    RooAddPdf pdf_sps_dps("pdf_sps_dps", "pdf_sps_dps", RooArgList(pdf_SPS, pdf_DPS), frac);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", pdf_sps_dps, n_sps_dps);
    // Template fit
    pdf_all.fitTo(dh_Diff);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1500, 1300);
    RooPlot* frame = evt_dY.frame();
    dh_Diff.plotOn(frame, Name("Data"));

    RooHistFunc func_SPS("func_SPS", "func_SPS", evt_dY, dh_SPS);
    RooHistFunc func_DPS("func_DPS", "func_DPS", evt_dY, dh_DPS);
    RooRealVar coef_SPS{"coef_SPS", "coef_SPS", frac.getVal() * n_sps_dps.getVal()};
    RooRealVar coef_DPS{"coef_DPS", "coef_DPS", (1 - frac.getVal()) * n_sps_dps.getVal()};
    RooRealSumFunc func("func", "func", {func_SPS, func_DPS}, {coef_SPS, coef_DPS});
    func.plotOn(frame, LineColor(kBlack), Name("All")); 
    func_SPS.plotOn(frame, Normalization(coef_SPS.getVal()), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    func_DPS.plotOn(frame, Normalization(coef_DPS.getVal()), LineColor(kRed), LineStyle(kDotted), Name("DPS"));

    // pdf_all.plotOn(frame, LineColor(kBlack), Name("All"));
    // pdf_all.plotOn(frame, Components(pdf_SPS), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    // pdf_all.plotOn(frame, Components(pdf_DPS), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "RunII 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend->AddEntry(frame->findObject("SPS"), "SPS p.d.f.", "L");
    legend->AddEntry(frame->findObject("DPS"), "DPS p.d.f.", "L");
    frame->SetXTitle("#Deltay");
    frame->SetYTitle("#frac{#partial#sigma}{#partial#Deltay}");
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs("fits/Diff_Xsec_DY.png");
    pdf_all.getVariables()->Print("v");
    return;
}
