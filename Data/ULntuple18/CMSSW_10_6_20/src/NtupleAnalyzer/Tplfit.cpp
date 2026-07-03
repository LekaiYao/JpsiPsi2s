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
#include "RooCategory.h"
#include "RooSimultaneous.h"

using namespace std;
using namespace RooFit;

constexpr long double operator"" _PI(long double f) {
    return 3.14159265359 * f;
}

void Tplfit() {
    // Differential cross-section results
    // DY
    int nDY = 6;
    double xDY[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0};
    TH1D *hDY = new TH1D("DY Xsection", "DY Xsection", nDY, xDY);
    for(int i = 0; i < nDY; i++) hDY->Fill(xDY[i] + 0.1);
    RooRealVar evt_dY("evt_dY", "evt_dY", 0., 4.);
    RooDataHist dh_dY("dh_dY", "dh_dY", evt_dY, hDY);
    double xSec_dY[] = {55.05, 26.73, 13.33, 12.72, 7.18, 2.23}, xSecErr_dY[] = {8.89, 5.05, 2.94, 2.82, 1.88, 0.6};
    for(int i = 0; i < nDY; i++) {
        dh_dY.get(i);
        dh_dY.set(xSec_dY[i], xSecErr_dY[i]);
    }
    // DPhi
    int nDPhi = 5;
    double xDPhi[] = {0, 0.2_PI, 0.4_PI, 0.6_PI, 0.8_PI, 1.0_PI};
    TH1D *hDPhi = new TH1D("DPhi Xsection", "DPhi Xsection", nDPhi, xDPhi);
    for(int i = 0; i < nDPhi; i++) hDPhi->Fill(xDPhi[i] + 0.1);
    RooRealVar evt_dPhi("evt_dPhi", "evt_dPhi", 0., 1.0_PI);
    RooDataHist dh_dPhi("dh_dPhi", "dh_dPhi", evt_dPhi, hDPhi);
    double xSec_dPhi[] = {34.75, 9.89, 16.09, 8.13, 15.51}, xSecErr_dPhi[] = {5.63, 2.13, 2.94, 1.94, 3.07};
    for(int i = 0; i < nDPhi; i++) {
        dh_dPhi.get(i);
        dh_dPhi.set(xSec_dPhi[i], xSecErr_dPhi[i]);
    }
    // SPS event yield
    TH1D *hSPS_dY = new TH1D("DY SPS histo", "DY SPS histo", nDY, xDY);
    TH1D *hSPS_dPhi = new TH1D("DPhi SPS histo", "DPhi SPS histo", nDPhi, xDPhi);
    TFile spsFile("UnweightSPS.root", "READ");
    TTree *spsTree = (TTree *)spsFile.Get("data");
    Double_t sps_dY, sps_dPhi, sps_weight;
    spsTree->SetBranchAddress("delta_y", &sps_dY);
    spsTree->SetBranchAddress("delta_phi", &sps_dPhi);
    spsTree->SetBranchAddress("evt_weight", &sps_weight);
    int nSPS = spsTree->GetEntries();
    Double_t sumSPS_dY = 0, sumSPS_dPhi = 0;
    for(int i = 0; i < nSPS; i++) {
        spsTree->GetEntry(i);
        if(sps_dY < 2) {
            hSPS_dY->Fill(sps_dY, sps_weight);
            sumSPS_dY += sps_weight;
        }else {
            hSPS_dY->Fill(sps_dY, sps_weight / 5);
            sumSPS_dY += sps_weight / 5;
        }
        hSPS_dPhi->Fill(sps_dPhi, sps_weight);
        sumSPS_dPhi += sps_weight;
    }
    hSPS_dY->Scale(1.0 / sumSPS_dY);
    hSPS_dPhi->Scale(1.0 / sumSPS_dPhi);
    RooDataHist dh_SPS_dY("dh_SPS_dY", "dh_SPS_dY", evt_dY, hSPS_dY);
    RooHistPdf pdf_SPS_dY("pdf_SPS_dY", "pdf_SPS_dY", evt_dY, dh_SPS_dY);
    RooDataHist dh_SPS_dPhi("dh_SPS_dPhi", "dh_SPS_dPhi", evt_dPhi, hSPS_dPhi);
    RooHistPdf pdf_SPS_dPhi("pdf_SPS_dPhi", "pdf_SPS_dPhi", evt_dPhi, dh_SPS_dPhi);
    // DPS event yield
    TH1D *hDPS_dY = new TH1D("DY DPS histo", "DY DPS histo", nDY, xDY);
    TH1D *hDPS_dPhi = new TH1D("DPhi DPS histo", "DPhi DPS histo", nDPhi, xDPhi);
    TFile dpsFile("UnweightDPS.root", "READ");
    TTree *dpsTree = (TTree *)dpsFile.Get("data");
    Double_t dps_dY, dps_dPhi, dps_weight;
    dpsTree->SetBranchAddress("delta_y", &dps_dY);
    dpsTree->SetBranchAddress("delta_phi", &dps_dPhi);
    dpsTree->SetBranchAddress("evt_weight", &dps_weight);
    int nDPS = dpsTree->GetEntries();
    Double_t sumDPS_dY = 0, sumDPS_dPhi = 0;
    for(int i = 0; i < nDPS; i++) {
        dpsTree->GetEntry(i);
        if(dps_dY < 2) {
            hDPS_dY->Fill(dps_dY, dps_weight);
            sumDPS_dY += dps_weight;
        }else {
            hDPS_dY->Fill(dps_dY, dps_weight / 5);
            sumDPS_dY += dps_weight / 5;
        }
        hDPS_dPhi->Fill(dps_dPhi, dps_weight);
        sumDPS_dPhi += dps_weight;
        
    }
    hDPS_dY->Scale(1.0 / sumDPS_dY);
    hDPS_dPhi->Scale(1.0 / sumDPS_dPhi);
    RooDataHist dh_DPS_dY("dh_DPS_dY", "dh_DPS_dY", evt_dY, hDPS_dY);
    RooHistPdf pdf_DPS_dY("pdf_DPS_dY", "pdf_DPS_dY", evt_dY, dh_DPS_dY);
    RooDataHist dh_DPS_dPhi("dh_DPS_dPhi", "dh_DPS_dPhi", evt_dPhi, hDPS_dPhi);
    RooHistPdf pdf_DPS_dPhi("pdf_DPS_dPhi", "pdf_DPS_dPhi", evt_dPhi, dh_DPS_dPhi);
    // Create category
    RooCategory cate("cate", "cate");
    cate.defineType("DY");
    cate.defineType("DPhi");
    RooDataHist combData("combData", "combData", RooArgList(evt_dY, evt_dPhi), Index(cate), Import("DY", dh_dY), Import("DPhi", dh_dPhi));
    // Construct p.d.f.
    RooRealVar frac("frac", "frac", 0.5, 0, 1);
    RooRealVar n_dY("n_dY", "n_dY", 10, 1, 1e4);
    RooRealVar n_dPhi("n_dPhi", "n_dPhi", 10, 1, 1e4);
    RooAddPdf pdf_dY("pdf_dY", "pdf_dY", RooArgList(pdf_SPS_dY, pdf_DPS_dY), frac);
    RooAddPdf pdf_dPhi("pdf_dPhi", "pdf_dPhi", RooArgList(pdf_SPS_dPhi, pdf_DPS_dPhi), frac);
    RooExtendPdf d_dY("d_dY", "d_dY", pdf_dY, n_dY);
    RooExtendPdf d_dPhi("d_dPhi", "d_dPhi", pdf_dPhi, n_dPhi);
    RooSimultaneous simPdf("simPdf", "simPdf", cate);
    simPdf.addPdf(d_dY, "DY");
    simPdf.addPdf(d_dPhi, "DPhi");
    // Template fit
    simPdf.fitTo(combData);
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 1500, 1300);
    RooPlot* frame1 = evt_dY.frame();
    dh_dY.plotOn(frame1, Name("Data"));
    // pdf_dY.plotOn(frame1, LineColor(kBlack), Name("All"));
    // pdf_dY.plotOn(frame1, Components(pdf_SPS_dY), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    // pdf_dY.plotOn(frame1, Components(pdf_DPS_dY), LineColor(kRed), LineStyle(kDotted), Name("DPS"));

    // Special plotting for variable binning
    RooHistFunc func_SPS("func_SPS", "func_SPS", evt_dY, dh_SPS_dY);
    RooHistFunc func_DPS("func_DPS", "func_DPS", evt_dY, dh_DPS_dY);
    RooRealVar coef_SPS{"coef_SPS", "coef_SPS", frac.getVal() * n_dY.getVal()};
    RooRealVar coef_DPS{"coef_DPS", "coef_DPS", (1 - frac.getVal()) * n_dY.getVal()};
    RooRealSumFunc func("func", "func", {func_SPS, func_DPS}, {coef_SPS, coef_DPS});
    func.plotOn(frame1, LineColor(kBlack), Name("All")); 
    func_SPS.plotOn(frame1, Normalization(coef_SPS.getVal()), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    func_DPS.plotOn(frame1, Normalization(coef_DPS.getVal()), LineColor(kRed), LineStyle(kDotted), Name("DPS"));

    TLegend *legend1 = new TLegend(.65, .60, .85, .85);
    legend1->AddEntry(frame1->findObject("Data"), "RunII 2018", "L");
    legend1->AddEntry(frame1->findObject("All"), "Total p.d.f.", "L");
    legend1->AddEntry(frame1->findObject("SPS"), "SPS p.d.f.", "L");
    legend1->AddEntry(frame1->findObject("DPS"), "DPS p.d.f.", "L");
    frame1->SetXTitle("#Deltay");
    frame1->SetYTitle("#frac{#partial#sigma}{#partial#Deltay}");
    frame1->Draw();
    legend1->DrawClone();
    canvas1->SaveAs("fits/Diff_Xsec_DY.png");
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 1500, 1300);
    RooPlot* frame2 = evt_dPhi.frame();
    dh_dPhi.plotOn(frame2, Name("Data"));
    pdf_dPhi.plotOn(frame2, LineColor(kBlack), Name("All"));
    pdf_dPhi.plotOn(frame2, Components(pdf_SPS_dPhi), LineColor(kBlue), LineStyle(kDashed), Name("SPS"));
    pdf_dPhi.plotOn(frame2, Components(pdf_DPS_dPhi), LineColor(kRed), LineStyle(kDotted), Name("DPS"));
    TLegend *legend2 = new TLegend(.65, .60, .85, .85);
    legend2->AddEntry(frame2->findObject("Data"), "RunII 2018", "L");
    legend2->AddEntry(frame2->findObject("All"), "Total p.d.f.", "L");
    legend2->AddEntry(frame2->findObject("SPS"), "SPS p.d.f.", "L");
    legend2->AddEntry(frame2->findObject("DPS"), "DPS p.d.f.", "L");
    frame2->SetXTitle("#Delta#phi");
    frame2->SetYTitle("#frac{#partial#sigma}{#partial#Delta#phi}");
    frame2->Draw();
    legend2->DrawClone();
    canvas2->SaveAs("fits/Diff_Xsec_DPhi.png");
    simPdf.getVariables()->Print("v");
    return;
}