#include "TCanvas.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGExpModel.h"
#include "RooCBShape.h"
// #include "RooCrystalBall.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
// #include "RooFFTConvPdf.h"

// #define JPSI
#define PSI2S
using namespace std;
using namespace RooFit;

void prefit() {
    Int_t N = 5000000, BinNum = 200;
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.55, 3.85);
	RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    // RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.02, 0.08);
	// RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.02, 0.02);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.03, 0.25);
	RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.03, 0.25);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
    // RooRealVar evt_vtxProb("evt_vtxProb", "evt_vtxProb", 0, 1);
    RooArgSet variables;
    variables.add(psi2S_mass);
    variables.add(Jpsi_mass);
    variables.add(psi2S_ctau);
    variables.add(Jpsi_ctau);
    variables.add(evt_weight);
    // variables.add(evt_vtxProb);
    // TFile *montFile = new TFile("WeightSPS.root", "READ");
    // TFile *montFile = new TFile("WeightDPS.root", "READ");
    // TFile *montFile = new TFile("WeightDPS_vali.root", "READ");
    TFile *montFile = new TFile("WeightMix.root", "READ");
    // TFile *montFile = new TFile("WeightBdecay.root", "READ");
    TTree *montTree = (TTree*)montFile->Get("data");
    RooDataSet *mont = new RooDataSet("mont", "mont", montTree, variables, "", "evt_weight");//evt_vtxProb != 0
    // RooDataSet *mont = new RooDataSet("mont", "mont", variables, Import(*montTree), WeightVar("evt_weight"));
    // Define J/psi Ctau p.d.f.
    #ifdef JPSI
    // prompt
    RooRealVar Jpsi_mu("Jpsi_mu", "Jpsi_mu", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.003, 0.003, 0.008);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", -1, -5, 0, "");
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 20);
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu, Jpsi_sigma2);
    // RooCBShape Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu, Jpsi_sigma2, Jpsi_alpha2, Jpsi_nx2);
    RooRealVar Jpsi_prop("Jpsi_prop", "Jpsi_prop", 0.5, 0, 1);
    RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop);
    // non-prompt
    // RooRealVar Jpsi_mu("Jpsi_mu", "Jpsi_mu", 0, -0.005, 0.005);
    // RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.001, 0, 0.003);
    // RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.003);
    // RooRealVar Jpsi_sigma4("Jpsi_sigma4", "Jpsi_sigma4", 0.001, 0.002, 0.008);
    // RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1, 0, 5);
    // RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 3, 1, 10);
    // RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.03, 0.01, 0.1);
    // RooRealVar Jpsi_coef2("Jpsi_coef2", "Jpsi_coef2", 0.03, 0.01, 0.1);
    // RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    // RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    // // RooCBShape Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu, Jpsi_sigma2, Jpsi_alpha2, Jpsi_nx2);
    // RooGExpModel Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // RooGExpModel Jpsi_gauss3("Jpsi_gauss3", "Jpsi_gauss3", Jpsi_ctau, Jpsi_sigma4, Jpsi_coef2, false, RooGExpModel::Type::Flipped);
    // RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss2, Jpsi_gauss3), Jpsi_prop1);
    // // RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2, Jpsi_gauss3), RooArgList(Jpsi_prop1, Jpsi_prop2));
    // // RooGExpModel JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    #endif

    // Define psi(2S) Ctau p.d.f.
    #ifdef PSI2S
    // prompt
    RooRealVar psi2S_mu("psi2S_mu", "psi2S_mu", 0, -1e-4, 1e-4);
    RooRealVar psi2S_sigma1("psi2S_sigma1", "psi2S_sigma1", 0.001, 0, 0.003);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.005, 0.003, 0.02);
    RooGaussian psi2S_gauss1("psi2S_gauss1", "psi2S_gauss1", psi2S_ctau, psi2S_mu, psi2S_sigma1);
    // RooGaussian psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_mu, psi2S_sigma2);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -1, -5, 0, "");
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 10);
    RooCBShape psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_mu, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    // RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1, 0, 20);
    // RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 20);
    // RooCrystalBall psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_mu, psi2S_sigma2, psi2S_sigma2, psi2S_alpha1, psi2S_nx1, psi2S_alpha2, psi2S_nx2);
    // RooRealVar psi2S_coef("psi2S_coef", "psi2S_coef", 0.06, 0.001, 0.1);
    // RooGExpModel psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_sigma2, psi2S_coef, false, RooGExpModel::Type::Flipped);
    RooRealVar psi2S_prop("psi2S_prop", "psi2S_prop", 0.5, 0, 1);
    RooAddPdf psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", RooArgList(psi2S_gauss1, psi2S_gauss2), psi2S_prop);
    // non-prompt
    // RooRealVar psi2S_mu("psi2S_mu", "psi2S_mu", 0, -0.005, 0.005);
    // RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.001, 0, 0.003);
    // RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.001, 0, 0.003);
    // RooRealVar psi2S_sigma4("psi2S_sigma4", "psi2S_sigma4", 0.001, 0.002, 0.008);
    // RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1, 0, 5);
    // RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 3, 1, 10);
    // RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.03, 0.01, 0.1);
    // RooRealVar psi2S_coef2("psi2S_coef2", "psi2S_coef2", 0.03, 0.01, 0.1);
    // RooRealVar psi2S_prop1("psi2S_prop1", "psi2S_prop1", 0.5, 0, 1);
    // RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.5, 0, 1);
    // RooGExpModel psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    // RooGExpModel psi2S_gauss3("psi2S_gauss3", "psi2S_gauss3", psi2S_ctau, psi2S_sigma4, psi2S_coef2, false, RooGExpModel::Type::Flipped);
    // RooAddPdf psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", RooArgList(psi2S_gauss2, psi2S_gauss3), psi2S_prop1);
    #endif

    #ifdef JPSI
    RooRealVar nSig("nSig", "nSig", 1, 0, N);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", JpsiCtauSig, nSig);
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*mont->reduce(Jpsi_ctau), Save());
        if(!res->status()) break;
    }
    // Draw data point and p.d.f. curve
    TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 2400, 1800);
    RooPlot *frame = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau 1D"), Bins(BinNum));
    mont->plotOn(frame, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame, LineColor(1), LineWidth(2), Name("All"));
    // pdf_all.plotOn(frame, Components(Jpsi_gauss1), LineColor(kBlue), LineWidth(2), Name("GS1"));
    // pdf_all.plotOn(frame, Components(Jpsi_gauss2), LineColor(kMagenta), LineWidth(2), Name("GS2"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "SPS MC 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    // legend->AddEntry(frame->findObject("GS1"), "Gaussian 1", "L");
    // legend->AddEntry(frame->findObject("GS2"), "Gaussian 2", "L");
    frame->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame->Draw();
    legend->DrawClone();
    canvas->SaveAs("fits/Mix_JpsiCtau.pdf");
    // canvas->SaveAs("fits/B_JpsiCtau.pdf");
    // canvas->SaveAs("fits/DPS_JpsiCtau.pdf");
    #endif

    #ifdef PSI2S
    RooRealVar nSig("nSig", "nSig", 1, 0, N);
    RooExtendPdf pdf_all("pdf_all", "pdf_all", psi2SCtauSig, nSig);
    RooFitResult *res;
    // res = pdf_all.fitTo(*mont->reduce(psi2S_ctau), Save());
    int x = 0;
    while(true && x++ < 20) {
        res = pdf_all.fitTo(*mont->reduce(psi2S_ctau), Save());
        if(!res->status()) break;
    }
    // Draw data point and p.d.f. curve
    TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 2400, 1800);
    RooPlot *frame = psi2S_ctau.frame(RooFit::Title("psi2S Ctau 1D"), Bins(BinNum));
    mont->plotOn(frame, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame, LineColor(1), LineWidth(2), Name("All"));
    // pdf_all.plotOn(frame, Components(psi2S_gauss1), LineColor(kBlue), LineWidth(2), Name("GS1"));
    // pdf_all.plotOn(frame, Components(psi2S_gauss2), LineColor(kMagenta), LineWidth(2), Name("GS2"));
    TLegend *legend = new TLegend(.65, .60, .85, .85);
    legend->AddEntry(frame->findObject("Data"), "SPS MC 2018", "L");
    legend->AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    // legend->AddEntry(frame->findObject("GS1"), "Gaussian 1", "L");
    // legend->AddEntry(frame->findObject("GS2"), "Gaussian 2", "L");
    frame->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame->Draw();
    legend->DrawClone();
    // canvas->SaveAs("fits/B_psi2SCtau.pdf");
    canvas->SaveAs("fits/Mix_psi2SCtau.pdf");
    // canvas->SaveAs("fits/DPS_psi2SCtau.pdf");
    #endif

    RooArgSet* params = pdf_all.getVariables();
    params->Print("v");
    cout<<endl<<res->status()<<endl;

    return;
}