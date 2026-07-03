#include "TCanvas.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
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
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooFitResult.h"

using namespace std;
using namespace RooFit;
// Switches for fitting process
// #define JPSI_MASS_1D
// #define PSI2S_MASS_1D
// #define JPSI_CTAU_SELF_SIDE
// #define PSI2S_CTAU_SELF_SIDE
// #define JPSI_CTAU_ELSE_SIDE
// #define PSI2S_CTAU_ELSE_SIDE
// #define ALL_4D

void fit() {
    Int_t N = 1500000, BinNum = 100;
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
    // RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.55, 3.85);
	RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.02, 0.16);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.02, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_ctau);
    variables.add(psi2S_ctau);
    variables.add(evt_weight);
    TFile *dataFile = new TFile("WeightData_prime.root", "READ");
    // TFile *dataFile = new TFile("WeightDPS.root", "READ");
    // TFile *dataFile = new TFile("WeightSPS.root", "READ");
    // TFile *dataFile = new TFile("WeightBdecaySmall.root", "READ");
    // TFile *dataFile = new TFile("WeightMix.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
    string selection = "1";
    RooFitResult *res;
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    #ifdef JPSI_MASS_1D
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.02);
    RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 1, 0.1, 5);
    RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 1, 0, 15);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.06);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 50);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    #else
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.09254);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.0236787);
    RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 3.63886);
    RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 8.44979);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.0695099);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.31692);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1.49545);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.644038);
    #endif
    RooCBShape Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    // RooGaussian Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    // RooGaussian Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    // RooCBShape JpsiMassSig("JpsiMassSig", "JpsiMassSig", Jpsi_mass, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    // Background p.d.f.
    #ifdef JPSI_MASS_1D
    // RooRealVar Jpsi_a("Jpsi_a", "Jpsi_a", 0, -5, 5);
    // RooRealVar Jpsi_b("Jpsi_b", "Jpsi_b", 0, -1.5, 1);
    #else
    // RooRealVar Jpsi_a("Jpsi_a", "Jpsi_a", 0.055558);
    // RooRealVar Jpsi_b("Jpsi_b", "Jpsi_b", -1.29935);
    #endif
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());
    // Pre-fit J/psi Mass
    #ifdef JPSI_MASS_1D
    RooRealVar nJpsiMassSig("nJpsiMassSig", "nJpsiMassSig", 1, 0, N);
    RooRealVar nJpsiMassComb("nJpsiMassComb", "nJpsiMassComb", 1, 0, N);
    RooExtendPdf pdf_JpsiMassSig("pdf_JpsiMassSig", "pdf_JpsiMassSig", JpsiMassSig, nJpsiMassSig);
    RooExtendPdf pdf_JpsiMassComb("pdf_JpsiMassComb", "pdf_JpsiMassComb", JpsiMassComb, nJpsiMassComb);
    RooAddPdf pdf_JpsiMass("pdf_JpsiMass", "pdf_JpsiMass", RooArgList(pdf_JpsiMassSig, pdf_JpsiMassComb));
    // RooAddPdf pdf_JpsiMass("pdf_JpsiMass", "pdf_JpsiMass", RooArgList(pdf_JpsiMassSig));
    while(true) {
        res = pdf_JpsiMass.fitTo(*data->reduce(Jpsi_mass, selection.c_str()), Save());
        if(!res->status()) break;
    }
    cout<<"Status: "<<res->status()<<endl;
    // Draw data point and p.d.f. curve
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 10, 10, 2400, 1800);
    RooPlot *frame1 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 1D"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame1, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiMass.plotOn(frame1, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiMass.plotOn(frame1, Components(Jpsi_crysBall1), LineColor(kBlue), LineWidth(2), Name("GS"));
    pdf_JpsiMass.plotOn(frame1, Components(Jpsi_crysBall2), LineColor(kRed), LineWidth(2), Name("CB"));
    pdf_JpsiMass.plotOn(frame1, Components(JpsiMassComb), LineColor(kMagenta), LineWidth(2), Name("U"));
    TLegend *legend1 = new TLegend(.65, .60, .85, .85);
    legend1->AddEntry(frame1->findObject("Data"), "RunII 2018", "L");
    legend1->AddEntry(frame1->findObject("All"), "Total p.d.f.", "L");
    legend1->AddEntry(frame1->findObject("GS"), "Gaussian", "L");
    legend1->AddEntry(frame1->findObject("CB"), "Crystall Ball", "L");
    legend1->AddEntry(frame1->findObject("U"), "Uniform", "L");
    frame1->Draw();
    legend1->DrawClone();
    canvas1->SaveAs("fits/1D_JpsiMass.pdf");
    pdf_JpsiMass.getVariables()->Print("v");
    #endif
    // Define psi(2S) Mass p.d.f.
    // Signal p.d.f.
    #ifdef PSI2S_MASS_1D
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.686, 3.63, 3.73);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.02, 0.01, 0.05);
    // RooRealVar psi2S_alpha1("psi2S_alpha1", "psi2S_alpha1", 1.5, 0.1, 15);
    // RooRealVar psi2S_nx1("psi2S_nx1", "psi2S_nx1", 3, 0.1, 15);
    // RooRealVar psi2S_devia2("psi2S_devia2", "psi2S_devia2", 0.07, 0.04, 0.2);
    // RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1.5, 0.1, 15);
    // RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 3, 0.1, 15);
    // RooRealVar psi2S_ratio("psi2S_ratio", "psi2S_ratio", 0.5, 0, 1);
    #else
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.68307);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.0330251);
    RooRealVar psi2S_alpha1("psi2S_alpha1", "psi2S_alpha1", 4.57881);
    RooRealVar psi2S_nx1("psi2S_nx1", "psi2S_nx1", 9.97682);
    RooRealVar psi2S_devia2("psi2S_devia2", "psi2S_devia2", 0.0799988);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1.86426);
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 0.500042);
    RooRealVar psi2S_ratio("psi2S_ratio", "psi2S_ratio", 0.830609);
    #endif
    // RooCBShape psi2S_crysBall1("psi2S_crysBall1", "psi2S_crysBall1", psi2S_mass, psi2S_mean, psi2S_devia1, psi2S_alpha1, psi2S_nx1);
    // RooGaussian psi2S_crysBall1("psi2S_crysBall1", "psi2S_crysBall1", psi2S_mass, psi2S_mean, psi2S_devia1);
    // RooCBShape psi2S_crysBall2("psi2S_crysBall2", "psi2S_crysBall2", psi2S_mass, psi2S_mean, psi2S_devia2, psi2S_alpha2, psi2S_nx2);
    // RooGaussian psi2S_crysBall2("psi2S_crysBall2", "psi2S_crysBall2", psi2S_mass, psi2S_mean, psi2S_devia2);
    // RooAddPdf psi2SMassSig("psi2SMassSig", "psi2SMassSig", RooArgList(psi2S_crysBall1, psi2S_crysBall2), psi2S_ratio);
    // RooCBShape psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1, psi2S_alpha1, psi2S_nx1);
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    // Background p.d.f.
    #ifdef PSI2S_MASS_1D
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    #else
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.450761);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0.0558413);
    #endif
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));
    // Pre-fit psi(2S) Mass
    #ifdef PSI2S_MASS_1D
    RooRealVar npsi2SMassSig("npsi2SMassSig", "npsi2SMassSig", 1, 0, N);
    RooRealVar npsi2SMassComb("npsi2SMassComb", "npsi2SMassComb", 1, 0, N);
    RooExtendPdf pdf_psi2SMassSig("pdf_psi2SMassSig", "pdf_psi2SMassSig", psi2SMassSig, npsi2SMassSig);
    RooExtendPdf pdf_psi2SMassComb("pdf_psi2SMassComb", "pdf_psi2SMassComb", psi2SMassComb, npsi2SMassComb);
    RooAddPdf pdf_psi2SMass("pdf_psi2SMass", "pdf_psi2SMass", RooArgList(pdf_psi2SMassSig, pdf_psi2SMassComb));
    // RooAddPdf pdf_psi2SMass("pdf_psi2SMass", "pdf_psi2SMass", RooArgList(pdf_psi2SMassSig));
    while(true) {
        res = pdf_psi2SMass.fitTo(*data->reduce(psi2S_mass, selection.c_str()), Save());
        if(!res->status()) break;
    }
    cout<<"Status: "<<res->status()<<endl;
    // Draw data point and p.d.f. curve
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 10, 10, 2400, 1800);
    RooPlot *frame2 = psi2S_mass.frame(RooFit::Title("psi2S Mass 1D"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame2, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SMass.plotOn(frame2, LineColor(1), LineWidth(2), Name("All"));
    // pdf_psi2SMass.plotOn(frame2, Components(psi2S_crysBall1), LineColor(kBlue), LineWidth(2), Name("CB"));
    // pdf_psi2SMass.plotOn(frame2, Components(psi2S_crysBall2), LineColor(kMagenta), LineWidth(2), Name("GS"));
    pdf_psi2SMass.plotOn(frame2, Components(psi2SMassSig), LineColor(kBlue), LineWidth(2), Name("CB"));
    pdf_psi2SMass.plotOn(frame2, Components(psi2SMassComb), LineColor(kMagenta), LineWidth(2), Name("CH"));
    TLegend *legend2 = new TLegend(.65, .60, .85, .85);
    legend2->AddEntry(frame2->findObject("Data"), "RunII 2018", "L");
    legend2->AddEntry(frame2->findObject("All"), "Total p.d.f.", "L");
    // legend2->AddEntry(frame2->findObject("CB"), "Crystall Ball 1", "L");
    // legend2->AddEntry(frame2->findObject("GS"), "Crystall Ball 2", "L");
    legend2->AddEntry(frame2->findObject("CB"), "Gaussian", "L");
    legend2->AddEntry(frame2->findObject("CH"), "Chebyshev", "L");
    frame2->Draw();
    legend2->DrawClone();
    canvas2->SaveAs("fits/1D_psi2SMass.pdf");
    pdf_psi2SMass.getVariables()->Print("v");
    #endif
    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 2.44665e-05);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.00144558);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.00354066);
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.817505);
    RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooGExpModel JpsiCtauBkg("JpsiCtauBkg", "JpsiCtauBkg", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar Jpsi_mu2("Jpsi_mu2", "Jpsi_mu2", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma4("Jpsi_sigma4", "Jpsi_sigma4", 0.005, 0, 0.02);
    RooRealVar Jpsi_sigma5("Jpsi_sigma5", "Jpsi_sigma5", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef2("Jpsi_coef2", "Jpsi_coef2", 0.06, 0.01, 0.1);
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    RooGaussian Jpsi_gauss3("Jpsi_gauss3", "Jpsi_gauss3", Jpsi_ctau, Jpsi_mu2, Jpsi_sigma4);
    RooGExpModel Jpsi_expGau2("Jpsi_expGau2", "Jpsi_expGau2", Jpsi_ctau, Jpsi_sigma5, Jpsi_coef2, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCombSig("JpsiCtauCombSig", "JpsiCtauCombSig", RooArgList(Jpsi_gauss3, Jpsi_expGau2), Jpsi_prop2);
    // Combinatorial background p.d.f.
    RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.1);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau, Jpsi_mu3, Jpsi_sigma6);
    RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_ctau, Jpsi_sigma6, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCombBkg("JpsiCtauCombBkg", "JpsiCtauCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);
    
    #ifdef JPSI_CTAU_SELF_SIDE
    // Fit J/psi Ctau in J/psi sideband
    selection = "Jpsi_mass > 3.19 || Jpsi_mass < 2.99";
    RooRealVar nJpsiCtauCombBkg("nJpsiCtauCombBkg", "nJpsiCtauCombBkg", 1, 0, N);
    RooExtendPdf pdf_JpsiCtauCombBkg("pdf_JpsiCtauCombBkg", "pdf_JpsiCtauCombBkg", JpsiCtauCombBkg, nJpsiCtauCombBkg);
    pdf_JpsiCtauCombBkg.fitTo(*data->reduce(Jpsi_ctau, selection.c_str()), Save());
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 10, 10, 2400, 1800);
    RooPlot *frame3 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau Self-sideband"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame3, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiCtauCombBkg.plotOn(frame3, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiCtauCombBkg.plotOn(frame3, Components(Jpsi_gauss4), LineColor(kBlue), LineWidth(2), Name("GS"));
    pdf_JpsiCtauCombBkg.plotOn(frame3, Components(Jpsi_expGau3), LineColor(kMagenta), LineWidth(2), Name("EX"));
    TLegend *legend3 = new TLegend(.65, .60, .85, .85);
    legend3->AddEntry(frame3->findObject("Data"), "RunII 2018", "L");
    legend3->AddEntry(frame3->findObject("All"), "Total p.d.f.", "L");
    legend3->AddEntry(frame3->findObject("GS"), "Gaussian", "L");
    legend3->AddEntry(frame3->findObject("EX"), "Exp#otimes Gauss", "L");
    frame3->SetAxisRange(1, 5e2, "Y");
    gPad->SetLogy();
    frame3->Draw();
    legend3->DrawClone();
    canvas3->SaveAs("fits/Side0_JpsiCtau.pdf");
    pdf_JpsiCtauCombBkg.getVariables()->Print("v");
    // Fix parameters as constant
    Jpsi_mu3.setConstant();
    Jpsi_sigma6.setConstant();
    Jpsi_coef3.setConstant();
    Jpsi_prop3.setConstant();
    #endif
    #ifdef JPSI_CTAU_ELSE_SIDE
    // Fit J/psi Ctau in psi(2S) sideband
    selection = "psi2S_mass > 3.85 || psi2S_mass < 3.55";
    RooRealVar nJpsiCtauCombSig("nJpsiCtauCombSig", "nJpsiCtauCombSig", 1, 0, N);
    RooExtendPdf pdf_JpsiCtauCombSig("pdf_JpsiCtauCombSig", "pdf_JpsiCtauCombSig", JpsiCtauCombSig, nJpsiCtauCombSig);
    RooAddPdf pdf_JpsiCtauComb("pdf_JpsiCtauComb", "pdf_JpsiCtauComb", RooArgList(pdf_JpsiCtauCombSig));
    pdf_JpsiCtauComb.fitTo(*data->reduce(Jpsi_ctau, selection.c_str()), Save());
    // Draw data point and p.d.f. curve
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 10, 10, 800, 600);
    RooPlot *frame4 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau Else-sideband"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame4, DataError(RooAbsData::SumW2));
    pdf_JpsiCtauComb.plotOn(frame4);
    pdf_JpsiCtauComb.plotOn(frame4, Components(Jpsi_gauss3), LineColor(1));
    pdf_JpsiCtauComb.plotOn(frame4, Components(Jpsi_expGau2), LineColor(3));
    // pdf_JpsiCtauComb.plotOn(frame4, Components(JpsiCtauCombBkg), LineColor(5));
    // frame4->SetAxisRange(1, 2e3, "Y");
    // gPad->SetLogy();
    frame4->Draw();
    canvas4->SaveAs("fits/Side1_JpsiCtau.pdf");
    pdf_JpsiCtauComb.getVariables()->Print("v");
    // Fix parameters as constant
    Jpsi_mu2.setConstant();
    Jpsi_sigma4.setConstant();
    Jpsi_sigma5.setConstant();
    Jpsi_coef2.setConstant();
    Jpsi_prop2.setConstant();
    // Jpsi_mu3.setConstant(kFALSE);
    // Jpsi_sigma6.setConstant(kFALSE);
    // Jpsi_coef3.setConstant(kFALSE);
    // Jpsi_prop3.setConstant(kFALSE);
    #endif
    // Define psi(2S) Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", -2.42944e-05);
    RooRealVar psi2S_sigma1("psi2S_sigma1", "psi2S_sigma1", 0.00137595);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.00316465);
    RooGaussian psi2S_gauss1("psi2S_gauss1", "psi2S_gauss1", psi2S_ctau, psi2S_mu1, psi2S_sigma1);
    RooGaussian psi2S_gauss2("psi2S_gauss2", "psi2S_gauss2", psi2S_ctau, psi2S_mu1, psi2S_sigma2);
    RooRealVar psi2S_prop1("psi2S_prop1", "psi2S_prop1", 0.726501);
    RooAddPdf psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", RooArgList(psi2S_gauss1, psi2S_gauss2), psi2S_prop1);
    // Background p.d.f.
    RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.001, 0, 0.02);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.06, 0.01, 0.1);
    RooGExpModel psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    // RooRealVar psi2S_mu2("psi2S_mu2", "psi2S_mu2", 0, -0.005, 0.005);
    // RooRealVar psi2S_sigma4("psi2S_sigma4", "psi2S_sigma4", 0.003, 0.001, 0.02);
    // RooRealVar psi2S_sigma5("psi2S_sigma5", "psi2S_sigma5", 0.01, 0.001, 0.03);
    // RooRealVar psi2S_coef2("psi2S_coef2", "psi2S_coef2", 0.06, 0.01, 0.1);
    // RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.5, 0, 1);
    // RooGaussian psi2S_gauss3("psi2S_gauss3", "psi2S_gauss3", psi2S_ctau, psi2S_mu2, psi2S_sigma4);
    // RooGExpModel psi2S_expGau2("psi2S_expGau2", "psi2S_expGau2", psi2S_ctau, psi2S_sigma5, psi2S_coef2, false, RooGExpModel::Type::Flipped);
    // RooAddPdf psi2SCtauCombSig("psi2SCtauCombSig", "psi2SCtauCombSig", RooArgList(psi2S_gauss3, psi2S_expGau2), psi2S_prop2);
    // Combinatorial background p.d.f.
    RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 0, -0.001, 0.001);
    RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.5, 0.1, 10);
    RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 3, 0.5, 20);
    RooRealVar psi2S_sigma6("psi2S_sigma6", "psi2S_sigma6", 0.001, 0, 0.007);
    RooRealVar psi2S_sigma7("psi2S_sigma7", "psi2S_sigma7", 0.01, 0.007, 0.02);
    RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 0.06, 0.01, 0.1);
    // RooRealVar psi2S_coef4("psi2S_coef4", "psi2S_coef4", 0.02, 0.01, 0.05);
    RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.5, 0, 1);
    // RooGaussian psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_mu3, psi2S_sigma6);
    RooCBShape psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_mu3, psi2S_sigma6, psi2S_alpha3, psi2S_nx3);
    // RooGExpModel psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_sigma6, psi2S_coef4, false, RooGExpModel::Type::Flipped);
    RooGExpModel psi2S_expGau3("psi2S_expGau3", "psi2S_expGau3", psi2S_ctau, psi2S_sigma7, psi2S_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCombBkg("psi2SCtauCombBkg", "psi2SCtauCombBkg", RooArgList(psi2S_gauss4, psi2S_expGau3), psi2S_prop3);
    
    #ifdef PSI2S_CTAU_SELF_SIDE
    // Fit psi(2S) Ctau in psi(2S) sideband
    selection = "psi2S_mass > 3.85 || psi2S_mass < 3.55";
    RooRealVar npsi2SCtauCombBkg("npsi2SCtauCombBkg", "npsi2SCtauCombBkg", 1, 0, N);
    RooExtendPdf pdf_psi2SCtauCombBkg("pdf_psi2SCtauCombBkg", "pdf_psi2SCtauCombBkg", psi2SCtauCombBkg, npsi2SCtauCombBkg);
    while(true) {
        res = pdf_psi2SCtauCombBkg.fitTo(*data->reduce(psi2S_ctau, selection.c_str()), Save());
        if(!res->status()) break;
    }
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 10, 10, 2400, 1800);
    RooPlot *frame5 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau Self-sideband"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame5, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SCtauCombBkg.plotOn(frame5, LineColor(1), LineWidth(2), Name("All"));
    pdf_psi2SCtauCombBkg.plotOn(frame5, Components(psi2S_gauss4), LineColor(kBlue), LineWidth(2), Name("GS"));
    pdf_psi2SCtauCombBkg.plotOn(frame5, Components(psi2S_expGau3), LineColor(kMagenta), LineWidth(2), Name("EX"));
    TLegend *legend5 = new TLegend(.65, .60, .85, .85);
    legend5->AddEntry(frame5->findObject("Data"), "RunII 2018", "L");
    legend5->AddEntry(frame5->findObject("All"), "Total p.d.f.", "L");
    legend5->AddEntry(frame5->findObject("GS"), "Gaussian", "L");
    legend5->AddEntry(frame5->findObject("EX"), "Exp#otimes Gauss", "L");
    frame5->SetAxisRange(1, 1e3, "Y");
    gPad->SetLogy();
    frame5->Draw();
    legend5->DrawClone();
    canvas5->SaveAs("fits/Side0_psi2SCtau.pdf");
    pdf_psi2SCtauCombBkg.getVariables()->Print("v");
    // Fix parameters as constant
    psi2S_mu3.setConstant();
    psi2S_sigma6.setConstant();
    psi2S_sigma7.setConstant();
    psi2S_coef3.setConstant();
    psi2S_prop3.setConstant();
    #endif
    #ifdef PSI2S_CTAU_ELSE_SIDE
    // Fit psi(2S) Ctau in J/psi sideband
    selection = "Jpsi_mass > 3.19 || Jpsi_mass < 2.99";
    RooRealVar npsi2SCtauCombSig("npsi2SCtauCombSig", "npsi2SCtauCombSig", 1, 0, N);
    RooExtendPdf pdf_psi2SCtauCombSig("pdf_psi2SCtauCombSig", "pdf_psi2SCtauCombSig", psi2SCtauCombSig, npsi2SCtauCombSig);
    RooAddPdf pdf_psi2SCtauComb("pdf_psi2SCtauComb", "pdf_psi2SCtauComb", RooArgList(pdf_psi2SCtauCombSig, pdf_psi2SCtauCombBkg));
    pdf_psi2SCtauComb.fitTo(*data->reduce(psi2S_ctau, selection.c_str()), Save());
    // Draw data point and p.d.f. curve
    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 10, 10, 800, 600);
    RooPlot *frame6 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau Else-sideband"), Bins(BinNum));
    data->reduce(selection.c_str())->plotOn(frame6, DataError(RooAbsData::SumW2));
    pdf_psi2SCtauComb.plotOn(frame6);
    pdf_psi2SCtauComb.plotOn(frame6, Components(psi2S_gauss3), LineColor(1));
    pdf_psi2SCtauComb.plotOn(frame6, Components(psi2S_expGau2), LineColor(3));
    pdf_psi2SCtauComb.plotOn(frame6, Components(psi2SCtauCombBkg), LineColor(5));
    frame6->SetAxisRange(1, 1e3, "Y");
    gPad->SetLogy();
    frame6->Draw();
    canvas6->SaveAs("fits/Side1_psi2SCtau.pdf");
    pdf_psi2SCtauComb.getVariables()->Print("v");
    // Fix parameters as constant
    // psi2S_mu2.setConstant();
    // psi2S_sigma4.setConstant();
    // psi2S_sigma5.setConstant();
    // psi2S_coef2.setConstant();
    // psi2S_prop2.setConstant();
    #endif
    // Form 4-dim p.d.f.
    #ifdef ALL_4D
    // J/psi Mass + psi(2S) Mass dimension
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    // RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    // RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    // J/psi Mass + psi(2S) Ctau dimension
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig, psi2SCtauSig);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig, psi2SCtauBkg);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg, psi2SCtauSig);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg, psi2SCtauBkg);
    RooProdPdf pdf_ctau_SigBkg("pdf_ctau_SigBkg", "pdf_ctau_SigBkg", JpsiCtauCombSig, psi2SCtauCombBkg);
    // RooProdPdf pdf_ctau_BkgSig("pdf_ctau_BkgSig", "pdf_ctau_BkgSig", JpsiCtauCombBkg, psi2SCtauCombSig);
    // RooProdPdf pdf_ctau_BkgBkg("pdf_ctau_BkgBkg", "pdf_ctau_BkgBkg", JpsiCtauCombBkg, psi2SCtauCombBkg);
    // Combine Mass and Ctau dimensions
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_ctau_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_ctau_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_ctau_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_ctau_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_ctau_SigBkg);
    // RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_ctau_BkgSig);
    // RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_ctau_BkgBkg);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e4, 1e2, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1e2, N);
    RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 1e2, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e2, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e3, 1, N);
    // RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e3, 1, N);
    // RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, 1e4);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb),
        RooArgList(n_P_P, n_P_NP, n_NP_P, n_NP_NP, n_Sig_Comb)
    );
    pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau)), Save());
    // Draw data point and p.d.f. curve
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 3000, 2000);
    RooPlot *frame7 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 4D"), Bins(BinNum));
    data->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame7, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame7, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame7, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame7, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame7, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame7, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    // pdf_all.plotOn(frame7, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    // pdf_all.plotOn(frame7, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    legend7->AddEntry(frame7->findObject("All"), "Total p.d.f.", "L");
    legend7->AddEntry(frame7->findObject("P_P"), "prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("P_NP"), "prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_P"), "non-prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    // legend7->AddEntry(frame7->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    // legend7->AddEntry(frame7->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame7->Draw();
    legend7->DrawClone();
    canvas7->SaveAs("fits/4D_JpsiMass.pdf");
    // canvas7->SaveAs("fits/4D_JpsiMass.png");
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 3000, 2000);
    RooPlot *frame8 = psi2S_mass.frame(RooFit::Title("psi2S Mass 4D"), Bins(BinNum));
    data->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame8, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame8, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame8, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame8, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame8, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame8, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    // pdf_all.plotOn(frame8, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    // pdf_all.plotOn(frame8, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    legend8->AddEntry(frame8->findObject("All"), "Total p.d.f.", "L");
    legend8->AddEntry(frame8->findObject("P_P"), "prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("P_NP"), "prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_P"), "non-prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    // legend8->AddEntry(frame8->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    // legend8->AddEntry(frame8->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame8->Draw();
    legend8->DrawClone();
    canvas8->SaveAs("fits/4D_psi2SMass.pdf");
    // canvas8->SaveAs("fits/4D_psi2SMass.png");
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 3000, 2000);
    RooPlot *frame9 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau 4D"), Bins(BinNum));
    data->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame9, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame9, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame9, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame9, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame9, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame9, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    // pdf_all.plotOn(frame9, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    // pdf_all.plotOn(frame9, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend9 = new TLegend(.65, .60, .85, .85);
    legend9->AddEntry(frame9->findObject("Data"), "RunII 2018", "L");
    legend9->AddEntry(frame9->findObject("All"), "Total p.d.f.", "L");
    legend9->AddEntry(frame9->findObject("P_P"), "prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("P_NP"), "prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_P"), "non-prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    // legend9->AddEntry(frame9->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    // legend9->AddEntry(frame9->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame9->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame9->Draw();
    legend9->DrawClone();
    canvas9->SaveAs("fits/4D_JpsiCtau.pdf");
    // canvas9->SaveAs("fits/4D_JpsiCtau.png");
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 3000, 2000);
    RooPlot *frame0 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau 4D"), Bins(BinNum));
    data->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_all.plotOn(frame0, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_all.plotOn(frame0, Components(pdf_P_P), LineColor(kBlue), LineStyle(kSolid), LineWidth(2), Name("P_P"));
    pdf_all.plotOn(frame0, Components(pdf_P_NP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_all.plotOn(frame0, Components(pdf_NP_P), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_all.plotOn(frame0, Components(pdf_NP_NP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    pdf_all.plotOn(frame0, Components(pdf_Sig_Comb), LineColor(kRed), LineWidth(1), Name("Sig_Comb"));
    // pdf_all.plotOn(frame0, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    // pdf_all.plotOn(frame0, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend0 = new TLegend(.65, .60, .85, .85);
    legend0->AddEntry(frame0->findObject("Data"), "RunII 2018", "L");
    legend0->AddEntry(frame0->findObject("All"), "Total p.d.f.", "L");
    legend0->AddEntry(frame0->findObject("P_P"), "prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("P_NP"), "prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_P"), "non-prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    // legend0->AddEntry(frame0->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    // legend0->AddEntry(frame0->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame0->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame0->Draw();
    legend0->DrawClone();
    canvas0->SaveAs("fits/4D_psi2SCtau.pdf");
    // canvas0->SaveAs("fits/4D_psi2SCtau.png");
    pdf_all.getVariables()->Print("v");
    // // Integral check
    // RooAbsReal* intpsi2SCtau_1 = pdf_all.createIntegral(RooArgList(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau));
    // cout<<"integral over cut region: "<<1 - intpsi2SCtau_1->getVal()<<endl;
    // Jpsi_mass.setRange(0, 5);
    // psi2S_mass.setRange(0, 5);
    // Jpsi_ctau.setRange(-1, 1);
    // psi2S_ctau.setRange(-1, 1);
    // RooAbsReal* intpsi2SCtau_2 = pdf_all.createIntegral(RooArgList(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau));
    // cout<<"integral outside cut region: "<<1 - intpsi2SCtau_2->getVal()<<endl;
    #endif
    return;
}