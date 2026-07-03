#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
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
#include "RooWorkspace.h"

using namespace std;
using namespace RooFit;

void Fit_1D(string var, double vmin, double vmax) {
    // Define variables
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
	RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.03, 0.25);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.03, 0.25);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
    RooRealVar kine(var.c_str(), var.c_str(), vmin, vmax);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_ctau);
    variables.add(psi2S_ctau);
    variables.add(evt_weight);
    variables.add(kine);
    // Read file and apply cut
    string sel = var + " > " + to_string(vmin) + " && " + var + " < " + to_string(vmax);
    TFile *BdecayFile = new TFile("WeightBdecaySmall.root", "READ");
    TFile *MixFile = new TFile("WeightMix.root", "READ");
    TTree *BdecayTree = (TTree*)BdecayFile->Get("data");
    TTree *MixTree = (TTree*)MixFile->Get("data");
    RooDataSet *Bdecay = new RooDataSet("Bdecay", "Bdecay", BdecayTree, variables, sel.c_str(), "evt_weight");
    RooDataSet *Mix = new RooDataSet("Mix", "Mix", MixTree, variables, sel.c_str(), "evt_weight");
    // Read parameters from differential 4D fit
    TFile *f = new TFile("Model_4D_diff.root");
    RooWorkspace *wsp = (RooWorkspace *)f->Get("wsp");
    
    // Define p.d.f.
    //// J/psi Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar Jpsi_m1("Jpsi_m1", "Jpsi_m1", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_d1("Jpsi_d1", "Jpsi_d1", 0.01, 0, 0.03);
    RooRealVar Jpsi_a1("Jpsi_a1", "Jpsi_a1", 1, 0.1, 5);
    RooRealVar Jpsi_n1("Jpsi_n1", "Jpsi_n1", 1, 0, 15);
    RooRealVar Jpsi_d2("Jpsi_d2", "Jpsi_d2", 0.05, 0.03, 0.08);
    RooRealVar Jpsi_r1("Jpsi_r1", "Jpsi_r1", 0.6, 0, 1);
    RooCBShape Jpsi_c1("Jpsi_c1", "Jpsi_c1", Jpsi_mass, Jpsi_m1, Jpsi_d1, Jpsi_a1, Jpsi_n1);
    RooGaussian Jpsi_g1("Jpsi_g1", "Jpsi_g1", Jpsi_mass, Jpsi_m1, Jpsi_d2);
    RooAddPdf JpsiMassMix("JpsiMassMix", "JpsiMassMix", RooArgList(Jpsi_c1, Jpsi_g1), Jpsi_r1);
    ////// Bdecay
    RooRealVar Jpsi_m2("Jpsi_m2", "Jpsi_m2", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_d3("Jpsi_d3", "Jpsi_d3", 0.01, 0, 0.03);
    RooRealVar Jpsi_a2("Jpsi_a2", "Jpsi_a2", 1, 0.1, 5);
    RooRealVar Jpsi_n2("Jpsi_n2", "Jpsi_n2", 1, 0, 15);
    RooRealVar Jpsi_d4("Jpsi_d4", "Jpsi_d4", 0.05, 0.03, 0.08);
    RooRealVar Jpsi_r2("Jpsi_r2", "Jpsi_r2", 0.6, 0, 1);
    RooCBShape Jpsi_c2("Jpsi_c2", "Jpsi_c2", Jpsi_mass, Jpsi_m2, Jpsi_d3, Jpsi_a2, Jpsi_n2);
    RooGaussian Jpsi_g2("Jpsi_g2", "Jpsi_g2", Jpsi_mass, Jpsi_m2, Jpsi_d4);
    RooAddPdf JpsiMassB("JpsiMassB", "JpsiMassB", RooArgList(Jpsi_c2, Jpsi_g2), Jpsi_r2);
    ////// Combinatorial signal
    RooRealVar Jpsi_r3("Jpsi_r3", "Jpsi_r3", wsp->var("Jpsi_prop2")->getVal());
    RooAddPdf JpsiMassCSig("JpsiMassCSig", "JpsiMassCSig", RooArgList(JpsiMassMix, JpsiMassB), Jpsi_r3);
    ////// Combinatorial background
    RooChebychev JpsiMassCBkg("JpsiMassCBkg", "JpsiMassCBkg", Jpsi_mass, RooArgList());

    //// psi(2S) Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar psi2S_m1("psi2S_m1", "psi2S_m1", 3.686, 3.63, 3.73);
    RooRealVar psi2S_d1("psi2S_d1", "psi2S_d1", 0.02, 0.01, 0.05);
    RooRealVar psi2S_a1("psi2S_a1", "psi2S_a1", 1.5, 0.1, 1);
    RooRealVar psi2S_n1("psi2S_n1", "psi2S_n1", 3, 0.1, 15);
    RooRealVar psi2S_d2("psi2S_d2", "psi2S_d2", 0.07, 0.04, 0.2);
    RooRealVar psi2S_r1("psi2S_r1", "psi2S_r1", 0.5, 0, 1);
    RooCBShape psi2S_c1("psi2S_c1", "psi2S_c1", psi2S_mass, psi2S_m1, psi2S_d1, psi2S_a1, psi2S_n1);
    RooGaussian psi2S_g1("psi2S_g1", "psi2S_g1", psi2S_mass, psi2S_m1, psi2S_d2);
    RooAddPdf psi2SMassMix("psi2SMassMix", "psi2SMassMix", RooArgList(psi2S_c1, psi2S_g1), psi2S_r1);
    ////// Bdecay
    RooRealVar psi2S_m2("psi2S_m2", "psi2S_m2", 3.686, 3.63, 3.73);
    RooRealVar psi2S_d3("psi2S_d3", "psi2S_d3", 0.02, 0.01, 0.05);
    RooRealVar psi2S_a2("psi2S_a2", "psi2S_a2", 1.5, 0.1, 1);
    RooRealVar psi2S_n2("psi2S_n2", "psi2S_n2", 3, 0.1, 15);
    RooRealVar psi2S_d4("psi2S_d4", "psi2S_d4", 0.07, 0.04, 0.2);
    RooRealVar psi2S_r2("psi2S_r2", "psi2S_r2", 0.5, 0, 1);
    RooCBShape psi2S_c2("psi2S_c2", "psi2S_c2", psi2S_mass, psi2S_m2, psi2S_d3, psi2S_a2, psi2S_n2);
    RooGaussian psi2S_g2("psi2S_g2", "psi2S_g2", psi2S_mass, psi2S_m2, psi2S_d4);
    RooAddPdf psi2SMassB("psi2SMassB", "psi2SMassB", RooArgList(psi2S_c2, psi2S_g2), psi2S_r2);
    ////// Combinatorial signal
    RooRealVar psi2S_r3("psi2S_r3", "psi2S_r3", wsp->var("psi2S_prop2")->getVal());
    RooAddPdf psi2SMassCSig("psi2SMassCSig", "psi2SMassCSig", RooArgList(psi2SMassMix, psi2SMassB), psi2S_r3);
    ////// Combinatorial background
    RooRealVar psi2S_a3("psi2S_a3", "psi2S_a3", wsp->var("psi2S_a")->getVal());
    RooRealVar psi2S_b1("psi2S_b1", "psi2S_b1", wsp->var("psi2S_b")->getVal());
    RooChebychev psi2SMassCBkg("psi2SMassCBkg", "psi2SMassCBkg", psi2S_mass, RooArgList(psi2S_a3, psi2S_b1));

    //// J/psi Ctau
    ////// Mixture of SPS/DPS=1:2
    RooRealVar Jpsi_u1("Jpsi_u1", "Jpsi_u1", 0, -0.005, 0.005);
    RooRealVar Jpsi_s1("Jpsi_s1", "Jpsi_s1", 0.001, 0, 0.003);
    RooRealVar Jpsi_s2("Jpsi_s2", "Jpsi_s2", 0.003, 0.003, 0.008);
    RooGaussian Jpsi_g3("Jpsi_g3", "Jpsi_g3", Jpsi_ctau, Jpsi_u1, Jpsi_s1);
    RooGaussian Jpsi_g4("Jpsi_g4", "Jpsi_g4", Jpsi_ctau, Jpsi_u1, Jpsi_s2);
    RooRealVar Jpsi_p1("Jpsi_p1", "Jpsi_p1", 0.5, 0, 1);
    RooAddPdf JpsiCtauMix("JpsiCtauMix", "JpsiCtauMix", RooArgList(Jpsi_g3, Jpsi_g4), Jpsi_p1);
    ////// Bdecay
    RooRealVar Jpsi_s3("Jpsi_s3", "Jpsi_s3", 0.001, 0, 0.003);
    RooRealVar Jpsi_s4("Jpsi_s4", "Jpsi_s4", 0.001, 0.002, 0.008);
    RooRealVar Jpsi_e1("Jpsi_e1", "Jpsi_e1", 0.03, 0.01, 0.1);
    RooRealVar Jpsi_e2("Jpsi_e2", "Jpsi_e2", 0.03, 0.01, 0.1);
    RooRealVar Jpsi_p2("Jpsi_p2", "Jpsi_p2", 0.5, 0, 1);
    RooGExpModel Jpsi_ge1("Jpsi_ge1", "Jpsi_ge1", Jpsi_ctau, Jpsi_s3, Jpsi_e1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_ge2("Jpsi_ge2", "Jpsi_ge2", Jpsi_ctau, Jpsi_s4, Jpsi_e2, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauB("JpsiCtauB", "JpsiCtauB", RooArgList(Jpsi_ge1, Jpsi_ge2), Jpsi_p2);
    ////// Combinatorial signal
    RooRealVar Jpsi_p3("Jpsi_p3", "Jpsi_p3", wsp->var("Jpsi_prop2")->getVal());
    RooAddPdf JpsiCtauCSig("JpsiCtauCSig", "JpsiCtauCSig", RooArgList(JpsiCtauMix, JpsiCtauB), Jpsi_p3);
    ////// Combinatorial background
    RooRealVar Jpsi_u2("Jpsi_u2", "Jpsi_u2", wsp->var("Jpsi_mu3")->getVal());
    RooRealVar Jpsi_s5("Jpsi_s5", "Jpsi_s5", wsp->var("Jpsi_sigma6")->getVal());
    RooRealVar Jpsi_e3("Jpsi_e3", "Jpsi_e3", wsp->var("Jpsi_coef3")->getVal());
    RooRealVar Jpsi_p4("Jpsi_p4", "Jpsi_p4", wsp->var("Jpsi_prop3")->getVal());
    RooGaussian Jpsi_g5("Jpsi_g5", "Jpsi_g5", Jpsi_ctau, Jpsi_u2, Jpsi_s5);
    RooGExpModel Jpsi_ge3("Jpsi_ge3", "Jpsi_ge3", Jpsi_ctau, Jpsi_s5, Jpsi_e3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCBkg("JpsiCtauCBkg", "JpsiCtauCBkg", RooArgList(Jpsi_g5, Jpsi_ge3), Jpsi_p4);

    //// psi(2S) Ctau
    ////// Signal
    RooRealVar psi2S_u1("psi2S_u1", "psi2S_u1", 0, -1e-4, 1e-4);
    RooRealVar psi2S_s1("psi2S_s1", "psi2S_s1", 0.001, 0, 0.003);
    RooRealVar psi2S_s2("psi2S_s2", "psi2S_s2", 0.005, 0.003, 0.02);
    RooRealVar psi2S_l2("psi2S_l2", "psi2S_l2", -1, -5, 0, "");
    RooRealVar psi2S_x2("psi2S_x2", "psi2S_x2", 1, 0, 10);
    RooRealVar psi2S_p1("psi2S_p1", "psi2S_p1", 0.5, 0, 1);
    RooGaussian psi2S_g3("psi2S_g3", "psi2S_g3", psi2S_ctau, psi2S_u1, psi2S_s1);
    RooCBShape psi2S_c3("psi2S_c3", "psi2S_c3", psi2S_ctau, psi2S_u1, psi2S_s2, psi2S_l2, psi2S_x2);
    RooAddPdf psi2SCtauMix("psi2SCtauMix", "psi2SCtauMix", RooArgList(psi2S_g3, psi2S_c3), psi2S_p1);
    ////// Bdecay
    RooRealVar psi2S_s3("psi2S_s3", "psi2S_s3", 0.001, 0, 0.003);
    RooRealVar psi2S_s4("psi2S_s4", "psi2S_s4", 0.001, 0.002, 0.008);
    RooRealVar psi2S_e1("psi2S_e1", "psi2S_e1", 0.03, 0.01, 0.1);
    RooRealVar psi2S_e2("psi2S_e2", "psi2S_e2", 0.03, 0.01, 0.1);
    RooRealVar psi2S_p2("psi2S_p2", "psi2S_p2", 0.5, 0, 1);
    RooGExpModel psi2S_ge1("psi2S_ge1", "psi2S_ge1", psi2S_ctau, psi2S_s3, psi2S_e1, false, RooGExpModel::Type::Flipped);
    RooGExpModel psi2S_ge2("psi2S_ge2", "psi2S_ge2", psi2S_ctau, psi2S_s4, psi2S_e2, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauB("psi2SCtauB", "psi2SCtauB", RooArgList(psi2S_ge1, psi2S_ge2), psi2S_p2);
    ////// Combinatorial signal
    RooRealVar psi2S_p3("psi2S_p3", "psi2S_p3", wsp->var("Jpsi_prop2")->getVal());
    RooAddPdf psi2SCtauCSig("psi2SCtauCSig", "psi2SCtauCSig", RooArgList(psi2SCtauMix, psi2SCtauB), psi2S_p3);
    ////// Combinatorial background
    RooRealVar psi2S_u2("psi2S_u2", "psi2S_u2", wsp->var("psi2S_mu3")->getVal());
    RooRealVar psi2S_l3("psi2S_l3", "psi2S_l3", wsp->var("psi2S_alpha3")->getVal());
    RooRealVar psi2S_x3("psi2S_x3", "psi2S_x3", wsp->var("psi2S_nx3")->getVal());
    RooRealVar psi2S_s5("psi2S_s5", "psi2S_s5", wsp->var("psi2S_sigma6")->getVal());
    RooRealVar psi2S_s6("psi2S_s6", "psi2S_s6", wsp->var("psi2S_sigma7")->getVal());
    RooRealVar psi2S_e3("psi2S_e3", "psi2S_e3", wsp->var("psi2S_coef3")->getVal());
    RooRealVar psi2S_p4("psi2S_p4", "psi2S_p4", wsp->var("psi2S_prop3")->getVal());
    RooCBShape psi2S_g4("psi2S_g4", "psi2S_g4", psi2S_ctau, psi2S_u2, psi2S_s5, psi2S_l3, psi2S_x3);
    RooGExpModel psi2S_ge3("psi2S_ge3", "psi2S_ge3", psi2S_ctau, psi2S_s6, psi2S_e3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCBkg("psi2SCtauCBkg", "psi2SCtauCBkg", RooArgList(psi2S_g4, psi2S_ge3), psi2S_p4);

    // 1D fit
    Int_t N = 1500000, BinNum = 100;
    //// J/psi Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar nJpsiMassMix("nJpsiMassMix", "nJpsiMassMix", 1, 0, N);
    RooExtendPdf pdf_JpsiMassMix("pdf_JpsiMassMix", "pdf_JpsiMassMix", JpsiMassMix, nJpsiMassMix);
    RooFitResult *res1;
    while(true) {
        res1 = pdf_JpsiMassMix.fitTo(*Mix->reduce(Jpsi_mass), Save());
        if(!res1->status()) break;
    }
    TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 10, 10, 2400, 1800);
    RooPlot *frame1 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 1D Mix"), Bins(BinNum));
    Mix->plotOn(frame1, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiMassMix.plotOn(frame1, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiMassMix.plotOn(frame1, Components(Jpsi_c1), LineColor(kBlue), LineWidth(2), Name("CB"));
    pdf_JpsiMassMix.plotOn(frame1, Components(Jpsi_g1), LineColor(kRed), LineWidth(2), Name("GS"));
    TLegend *legend1 = new TLegend(.65, .60, .85, .85);
    legend1->AddEntry(frame1->findObject("Data"), "RunII 2018", "L");
    legend1->AddEntry(frame1->findObject("All"), "Total p.d.f.", "L");
    legend1->AddEntry(frame1->findObject("CB"), "Crystall Ball", "L");
    legend1->AddEntry(frame1->findObject("GS"), "Gaussian", "L");
    frame1->Draw();
    legend1->DrawClone();
    canvas1->SaveAs("fig/1D_JpsiMassMix.png");
    pdf_JpsiMassMix.getVariables()->Print("v");
    ////// Bdecay
    RooRealVar nJpsiMassB("nJpsiMassB", "nJpsiMassB", 1, 0, N);
    RooExtendPdf pdf_JpsiMassB("pdf_JpsiMassB", "pdf_JpsiMassB", JpsiMassB, nJpsiMassB);
    RooFitResult *res2;
    while(true) {
        res2 = pdf_JpsiMassB.fitTo(*Bdecay->reduce(Jpsi_mass), Save());
        if(!res2->status()) break;
    }
    TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 10, 10, 2400, 1800);
    RooPlot *frame2 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 1D Bdecay"), Bins(BinNum));
    Bdecay->plotOn(frame2, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiMassB.plotOn(frame2, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiMassB.plotOn(frame2, Components(Jpsi_c2), LineColor(kBlue), LineWidth(2), Name("CB"));
    pdf_JpsiMassB.plotOn(frame2, Components(Jpsi_g2), LineColor(kRed), LineWidth(2), Name("GS"));
    TLegend *legend2 = new TLegend(.65, .60, .85, .85);
    legend2->AddEntry(frame2->findObject("Data"), "RunII 2018", "L");
    legend2->AddEntry(frame2->findObject("All"), "Total p.d.f.", "L");
    legend2->AddEntry(frame2->findObject("CB"), "Crystall Ball", "L");
    legend2->AddEntry(frame2->findObject("GS"), "Gaussian", "L");
    frame2->Draw();
    legend2->DrawClone();
    canvas2->SaveAs("fig/1D_JpsiMassBdecay.png");
    pdf_JpsiMassB.getVariables()->Print("v");

    //// psi(2S) Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar npsi2SMassMix("npsi2SMassMix", "npsi2SMassMix", 1, 0, N);
    RooExtendPdf pdf_psi2SMassMix("pdf_psi2SMassMix", "pdf_psi2SMassMix", psi2SMassMix, npsi2SMassMix);
    RooFitResult *res3;
    while(true) {
        res3 = pdf_psi2SMassMix.fitTo(*Mix->reduce(psi2S_mass, "psi2S_mass > 3.55 && psi2S_mass < 3.85"), Save());
        if(!res3->status()) break;
    }
    TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 10, 10, 2400, 1800);
    RooPlot *frame3 = new RooPlot(psi2S_mass, 3.55, 3.85, BinNum);//psi2S_mass.frame(RooFit::Title("psi2S Mass 1D Mix"), Bins(BinNum));
    frame3->SetTitle("psi2S Mass 1D Mix");
    Mix->plotOn(frame3, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SMassMix.plotOn(frame3, LineColor(1), LineWidth(2), Name("All"));
    pdf_psi2SMassMix.plotOn(frame3, Components(psi2S_c1), LineColor(kBlue), LineWidth(2), Name("CB"));
    pdf_psi2SMassMix.plotOn(frame3, Components(psi2S_g1), LineColor(kRed), LineWidth(2), Name("GS"));
    TLegend *legend3 = new TLegend(.65, .60, .85, .85);
    legend3->AddEntry(frame3->findObject("Data"), "RunII 2018", "L");
    legend3->AddEntry(frame3->findObject("All"), "Total p.d.f.", "L");
    legend3->AddEntry(frame3->findObject("CB"), "Crystall Ball", "L");
    legend3->AddEntry(frame3->findObject("GS"), "Gaussian", "L");
    frame3->Draw();
    legend3->DrawClone();
    canvas3->SaveAs("fig/1D_psi2SMassMix.png");
    pdf_psi2SMassMix.getVariables()->Print("v");
    ////// Bdecay
    RooRealVar npsi2SMassB("npsi2SMassB", "npsi2SMassB", 1, 0, N);
    RooExtendPdf pdf_psi2SMassB("pdf_psi2SMassB", "pdf_psi2SMassB", psi2SMassB, npsi2SMassB);
    RooFitResult *res4;
    while(true) {
        res4 = pdf_psi2SMassB.fitTo(*Bdecay->reduce(psi2S_mass, "psi2S_mass > 3.55 && psi2S_mass < 3.85"), Save());
        if(!res4->status()) break;
    }
    TCanvas *canvas4 = new TCanvas("canvas4", "canvas4", 10, 10, 2400, 1800);
    RooPlot *frame4 = new RooPlot(psi2S_mass, 3.55, 3.85, BinNum);
    frame4->SetTitle("psi2S Mass 1D Bdecay");
    Bdecay->plotOn(frame4, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SMassB.plotOn(frame4, LineColor(1), LineWidth(2), Name("All"));
    pdf_psi2SMassB.plotOn(frame4, Components(psi2S_c2), LineColor(kBlue), LineWidth(2), Name("CB"));
    pdf_psi2SMassB.plotOn(frame4, Components(psi2S_g2), LineColor(kRed), LineWidth(2), Name("GS"));
    TLegend *legend4 = new TLegend(.65, .60, .85, .85);
    legend4->AddEntry(frame4->findObject("Data"), "RunII 2018", "L");
    legend4->AddEntry(frame4->findObject("All"), "Total p.d.f.", "L");
    legend4->AddEntry(frame4->findObject("CB"), "Crystall Ball", "L");
    legend4->AddEntry(frame4->findObject("GS"), "Gaussian", "L");
    frame4->Draw();
    legend4->DrawClone();
    canvas4->SaveAs("fig/1D_psi2SMassBdecay.png");
    pdf_psi2SMassB.getVariables()->Print("v");

    //// J/psi Ctau
    ////// Mixture of SPS/DPS=1:2
    RooRealVar nJpsiCtauMix("nJpsiCtauMix", "nJpsiCtauMix", 1, 0, N);
    RooExtendPdf pdf_JpsiCtauMix("pdf_JpsiCtauMix", "pdf_JpsiCtauMix", JpsiCtauMix, nJpsiCtauMix);
    RooFitResult *res5;
    while(true) {
        res5 = pdf_JpsiCtauMix.fitTo(*Mix->reduce(Jpsi_ctau, "Jpsi_ctau > -0.02 && Jpsi_ctau < 0.02"), Save());
        if(!res5->status()) break;
    }
    TCanvas *canvas5 = new TCanvas("canvas5", "canvas5", 10, 10, 2400, 1800);
    RooPlot *frame5 = new RooPlot(Jpsi_ctau, -0.02, 0.02, BinNum);
    frame5->SetTitle("Jpsi Ctau 1D Mix");
    Mix->plotOn(frame5, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiCtauMix.plotOn(frame5, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiCtauMix.plotOn(frame5, Components(Jpsi_g3), LineColor(kBlue), LineWidth(2), Name("GS1"));
    pdf_JpsiCtauMix.plotOn(frame5, Components(Jpsi_g4), LineColor(kRed), LineWidth(2), Name("GS2"));
    TLegend *legend5 = new TLegend(.65, .60, .85, .85);
    legend5->AddEntry(frame5->findObject("Data"), "RunII 2018", "L");
    legend5->AddEntry(frame5->findObject("All"), "Total p.d.f.", "L");
    legend5->AddEntry(frame5->findObject("GS1"), "Gaussian 1", "L");
    legend5->AddEntry(frame5->findObject("GS2"), "Gaussian 2", "L");
    frame5->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame5->Draw();
    legend5->DrawClone();
    canvas5->SaveAs("fig/1D_JpsiCtauMix.png");
    pdf_JpsiCtauMix.getVariables()->Print("v");
    ////// Bdecay
    RooRealVar nJpsiCtauB("nJpsiCtauB", "nJpsiCtauB", 1, 0, N);
    RooExtendPdf pdf_JpsiCtauB("pdf_JpsiCtauB", "pdf_JpsiCtauB", JpsiCtauB, nJpsiCtauB);
    RooFitResult *res6;
    while(true) {
        res6 = pdf_JpsiCtauB.fitTo(*Bdecay->reduce(Jpsi_ctau, "Jpsi_ctau > -0.03 && Jpsi_ctau < 0.25"), Save());
        if(!res6->status()) break;
    }
    TCanvas *canvas6 = new TCanvas("canvas6", "canvas6", 10, 10, 2400, 1800);
    RooPlot *frame6 = new RooPlot(Jpsi_ctau, -0.03, 0.25, BinNum);
    frame6->SetTitle("Jpsi Ctau 1D Bdecay");
    Bdecay->plotOn(frame6, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_JpsiCtauB.plotOn(frame6, LineColor(1), LineWidth(2), Name("All"));
    pdf_JpsiCtauB.plotOn(frame6, Components(Jpsi_ge1), LineColor(kBlue), LineWidth(2), Name("GExp1"));
    pdf_JpsiCtauB.plotOn(frame6, Components(Jpsi_ge2), LineColor(kRed), LineWidth(2), Name("GExp2"));
    TLegend *legend6 = new TLegend(.65, .60, .85, .85);
    legend6->AddEntry(frame6->findObject("Data"), "RunII 2018", "L");
    legend6->AddEntry(frame6->findObject("All"), "Total p.d.f.", "L");
    legend6->AddEntry(frame6->findObject("GExp1"), "Gaussian Exp 1", "L");
    legend6->AddEntry(frame6->findObject("GExp2"), "Gaussian Exp 2", "L");
    frame6->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame6->Draw();
    legend6->DrawClone();
    canvas6->SaveAs("fig/1D_JpsiCtauBdecay.png");
    pdf_JpsiCtauB.getVariables()->Print("v");

    //// psi(2S) Ctau
    ////// Mixture of SPS/DPS=1:2
    RooRealVar npsi2SCtauMix("npsi2SCtauMix", "npsi2SCtauMix", 1, 0, N);
    RooExtendPdf pdf_psi2SCtauMix("pdf_psi2SCtauMix", "pdf_psi2SCtauMix", psi2SCtauMix, npsi2SCtauMix);
    RooFitResult *res7;
    while(true) {
        res7 = pdf_psi2SCtauMix.fitTo(*Mix->reduce(psi2S_ctau, "psi2S_ctau > -0.02 && psi2S_ctau < 0.08"), Save());
        if(!res7->status()) break;
    }
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 10, 10, 2400, 1800);
    RooPlot *frame7 = new RooPlot(psi2S_ctau, -0.02, 0.02, BinNum);
    frame7->SetTitle("psi2S Ctau 1D Mix");
    Mix->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SCtauMix.plotOn(frame7, LineColor(1), LineWidth(2), Name("All"));
    pdf_psi2SCtauMix.plotOn(frame7, Components(psi2S_g3), LineColor(kBlue), LineWidth(2), Name("GS"));
    pdf_psi2SCtauMix.plotOn(frame7, Components(psi2S_c3), LineColor(kRed), LineWidth(2), Name("CB"));
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    legend7->AddEntry(frame7->findObject("All"), "Total p.d.f.", "L");
    legend7->AddEntry(frame7->findObject("GS"), "Gaussian", "L");
    legend7->AddEntry(frame7->findObject("CB"), "Crystall Ball", "L");
    frame7->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame7->Draw();
    legend7->DrawClone();
    canvas7->SaveAs("fig/1D_psi2SCtauMix.png");
    pdf_psi2SCtauMix.getVariables()->Print("v");
    ////// Bdecay
    RooRealVar npsi2SCtauB("npsi2SCtauB", "npsi2SCtauB", 1, 0, N);
    RooExtendPdf pdf_psi2SCtauB("pdf_psi2SCtauB", "pdf_psi2SCtauB", psi2SCtauB, npsi2SCtauB);
    RooFitResult *res8;
    while(true) {
        res8 = pdf_psi2SCtauB.fitTo(*Bdecay->reduce(psi2S_ctau, "psi2S_ctau > -0.03 && psi2S_ctau < 0.25"), Save());
        if(!res8->status()) break;
    }
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 10, 10, 2400, 1800);
    RooPlot *frame8 = new RooPlot(psi2S_ctau, -0.03, 0.25, BinNum);
    frame8->SetTitle("psi2S Ctau 1D Bdecay");
    Bdecay->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_psi2SCtauB.plotOn(frame8, LineColor(1), LineWidth(2), Name("All"));
    pdf_psi2SCtauB.plotOn(frame8, Components(psi2S_ge1), LineColor(kBlue), LineWidth(2), Name("GExp1"));
    pdf_psi2SCtauB.plotOn(frame8, Components(psi2S_ge2), LineColor(kRed), LineWidth(2), Name("GExp2"));
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    legend8->AddEntry(frame8->findObject("All"), "Total p.d.f.", "L");
    legend8->AddEntry(frame8->findObject("GExp1"), "Gaussian Exp 1", "L");
    legend8->AddEntry(frame8->findObject("GExp2"), "Gaussian Exp 2", "L");
    frame8->SetAxisRange(1, 5e4, "Y");
    gPad->SetLogy();
    frame8->Draw();
    legend8->DrawClone();
    canvas8->SaveAs("fig/1D_psi2SCtauBdecay.png");
    pdf_psi2SCtauB.getVariables()->Print("v");

    // Form 4-dim p.d.f.
    double scale = 1500 / wsp->var("n_P_P")->getVal();
    RooRealVar vScale("vScale", "vScale", scale);
    RooProdPdf model_P_P("model_P_P", "model_P_P", RooArgList(JpsiMassMix, psi2SMassMix, JpsiCtauMix, psi2SCtauMix));
    RooProdPdf model_P_NP("model_P_NP", "model_P_NP", RooArgList(JpsiMassMix, psi2SMassB, JpsiCtauMix, psi2SCtauB));
    RooProdPdf model_NP_P("model_NP_P", "model_NP_P", RooArgList(JpsiMassB, psi2SMassMix, JpsiCtauB, psi2SCtauMix));
    RooProdPdf model_NP_NP("model_NP_NP", "model_NP_NP", RooArgList(JpsiMassB, psi2SMassB, JpsiCtauB, psi2SCtauB));
    RooProdPdf model_Sig_Comb("model_Sig_Comb", "model_Sig_Comb", RooArgList(JpsiMassCSig, psi2SMassCBkg, JpsiCtauCSig, psi2SCtauCBkg));
    RooProdPdf model_Comb_Sig("model_Comb_Sig", "model_Comb_Sig", RooArgList(JpsiMassCBkg, psi2SMassCSig, JpsiCtauCBkg, psi2SCtauCSig));
    RooProdPdf model_Comb_Comb("model_Comb_Comb", "model_Comb_Comb", RooArgList(JpsiMassCBkg, psi2SMassCBkg, JpsiCtauCBkg, psi2SCtauCBkg));
    RooRealVar y_P_P("y_P_P", "y_P_P", 1500);
    RooRealVar y_P_NP("y_P_NP", "y_P_NP", wsp->var("n_P_NP")->getVal() * scale);
    RooRealVar y_NP_P("y_NP_P", "y_NP_P", wsp->var("n_NP_P")->getVal() * scale);
    RooRealVar y_NP_NP("y_NP_NP", "y_NP_NP", wsp->var("n_NP_NP")->getVal() * scale);
    RooRealVar y_Sig_Comb("y_Sig_Comb", "y_Sig_Comb", wsp->var("n_Sig_Comb")->getVal() * scale);
    RooRealVar y_Comb_Sig("y_Comb_Sig", "y_Comb_Sig", wsp->var("n_Comb_Sig")->getVal() * scale);
    RooRealVar y_Comb_Comb("y_Comb_Comb", "y_Comb_Comb", wsp->var("n_Comb_Comb")->getVal() * scale);
    RooAddPdf model_all("model_all", "model_all",
        RooArgList(model_P_P, model_P_NP, model_NP_P, model_NP_NP, model_Sig_Comb, model_Comb_Sig, model_Comb_Comb),
        RooArgList(y_P_P, y_P_NP, y_NP_P, y_NP_NP, y_Sig_Comb, y_Comb_Sig, y_Comb_Comb)
    );
    // Generate data
    Jpsi_m1.setConstant(kTRUE);
    Jpsi_d1.setConstant(kTRUE);
    Jpsi_a1.setConstant(kTRUE);
    Jpsi_n1.setConstant(kTRUE);
    Jpsi_d2.setConstant(kTRUE);
    Jpsi_r1.setConstant(kTRUE);

    Jpsi_m2.setConstant(kTRUE);
    Jpsi_d3.setConstant(kTRUE);
    Jpsi_a2.setConstant(kTRUE);
    Jpsi_n2.setConstant(kTRUE);
    Jpsi_d4.setConstant(kTRUE);
    Jpsi_r2.setConstant(kTRUE);

    Jpsi_r3.setConstant(kTRUE);

    psi2S_m1.setConstant(kTRUE);
    psi2S_d1.setConstant(kTRUE);
    psi2S_a1.setConstant(kTRUE);
    psi2S_n1.setConstant(kTRUE);
    psi2S_d2.setConstant(kTRUE);
    psi2S_r1.setConstant(kTRUE);

    psi2S_m2.setConstant(kTRUE);
    psi2S_d3.setConstant(kTRUE);
    psi2S_a2.setConstant(kTRUE);
    psi2S_n2.setConstant(kTRUE);
    psi2S_d4.setConstant(kTRUE);
    psi2S_r2.setConstant(kTRUE);

    psi2S_r3.setConstant(kTRUE);

    psi2S_a3.setConstant(kTRUE);
    psi2S_b1.setConstant(kTRUE);

    Jpsi_u1.setConstant(kTRUE);
    Jpsi_s1.setConstant(kTRUE);
    Jpsi_s2.setConstant(kTRUE);
    Jpsi_p1.setConstant(kTRUE);

    Jpsi_s3.setConstant(kTRUE);
    Jpsi_s4.setConstant(kTRUE);
    Jpsi_e1.setConstant(kTRUE);
    Jpsi_e2.setConstant(kTRUE);
    Jpsi_p2.setConstant(kTRUE);

    Jpsi_p3.setConstant(kTRUE);

    Jpsi_u2.setConstant(kTRUE);
    Jpsi_s5.setConstant(kTRUE);
    Jpsi_e3.setConstant(kTRUE);
    Jpsi_p4.setConstant(kTRUE);

    psi2S_u1.setConstant(kTRUE);
    psi2S_s1.setConstant(kTRUE);
    psi2S_s2.setConstant(kTRUE);
    psi2S_l2.setConstant(kTRUE);
    psi2S_x2.setConstant(kTRUE);
    psi2S_p1.setConstant(kTRUE);

    psi2S_s3.setConstant(kTRUE);
    psi2S_s4.setConstant(kTRUE);
    psi2S_e1.setConstant(kTRUE);
    psi2S_e2.setConstant(kTRUE);
    psi2S_p2.setConstant(kTRUE);

    psi2S_p3.setConstant(kTRUE);

    psi2S_u2.setConstant(kTRUE);
    psi2S_l3.setConstant(kTRUE);
    psi2S_x3.setConstant(kTRUE);
    psi2S_s5.setConstant(kTRUE);
    psi2S_s6.setConstant(kTRUE);
    psi2S_e3.setConstant(kTRUE);
    psi2S_p4.setConstant(kTRUE);

    RooWorkspace *wk = new RooWorkspace("wk", "wk");
    wk->import(vScale);
    wk->import(model_all);
    wk->writeToFile("Data_4D.root");
    // Display fitting status
    cout<<"Status of all 1D fits:"<<endl;
    cout<<res1->status()<<' '<<res2->status()<<' '<<res3->status()<<' '<<res4->status()<<endl;
    cout<<res5->status()<<' '<<res6->status()<<' '<<res7->status()<<' '<<res8->status()<<endl;
    return;
}