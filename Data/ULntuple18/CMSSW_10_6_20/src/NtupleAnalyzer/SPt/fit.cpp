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

void fit() {
    Int_t N = 60000, BinNum = 120;
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
    RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.03, 0.16);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.03, 0.16);
    RooRealVar evt_pt("evt_pt", "evt_pt", 40, 80);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_ctau);
    variables.add(psi2S_ctau);
    variables.add(evt_pt);
    variables.add(evt_weight);
    string bin = "40_80";
    TFile *dataFile = new TFile("../WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
    
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.09264);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.0190165);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.0413783);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.8017);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 2.06686);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.422988);
    RooGaussian Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    // Background p.d.f.
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());

    // Define psi(2S) Mass p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.68232);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.0379327);// 0.05
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    // Background p.d.f.
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));
    
    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", -9.26421e-5);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.00100054);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.00297797);//0.005
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.270188);
    RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.0026519);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.0375823);
    RooGExpModel JpsiCtauBkg("JpsiCtauBkg", "JpsiCtauBkg", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.299304);
    RooAddPdf JpsiCtauCombSig("JpsiCtauCombSig", "JpsiCtauCombSig", RooArgList(JpsiCtauSig, JpsiCtauBkg), Jpsi_prop2);
    // Combinatorial background p.d.f.
    RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0.00253055);
    RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.00684072);//0.02
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.0655437);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.670875);
    RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau, Jpsi_mu3, Jpsi_sigma6);
    RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_ctau, Jpsi_sigma6, Jpsi_coef3, false, RooGExpModel::Type::Flipped);//Jpsi_sigma7
    RooAddPdf JpsiCtauCombBkg("JpsiCtauCombBkg", "JpsiCtauCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);
    
    // Define psi(2S) Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", -0.000104981);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.00228342);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -3.36064);
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 0.0239843);
    RooCBShape psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", psi2S_ctau, psi2S_mu1, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    // Background p.d.f.
    RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.00265165);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.0398886);
    RooGExpModel psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.887938);
    RooAddPdf psi2SCtauCombSig("psi2SCtauCombSig", "psi2SCtauCombSig", RooArgList(psi2SCtauSig, psi2SCtauBkg), psi2S_prop2);
    // Combinatorial background p.d.f.
    RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 4.54293e-5);
    RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.1051);
    RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 2.53908);
    RooRealVar psi2S_sigma6("psi2S_sigma6", "psi2S_sigma6", 0.00364727);
    RooRealVar psi2S_sigma7("psi2S_sigma7", "psi2S_sigma7", 0.00305531);
    RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 0.0319934);
    RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.496629);
    RooCBShape psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_mu3, psi2S_sigma6, psi2S_alpha3, psi2S_nx3);
    RooGExpModel psi2S_expGau3("psi2S_expGau3", "psi2S_expGau3", psi2S_ctau, psi2S_sigma7, psi2S_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCombBkg("psi2SCtauCombBkg", "psi2SCtauCombBkg", RooArgList(psi2S_gauss4, psi2S_expGau3), psi2S_prop3);
    
    // Form 4-dim p.d.f.
    // J/psi Mass + psi(2S) Mass dimension
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    // J/psi Mass + psi(2S) Ctau dimension
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig, psi2SCtauSig);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig, psi2SCtauBkg);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg, psi2SCtauSig);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg, psi2SCtauBkg);
    RooProdPdf pdf_ctau_SigBkg("pdf_ctau_SigBkg", "pdf_ctau_SigBkg", JpsiCtauCombSig, psi2SCtauCombBkg);
    RooProdPdf pdf_ctau_BkgSig("pdf_ctau_BkgSig", "pdf_ctau_BkgSig", JpsiCtauCombBkg, psi2SCtauCombSig);
    RooProdPdf pdf_ctau_BkgBkg("pdf_ctau_BkgBkg", "pdf_ctau_BkgBkg", JpsiCtauCombBkg, psi2SCtauCombBkg);
    // Combine Mass and Ctau dimensions
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_ctau_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_ctau_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_ctau_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_ctau_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_ctau_SigBkg);
    RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_ctau_BkgSig);
    RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_ctau_BkgBkg);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 0, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 0, N);
    RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 0, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e2, 0, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e3, 0, N);
    RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e2, 0, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 0, N);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),
        RooArgList(n_P_P, n_P_NP, n_NP_P, n_NP_NP, n_Sig_Comb, n_Comb_Sig, n_Comb_Comb)
    );
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau)), Save());
        if(!res->status()) break;
    }
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
    pdf_all.plotOn(frame7, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame7, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    legend7->AddEntry(frame7->findObject("All"), "Total p.d.f.", "L");
    legend7->AddEntry(frame7->findObject("P_P"), "prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("P_NP"), "prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_P"), "non-prompt, prompt", "L");
    legend7->AddEntry(frame7->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend7->AddEntry(frame7->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend7->AddEntry(frame7->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend7->AddEntry(frame7->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame7->Draw();
    legend7->DrawClone();
    // canvas7->SaveAs("fits/4D_JpsiMass.pdf");
    canvas7->SaveAs(("fits_" + bin + "/4D_JpsiMass.png").c_str());
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
    pdf_all.plotOn(frame8, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame8, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    legend8->AddEntry(frame8->findObject("All"), "Total p.d.f.", "L");
    legend8->AddEntry(frame8->findObject("P_P"), "prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("P_NP"), "prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_P"), "non-prompt, prompt", "L");
    legend8->AddEntry(frame8->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend8->AddEntry(frame8->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend8->AddEntry(frame8->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend8->AddEntry(frame8->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame8->Draw();
    legend8->DrawClone();
    // canvas8->SaveAs("fits/4D_psi2SMass.pdf");
    canvas8->SaveAs(("fits_" + bin + "/4D_psi2SMass.png").c_str());
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
    pdf_all.plotOn(frame9, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame9, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend9 = new TLegend(.65, .60, .85, .85);
    legend9->AddEntry(frame9->findObject("Data"), "RunII 2018", "L");
    legend9->AddEntry(frame9->findObject("All"), "Total p.d.f.", "L");
    legend9->AddEntry(frame9->findObject("P_P"), "prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("P_NP"), "prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_P"), "non-prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend9->AddEntry(frame9->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend9->AddEntry(frame9->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame9->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame9->Draw();
    legend9->DrawClone();
    // canvas9->SaveAs("fits/4D_JpsiCtau.pdf");
    canvas9->SaveAs(("fits_" + bin + "/4D_JpsiCtau.png").c_str());
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
    pdf_all.plotOn(frame0, Components(pdf_Comb_Sig), LineColor(kGreen), LineWidth(1), Name("Comb_Sig"));
    pdf_all.plotOn(frame0, Components(pdf_Comb_Comb), LineColor(kMagenta), LineWidth(1), Name("Comb_Comb"));
    TLegend *legend0 = new TLegend(.65, .60, .85, .85);
    legend0->AddEntry(frame0->findObject("Data"), "RunII 2018", "L");
    legend0->AddEntry(frame0->findObject("All"), "Total p.d.f.", "L");
    legend0->AddEntry(frame0->findObject("P_P"), "prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("P_NP"), "prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_P"), "non-prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend0->AddEntry(frame0->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend0->AddEntry(frame0->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame0->SetAxisRange(1, 5e3, "Y");
    gPad->SetLogy();
    frame0->Draw();
    legend0->DrawClone();
    // canvas0->SaveAs("fits/4D_psi2SCtau.pdf");
    canvas0->SaveAs(("fits_" + bin + "/4D_psi2SCtau.png").c_str());
    // canvas0->SaveAs("fits/4D_psi2SCtau.png");
    pdf_all.getVariables()->Print("v");
    return;
}