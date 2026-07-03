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

void sPlot() {
    Int_t N = 60000, BinNum = 100;
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
	RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.03, 0.16);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_ctau);
    variables.add(psi2S_ctau);
    variables.add(evt_weight);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables);
    // TFile *fData = new TFile("Data_4D.root");
    // RooWorkspace *wk = (RooWorkspace *)fData->Get("wk");
    // RooAbsPdf *model_all = wk->pdf("model_all");
    // srand(time(0));
    // RooRandom::randomGenerator()->SetSeed(rand());
    // RooDataSet *data = model_all->generate(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau), Extended(kTRUE));
    
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.03, 0.08);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 10);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooGaussian Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1);
    RooCBShape Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    // Background p.d.f.
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());
    // Define psi(2S) Mass p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.686, 3.63, 3.73);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.02, 0.01, 0.05);// 0.05
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    // Background p.d.f.
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));

    // 2D fit to J/psi Mass + psi(2S) Mass
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    RooRealVar n_Sig_Sig("n_Sig_Sig", "n_Sig_Sig", 1e4, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e2, 1, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);
    RooAddPdf pdf_mass("pdf_mass", "pdf_mass",
        RooArgList(pdf_mass_SigSig, pdf_mass_SigComb, pdf_mass_CombSig, pdf_mass_CombComb),
        RooArgList(n_Sig_Sig, n_Sig_Comb, n_Comb_Sig, n_Comb_Comb)
    );
    RooFitResult *res1;
    while(true) {
        res1 = pdf_mass.fitTo(*data->reduce(RooArgSet(Jpsi_mass, psi2S_mass)), Save());
        if(!res1->status() && res1->edm()<0.01) break;
    }

    // Draw data point and p.d.f. curve
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 1500, 1000);
    RooPlot *frame7 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 2D"), Bins(BinNum));
    data->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_mass.plotOn(frame7, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_mass.plotOn(frame7, Components(pdf_mass_SigSig), LineColor(kBlue), LineStyle(kSolid), LineWidth(3), Name("Sig_Sig"));
    pdf_mass.plotOn(frame7, Components(pdf_mass_SigComb), LineColor(kRed), LineWidth(2), Name("Sig_Comb"));
    pdf_mass.plotOn(frame7, Components(pdf_mass_CombSig), LineColor(kGreen), LineWidth(2), Name("Comb_Sig"));
    pdf_mass.plotOn(frame7, Components(pdf_mass_CombComb), LineColor(kMagenta), LineWidth(2), Name("Comb_Comb"));
    TLegend *legend7 = new TLegend(.65, .60, .85, .85);
    legend7->AddEntry(frame7->findObject("Data"), "RunII 2018", "L");
    legend7->AddEntry(frame7->findObject("All"), "Total p.d.f.", "L");
    legend7->AddEntry(frame7->findObject("Sig_Sig"), "J/#psi, #psi(2S)", "L");
    legend7->AddEntry(frame7->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend7->AddEntry(frame7->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend7->AddEntry(frame7->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame7->Draw();
    legend7->DrawClone();
    canvas7->SaveAs("splot/2D_JpsiMass.png");
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 1500, 1000);
    RooPlot *frame8 = psi2S_mass.frame(RooFit::Title("psi2S Mass 2D"), Bins(BinNum));
    data->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_mass.plotOn(frame8, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_mass.plotOn(frame8, Components(pdf_mass_SigSig), LineColor(kBlue), LineStyle(kSolid), LineWidth(3), Name("Sig_Sig"));
    pdf_mass.plotOn(frame8, Components(pdf_mass_SigComb), LineColor(kRed), LineWidth(2), Name("Sig_Comb"));
    pdf_mass.plotOn(frame8, Components(pdf_mass_CombSig), LineColor(kGreen), LineWidth(2), Name("Comb_Sig"));
    pdf_mass.plotOn(frame8, Components(pdf_mass_CombComb), LineColor(kMagenta), LineWidth(2), Name("Comb_Comb"));
    TLegend *legend8 = new TLegend(.65, .60, .85, .85);
    legend8->AddEntry(frame8->findObject("Data"), "RunII 2018", "L");
    legend8->AddEntry(frame8->findObject("All"), "Total p.d.f.", "L");
    legend8->AddEntry(frame8->findObject("Sig_Sig"), "J/#psi, #psi(2S)", "L");
    legend8->AddEntry(frame8->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend8->AddEntry(frame8->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, #psi(2S)", "L");
    legend8->AddEntry(frame8->findObject("Comb_Comb"), "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    frame8->Draw();
    legend8->DrawClone();
    canvas8->SaveAs("splot/2D_psi2SMass.png");
    pdf_mass.getVariables()->Print("v");

    // Create sPlot from 2D mass fit
    RooStats::SPlot splot_mass("splot_mass", "splot_mass", *data, &pdf_mass, RooArgSet(n_Sig_Sig, n_Sig_Comb, n_Comb_Sig, n_Comb_Comb));
    RooDataSet *data_sw = new RooDataSet("data_sw", "data_sw", data, *data->get(), 0, "n_Sig_Sig_sw");

    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);//, -0.005, 0.005
    // RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0.002, 0.006);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.005, 0, 0.006);//0.005
    // RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma1);
    // RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    // RooRealVar Jpsi_a2("Jpsi_a2", "Jpsi_a2", 1, 0, 5);
    // RooRealVar Jpsi_x2("Jpsi_x2", "Jpsi_x2", 1, 0, 5);
    // RooCBShape Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2, Jpsi_a2, Jpsi_x2);
    // RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    // RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    RooGaussian JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.003, 0.001, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooGExpModel JpsiCtauBkg("JpsiCtauBkg", "JpsiCtauBkg", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Wrong PV p.d.f
    // RooRealVar Jpsi_mu2("Jpsi_mu2", "Jpsi_mu2", 0);
    // RooRealVar Jpsi_sigma4("Jpsi_sigma4", "Jpsi_sigma4", 0.01, 0.005, 0.2);
    // RooGaussian JpsiWrPV("JpsiWrPV", "JpsiWrPV", Jpsi_ctau, Jpsi_mu2, Jpsi_sigma4);

    // Define psi(2S) Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", 0, -0.005, 0.005);
    RooRealVar psi2S_s1("psi2S_s1", "psi2S_s1", 0.001, 0, 0.003);
    // RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.005, 0.002, 0.01);
    // RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1, 0, 5);
    // RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0.1, 5);
    // RooRealVar psi2S_p1("psi2S_p1", "psi2S_p1", 0.5, 0, 1);
    // RooGaussian psi2S_g3("psi2S_g3", "psi2S_g3", psi2S_ctau, psi2S_mu1, psi2S_s1);
    // RooCBShape psi2S_c3("psi2S_c3", "psi2S_c3", psi2S_ctau, psi2S_mu1, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    // RooAddPdf psi2SCtauSig("psi2SCtauMix", "psi2SCtauMix", RooArgList(psi2S_g3, psi2S_c3), psi2S_p1);
    // RooCBShape psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", psi2S_ctau, psi2S_mu1, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    RooGaussian psi2SCtauSig("psi2S_g3", "psi2S_g3", psi2S_ctau, psi2S_mu1, psi2S_s1);
    // Background p.d.f.
    RooRealVar psi2S_mu2("psi2S_mu2", "psi2S_mu2", 0, -0.005, 0.005);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.005, 0, 0.012);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", 1, 0.1, 5);
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 5);
    RooCBShape psi2S_c3("psi2S_c3", "psi2S_c3", psi2S_ctau, psi2S_mu2, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    RooRealVar psi2S_p1("psi2S_p1", "psi2S_p1", 0.5, 0, 1);
    RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.003, 0.001, 0.01);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.06, 0.01, 0.1);
    // RooGExpModel psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel psi2S_e1("psi2S_e1", "psi2S_e1", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", RooArgList(psi2S_c3, psi2S_e1), psi2S_p1);
    // Wrong PV p.d.f
    // RooRealVar psi2S_mu2("psi2S_mu2", "psi2S_mu2", 0);
    // RooRealVar psi2S_sigma4("psi2S_sigma4", "psi2S_sigma4", 0.1, 0.05, 0.5);
    // RooGaussian psi2SWrPV("psi2SWrPV", "psi2SWrPV", psi2S_ctau, psi2S_mu2, psi2S_sigma4);

    // 2D fit to J/psi Ctau + psi(2S) Ctau
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig, psi2SCtauSig);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig, psi2SCtauBkg);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg, psi2SCtauSig);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg, psi2SCtauBkg);
    // RooProdPdf pdf_ctau_WrPV("pdf_ctau_WrPV", "pdf_ctau_WrPV", JpsiWrPV, psi2SWrPV);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1, N);
    RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 1, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    // RooRealVar n_WrPV("n_WrPV", "n_WrPV", 1, 0.1, N);
    RooAddPdf pdf_ctau("pdf_ctau", "pdf_ctau",
        RooArgList(pdf_ctau_PP, pdf_ctau_PNP, pdf_ctau_NPP, pdf_ctau_NPNP),//, pdf_ctau_WrPV
        RooArgList(n_P_P, n_P_NP, n_NP_P, n_NP_NP)//, n_WrPV
    );
    RooFitResult *res2;
    while(true) {
        res2 = pdf_ctau.fitTo(*data_sw->reduce(RooArgSet(Jpsi_ctau, psi2S_ctau)), Save());
        if(!res2->status() && res2->edm()<0.01) break;
    }

    // Draw data point and p.d.f. curve
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 1500, 1000);
    RooPlot *frame9 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau 2D sPlot"), Bins(BinNum));
    data_sw->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_ctau.plotOn(frame9, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_ctau.plotOn(frame9, Components(pdf_ctau_PP), LineColor(kBlue), LineStyle(kSolid), LineWidth(3), Name("P_P"));
    pdf_ctau.plotOn(frame9, Components(pdf_ctau_PNP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_ctau.plotOn(frame9, Components(pdf_ctau_NPP), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_ctau.plotOn(frame9, Components(pdf_ctau_NPNP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    // pdf_ctau.plotOn(frame9, Components(pdf_ctau_WrPV), LineColor(kGreen), LineStyle(kDashed), LineWidth(2), Name("WrPV"));
    TLegend *legend9 = new TLegend(.65, .60, .85, .85);
    legend9->AddEntry(frame9->findObject("Data"), "RunII 2018", "L");
    legend9->AddEntry(frame9->findObject("All"), "Total p.d.f.", "L");
    legend9->AddEntry(frame9->findObject("P_P"), "prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("P_NP"), "prompt, non-prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_P"), "non-prompt, prompt", "L");
    legend9->AddEntry(frame9->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    frame9->SetAxisRange(0.01, 1e2, "Y");//2e2
    gPad->SetLogy();
    frame9->Draw();
    legend9->DrawClone();
    canvas9->SaveAs("splot/2D_JpsiCtau.png");
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 1500, 1000);
    RooPlot *frame0 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau 4D sPlot"), Bins(BinNum));
    data_sw->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
    pdf_ctau.plotOn(frame0, LineColor(kBlack), LineWidth(2), Name("All"));
    pdf_ctau.plotOn(frame0, Components(pdf_ctau_PP), LineColor(kBlue), LineStyle(kSolid), LineWidth(3), Name("P_P"));
    pdf_ctau.plotOn(frame0, Components(pdf_ctau_PNP), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("P_NP"));
    pdf_ctau.plotOn(frame0, Components(pdf_ctau_NPP), LineColor(kBlue), LineStyle(kDotted), LineWidth(2), Name("NP_P"));
    pdf_ctau.plotOn(frame0, Components(pdf_ctau_NPNP), LineColor(kBlue), LineStyle(kDashDotted), LineWidth(2), Name("NP_NP"));
    // pdf_ctau.plotOn(frame0, Components(pdf_ctau_WrPV), LineColor(kGreen), LineStyle(kDashed), LineWidth(2), Name("WrPV"));
    TLegend *legend0 = new TLegend(.65, .60, .85, .85);
    legend0->AddEntry(frame0->findObject("Data"), "RunII 2018", "L");
    legend0->AddEntry(frame0->findObject("All"), "Total p.d.f.", "L");
    legend0->AddEntry(frame0->findObject("P_P"), "prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("P_NP"), "prompt, non-prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_P"), "non-prompt, prompt", "L");
    legend0->AddEntry(frame0->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    frame0->SetAxisRange(0.01, 1e2, "Y");//2e2
    gPad->SetLogy();
    frame0->Draw();
    legend0->DrawClone();
    canvas0->SaveAs("splot/2D_psi2SCtau.png");

    // Create sPlot from 2D ctau fit
    RooStats::SPlot splot_ctau("splot_ctau", "splot_ctau", *data_sw, &pdf_ctau, RooArgSet(n_P_P, n_P_NP, n_NP_P, n_NP_NP));//, n_WrPV
    int nEntry = data_sw->numEntries();
    double Ncorr = 0, Nyield = 0;
    for(int i = 0; i < nEntry; i++) {
        Nyield += ((RooRealVar&)((*data_sw->get(i))["n_P_P_sw"])).getVal();
        Ncorr += ((RooRealVar&)((*data_sw->get(i))["n_P_P_sw"])).getVal() * ((RooRealVar&)((*data_sw->get(i))["evt_weight"])).getVal();
    }
    printf("\nSum of sWeights: %.2f\nSum of sWeight-times-correction-factors: %.2f\n", Nyield, Ncorr);
    // cout<<"Status: "<<res->status()<<endl;
    return;
}