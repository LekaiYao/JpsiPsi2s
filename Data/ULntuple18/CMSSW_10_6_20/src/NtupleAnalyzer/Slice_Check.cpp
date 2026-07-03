#include "Plot_4D.hpp"

void Slice_Check(string var, double vmin, double vmax) {
    // Define variables
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
	RooRealVar Jpsi_ctau("Jpsi_ctau", "Jpsi_ctau", -0.03, 0.16);
    RooRealVar psi2S_ctau("psi2S_ctau", "psi2S_ctau", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 100);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_ctau);
    variables.add(psi2S_ctau);
    variables.add(evt_weight);
    // Read file and apply cut
    string sel = var + " > " + to_string(vmin) + " && " + var + " < " + to_string(vmax);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
    
    // Parameters below are imported from total cross section 4D fit
    TFile *f = new TFile("Model_4D_tot.root");
    RooWorkspace *wsp = (RooWorkspace *)f->Get("wsp");
    RooAddPdf &pdf_all = dynamic_cast<RooAddPdf &>(wsp->allPdfs()["pdf_all"]);
    RooAbsPdf &pdf_P_P = *(wsp->pdf("pdf_P_P"));
    RooAbsPdf &pdf_P_NP = *(wsp->pdf("pdf_P_NP"));
    RooAbsPdf &pdf_NP_P = *(wsp->pdf("pdf_NP_P"));
    RooAbsPdf &pdf_NP_NP = *(wsp->pdf("pdf_NP_NP"));
    RooAbsPdf &pdf_Sig_Comb = *(wsp->pdf("pdf_Sig_Comb"));
    RooAbsPdf &pdf_Comb_Sig = *(wsp->pdf("pdf_Comb_Sig"));
    RooAbsPdf &pdf_Comb_Comb = *(wsp->pdf("pdf_Comb_Comb"));
    RooRealVar *cutVar = wsp->var(var.c_str());

    // Draw data point and p.d.f. curve
    string prefix = "fig/slice/";
    cutVar->setRange("central", vmin, vmax);
    Plot_4D(data->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau), sel.c_str()),
        pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau, "central");
    return;
}