// Lifetime systematic: replace ctau with Sig_Lxy in 4D fit.
// Compare n_P_P yield with Fit_4D_tot.cpp (ctau version) to estimate
// the lifetime-variable systematic uncertainty.
// Reference: JJ AN-22-160 §7.1

#include "Plot_4D.hpp"
#include <fstream>

void Fit_4D_tot_sigLxy() {
    // Define variables and read input file
    RooRealVar Jpsi_mass("Jpsi_mass", "Jpsi_mass", 2.95, 3.25);
    RooRealVar psi2S_mass("psi2S_mass", "psi2S_mass", 3.35, 4.05);
    RooRealVar Jpsi_sigLxy("Jpsi_sigLxy", "Jpsi_sigLxy", 0, 20);       // Sig_Lxy replaces ctau
    RooRealVar psi2S_sigLxy("psi2S_sigLxy", "psi2S_sigLxy", 0, 20);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooArgSet variables;
    variables.add(Jpsi_mass);
    variables.add(psi2S_mass);
    variables.add(Jpsi_sigLxy);
    variables.add(psi2S_sigLxy);
    variables.add(evt_weight);
    TFile *dataFile = new TFile("WeightData_prime.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");

    // Define J/psi Mass p.d.f.  (same as baseline)
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 1, 0.1, 5);
    RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 1, 0, 15);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooCBShape Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooGaussian Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());

    // Define psi(2S) Mass p.d.f.  (same as baseline)
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.686, 3.63, 3.73);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.02, 0.01, 0.05);
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));

    // ── J/psi Sig_Lxy p.d.f.  (Sig_Lxy ~ [0, 20], different scale from ctau) ──
    // Signal: double Gaussian
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.5, 1);
    RooRealVar Jpsi_sigL1("Jpsi_sigL1", "Jpsi_sigL1", 2, 0.5, 6);
    RooRealVar Jpsi_sigL2("Jpsi_sigL2", "Jpsi_sigL2", 5, 2, 15);
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_sigLxy, Jpsi_mu1, Jpsi_sigL1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_sigLxy, Jpsi_mu1, Jpsi_sigL2);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooAddPdf JpsiSigLxySig("JpsiSigLxySig", "JpsiSigLxySig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    // Background: exponential (B → long tail)
    RooRealVar Jpsi_sigL3("Jpsi_sigL3", "Jpsi_sigL3", 1, 0, 8);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 1, 0.1, 10);
    RooGExpModel JpsiSigLxyBkg("JpsiSigLxyBkg", "JpsiSigLxyBkg", Jpsi_sigLxy, Jpsi_sigL3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    RooAddPdf JpsiSigLxyCombSig("JpsiSigLxyCombSig", "JpsiSigLxyCombSig", RooArgList(JpsiSigLxySig, JpsiSigLxyBkg), Jpsi_prop2);
    // Combinatorial background
    RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0, -0.5, 1);
    RooRealVar Jpsi_sigL6("Jpsi_sigL6", "Jpsi_sigL6", 3, 1, 10);
    RooRealVar Jpsi_sigL7("Jpsi_sigL7", "Jpsi_sigL7", 6, 2, 20);
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 1, 0.1, 10);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_sigLxy, Jpsi_mu3, Jpsi_sigL6);
    RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_sigLxy, Jpsi_sigL7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiSigLxyCombBkg("JpsiSigLxyCombBkg", "JpsiSigLxyCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);

    // ── psi(2S) Sig_Lxy p.d.f. ──
    // Signal: Crystal Ball
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", 0, -0.5, 1);
    RooRealVar psi2S_sigL2("psi2S_sigL2", "psi2S_sigL2", 2, 0.5, 6);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -1, -5, 0);
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 5);
    RooCBShape psi2SSigLxySig("psi2SSigLxySig", "psi2SSigLxySig", psi2S_sigLxy, psi2S_mu1, psi2S_sigL2, psi2S_alpha2, psi2S_nx2);
    // Background: exponential
    RooRealVar psi2S_sigL3("psi2S_sigL3", "psi2S_sigL3", 1, 0, 8);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 1, 0.1, 10);
    RooGExpModel psi2SSigLxyBkg("psi2SSigLxyBkg", "psi2SSigLxyBkg", psi2S_sigLxy, psi2S_sigL3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal
    RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.5, 0, 1);
    RooAddPdf psi2SSigLxyCombSig("psi2SSigLxyCombSig", "psi2SSigLxyCombSig", RooArgList(psi2SSigLxySig, psi2SSigLxyBkg), psi2S_prop2);
    // Combinatorial background
    RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 0, -0.5, 1);
    RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.5, 0.1, 10);
    RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 3, 0.5, 20);
    RooRealVar psi2S_sigL6("psi2S_sigL6", "psi2S_sigL6", 3, 1, 10);
    RooRealVar psi2S_sigL7("psi2S_sigL7", "psi2S_sigL7", 6, 2, 20);
    RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 1, 0.1, 10);
    RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.5, 0, 1);
    RooCBShape psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_sigLxy, psi2S_mu3, psi2S_sigL6, psi2S_alpha3, psi2S_nx3);
    RooGExpModel psi2S_expGau3("psi2S_expGau3", "psi2S_expGau3", psi2S_sigLxy, psi2S_sigL7, psi2S_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SSigLxyCombBkg("psi2SSigLxyCombBkg", "psi2SSigLxyCombBkg", RooArgList(psi2S_gauss4, psi2S_expGau3), psi2S_prop3);

    // Form 4-dim p.d.f.
    Int_t N = 60000;
    // Mass dimension (same as baseline)
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    // Sig_Lxy dimension (replaces ctau)
    RooProdPdf pdf_sigL_PP("pdf_sigL_PP", "pdf_sigL_PP", JpsiSigLxySig, psi2SSigLxySig);
    RooProdPdf pdf_sigL_PNP("pdf_sigL_PNP", "pdf_sigL_PNP", JpsiSigLxySig, psi2SSigLxyBkg);
    RooProdPdf pdf_sigL_NPP("pdf_sigL_NPP", "pdf_sigL_NPP", JpsiSigLxyBkg, psi2SSigLxySig);
    RooProdPdf pdf_sigL_NPNP("pdf_sigL_NPNP", "pdf_sigL_NPNP", JpsiSigLxyBkg, psi2SSigLxyBkg);
    RooProdPdf pdf_sigL_SigBkg("pdf_sigL_SigBkg", "pdf_sigL_SigBkg", JpsiSigLxyCombSig, psi2SSigLxyCombBkg);
    RooProdPdf pdf_sigL_BkgSig("pdf_sigL_BkgSig", "pdf_sigL_BkgSig", JpsiSigLxyCombBkg, psi2SSigLxyCombSig);
    RooProdPdf pdf_sigL_BkgBkg("pdf_sigL_BkgBkg", "pdf_sigL_BkgBkg", JpsiSigLxyCombBkg, psi2SSigLxyCombBkg);
    // Combine Mass and Sig_Lxy dimensions
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_sigL_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_sigL_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_sigL_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_sigL_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_sigL_SigBkg);
    RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_sigL_BkgSig);
    RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_sigL_BkgBkg);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1e1, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1, N);
    RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 0, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e2, 0, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),
        RooArgList(n_P_P, n_P_NP, n_NP_P, n_NP_NP, n_Sig_Comb, n_Comb_Sig, n_Comb_Comb)
    );
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_sigLxy, psi2S_sigLxy)), Save());
        if(!res->status() && res->edm()<0.01) break;
    }

    // Draw
    string prefix = "fig/";
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass, psi2S_mass, Jpsi_sigLxy, psi2S_sigLxy);
    pdf_all.getVariables()->Print("v");

    // Set parameters constant and save workspace
    Jpsi_mean.setConstant(kTRUE);
    Jpsi_devia1.setConstant(kTRUE);
    Jpsi_alpha1.setConstant(kTRUE);
    Jpsi_nx1.setConstant(kTRUE);
    Jpsi_devia2.setConstant(kTRUE);
    Jpsi_ratio.setConstant(kTRUE);

    psi2S_mean.setConstant(kTRUE);
    psi2S_devia1.setConstant(kTRUE);

    Jpsi_mu1.setConstant(kTRUE);
    Jpsi_sigL1.setConstant(kTRUE);
    Jpsi_sigL2.setConstant(kTRUE);
    Jpsi_prop1.setConstant(kTRUE);
    Jpsi_sigL3.setConstant(kTRUE);
    Jpsi_coef1.setConstant(kTRUE);
    Jpsi_prop2.setConstant(kTRUE);
    Jpsi_mu3.setConstant(kTRUE);
    Jpsi_sigL6.setConstant(kTRUE);
    Jpsi_sigL7.setConstant(kTRUE);
    Jpsi_coef3.setConstant(kTRUE);
    Jpsi_prop3.setConstant(kTRUE);

    psi2S_mu1.setConstant(kTRUE);
    psi2S_sigL2.setConstant(kTRUE);
    psi2S_alpha2.setConstant(kTRUE);
    psi2S_nx2.setConstant(kTRUE);
    psi2S_sigL3.setConstant(kTRUE);
    psi2S_coef1.setConstant(kTRUE);
    psi2S_prop2.setConstant(kTRUE);
    psi2S_mu3.setConstant(kTRUE);
    psi2S_alpha3.setConstant(kTRUE);
    psi2S_nx3.setConstant(kTRUE);
    psi2S_sigL6.setConstant(kTRUE);
    psi2S_sigL7.setConstant(kTRUE);
    psi2S_coef3.setConstant(kTRUE);
    psi2S_prop3.setConstant(kTRUE);

    RooWorkspace *wsp = new RooWorkspace("wsp", "wsp");
    wsp->import(pdf_all);
    wsp->writeToFile("Model_4D_tot_sigLxy.root");

    // Compute total cross section
    const double L_int = 61.31;
    const double B_Jpsi = 0.05961;
    const double B_psi2S = 0.00794;
    double norm = 1e-3 / (L_int * B_Jpsi * B_psi2S);
    double sigma_tot = n_P_P.getVal() * norm;
    double sigma_tot_err = n_P_P.getError() * norm;

    // Write results
    ofstream fout("fit_results_4D_tot_sigLxy.txt");
    fout << "# 4D total fit results (Sig_Lxy variant — lifetime systematic)" << endl;
    fout << "# L_int = " << L_int << " fb^-1" << endl;
    fout << "# B(J/psi -> mu mu) = " << B_Jpsi << endl;
    fout << "# B(psi(2S) -> ee) = " << B_psi2S << endl;
    fout << "fit_status " << res->status() << endl;
    fout << "n_P_P " << n_P_P.getVal() << " " << n_P_P.getError() << endl;
    fout << "sigma_tot_pb " << sigma_tot << " " << sigma_tot_err << endl;
    fout.close();

    cout << "=== Fit_4D_tot_sigLxy (lifetime systematic) ===" << endl;
    cout << "Status: " << res->status() << endl;
    cout << "Event yield: " << n_P_P.getVal() << " +/- " << n_P_P.getError() << endl;
    cout << "sigma_tot [pb]: " << sigma_tot << " +/- " << sigma_tot_err << " (stat only)" << endl;
    cout << "Compare with fit_results_4D_tot.txt for lifetime systematic." << endl;
    return;
}
