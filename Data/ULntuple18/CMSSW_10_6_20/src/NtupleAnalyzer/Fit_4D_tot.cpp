#include "Plot_4D.hpp"
#include <fstream>

void Fit_4D_tot() {
    // Define variables and read input file
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
    TFile *dataFile = new TFile("WeightData_prime.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
    
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 1, 0.1, 5);
    RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 1, 0, 15);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooCBShape Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooGaussian Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    // Background p.d.f.
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());

    // Define psi(2S) Mass p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.686, 3.63, 3.73);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.02, 0.01, 0.05);
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    // Background p.d.f.
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));
    
    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.004, 0.002, 0.008);
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooGExpModel JpsiCtauBkg("JpsiCtauBkg", "JpsiCtauBkg", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    RooAddPdf JpsiCtauCombSig("JpsiCtauCombSig", "JpsiCtauCombSig", RooArgList(JpsiCtauSig, JpsiCtauBkg), Jpsi_prop2);
    // Combinatorial background p.d.f.
    RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.003, 0.001, 0.01);
    RooRealVar Jpsi_sigma7("Jpsi_sigma7", "Jpsi_sigma7", 0.003, 0.001, 0.02);
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.2);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau, Jpsi_mu3, Jpsi_sigma6);
    RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_ctau, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCombBkg("JpsiCtauCombBkg", "JpsiCtauCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);
    
    // Define psi(2S) Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", 0, -5e-4, 5e-4);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.005, 0.001, 0.012);
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -1, -5, 0, "");
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 5);
    RooCBShape psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", psi2S_ctau, psi2S_mu1, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    // Background p.d.f.
    RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.001, 0, 0.02);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.06, 0.01, 0.1);
    RooGExpModel psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.5, 0, 1);
    RooAddPdf psi2SCtauCombSig("psi2SCtauCombSig", "psi2SCtauCombSig", RooArgList(psi2SCtauSig, psi2SCtauBkg), psi2S_prop2);
    // Combinatorial background p.d.f.
    RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 0, -0.001, 0.001);
    RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.5, 0.1, 10);
    RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 3, 0.5, 20);
    RooRealVar psi2S_sigma6("psi2S_sigma6", "psi2S_sigma6", 0.001, 0, 0.007);
    RooRealVar psi2S_sigma7("psi2S_sigma7", "psi2S_sigma7", 0.001, 0, 0.02);
    RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 0.06, 0.01, 0.1);
    RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.5, 0, 1);
    RooCBShape psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_mu3, psi2S_sigma6, psi2S_alpha3, psi2S_nx3);
    RooGExpModel psi2S_expGau3("psi2S_expGau3", "psi2S_expGau3", psi2S_ctau, psi2S_sigma7, psi2S_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCombBkg("psi2SCtauCombBkg", "psi2SCtauCombBkg", RooArgList(psi2S_gauss4, psi2S_expGau3), psi2S_prop3);
    
    // Form 4-dim p.d.f.
    Int_t N = 60000;
    // J/psi Mass + psi(2S) Mass dimension
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    // J/psi Ctau + psi(2S) Ctau dimension
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
        res = pdf_all.fitTo(*data->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau)), Save());
        if(!res->status() && res->edm()<0.01) break;//
    }

    // Draw data point and p.d.f. curve
    string prefix = "fig/";
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau);
    pdf_all.getVariables()->Print("v");

    // Set parameters to constant save model to file
    Jpsi_mean.setConstant(kTRUE);
    Jpsi_devia1.setConstant(kTRUE);
    Jpsi_alpha1.setConstant(kTRUE);
    Jpsi_nx1.setConstant(kTRUE);
    Jpsi_devia2.setConstant(kTRUE);
    Jpsi_ratio.setConstant(kTRUE);

    psi2S_mean.setConstant(kTRUE);
    psi2S_devia1.setConstant(kTRUE);

    Jpsi_mu1.setConstant(kTRUE);
    Jpsi_sigma1.setConstant(kTRUE);
    Jpsi_sigma2.setConstant(kTRUE);
    Jpsi_prop1.setConstant(kTRUE);
    Jpsi_sigma3.setConstant(kTRUE);
    Jpsi_coef1.setConstant(kTRUE);
    Jpsi_prop2.setConstant(kTRUE);
    Jpsi_mu3.setConstant(kTRUE);
    Jpsi_sigma6.setConstant(kTRUE);
    Jpsi_sigma7.setConstant(kTRUE);
    Jpsi_coef3.setConstant(kTRUE);
    Jpsi_prop3.setConstant(kTRUE);

    psi2S_mu1.setConstant(kTRUE);
    psi2S_sigma2.setConstant(kTRUE);
    psi2S_alpha2.setConstant(kTRUE);
    psi2S_nx2.setConstant(kTRUE);
    psi2S_sigma3.setConstant(kTRUE);
    psi2S_coef1.setConstant(kTRUE);
    psi2S_prop2.setConstant(kTRUE);
    psi2S_mu3.setConstant(kTRUE);
    psi2S_alpha3.setConstant(kTRUE);
    psi2S_nx3.setConstant(kTRUE);
    psi2S_sigma6.setConstant(kTRUE);
    psi2S_sigma7.setConstant(kTRUE);
    psi2S_coef3.setConstant(kTRUE);
    psi2S_prop3.setConstant(kTRUE);
    RooWorkspace *wsp = new RooWorkspace("wsp", "wsp");
    wsp->import(pdf_all);
    wsp->writeToFile("Model_4D_tot.root");
    // Compute total cross section
    const double L_int = 61.31;          // fb^-1
    const double B_Jpsi = 0.05961;       // B(J/psi -> mu mu)
    const double B_psi2S = 0.00794;      // B(psi(2S) -> ee)
    double norm = 1e-3 / (L_int * B_Jpsi * B_psi2S);
    double sigma_tot = n_P_P.getVal() * norm;
    double sigma_tot_err = n_P_P.getError() * norm;

    // Write results to file
    ofstream fout("fit_results_4D_tot.txt");
    fout << "# 4D total fit results" << endl;
    fout << "# L_int = " << L_int << " fb^-1" << endl;
    fout << "# B(J/psi -> mu mu) = " << B_Jpsi << endl;
    fout << "# B(psi(2S) -> ee) = " << B_psi2S << endl;
    fout << "fit_status " << res->status() << endl;
    fout << "n_P_P " << n_P_P.getVal() << " " << n_P_P.getError() << endl;
    fout << "sigma_tot_pb " << sigma_tot << " " << sigma_tot_err << endl;
    fout.close();

    // Log message to screen
    cout << "Status: " << res->status() << endl;
    cout << "Event yield: " << n_P_P.getVal() << " +/- " << n_P_P.getError() << endl;
    cout << "sigma_tot [pb]: " << sigma_tot << " +/- " << sigma_tot_err << " (stat only)" << endl;
    cout << "Results written to fit_results_4D_tot.txt" << endl;
    return;
}