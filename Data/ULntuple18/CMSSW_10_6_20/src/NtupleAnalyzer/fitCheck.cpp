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
#include "RooMCStudy.h"
#include "RooFitResult.h"

using namespace std;
using namespace RooFit;

void fitCheck() {
    double scale = 1;
    // Define fitting variables
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
    // Read file: SPS, DPS, Bdecay
    // TFile *mixFile = new TFile("WeightMix.root", "READ");
    // TTree *mixTree = (TTree*)mixFile->Get("data");
    // RooDataSet *mix = new RooDataSet("mix", "mix", mixTree, variables, "", "evt_weight");
    
    ////////////////////////////////
    // Model generation
    ////////////////////////////////

    // Define p.d.f.
    //// J/psi Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar Jpsi_m1("Jpsi_m1", "Jpsi_m1", 3.09708);
    RooRealVar Jpsi_d1("Jpsi_d1", "Jpsi_d1", 0.0226933);
    RooRealVar Jpsi_a1("Jpsi_a1", "Jpsi_a1", 1.72006);
    RooRealVar Jpsi_n1("Jpsi_n1", "Jpsi_n1", 2.01234);
    RooRealVar Jpsi_d2("Jpsi_d2", "Jpsi_d2", 0.0454891);
    RooRealVar Jpsi_r1("Jpsi_r1", "Jpsi_r1", 0.729023);
    RooCBShape Jpsi_c1("Jpsi_c1", "Jpsi_c1", Jpsi_mass, Jpsi_m1, Jpsi_d1, Jpsi_a1, Jpsi_n1);
    RooGaussian Jpsi_g1("Jpsi_g1", "Jpsi_g1", Jpsi_mass, Jpsi_m1, Jpsi_d2);
    RooAddPdf JpsiMassMix("JpsiMassMix", "JpsiMassMix", RooArgList(Jpsi_c1, Jpsi_g1), Jpsi_r1);
    ////// Bdecay
    RooRealVar Jpsi_m2("Jpsi_m2", "Jpsi_m2", 3.09628);
    RooRealVar Jpsi_d3("Jpsi_d3", "Jpsi_d3", 0.0223427);
    RooRealVar Jpsi_a2("Jpsi_a2", "Jpsi_a2", 1.82166);
    RooRealVar Jpsi_n2("Jpsi_n2", "Jpsi_n2", 1.82579);
    RooRealVar Jpsi_d4("Jpsi_d4", "Jpsi_d4", 0.0428065);
    RooRealVar Jpsi_r2("Jpsi_r2", "Jpsi_r2", 0.668373);
    RooCBShape Jpsi_c2("Jpsi_c2", "Jpsi_c2", Jpsi_mass, Jpsi_m2, Jpsi_d3, Jpsi_a2, Jpsi_n2);
    RooGaussian Jpsi_g2("Jpsi_g2", "Jpsi_g2", Jpsi_mass, Jpsi_m2, Jpsi_d4);
    RooAddPdf JpsiMassB("JpsiMassB", "JpsiMassB", RooArgList(Jpsi_c2, Jpsi_g2), Jpsi_r2);
    ////// Combinatorial signal
    RooRealVar Jpsi_r3("Jpsi_r3", "Jpsi_r3", 0.299304);
    RooAddPdf JpsiMassCSig("JpsiMassCSig", "JpsiMassCSig", RooArgList(JpsiMassMix, JpsiMassB), Jpsi_r3);
    ////// Combinatorial background
    RooChebychev JpsiMassCBkg("JpsiMassCBkg", "JpsiMassCBkg", Jpsi_mass, RooArgList());
    //// psi(2S) Mass
    ////// Mixture of SPS/DPS=1:2
    RooRealVar psi2S_m1("psi2S_m1", "psi2S_m1", 3.68653);
    RooRealVar psi2S_d1("psi2S_d1", "psi2S_d1", 0.0263221);
    RooRealVar psi2S_a1("psi2S_a1", "psi2S_a1", 1.69063);
    RooRealVar psi2S_n1("psi2S_n1", "psi2S_n1", 2.42646);
    RooRealVar psi2S_d2("psi2S_d2", "psi2S_d2", 0.0539009);
    RooRealVar psi2S_r1("psi2S_r1", "psi2S_r1", 0.696941);
    RooCBShape psi2S_c1("psi2S_c1", "psi2S_c1", psi2S_mass, psi2S_m1, psi2S_d1, psi2S_a1, psi2S_n1);
    RooGaussian psi2S_g1("psi2S_g1", "psi2S_g1", psi2S_mass, psi2S_m1, psi2S_d2);
    RooAddPdf psi2SMassMix("psi2SMassMix", "psi2SMassMix", RooArgList(psi2S_c1, psi2S_g1), psi2S_r1);
    ////// Bdecay
    RooRealVar psi2S_m2("psi2S_m2", "psi2S_m2", 3.68564);
    RooRealVar psi2S_d3("psi2S_d3", "psi2S_d3", 0.0249261);
    RooRealVar psi2S_a2("psi2S_a2", "psi2S_a2", 1.82743);
    RooRealVar psi2S_n2("psi2S_n2", "psi2S_n2", 1.48718);
    RooRealVar psi2S_d4("psi2S_d4", "psi2S_d4", 0.0482346);
    RooRealVar psi2S_r2("psi2S_r2", "psi2S_r2", 0.630784);
    RooCBShape psi2S_c2("psi2S_c2", "psi2S_c2", psi2S_mass, psi2S_m2, psi2S_d3, psi2S_a2, psi2S_n2);
    RooGaussian psi2S_g2("psi2S_g2", "psi2S_g2", psi2S_mass, psi2S_m2, psi2S_d4);
    RooAddPdf psi2SMassB("psi2SMassB", "psi2SMassB", RooArgList(psi2S_c2, psi2S_g2), psi2S_r2);
    ////// Combinatorial signal
    RooRealVar psi2S_r3("psi2S_r3", "psi2S_r3", 0.887774);
    RooAddPdf psi2SMassCSig("psi2SMassCSig", "psi2SMassCSig", RooArgList(psi2SMassMix, psi2SMassB), psi2S_r3);
    ////// Combinatorial background
    RooRealVar psi2S_a3("psi2S_a3", "psi2S_a3", -0.441759);
    RooRealVar psi2S_b1("psi2S_b1", "psi2S_b1", 0.0671148);
    RooChebychev psi2SMassCBkg("psi2SMassCBkg", "psi2SMassCBkg", psi2S_mass, RooArgList(psi2S_a3, psi2S_b1));
    //// J/psi Ctau
    ////// Mixture of SPS/DPS=1:2
    RooRealVar Jpsi_u1("Jpsi_u1", "Jpsi_u1", 1.5898e-05);
    RooRealVar Jpsi_s1("Jpsi_s1", "Jpsi_s1", 0.00177877);
    RooRealVar Jpsi_s2("Jpsi_s2", "Jpsi_s2", 0.00518729);
    RooGaussian Jpsi_g3("Jpsi_g3", "Jpsi_g3", Jpsi_ctau, Jpsi_u1, Jpsi_s1);
    RooGaussian Jpsi_g4("Jpsi_g4", "Jpsi_g4", Jpsi_ctau, Jpsi_u1, Jpsi_s2);
    RooRealVar Jpsi_p1("Jpsi_p1", "Jpsi_p1", 0.890139);
    RooAddPdf JpsiCtauMix("JpsiCtauMix", "JpsiCtauMix", RooArgList(Jpsi_g3, Jpsi_g4), Jpsi_p1);
    ////// Bdecay
    RooRealVar Jpsi_s3("Jpsi_s3", "Jpsi_s3", 0.0018038);
    RooRealVar Jpsi_s4("Jpsi_s4", "Jpsi_s4", 0.00615102);
    RooRealVar Jpsi_e1("Jpsi_e1", "Jpsi_e1", 0.0320478);
    RooRealVar Jpsi_e2("Jpsi_e2", "Jpsi_e2", 0.061382);
    RooRealVar Jpsi_p2("Jpsi_p2", "Jpsi_p2", 0.887466);
    RooGExpModel Jpsi_ge1("Jpsi_ge1", "Jpsi_ge1", Jpsi_ctau, Jpsi_s3, Jpsi_e1, false, RooGExpModel::Type::Flipped);
    RooGExpModel Jpsi_ge2("Jpsi_ge2", "Jpsi_ge2", Jpsi_ctau, Jpsi_s4, Jpsi_e2, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauB("JpsiCtauB", "JpsiCtauB", RooArgList(Jpsi_ge1, Jpsi_ge2), Jpsi_p2);
    ////// Combinatorial signal
    RooRealVar Jpsi_p3("Jpsi_p3", "Jpsi_p3", 0.299304);
    RooAddPdf JpsiCtauCSig("JpsiCtauCSig", "JpsiCtauCSig", RooArgList(JpsiCtauMix, JpsiCtauB), Jpsi_p3);
    ////// Combinatorial background
    RooRealVar Jpsi_u2("Jpsi_u2", "Jpsi_u2", 0.00253027);
    RooRealVar Jpsi_s5("Jpsi_s5", "Jpsi_s5", 0.00684047);
    RooRealVar Jpsi_e3("Jpsi_e3", "Jpsi_e3", 0.0655259);
    RooRealVar Jpsi_p4("Jpsi_p4", "Jpsi_p4", 0.670822);
    RooGaussian Jpsi_g5("Jpsi_g5", "Jpsi_g5", Jpsi_ctau, Jpsi_u2, Jpsi_s5);
    RooGExpModel Jpsi_ge3("Jpsi_ge3", "Jpsi_ge3", Jpsi_ctau, Jpsi_s5, Jpsi_e3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCBkg("JpsiCtauCBkg", "JpsiCtauCBkg", RooArgList(Jpsi_g5, Jpsi_ge3), Jpsi_p4);
    //// psi(2S) Ctau
    ////// Signal
    RooRealVar psi2S_u1("psi2S_u1", "psi2S_u1", 1.17858e-05);
    RooRealVar psi2S_s1("psi2S_s1", "psi2S_s1", 0.00164839);
    RooRealVar psi2S_s2("psi2S_s2", "psi2S_s2", 0.00410447);
    RooRealVar psi2S_l2("psi2S_l2", "psi2S_l2", -2.15175);
    RooRealVar psi2S_x2("psi2S_x2", "psi2S_x2", 1.12331);
    RooRealVar psi2S_p1("psi2S_p1", "psi2S_p1", 0.804393);
    RooGaussian psi2S_g3("psi2S_g3", "psi2S_g3", psi2S_ctau, psi2S_u1, psi2S_s1);
    RooCBShape psi2S_c3("psi2S_c3", "psi2S_c3", psi2S_ctau, psi2S_u1, psi2S_s2, psi2S_l2, psi2S_x2);
    RooAddPdf psi2SCtauMix("psi2SCtauMix", "psi2SCtauMix", RooArgList(psi2S_g3, psi2S_c3), psi2S_p1);
    ////// Bdecay
    RooRealVar psi2S_s3("psi2S_s3", "psi2S_s3", 0.00164627);
    RooRealVar psi2S_s4("psi2S_s4", "psi2S_s4", 0.00522589);
    RooRealVar psi2S_e1("psi2S_e1", "psi2S_e1", 0.035629);
    RooRealVar psi2S_e2("psi2S_e2", "psi2S_e2", 0.055781);
    RooRealVar psi2S_p2("psi2S_p2", "psi2S_p2", 0.810501);
    RooGExpModel psi2S_ge1("psi2S_ge1", "psi2S_ge1", psi2S_ctau, psi2S_s3, psi2S_e1, false, RooGExpModel::Type::Flipped);
    RooGExpModel psi2S_ge2("psi2S_ge2", "psi2S_ge2", psi2S_ctau, psi2S_s4, psi2S_e2, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauB("psi2SCtauB", "psi2SCtauB", RooArgList(psi2S_ge1, psi2S_ge2), psi2S_p2);
    ////// Combinatorial signal
    RooRealVar psi2S_p3("psi2S_p3", "psi2S_p3", 0.887774);
    RooAddPdf psi2SCtauCSig("psi2SCtauCSig", "psi2SCtauCSig", RooArgList(psi2SCtauMix, psi2SCtauB), psi2S_p3);
    ////// Combinatorial background
    RooRealVar psi2S_u2("psi2S_u2", "psi2S_u2", 4.5446e-5);
    RooRealVar psi2S_l3("psi2S_l3", "psi2S_l3", 1.1051);
    RooRealVar psi2S_x3("psi2S_x3", "psi2S_x3", 2.53902);
    RooRealVar psi2S_s5("psi2S_s5", "psi2S_s5", 0.00364726);
    RooRealVar psi2S_s6("psi2S_s6", "psi2S_s6", 0.0030555);
    RooRealVar psi2S_e3("psi2S_e3", "psi2S_e3", 0.0319934);
    RooRealVar psi2S_p4("psi2S_p4", "psi2S_p4", 0.496629);
    RooCBShape psi2S_g4("psi2S_g4", "psi2S_g4", psi2S_ctau, psi2S_u2, psi2S_s5, psi2S_l3, psi2S_x3);
    RooGExpModel psi2S_ge3("psi2S_ge3", "psi2S_ge3", psi2S_ctau, psi2S_s6, psi2S_e3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCBkg("psi2SCtauCBkg", "psi2SCtauCBkg", RooArgList(psi2S_g4, psi2S_ge3), psi2S_p4);
    // Form 4-dim p.d.f.
    RooProdPdf model_P_P("model_P_P", "model_P_P", RooArgList(JpsiMassMix, psi2SMassMix, JpsiCtauMix, psi2SCtauMix));
    RooProdPdf model_P_NP("model_P_NP", "model_P_NP", RooArgList(JpsiMassMix, psi2SMassB, JpsiCtauMix, psi2SCtauB));
    RooProdPdf model_NP_P("model_NP_P", "model_NP_P", RooArgList(JpsiMassB, psi2SMassMix, JpsiCtauB, psi2SCtauMix));
    RooProdPdf model_NP_NP("model_NP_NP", "model_NP_NP", RooArgList(JpsiMassB, psi2SMassB, JpsiCtauB, psi2SCtauB));
    RooProdPdf model_Sig_Comb("model_Sig_Comb", "model_Sig_Comb", RooArgList(JpsiMassCSig, psi2SMassCBkg, JpsiCtauCSig, psi2SCtauCBkg));
    RooProdPdf model_Comb_Sig("model_Comb_Sig", "model_Comb_Sig", RooArgList(JpsiMassCBkg, psi2SMassCSig, JpsiCtauCBkg, psi2SCtauCSig));
    RooProdPdf model_Comb_Comb("model_Comb_Comb", "model_Comb_Comb", RooArgList(JpsiMassCBkg, psi2SMassCBkg, JpsiCtauCBkg, psi2SCtauCBkg));
    RooRealVar y_P_P("y_P_P", "y_P_P", 1472.15 * scale);
    RooRealVar y_P_NP("y_P_NP", "y_P_NP", 376.853 * scale);
    RooRealVar y_NP_P("y_NP_P", "y_NP_P", 219.528 * scale);
    RooRealVar y_NP_NP("y_NP_NP", "y_NP_NP", 3066.57 * scale);
    RooRealVar y_Sig_Comb("y_Sig_Comb", "y_Sig_Comb", 15554.5 * scale);
    RooRealVar y_Comb_Sig("y_Comb_Sig", "y_Comb_Sig", 404.173 * scale);
    RooRealVar y_Comb_Comb("y_Comb_Comb", "y_Comb_Comb", 1807.7 * scale);
    RooAddPdf model_all("model_all", "model_all",
        RooArgList(model_P_P, model_P_NP, model_NP_P, model_NP_NP, model_Sig_Comb, model_Comb_Sig, model_Comb_Comb),
        RooArgList(y_P_P, y_P_NP, y_NP_P, y_NP_NP, y_Sig_Comb, y_Comb_Sig, y_Comb_Comb)
    );
    // Generate model: P+NP, combinatorial
    srand(time(0));
    RooRandom::randomGenerator()->SetSeed(rand());
    RooDataSet *mix = model_all.generate(variables, Extended(kTRUE));

    ////////////////////////////////
    // PDF and fitting
    ////////////////////////////////

    // Define p.d.f.
    //// J/psi Mass
    ////// Signal
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.03, 0.08);
    RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.5, 0.1, 3.5);
    RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 1, 0, 10);//!
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    // RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.09264);
    // RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.0190156);
    // RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.0413757);
    // RooRealVar Jpsi_alpha2("Jpsi_alpha2", "Jpsi_alpha2", 1.80134);
    // RooRealVar Jpsi_nx2("Jpsi_nx2", "Jpsi_nx2", 2.07041);
    // RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.422938);
    RooGaussian Jpsi_crysBall1("Jpsi_crysBall1", "Jpsi_crysBall1", Jpsi_mass, Jpsi_mean, Jpsi_devia1);
    RooGaussian Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2);
    // RooCBShape Jpsi_crysBall2("Jpsi_crysBall2", "Jpsi_crysBall2", Jpsi_mass, Jpsi_mean, Jpsi_devia2, Jpsi_alpha2, Jpsi_nx2);
    RooAddPdf JpsiMassSig("JpsiMassSig", "JpsiMassSig", RooArgList(Jpsi_crysBall1, Jpsi_crysBall2), Jpsi_ratio);
    ////// Background
    RooChebychev JpsiMassComb("JpsiMassComb", "JpsiMassComb", Jpsi_mass, RooArgList());
    //// psi(2S) Mass
    ////// Signal
    RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.686, 3.63, 3.73);
    RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.02, 0.01, 0.05);
    // RooRealVar psi2S_mean("psi2S_mean", "psi2S_mean", 3.68232);
    // RooRealVar psi2S_devia1("psi2S_devia1", "psi2S_devia1", 0.0379327);
    RooGaussian psi2SMassSig("psi2SMassSig", "psi2SMassSig", psi2S_mass, psi2S_mean, psi2S_devia1);
    ////// Background
    RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.5, -1, 0.1);
    RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0, -0.5, 0.5);
    // RooRealVar psi2S_a("psi2S_a", "psi2S_a", -0.441759);//c
    // RooRealVar psi2S_b("psi2S_b", "psi2S_b", 0.0671148);//c
    RooChebychev psi2SMassComb("psi2SMassComb", "psi2SMassComb", psi2S_mass, RooArgList(psi2S_a, psi2S_b));
    //// J/psi Ctau
    ////// Signal
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.002);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.003, 0.002, 0.006);//0.005
    // RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", -9.26824e-5);
    // RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.00100051);
    // RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.00297796);
    RooGaussian Jpsi_gauss1("Jpsi_gauss1", "Jpsi_gauss1", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2("Jpsi_gauss2", "Jpsi_gauss2", Jpsi_ctau, Jpsi_mu1, Jpsi_sigma2);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    // RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.270181);
    RooAddPdf JpsiCtauSig("JpsiCtauSig", "JpsiCtauSig", RooArgList(Jpsi_gauss1, Jpsi_gauss2), Jpsi_prop1);
    ////// Background
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    // RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.0026517);
    // RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.0375823);
    RooGExpModel JpsiCtauBkg("JpsiCtauBkg", "JpsiCtauBkg", Jpsi_ctau, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    ////// Combinatorial signal
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    // RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.299304);
    RooAddPdf JpsiCtauCombSig("JpsiCtauCombSig", "JpsiCtauCombSig", RooArgList(JpsiCtauSig, JpsiCtauBkg), Jpsi_prop2);
    ////// Combinatorial background
    RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.003, 0.001, 0.01);//0.02
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.1);
    RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.5, 0, 1);
    // RooRealVar Jpsi_mu3("Jpsi_mu3", "Jpsi_mu3", 0.00253027);//c
    // RooRealVar Jpsi_sigma6("Jpsi_sigma6", "Jpsi_sigma6", 0.00684047);//c
    // RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.0655259);//c
    // RooRealVar Jpsi_prop3("Jpsi_prop3", "Jpsi_prop3", 0.670822);//c
    RooGaussian Jpsi_gauss4("Jpsi_gauss4", "Jpsi_gauss4", Jpsi_ctau, Jpsi_mu3, Jpsi_sigma6);
    RooGExpModel Jpsi_expGau3("Jpsi_expGau3", "Jpsi_expGau3", Jpsi_ctau, Jpsi_sigma6, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf JpsiCtauCombBkg("JpsiCtauCombBkg", "JpsiCtauCombBkg", RooArgList(Jpsi_gauss4, Jpsi_expGau3), Jpsi_prop3);
    //// psi(2S) Ctau
    ////// Signal
    RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", 0, -5e-4, 5e-4);
    RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.005, 0.001, 0.012);//0.002
    RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -1, -5, 0, "");
    RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 1, 0, 5);
    // RooRealVar psi2S_mu1("psi2S_mu1", "psi2S_mu1", -0.000104962);
    // RooRealVar psi2S_sigma2("psi2S_sigma2", "psi2S_sigma2", 0.00228344);
    // RooRealVar psi2S_alpha2("psi2S_alpha2", "psi2S_alpha2", -3.36028);
    // RooRealVar psi2S_nx2("psi2S_nx2", "psi2S_nx2", 0.239865);
    RooCBShape psi2SCtauSig("psi2SCtauSig", "psi2SCtauSig", psi2S_ctau, psi2S_mu1, psi2S_sigma2, psi2S_alpha2, psi2S_nx2);
    ////// Background
    RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.001, 0, 0.02);
    RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.06, 0.01, 0.1);
    // RooRealVar psi2S_sigma3("psi2S_sigma3", "psi2S_sigma3", 0.00265165);
    // RooRealVar psi2S_coef1("psi2S_coef1", "psi2S_coef1", 0.0398885);
    RooGExpModel psi2SCtauBkg("psi2SCtauBkg", "psi2SCtauBkg", psi2S_ctau, psi2S_sigma3, psi2S_coef1, false, RooGExpModel::Type::Flipped);
    ////// Combinatorial signal
    RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.5, 0, 1);
    // RooRealVar psi2S_prop2("psi2S_prop2", "psi2S_prop2", 0.887774);
    RooAddPdf psi2SCtauCombSig("psi2SCtauCombSig", "psi2SCtauCombSig", RooArgList(psi2SCtauSig, psi2SCtauBkg), psi2S_prop2);
    ////// Combinatorial background
    RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 0, -0.001, 0.001);
    RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.5, 0.1, 10);
    RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 3, 0.5, 20);
    RooRealVar psi2S_sigma6("psi2S_sigma6", "psi2S_sigma6", 0.001, 0, 0.007);
    RooRealVar psi2S_sigma7("psi2S_sigma7", "psi2S_sigma7", 0.001, 0, 0.02);
    RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 0.06, 0.01, 0.1);
    RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.5, 0, 1);
    // RooRealVar psi2S_mu3("psi2S_mu3", "psi2S_mu3", 4.5446e-5);//c
    // RooRealVar psi2S_alpha3("psi2S_alpha3", "psi2S_alpha3", 1.1051);//c
    // RooRealVar psi2S_nx3("psi2S_nx3", "psi2S_nx3", 2.53902);//c
    // RooRealVar psi2S_sigma6("psi2S_sigma6", "psi2S_sigma6", 0.00364726);//c
    // RooRealVar psi2S_sigma7("psi2S_sigma7", "psi2S_sigma7", 0.0030555);//c
    // RooRealVar psi2S_coef3("psi2S_coef3", "psi2S_coef3", 0.0319934);//c
    // RooRealVar psi2S_prop3("psi2S_prop3", "psi2S_prop3", 0.496629);//c
    RooCBShape psi2S_gauss4("psi2S_gauss4", "psi2S_gauss4", psi2S_ctau, psi2S_mu3, psi2S_sigma6, psi2S_alpha3, psi2S_nx3);
    RooGExpModel psi2S_expGau3("psi2S_expGau3", "psi2S_expGau3", psi2S_ctau, psi2S_sigma7, psi2S_coef3, false, RooGExpModel::Type::Flipped);
    RooAddPdf psi2SCtauCombBkg("psi2SCtauCombBkg", "psi2SCtauCombBkg", RooArgList(psi2S_gauss4, psi2S_expGau3), psi2S_prop3);
    // Form 4-dim p.d.f.
    //// J/psi Mass + psi(2S) Mass
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig, psi2SMassSig);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig, psi2SMassComb);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb, psi2SMassSig);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb, psi2SMassComb);
    //// J/psi Ctau + psi(2S) Ctau
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig, psi2SCtauSig);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig, psi2SCtauBkg);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg, psi2SCtauSig);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg, psi2SCtauBkg);
    RooProdPdf pdf_ctau_SigBkg("pdf_ctau_SigBkg", "pdf_ctau_SigBkg", JpsiCtauCombSig, psi2SCtauCombBkg);
    RooProdPdf pdf_ctau_BkgSig("pdf_ctau_BkgSig", "pdf_ctau_BkgSig", JpsiCtauCombBkg, psi2SCtauCombSig);
    RooProdPdf pdf_ctau_BkgBkg("pdf_ctau_BkgBkg", "pdf_ctau_BkgBkg", JpsiCtauCombBkg, psi2SCtauCombBkg);
    //// Combine Mass and Ctau
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_ctau_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_ctau_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_ctau_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_ctau_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_ctau_SigBkg);
    RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_ctau_BkgSig);
    RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_ctau_BkgBkg);
    double N = 60000 * scale;
    Int_t BinNum = 200;
    // RooRealVar n_P_P("n_P_P", "n_P_P", 1472.15 * scale, 0, N);
    // RooRealVar n_P_NP("n_P_NP", "n_P_NP", 376.853 * scale, 0, N);
    // RooRealVar n_NP_P("n_NP_P", "n_NP_P", 219.528 * scale, 0, N);
    // RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 3066.57 * scale, 0, N);
    // RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 15554.5 * scale, 0, N);
    // RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 404.173 * scale, 0, N);
    // RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1807.7 * scale, 0, N);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1e2, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 1e3, 1e2, N);
    RooRealVar n_NP_P("n_NP_P", "n_NP_P", 1e3, 1e2, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    RooRealVar n_Comb_Sig("n_Comb_Sig", "n_Comb_Sig", 1e2, 1, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);

    // RooExtendPdf dis_P_P("dis_P_P", "dis_P_P", pdf_P_P, n_P_P);
    // RooExtendPdf dis_P_NP("dis_P_NP", "dis_P_NP", pdf_P_NP, n_P_NP);
    // RooExtendPdf dis_NP_P("dis_NP_P", "dis_NP_P", pdf_NP_P, n_NP_P);
    // RooExtendPdf dis_NP_NP("dis_NP_NP", "dis_NP_NP", pdf_NP_NP, n_NP_NP);
    // RooExtendPdf dis_Sig_Comb("dis_Sig_Comb", "dis_Sig_Comb", pdf_Sig_Comb, n_Sig_Comb);
    // RooExtendPdf dis_Comb_Sig("dis_Comb_Sig", "dis_Comb_Sig", pdf_Comb_Sig, n_Comb_Sig);
    // RooExtendPdf dis_Comb_Comb("dis_Comb_Comb", "dis_Comb_Comb", pdf_Comb_Comb, n_Comb_Comb);

    // RooAddPdf pdf_all("pdf_all", "pdf_all",
    //     RooArgList(dis_P_P, dis_P_NP, dis_NP_P, dis_NP_NP, dis_Sig_Comb, dis_Comb_Sig, dis_Comb_Comb)
    // );
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),
        RooArgList(n_P_P, n_P_NP, n_NP_P, n_NP_NP, n_Sig_Comb, n_Comb_Sig, n_Comb_Comb)
    );
    // Generate model: P+NP, combinatorial
    // srand(time(0));
    // RooRandom::randomGenerator()->SetSeed(rand());
    // RooDataSet *d_P_NP = dis_P_NP.generate(variables, Extended(kTRUE));
    // mix->append(*d_P_NP);
    // RooDataSet *d_NP_P = dis_NP_P.generate(variables, Extended(kTRUE));
    // mix->append(*d_NP_P);
    // RooDataSet *d_Sig_Comb = dis_Sig_Comb.generate(variables, Extended(kTRUE));
    // mix->append(*d_Sig_Comb);
    // RooDataSet *d_Comb_Sig = dis_Comb_Sig.generate(variables, Extended(kTRUE));
    // mix->append(*d_Comb_Sig);
    // RooDataSet *d_Comb_Comb = dis_Comb_Comb.generate(variables, Extended(kTRUE));
    // mix->append(*d_Comb_Comb);
    // 4D fitting
    RooFitResult *res;
    while(true) {
        res = pdf_all.fitTo(*mix->reduce(RooArgSet(Jpsi_mass, psi2S_mass, Jpsi_ctau, psi2S_ctau)), Save());
        if(!res->status()) break;
    }
    // Draw data point and p.d.f. curve
    TCanvas *canvas7 = new TCanvas("canvas7", "canvas7", 3000, 2000);
    RooPlot *frame7 = Jpsi_mass.frame(RooFit::Title("Jpsi Mass 4D"), Bins(BinNum));
    mix->plotOn(frame7, DataError(RooAbsData::SumW2), Name("Data"));
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
    // canvas7->SaveAs("stabTest/4D_JpsiMass.pdf");
    canvas7->SaveAs("stabTest/4D_JpsiMass.png");
    TCanvas *canvas8 = new TCanvas("canvas8", "canvas8", 3000, 2000);
    RooPlot *frame8 = psi2S_mass.frame(RooFit::Title("psi2S Mass 4D"), Bins(BinNum));
    mix->plotOn(frame8, DataError(RooAbsData::SumW2), Name("Data"));
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
    // canvas8->SaveAs("stabTest/4D_psi2SMass.pdf");
    canvas8->SaveAs("stabTest/4D_psi2SMass.png");
    TCanvas *canvas9 = new TCanvas("canvas9", "canvas9", 3000, 2000);
    RooPlot *frame9 = Jpsi_ctau.frame(RooFit::Title("Jpsi Ctau 4D"), Bins(BinNum));
    mix->plotOn(frame9, DataError(RooAbsData::SumW2), Name("Data"));
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
    // canvas9->SaveAs("stabTest/4D_JpsiCtau.pdf");
    canvas9->SaveAs("stabTest/4D_JpsiCtau.png");
    TCanvas *canvas0 = new TCanvas("canvas0", "canvas0", 3000, 2000);
    RooPlot *frame0 = psi2S_ctau.frame(RooFit::Title("psi2S Ctau 4D"), Bins(BinNum));
    mix->plotOn(frame0, DataError(RooAbsData::SumW2), Name("Data"));
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
    // canvas0->SaveAs("stabTest/4D_psi2SCtau.pdf");
    canvas0->SaveAs("stabTest/4D_psi2SCtau.png");
    pdf_all.getVariables()->Print("v");
    cout<<"Status: "<<res->status()<<endl;
    return;
}