// plot_weightdata_compare.C
// Compare our WeightData.root vs chensh WeightData_prime.root on fit variables.
// Usage: root -l -q -b 'plot_weightdata_compare.C'

#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include <iostream>

void plot_weightdata_compare() {
    gStyle->SetOptStat(0);

    struct VarConfig {
        const char* name;
        double xmin, xmax;
        int nbins;
        const char* title;
    };

    VarConfig vars[] = {
        {"Jpsi_mass",    2.95, 3.25, 30, "M(J/#psi) [GeV]"},
        {"psi2S_mass",   3.35, 4.05, 30, "M(#psi(2S)) [GeV]"},
        {"Jpsi_ctau",  -0.03, 0.16, 30, "c#tau(J/#psi) [cm]"},
        {"psi2S_ctau", -0.03, 0.16, 30, "c#tau(#psi(2S)) [cm]"},
    };
    const int nVars = 4;

    TFile *f_ours   = TFile::Open("WeightData.root");
    TFile *f_chensh = TFile::Open("/eos/user/c/chensh/JPsiPsi2s/Data/ULntuple18/CMSSW_10_6_20/src/NtupleAnalyzer/WeightData_prime.root");
    TTree *t_ours   = (TTree*)f_ours->Get("data");
    TTree *t_chensh = (TTree*)f_chensh->Get("data");

    cout << "Our entries:     " << t_ours->GetEntries() << endl;
    cout << "Chensh entries:  " << t_chensh->GetEntries() << endl;

    // Branch variables
    Double_t v_ours, v_chensh, w_ours, w_chensh;

    for (int iv = 0; iv < nVars; iv++) {
        VarConfig &vc = vars[iv];

        t_ours->SetBranchAddress(vc.name, &v_ours);
        t_ours->SetBranchAddress("evt_weight", &w_ours);
        t_chensh->SetBranchAddress(vc.name, &v_chensh);
        t_chensh->SetBranchAddress("evt_weight", &w_chensh);

        // Build histograms (weighted)
        TH1D *h_ours   = new TH1D(Form("h_ours_%d", iv),   Form(";%s;Weighted events", vc.title),  vc.nbins, vc.xmin, vc.xmax);
        TH1D *h_chensh = new TH1D(Form("h_chensh_%d", iv), Form(";%s;Weighted events", vc.title),  vc.nbins, vc.xmin, vc.xmax);
        TH1D *h_ours_norm   = new TH1D(Form("h_ours_n_%d", iv),   Form(";%s;dN/d%s (norm.)", vc.title, vc.title), vc.nbins, vc.xmin, vc.xmax);
        TH1D *h_chensh_norm = new TH1D(Form("h_chensh_n_%d", iv), Form(";%s;dN/d%s (norm.)", vc.title, vc.title), vc.nbins, vc.xmin, vc.xmax);

        for (int i = 0; i < t_ours->GetEntries(); i++) {
            t_ours->GetEntry(i);
            h_ours->Fill(v_ours, w_ours);
            h_ours_norm->Fill(v_ours, w_ours);
        }
        for (int i = 0; i < t_chensh->GetEntries(); i++) {
            t_chensh->GetEntry(i);
            h_chensh->Fill(v_chensh, w_chensh);
            h_chensh_norm->Fill(v_chensh, w_chensh);
        }

        // Normalise to unit area
        if (h_ours_norm->Integral() > 0)   h_ours_norm->Scale(1.0 / h_ours_norm->Integral(), "width");
        if (h_chensh_norm->Integral() > 0) h_chensh_norm->Scale(1.0 / h_chensh_norm->Integral(), "width");

        // Ratio histogram
        TH1D *h_ratio = (TH1D*)h_ours_norm->Clone(Form("ratio_%d", iv));
        h_ratio->Divide(h_chensh_norm);
        double max_dev = 0;
        for (int j = 1; j <= vc.nbins; j++) {
            double c = h_ratio->GetBinContent(j);
            if (c > 0) { double d = fabs(c - 1.0); if (d > max_dev) max_dev = d; }
        }
        if (max_dev < 0.05) max_dev = 0.5;
        double ratio_lo = 1.0 - 1.1 * max_dev;
        double ratio_hi = 1.0 + 1.1 * max_dev;

        // Canvas
        TCanvas *c = new TCanvas(Form("c_%s", vc.name), vc.title, 900, 700);
        TPad *pad_top = new TPad("ptop", "top", 0, 0.30, 1, 1);
        TPad *pad_bot = new TPad("pbot", "bot", 0, 0, 1, 0.30);
        pad_top->SetBottomMargin(0.02); pad_top->SetTopMargin(0.08); pad_top->SetRightMargin(0.05);
        pad_bot->SetTopMargin(0.02); pad_bot->SetBottomMargin(0.30); pad_bot->SetRightMargin(0.05);
        pad_top->Draw(); pad_bot->Draw();

        // Top: unnormalised overlay
        pad_top->cd();
        h_ours->SetLineColor(kBlue+1);    h_ours->SetLineWidth(2);
        h_chensh->SetLineColor(kRed+1);   h_chensh->SetLineWidth(2);
        double peak = h_ours->GetMaximum();
        if (h_chensh->GetMaximum() > peak) peak = h_chensh->GetMaximum();
        h_ours->SetMinimum(0); h_ours->SetMaximum(1.2 * peak);
        h_ours->Draw("HIST");
        h_chensh->Draw("HIST SAME");
        TLegend *leg = new TLegend(0.68, 0.78, 0.88, 0.88);
        leg->AddEntry(h_ours,   Form("Ours (%d evt)",   t_ours->GetEntries()), "l");
        leg->AddEntry(h_chensh, Form("Chensh (%d evt)", t_chensh->GetEntries()), "l");
        leg->Draw();

        // Bottom: normalised ratio
        pad_bot->cd();
        h_ratio->SetTitle(Form(";%s;Ours / Chensh", vc.title));
        h_ratio->SetLineColor(kViolet+2); h_ratio->SetLineWidth(2);
        h_ratio->SetMarkerStyle(20); h_ratio->SetMarkerSize(0.5);
        h_ratio->SetMinimum(ratio_lo); h_ratio->SetMaximum(ratio_hi);
        h_ratio->GetYaxis()->SetNdivisions(505);
        h_ratio->GetYaxis()->SetTitleSize(0.12);
        h_ratio->GetYaxis()->SetTitleOffset(0.35);
        h_ratio->GetYaxis()->SetLabelSize(0.10);
        h_ratio->GetXaxis()->SetTitleSize(0.14);
        h_ratio->GetXaxis()->SetTitleOffset(0.80);
        h_ratio->GetXaxis()->SetLabelSize(0.10);
        h_ratio->Draw("HIST E");
        TLine *line = new TLine(vc.xmin, 1.0, vc.xmax, 1.0);
        line->SetLineColor(kGray+2); line->SetLineStyle(2);
        line->Draw("SAME");

        c->SaveAs(Form("fig/compare_%s.png", vc.name));
        cout << "  Saved fig/compare_" << vc.name << ".png" << endl;
        delete c;
    }

    f_ours->Close();
    f_chensh->Close();
    cout << "Done." << endl;
}
