#include "TCanvas.h"
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
#include <iostream>
#include <fstream>

using namespace std;
using namespace RooFit;

void test() {
    ifstream inFile("summary.txt");
    double x;
    TTree *t = new TTree("tree", "tree");
    t->Branch("x", &x);
    if(!inFile.is_open()) return;
    while(inFile>>x) {
        if(x == -1) continue;
        t->Fill();
    }
    inFile.close();
    
    RooRealVar n_P_P("x", "x", -4, 4);
    RooDataSet *data = new RooDataSet("data", "data", t, RooArgSet(n_P_P));

    RooRealVar mean("mean", "mean", 0, -1, 1);
    RooRealVar sigma("sigma", "sigma", 1, 0.2, 5);
    RooGaussian gauss("gauss", "gauss", n_P_P, mean, sigma);

    RooRealVar n("n", "n", 1e3, 1, 1e4);
    RooAddPdf pdf("pdf", "pdf", RooArgList(gauss), RooArgList(n));
    // RooDataSet *data = pdf.generate(x, Extended(kTRUE));
    RooFitResult *res;
    while(true) {
        res = pdf.fitTo(*data, Save());
        if(!res->status()) break;
    }
    
    TCanvas *c = new TCanvas("c", "c", 2000, 1800);
    RooPlot *f = n_P_P.frame(RooFit::Title("n_P_P Pull Distribution;N_{corr}"), Bins(40));
    data->plotOn(f, DataError(RooAbsData::SumW2), Name("Data"));//
    pdf.plotOn(f, LineColor(kBlue), LineWidth(2), Name("GS"));
    TLegend *l = new TLegend(.65, .60, .85, .85);
    l->AddEntry(f->findObject("Data"), "N_{corr} Pull", "L");
    l->AddEntry(f->findObject("GS"), "Gaussian", "L");
    f->Draw();
    l->DrawClone();
    c->SaveAs("pull.png");
    pdf.getVariables()->Print("v");
    return;
}
