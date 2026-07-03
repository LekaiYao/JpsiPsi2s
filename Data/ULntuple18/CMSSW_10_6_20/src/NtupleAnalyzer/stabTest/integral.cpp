#include "TCanvas.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"

using namespace std;
using namespace RooFit;

void integral() {
    RooRealVar x("x", "x", -2, 2);
    RooRealVar y("y", "y", -2, 2);

    RooRealVar mean1("mean1", "mean1", 0);
    RooRealVar sigma1("sigma1", "sigma1", 2);
    RooGaussian gauss1("gauss1", "gauss1", x, mean1, sigma1);
    RooRealVar mean2("mean2", "mean2", 0);
    RooRealVar sigma2("sigma2", "sigma2", 3);
    RooGaussian gauss2("gauss2", "gauss2", y, mean2, sigma2);

    RooRealVar n("n", "n", 10);
    RooProdPdf GsGs("GsGs", "GsGs", gauss1, gauss2);
    RooAddPdf pdf("pdf", "pdf", RooArgList(GsGs), RooArgList(n));
    x.setRange("window", -2, 2);
    RooAbsReal* intpdf = pdf.createIntegral(x, x, "window");
    
    cout<<intpdf->getVal()<<endl;

    return;
}
