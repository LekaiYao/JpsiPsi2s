#include "TCanvas.h"
#include "TH1D.h"

using namespace std;

void test() {
    int nDY = 6;
    double xDY[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0};
    TH1D *h = new TH1D("h", "h", nDY, xDY);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    for(int i = 0; i < nDY; i++) h->Fill(xDY[i] + 0.1);
    h->Draw();
    canvas->SaveAs("test.png");
    return;
}