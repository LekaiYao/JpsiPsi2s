// Author: Shiyang CHEN
// Description: 
// Implementation: 
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
using namespace std;

constexpr long double operator"" _PI(long double f) {
    return 3.14159265359 * f;
}

void diff_plot() {
    int nDY = 6;
    double xDY[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0, 4.0};
    TH1D *hDY1 = new TH1D("DY Xsection 1", "DY Xsection 1", nDY, xDY), *hDY2 = new TH1D("DY Xsection 2", "DY Xsection 2", nDY, xDY);
    double xSec_dY[] = {55.05, 26.73, 13.33, 12.72, 7.18, 2.23}, xSecErr_dY[] = {8.89, 5.05, 2.94, 2.82, 1.88, 0.6}, xSecSta_dY[] = {2.84, 2.11, 1.48, 1.42, 1.09, 0.35};
    for(int i = 0; i < nDY; i++) {
        hDY1->Fill(xDY[i] + 0.1, xSec_dY[i]);
        hDY1->SetBinError(i + 1, xSecSta_dY[i]);
        hDY2->Fill(xDY[i] + 0.1, xSec_dY[i]);
        hDY2->SetBinError(i + 1, xSecErr_dY[i]);
    }
    int nDPhi = 5;
    double xDPhi[] = {0, 0.2_PI, 0.4_PI, 0.6_PI, 0.8_PI, 1.0_PI};
    TH1D *hDPhi1 = new TH1D("DPhi Xsection 1", "DPhi Xsection 1", nDPhi, xDPhi), *hDPhi2 = new TH1D("DPhi Xsection 2", "DPhi Xsection 2", nDPhi, xDPhi);
    double xSec_dPhi[] = {34.75, 9.89, 16.09, 8.13, 15.51}, xSecErr_dPhi[] = {5.63, 2.13, 2.94, 1.94, 3.07}, xSecSta_dPhi[] = {1.81, 1.04, 1.17, 1.05, 1.37};
    for(int i = 0; i < nDPhi; i++) {
        hDPhi1->Fill(xDPhi[i] + 0.1, xSec_dPhi[i]);
        hDPhi1->SetBinError(i + 1, xSecSta_dPhi[i]);
        hDPhi2->Fill(xDPhi[i] + 0.1, xSec_dPhi[i]);
        hDPhi2->SetBinError(i + 1, xSecErr_dPhi[i]);
    }
    int nMinv = 5;
    double xMinv[] = {7.5, 17.5, 27.5, 37.5, 47.5, 107.5};
    TH1D *hMinv1 = new TH1D("Minv Xsection 1", "Minv Xsection 1", nMinv, xMinv), *hMinv2 = new TH1D("Minv Xsection 2", "Minv Xsection 2", nMinv, xMinv);
    double xSec_Minv[] = {2.77, 0.87, 0.87, 0.23, 0.043}, xSecErr_Minv[] = {0.3, 0.14, 0.14, 0.08, 0.01}, xSecSta_Minv[] = {0.12, 0.08, 0.08, 0.06, 0.0083};
    for(int i = 0; i < nMinv; i++) {
        hMinv1->Fill(xMinv[i] + 1, xSec_Minv[i]);
        hMinv1->SetBinError(i + 1, xSecSta_Minv[i]);
        hMinv2->Fill(xMinv[i] + 1, xSec_Minv[i]);
        hMinv2->SetBinError(i + 1, xSecErr_Minv[i]);
    }
    int nSY = 5;
    double xSY[] = {0, 0.4, 0.8, 1.2, 1.6, 2.0};
    TH1D *hSY1 = new TH1D("SY Xsection 1", "SY Xsection 1", nSY, xSY), *hSY2 = new TH1D("SY Xsection 2", "SY Xsection 2", nSY, xSY);
    double xSec_sY[] = {41.98, 26.64, 31.55, 19.91, 7.2}, xSecErr_sY[] = {7.17, 5.19, 5.74, 3.95, 1.95}, xSecSta_sY[] = {2.56, 2.26, 2.27, 1.76, 1.16};
    for(int i = 0; i < nSY; i++) {
        hSY1->Fill(xSY[i] + 0.1, xSec_sY[i]);
        hSY1->SetBinError(i + 1, xSecSta_sY[i]);
        hSY2->Fill(xSY[i] + 0.1, xSec_sY[i]);
        hSY2->SetBinError(i + 1, xSecErr_sY[i]);
    }
    int nSPt = 6;
    double xSPt[] = {0, 8, 16, 24, 32, 40, 80};
    TH1D *hSPt1 = new TH1D("SPt Xsection 1", "SPt Xsection 1", nSPt, xSPt), *hSPt2 = new TH1D("SPt Xsection 2", "SPt Xsection 2", nSPt, xSPt);
    double xSec_sPt[] = {1.11, 0.94, 1.71, 1.60, 0.61, 0.086}, xSecErr_sPt[] = {0.21, 0.2, 0.31, 0.3, 0.14, 0.02}, xSecSta_sPt[] = {0.09, 0.1, 0.12, 0.12, 0.07, 0.01};
    for(int i = 0; i < nSPt; i++) {
        hSPt1->Fill(xSPt[i] + 1, xSec_sPt[i]);
        hSPt1->SetBinError(i + 1, xSecSta_sPt[i]);
        hSPt2->Fill(xSPt[i] + 1, xSec_sPt[i]);
        hSPt2->SetBinError(i + 1, xSecErr_sPt[i]);
    }

    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    hDY1->SetMarkerSize(0.);
    hDY1->SetFillStyle(3003);
    hDY1->SetMarkerColor(kRed);
    hDY1->SetFillColor(kRed);
    hDY1->GetYaxis()->SetRangeUser(0, 1.2 * hDY1->GetMaximum());
    hDY1->GetXaxis()->SetTitle("#Deltay");
    hDY1->GetYaxis()->SetTitle("d#sigma/d#Deltay");
    hDY1->SetStats(0);
    hDY1->Draw("e2");
    hDY2->SetStats(0);
    hDY2->Draw("same");
    canvas->SaveAs("fits/DiffXsec_DY.png");
    delete canvas;
    canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    hDPhi1->SetMarkerSize(0.);
    hDPhi1->SetFillStyle(3003);
    hDPhi1->SetMarkerColor(kRed);
    hDPhi1->SetFillColor(kRed);
    hDPhi1->GetYaxis()->SetRangeUser(0, 1.2 * hDPhi1->GetMaximum());
    hDPhi1->GetXaxis()->SetTitle("#Delta#phi");
    hDPhi1->GetYaxis()->SetTitle("d#sigma/d#Delta#phi");
    hDPhi1->SetStats(0);
    hDPhi1->Draw("e2");
    hDPhi2->SetStats(0);
    hDPhi2->Draw("same");
    canvas->SaveAs("fits/DiffXsec_DPhi.png");
    delete canvas;
    canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    hMinv1->SetMarkerSize(0.);
    hMinv1->SetFillStyle(3003);
    hMinv1->SetMarkerColor(kRed);
    hMinv1->SetFillColor(kRed);
    hMinv1->GetYaxis()->SetRangeUser(0, 1.2 * hMinv1->GetMaximum());
    hMinv1->GetXaxis()->SetTitle("M_{inv}");
    hMinv1->GetYaxis()->SetTitle("d#sigma/dM_{inv}");
    hMinv1->SetStats(0);
    hMinv1->Draw("e2");
    hMinv2->SetStats(0);
    hMinv2->Draw("same");
    canvas->SaveAs("fits/DiffXsec_Minv.png");
    delete canvas;
    canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    hSY1->SetMarkerSize(0.);
    hSY1->SetFillStyle(3003);
    hSY1->SetMarkerColor(kRed);
    hSY1->SetFillColor(kRed);
    hSY1->GetYaxis()->SetRangeUser(0, 1.2 * hSY1->GetMaximum());
    hSY1->GetXaxis()->SetTitle("y");
    hSY1->GetYaxis()->SetTitle("d#sigma/dy");
    hSY1->SetStats(0);
    hSY1->Draw("e2");
    hSY2->SetStats(0);
    hSY2->Draw("same");
    canvas->SaveAs("fits/DiffXsec_SY.png");
    delete canvas;
    canvas = new TCanvas("canvas", "canvas", 1000, 1000);
    hSPt1->SetMarkerSize(0.);
    hSPt1->SetFillStyle(3003);
    hSPt1->SetMarkerColor(kRed);
    hSPt1->SetFillColor(kRed);
    hSPt1->GetYaxis()->SetRangeUser(0, 1.2 * hSPt1->GetMaximum());
    hSPt1->GetXaxis()->SetTitle("p_{T}");
    hSPt1->GetYaxis()->SetTitle("d#sigma/dp_{T}");
    hSPt1->SetStats(0);
    hSPt1->Draw("e2");
    hSPt2->SetStats(0);
    hSPt2->Draw("same");
    canvas->SaveAs("fits/DiffXsec_SPt.png");
    delete canvas;
}