#define EventMixing_cxx
#include "EventMixing.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//2024.9.30
// Event mixing method to extract another DPS sample 
// Canidates outside the "prompt peak" are excluded
// Candidates with four muon mass smaller than 7.5 are excluded

//2024.12.23
//Candidate will now get mixed with more other candidates 

float J1Mass_EventMixing;
float J2Mass_EventMixing;

float J1Pt_EventMixing;
float J2Pt_EventMixing;

float J1y_EventMixing;
float J2y_EventMixing;

float J1Phi_EventMixing;
float J2Phi_EventMixing;

float x1_EventMixing;
float x2_EventMixing;

float DeltaY_EventMixing;
float DeltaPhi_EventMixing;

//Vector for storage
vector<float> Vector_J1Mass;
vector<float> Vector_J2Mass;

vector<float> Vector_J1Pt;
vector<float> Vector_J2Pt;

vector<float> Vector_J1y;
vector<float> Vector_J2y;

vector<float> Vector_J1Phi;
vector<float> Vector_J2Phi;

vector<float> Vector_x1;
vector<float> Vector_x2;


void EventMixing::Loop()
{
    //Output file
    TFile* file = new TFile("EventMixing_2016.root", "RECREATE");

    //Tree
    TTree EventMixingTree("EventMixingTree", "EventMixingTree");

    //Branch of the Tree
    EventMixingTree.Branch("J1Mass_EventMixing", &J1Mass_EventMixing, "J1Mass_EventMixing/F");
    EventMixingTree.Branch("J2Mass_EventMixing", &J2Mass_EventMixing, "J2Mass_EventMixing/F");

    EventMixingTree.Branch("J1Pt_EventMixing", &J1Pt_EventMixing, "J1Pt_EventMixing/F");
    EventMixingTree.Branch("J2Pt_EventMixing", &J2Pt_EventMixing, "J2Pt_EventMixing/F");

    EventMixingTree.Branch("J1Phi_EventMixing", &J1Phi_EventMixing, "J1Phi_EventMixing/F");
    EventMixingTree.Branch("J2Phi_EventMixing", &J2Phi_EventMixing, "J2Phi_EventMixing/F");

    EventMixingTree.Branch("J1y_EventMixing", &J1y_EventMixing, "J1y_EventMixing/F");
    EventMixingTree.Branch("J2y_EventMixing", &J2y_EventMixing, "J2y_EventMixing/F");

    EventMixingTree.Branch("x1_EventMixing", &x1_EventMixing, "x1_EventMixing/F");
    EventMixingTree.Branch("x2_EventMixing", &x2_EventMixing, "x2_EventMixing/F");

    EventMixingTree.Branch("DeltaY_EventMixing", &DeltaY_EventMixing, "DeltaY_EventMixing/F");
    EventMixingTree.Branch("DeltaPhi_EventMixing", &DeltaPhi_EventMixing, "DeltaPhi_EventMixing/F");

//   In a ROOT session, you can do:
//      root> .L EventMixing.C
//      root> EventMixing t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if ((x1 < -0.02) || (x1 > 0.02) || (x2 < -0.02) || (x2 > 0.02))
          continue;

      if (FourMuonMass < 7.5)
          continue;

      Vector_J1Mass.push_back(J1Mass);
      Vector_J2Mass.push_back(J2Mass);

      Vector_J1Pt.push_back(J1Pt);
      Vector_J2Pt.push_back(J2Pt);

      Vector_J1y.push_back(J1y);
      Vector_J2y.push_back(J2y);

      Vector_J1Phi.push_back(J1Phi);
      Vector_J2Phi.push_back(J2Phi);

      Vector_x1.push_back(x1);
      Vector_x2.push_back(x2);

   }

   for (int i = 0; i < Vector_J1Mass.size(); i++)
   {
       int j = 0;

       for (int k = 0; k < 5; k++)
       {
           for (;;)
           {
               j = gRandom->Integer(Vector_J1Mass.size());
               if (i == j)
                   continue;
               else
                   break;
           }

           J1Mass_EventMixing = Vector_J1Mass[i];
           J2Mass_EventMixing = Vector_J2Mass[j];

           J1Pt_EventMixing = Vector_J1Pt[i];
           J2Pt_EventMixing = Vector_J2Pt[j];

           J1y_EventMixing = Vector_J1y[i];
           J2y_EventMixing = Vector_J2y[j];

           J1Phi_EventMixing = Vector_J1Phi[i];
           J2Phi_EventMixing = Vector_J2Phi[j];

           x1_EventMixing = Vector_x1[i];
           x2_EventMixing = Vector_x2[j];

           DeltaY_EventMixing = fabs(J1y_EventMixing - J2y_EventMixing);
           DeltaPhi_EventMixing = (fabs(J1Phi_EventMixing - J2Phi_EventMixing) < TMath::Pi()) ? fabs(J1Phi_EventMixing - J2Phi_EventMixing) : (2 * TMath::Pi() - fabs(J1Phi_EventMixing - J2Phi_EventMixing));

           EventMixingTree.Fill();

       }

   }

   file->Write();
}
