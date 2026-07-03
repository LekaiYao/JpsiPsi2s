//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 30 13:34:02 2024 by ROOT version 6.32.04
// from TTree TreeAfterCut/TreeAfterCut
// found on file: root://eoshome-j.cern.ch//eos/user/j/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/Data/2016/2016.root
//////////////////////////////////////////////////////////

#ifndef EventMixing_h
#define EventMixing_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class EventMixing {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       Run;
   ULong64_t       Lumi;
   ULong64_t       Event;
   Float_t         J1Mass;
   Float_t         J2Mass;
   Float_t         J1MassErr;
   Float_t         J2MassErr;
   Float_t         J1Pt;
   Float_t         J2Pt;
   Float_t         J1y;
   Float_t         J2y;
   Float_t         J1Phi;
   Float_t         J2Phi;
   Float_t         Mu1Pt;
   Float_t         Mu2Pt;
   Float_t         Mu3Pt;
   Float_t         Mu4Pt;
   Float_t         Mu1Eta;
   Float_t         Mu2Eta;
   Float_t         Mu3Eta;
   Float_t         Mu4Eta;
   Float_t         Mu1Phi;
   Float_t         Mu2Phi;
   Float_t         Mu3Phi;
   Float_t         Mu4Phi;
   Float_t         Mu1E;
   Float_t         Mu2E;
   Float_t         Mu3E;
   Float_t         Mu4E;
   Float_t         x1;
   Float_t         x2;
   Float_t         LxyPV1;
   Float_t         LxyPV2;
   Float_t         LxyPVSig1;
   Float_t         LxyPVSig2;
   Float_t         VertexDistanceSig;
   Float_t         VertexProbability;
   Float_t         FourMuonMass;
   Float_t         FourMuonPt;
   Float_t         FourMuonY;
   Float_t         DeltaY;
   Float_t         DeltaPhi;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_J1Mass;   //!
   TBranch        *b_J2Mass;   //!
   TBranch        *b_J1MassErr;   //!
   TBranch        *b_J2MassErr;   //!
   TBranch        *b_J1Pt;   //!
   TBranch        *b_J2Pt;   //!
   TBranch        *b_J1y;   //!
   TBranch        *b_J2y;   //!
   TBranch        *b_J1Phi;   //!
   TBranch        *b_J2Phi;   //!
   TBranch        *b_Mu1Pt;   //!
   TBranch        *b_Mu2Pt;   //!
   TBranch        *b_Mu3Pt;   //!
   TBranch        *b_Mu4Pt;   //!
   TBranch        *b_Mu1Eta;   //!
   TBranch        *b_Mu2Eta;   //!
   TBranch        *b_Mu3Eta;   //!
   TBranch        *b_Mu4Eta;   //!
   TBranch        *b_Mu1Phi;   //!
   TBranch        *b_Mu2Phi;   //!
   TBranch        *b_Mu3Phi;   //!
   TBranch        *b_Mu4Phi;   //!
   TBranch        *b_Mu1E;   //!
   TBranch        *b_Mu2E;   //!
   TBranch        *b_Mu3E;   //!
   TBranch        *b_Mu4E;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_x2;   //!
   TBranch        *b_LxyPV1;   //!
   TBranch        *b_LxyPV2;   //!
   TBranch        *b_LxyPV1Sig;   //!
   TBranch        *b_LxyPV2Sig;   //!
   TBranch        *b_VertexDistanceSig;   //!
   TBranch        *b_VertexProbability;   //!
   TBranch        *b_FourMuonMass;   //!
   TBranch        *b_FourMuonPt;   //!
   TBranch        *b_FourMuonY;   //!
   TBranch        *b_DeltaY;   //!
   TBranch        *b_DeltaPhi;   //!

   EventMixing(TTree *tree=0);
   virtual ~EventMixing();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventMixing_cxx
EventMixing::EventMixing(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root://eoshome-j.cern.ch//eos/user/j/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/Data/2016/2016.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root://eoshome-j.cern.ch//eos/user/j/jinfeng/CrossSection/CMSSW_10_2_5/src/PreSelection/Feb/PreSelection/Data/2016/2016.root");
      }
      f->GetObject("TreeAfterCut",tree);

   }
   Init(tree);
}

EventMixing::~EventMixing()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventMixing::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventMixing::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventMixing::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("J1Mass", &J1Mass, &b_J1Mass);
   fChain->SetBranchAddress("J2Mass", &J2Mass, &b_J2Mass);
   fChain->SetBranchAddress("J1MassErr", &J1MassErr, &b_J1MassErr);
   fChain->SetBranchAddress("J2MassErr", &J2MassErr, &b_J2MassErr);
   fChain->SetBranchAddress("J1Pt", &J1Pt, &b_J1Pt);
   fChain->SetBranchAddress("J2Pt", &J2Pt, &b_J2Pt);
   fChain->SetBranchAddress("J1y", &J1y, &b_J1y);
   fChain->SetBranchAddress("J2y", &J2y, &b_J2y);
   fChain->SetBranchAddress("J1Phi", &J1Phi, &b_J1Phi);
   fChain->SetBranchAddress("J2Phi", &J2Phi, &b_J2Phi);
   fChain->SetBranchAddress("Mu1Pt", &Mu1Pt, &b_Mu1Pt);
   fChain->SetBranchAddress("Mu2Pt", &Mu2Pt, &b_Mu2Pt);
   fChain->SetBranchAddress("Mu3Pt", &Mu3Pt, &b_Mu3Pt);
   fChain->SetBranchAddress("Mu4Pt", &Mu4Pt, &b_Mu4Pt);
   fChain->SetBranchAddress("Mu1Eta", &Mu1Eta, &b_Mu1Eta);
   fChain->SetBranchAddress("Mu2Eta", &Mu2Eta, &b_Mu2Eta);
   fChain->SetBranchAddress("Mu3Eta", &Mu3Eta, &b_Mu3Eta);
   fChain->SetBranchAddress("Mu4Eta", &Mu4Eta, &b_Mu4Eta);
   fChain->SetBranchAddress("Mu1Phi", &Mu1Phi, &b_Mu1Phi);
   fChain->SetBranchAddress("Mu2Phi", &Mu2Phi, &b_Mu2Phi);
   fChain->SetBranchAddress("Mu3Phi", &Mu3Phi, &b_Mu3Phi);
   fChain->SetBranchAddress("Mu4Phi", &Mu4Phi, &b_Mu4Phi);
   fChain->SetBranchAddress("Mu1E", &Mu1E, &b_Mu1E);
   fChain->SetBranchAddress("Mu2E", &Mu2E, &b_Mu2E);
   fChain->SetBranchAddress("Mu3E", &Mu3E, &b_Mu3E);
   fChain->SetBranchAddress("Mu4E", &Mu4E, &b_Mu4E);
   fChain->SetBranchAddress("x1", &x1, &b_x1);
   fChain->SetBranchAddress("x2", &x2, &b_x2);
   fChain->SetBranchAddress("LxyPV1", &LxyPV1, &b_LxyPV1);
   fChain->SetBranchAddress("LxyPV2", &LxyPV2, &b_LxyPV2);
   fChain->SetBranchAddress("LxyPVSig1", &LxyPVSig1, &b_LxyPV1Sig);
   fChain->SetBranchAddress("LxyPVSig2", &LxyPVSig2, &b_LxyPV2Sig);
   fChain->SetBranchAddress("VertexDistanceSig", &VertexDistanceSig, &b_VertexDistanceSig);
   fChain->SetBranchAddress("VertexProbability", &VertexProbability, &b_VertexProbability);
   fChain->SetBranchAddress("FourMuonMass", &FourMuonMass, &b_FourMuonMass);
   fChain->SetBranchAddress("FourMuonPt", &FourMuonPt, &b_FourMuonPt);
   fChain->SetBranchAddress("FourMuonY", &FourMuonY, &b_FourMuonY);
   fChain->SetBranchAddress("DeltaY", &DeltaY, &b_DeltaY);
   fChain->SetBranchAddress("DeltaPhi", &DeltaPhi, &b_DeltaPhi);
   Notify();
}

bool EventMixing::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void EventMixing::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventMixing::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventMixing_cxx
