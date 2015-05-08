//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 22 14:36:46 2014 by ROOT version 5.34/04
// from TTree Gen4Taus/Gen4Taus
// found on file: tautau_NMSSM_H2ToH1H1_H1ToTauTau_mH2-125_mH1-8_8TeV.root
//////////////////////////////////////////////////////////

#ifndef HiggsAnalysis_h
#define HiggsAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class HiggsAnalysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         H2pt;
   Float_t         H2px;
   Float_t         H2py;
   Float_t         H2pz;
   Float_t         H2phi;
   Float_t         H2eta;
   Float_t         H2mass;
   Float_t         PtTau1;
   Float_t         PxTau1;
   Float_t         PyTau1;
   Float_t         PzTau1;
   Float_t         PhiTau1;
   Float_t         EtaTau1;
   Int_t           Tau1FinalStates;
   Int_t           FinalStateIdTau1[4];   //[Tau1FinalStates]
   Float_t         PtTau2;
   Float_t         PxTau2;
   Float_t         PyTau2;
   Float_t         PzTau2;
   Float_t         PhiTau2;
   Float_t         EtaTau2;
   Int_t           Tau2FinalStates;
   Int_t           FinalStateIdTau2[1];   //[Tau2FinalStates]
   Float_t         PtTau3;
   Float_t         PxTau3;
   Float_t         PyTau3;
   Float_t         PzTau3;
   Float_t         PhiTau3;
   Float_t         EtaTau3;
   Int_t           Tau3FinalStates;
   Int_t           FinalStateIdTau3[4];   //[Tau3FinalStates]
   Float_t         PtTau4;
   Float_t         PxTau4;
   Float_t         PyTau4;
   Float_t         PzTau4;
   Float_t         PhiTau4;
   Float_t         EtaTau4;
   Int_t           Tau4FinalStates;
   Int_t           FinalStateIdTau4[4];   //[Tau4FinalStates]
   Float_t         PtFirstH1;
   Float_t         PxFirstH1;
   Float_t         PyFirstH1;
   Float_t         PzFirstH1;
   Float_t         PhiFirstH1;
   Float_t         EtaFirstH1;
   Float_t         MassFirstH1;
   Float_t         PtSecondH1;
   Float_t         PxSecondH1;
   Float_t         PySecondH1;
   Float_t         PzSecondH1;
   Float_t         PhiSecondH1;
   Float_t         EtaSecondH1;
   Float_t         MassSecondH1;
   Int_t           NPartFirstH1;
   Int_t           IdPartFirstH1[19];   //[NPartFirstH1]
   Int_t           QPartFirstH1[19];   //[NPartFirstH1]
   Float_t         PxPartFirstH1[19];   //[NPartFirstH1]
   Float_t         PyPartFirstH1[19];   //[NPartFirstH1]
   Float_t         PzPartFirstH1[19];   //[NPartFirstH1]
   Float_t         PtPartFirstH1[19];   //[NPartFirstH1]
   Float_t         EtaPartFirstH1[19];   //[NPartFirstH1]
   Float_t         PhiPartFirstH1[19];   //[NPartFirstH1]
   Int_t           NPartAroundMuFirstH1;
   Int_t           IdPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Int_t           QPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         PxPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         PyPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         PzPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         PtPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         EtaPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Float_t         PhiPartAroundMuFirstH1[15];   //[NPartAroundMuFirstH1]
   Int_t           NPartSecondH1;
   Int_t           IdPartSecondH1[19];   //[NPartSecondH1]
   Int_t           QPartSecondH1[19];   //[NPartSecondH1]
   Float_t         PxPartSecondH1[19];   //[NPartSecondH1]
   Float_t         PyPartSecondH1[19];   //[NPartSecondH1]
   Float_t         PzPartSecondH1[19];   //[NPartSecondH1]
   Float_t         PtPartSecondH1[19];   //[NPartSecondH1]
   Float_t         EtaPartSecondH1[19];   //[NPartSecondH1]
   Float_t         PhiPartSecondH1[19];   //[NPartSecondH1]
   Int_t           NPartAroundMuSecondH1;
   Int_t           IdPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Int_t           QPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         PxPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         PyPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         PzPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         PtPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         EtaPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]
   Float_t         PhiPartAroundMuSecondH1[13];   //[NPartAroundMuSecondH1]

   // List of branches
   TBranch        *b_H2pt;   //!
   TBranch        *b_H2px;   //!
   TBranch        *b_H2py;   //!
   TBranch        *b_H2pz;   //!
   TBranch        *b_H2phi;   //!
   TBranch        *b_H2eta;   //!
   TBranch        *b_H2mass;   //!
   TBranch        *b_PtTau1;   //!
   TBranch        *b_PxTau1;   //!
   TBranch        *b_PyTau1;   //!
   TBranch        *b_PzTau1;   //!
   TBranch        *b_PhiTau1;   //!
   TBranch        *b_EtaTau1;   //!
   TBranch        *b_Tau1FinalStates;   //!
   TBranch        *b_FinalStateIdTau1;   //!
   TBranch        *b_PtTau2;   //!
   TBranch        *b_PxTau2;   //!
   TBranch        *b_PyTau2;   //!
   TBranch        *b_PzTau2;   //!
   TBranch        *b_PhiTau2;   //!
   TBranch        *b_EtaTau2;   //!
   TBranch        *b_Tau2FinalStates;   //!
   TBranch        *b_FinalStateIdTau2;   //!
   TBranch        *b_PtTau3;   //!
   TBranch        *b_PxTau3;   //!
   TBranch        *b_PyTau3;   //!
   TBranch        *b_PzTau3;   //!
   TBranch        *b_PhiTau3;   //!
   TBranch        *b_EtaTau3;   //!
   TBranch        *b_Tau3FinalStates;   //!
   TBranch        *b_FinalStateIdTau3;   //!
   TBranch        *b_PtTau4;   //!
   TBranch        *b_PxTau4;   //!
   TBranch        *b_PyTau4;   //!
   TBranch        *b_PzTau4;   //!
   TBranch        *b_PhiTau4;   //!
   TBranch        *b_EtaTau4;   //!
   TBranch        *b_Tau4FinalStates;   //!
   TBranch        *b_FinalStateIdTau4;   //!
   TBranch        *b_PtFirstH1;   //!
   TBranch        *b_PxFirstH1;   //!
   TBranch        *b_PyFirstH1;   //!
   TBranch        *b_PzFirstH1;   //!
   TBranch        *b_PhiFirstH1;   //!
   TBranch        *b_EtaFirstH1;   //!
   TBranch        *b_MassFirstH1;   //!
   TBranch        *b_PtSecondH1;   //!
   TBranch        *b_PxSecondH1;   //!
   TBranch        *b_PySecondH1;   //!
   TBranch        *b_PzSecondH1;   //!
   TBranch        *b_PhiSecondH1;   //!
   TBranch        *b_EtaSecondH1;   //!
   TBranch        *b_MassSecondH1;   //!
   TBranch        *b_NPartFirstH1;   //!
   TBranch        *b_IdPartFirstH1;   //!
   TBranch        *b_QPartFirstH1;   //!
   TBranch        *b_PxPartFirstH1;   //!
   TBranch        *b_PyPartFirstH1;   //!
   TBranch        *b_PzPartFirstH1;   //!
   TBranch        *b_PtPartFirstH1;   //!
   TBranch        *b_EtaPartFirstH1;   //!
   TBranch        *b_PhiPartFirstH1;   //!
   TBranch        *b_NPartAroundMuFirstH1;   //!
   TBranch        *b_IdPartAroundMuFirstH1;   //!
   TBranch        *b_QPartAroundMuFirstH1;   //!
   TBranch        *b_PxPartAroundMuFirstH1;   //!
   TBranch        *b_PyPartAroundMuFirstH1;   //!
   TBranch        *b_PzPartAroundMuFirstH1;   //!
   TBranch        *b_PtPartAroundMuFirstH1;   //!
   TBranch        *b_EtaPartAroundMuFirstH1;   //!
   TBranch        *b_PhiPartAroundMuFirstH1;   //!
   TBranch        *b_NPartSecondH1;   //!
   TBranch        *b_IdPartSecondH1;   //!
   TBranch        *b_QPartSecondH1;   //!
   TBranch        *b_PxPartSecondH1;   //!
   TBranch        *b_PyPartSecondH1;   //!
   TBranch        *b_PzPartSecondH1;   //!
   TBranch        *b_PtPartSecondH1;   //!
   TBranch        *b_EtaPartSecondH1;   //!
   TBranch        *b_PhiPartSecondH1;   //!
   TBranch        *b_NPartAroundMuSecondH1;   //!
   TBranch        *b_IdPartAroundMuSecondH1;   //!
   TBranch        *b_QPartAroundMuSecondH1;   //!
   TBranch        *b_PxPartAroundMuSecondH1;   //!
   TBranch        *b_PyPartAroundMuSecondH1;   //!
   TBranch        *b_PzPartAroundMuSecondH1;   //!
   TBranch        *b_PtPartAroundMuSecondH1;   //!
   TBranch        *b_EtaPartAroundMuSecondH1;   //!
   TBranch        *b_PhiPartAroundMuSecondH1;   //!

   HiggsAnalysis(TTree *tree=0);
   virtual ~HiggsAnalysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef HiggsAnalysis_cxx
HiggsAnalysis::HiggsAnalysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("tautau_NMSSM_H2ToH1H1_H1ToTau_allfiles.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("tautau_NMSSM_H2ToH1H1_H1ToTau_allfiles.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("tautau_NMSSM_H2ToH1H1_H1ToTau_allfiles.root:/generatorNMSSM");
      dir->GetObject("Gen4Taus",tree);

   }
   Init(tree);
}

HiggsAnalysis::~HiggsAnalysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HiggsAnalysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HiggsAnalysis::LoadTree(Long64_t entry)
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

void HiggsAnalysis::Init(TTree *tree)
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

   fChain->SetBranchAddress("H2pt", &H2pt, &b_H2pt);
   fChain->SetBranchAddress("H2px", &H2px, &b_H2px);
   fChain->SetBranchAddress("H2py", &H2py, &b_H2py);
   fChain->SetBranchAddress("H2pz", &H2pz, &b_H2pz);
   fChain->SetBranchAddress("H2phi", &H2phi, &b_H2phi);
   fChain->SetBranchAddress("H2eta", &H2eta, &b_H2eta);
   fChain->SetBranchAddress("H2mass", &H2mass, &b_H2mass);
   fChain->SetBranchAddress("PtTau1", &PtTau1, &b_PtTau1);
   fChain->SetBranchAddress("PxTau1", &PxTau1, &b_PxTau1);
   fChain->SetBranchAddress("PyTau1", &PyTau1, &b_PyTau1);
   fChain->SetBranchAddress("PzTau1", &PzTau1, &b_PzTau1);
   fChain->SetBranchAddress("PhiTau1", &PhiTau1, &b_PhiTau1);
   fChain->SetBranchAddress("EtaTau1", &EtaTau1, &b_EtaTau1);
   fChain->SetBranchAddress("Tau1FinalStates", &Tau1FinalStates, &b_Tau1FinalStates);
   fChain->SetBranchAddress("FinalStateIdTau1", FinalStateIdTau1, &b_FinalStateIdTau1);
   fChain->SetBranchAddress("PtTau2", &PtTau2, &b_PtTau2);
   fChain->SetBranchAddress("PxTau2", &PxTau2, &b_PxTau2);
   fChain->SetBranchAddress("PyTau2", &PyTau2, &b_PyTau2);
   fChain->SetBranchAddress("PzTau2", &PzTau2, &b_PzTau2);
   fChain->SetBranchAddress("PhiTau2", &PhiTau2, &b_PhiTau2);
   fChain->SetBranchAddress("EtaTau2", &EtaTau2, &b_EtaTau2);
   fChain->SetBranchAddress("Tau2FinalStates", &Tau2FinalStates, &b_Tau2FinalStates);
   fChain->SetBranchAddress("FinalStateIdTau2", &FinalStateIdTau2, &b_FinalStateIdTau2);
   fChain->SetBranchAddress("PtTau3", &PtTau3, &b_PtTau3);
   fChain->SetBranchAddress("PxTau3", &PxTau3, &b_PxTau3);
   fChain->SetBranchAddress("PyTau3", &PyTau3, &b_PyTau3);
   fChain->SetBranchAddress("PzTau3", &PzTau3, &b_PzTau3);
   fChain->SetBranchAddress("PhiTau3", &PhiTau3, &b_PhiTau3);
   fChain->SetBranchAddress("EtaTau3", &EtaTau3, &b_EtaTau3);
   fChain->SetBranchAddress("Tau3FinalStates", &Tau3FinalStates, &b_Tau3FinalStates);
   fChain->SetBranchAddress("FinalStateIdTau3", FinalStateIdTau3, &b_FinalStateIdTau3);
   fChain->SetBranchAddress("PtTau4", &PtTau4, &b_PtTau4);
   fChain->SetBranchAddress("PxTau4", &PxTau4, &b_PxTau4);
   fChain->SetBranchAddress("PyTau4", &PyTau4, &b_PyTau4);
   fChain->SetBranchAddress("PzTau4", &PzTau4, &b_PzTau4);
   fChain->SetBranchAddress("PhiTau4", &PhiTau4, &b_PhiTau4);
   fChain->SetBranchAddress("EtaTau4", &EtaTau4, &b_EtaTau4);
   fChain->SetBranchAddress("Tau4FinalStates", &Tau4FinalStates, &b_Tau4FinalStates);
   fChain->SetBranchAddress("FinalStateIdTau4", FinalStateIdTau4, &b_FinalStateIdTau4);
   fChain->SetBranchAddress("PtFirstH1", &PtFirstH1, &b_PtFirstH1);
   fChain->SetBranchAddress("PxFirstH1", &PxFirstH1, &b_PxFirstH1);
   fChain->SetBranchAddress("PyFirstH1", &PyFirstH1, &b_PyFirstH1);
   fChain->SetBranchAddress("PzFirstH1", &PzFirstH1, &b_PzFirstH1);
   fChain->SetBranchAddress("PhiFirstH1", &PhiFirstH1, &b_PhiFirstH1);
   fChain->SetBranchAddress("EtaFirstH1", &EtaFirstH1, &b_EtaFirstH1);
   fChain->SetBranchAddress("MassFirstH1", &MassFirstH1, &b_MassFirstH1);
   fChain->SetBranchAddress("PtSecondH1", &PtSecondH1, &b_PtSecondH1);
   fChain->SetBranchAddress("PxSecondH1", &PxSecondH1, &b_PxSecondH1);
   fChain->SetBranchAddress("PySecondH1", &PySecondH1, &b_PySecondH1);
   fChain->SetBranchAddress("PzSecondH1", &PzSecondH1, &b_PzSecondH1);
   fChain->SetBranchAddress("PhiSecondH1", &PhiSecondH1, &b_PhiSecondH1);
   fChain->SetBranchAddress("EtaSecondH1", &EtaSecondH1, &b_EtaSecondH1);
   fChain->SetBranchAddress("MassSecondH1", &MassSecondH1, &b_MassSecondH1);
   fChain->SetBranchAddress("NPartFirstH1", &NPartFirstH1, &b_NPartFirstH1);
   fChain->SetBranchAddress("IdPartFirstH1", IdPartFirstH1, &b_IdPartFirstH1);
   fChain->SetBranchAddress("QPartFirstH1", QPartFirstH1, &b_QPartFirstH1);
   fChain->SetBranchAddress("PxPartFirstH1", PxPartFirstH1, &b_PxPartFirstH1);
   fChain->SetBranchAddress("PyPartFirstH1", PyPartFirstH1, &b_PyPartFirstH1);
   fChain->SetBranchAddress("PzPartFirstH1", PzPartFirstH1, &b_PzPartFirstH1);
   fChain->SetBranchAddress("PtPartFirstH1", PtPartFirstH1, &b_PtPartFirstH1);
   fChain->SetBranchAddress("EtaPartFirstH1", EtaPartFirstH1, &b_EtaPartFirstH1);
   fChain->SetBranchAddress("PhiPartFirstH1", PhiPartFirstH1, &b_PhiPartFirstH1);
   fChain->SetBranchAddress("NPartAroundMuFirstH1", &NPartAroundMuFirstH1, &b_NPartAroundMuFirstH1);
   fChain->SetBranchAddress("IdPartAroundMuFirstH1", IdPartAroundMuFirstH1, &b_IdPartAroundMuFirstH1);
   fChain->SetBranchAddress("QPartAroundMuFirstH1", QPartAroundMuFirstH1, &b_QPartAroundMuFirstH1);
   fChain->SetBranchAddress("PxPartAroundMuFirstH1", PxPartAroundMuFirstH1, &b_PxPartAroundMuFirstH1);
   fChain->SetBranchAddress("PyPartAroundMuFirstH1", PyPartAroundMuFirstH1, &b_PyPartAroundMuFirstH1);
   fChain->SetBranchAddress("PzPartAroundMuFirstH1", PzPartAroundMuFirstH1, &b_PzPartAroundMuFirstH1);
   fChain->SetBranchAddress("PtPartAroundMuFirstH1", PtPartAroundMuFirstH1, &b_PtPartAroundMuFirstH1);
   fChain->SetBranchAddress("EtaPartAroundMuFirstH1", EtaPartAroundMuFirstH1, &b_EtaPartAroundMuFirstH1);
   fChain->SetBranchAddress("PhiPartAroundMuFirstH1", PhiPartAroundMuFirstH1, &b_PhiPartAroundMuFirstH1);
   fChain->SetBranchAddress("NPartSecondH1", &NPartSecondH1, &b_NPartSecondH1);
   fChain->SetBranchAddress("IdPartSecondH1", IdPartSecondH1, &b_IdPartSecondH1);
   fChain->SetBranchAddress("QPartSecondH1", QPartSecondH1, &b_QPartSecondH1);
   fChain->SetBranchAddress("PxPartSecondH1", PxPartSecondH1, &b_PxPartSecondH1);
   fChain->SetBranchAddress("PyPartSecondH1", PyPartSecondH1, &b_PyPartSecondH1);
   fChain->SetBranchAddress("PzPartSecondH1", PzPartSecondH1, &b_PzPartSecondH1);
   fChain->SetBranchAddress("PtPartSecondH1", PtPartSecondH1, &b_PtPartSecondH1);
   fChain->SetBranchAddress("EtaPartSecondH1", EtaPartSecondH1, &b_EtaPartSecondH1);
   fChain->SetBranchAddress("PhiPartSecondH1", PhiPartSecondH1, &b_PhiPartSecondH1);
   fChain->SetBranchAddress("NPartAroundMuSecondH1", &NPartAroundMuSecondH1, &b_NPartAroundMuSecondH1);
   fChain->SetBranchAddress("IdPartAroundMuSecondH1", IdPartAroundMuSecondH1, &b_IdPartAroundMuSecondH1);
   fChain->SetBranchAddress("QPartAroundMuSecondH1", QPartAroundMuSecondH1, &b_QPartAroundMuSecondH1);
   fChain->SetBranchAddress("PxPartAroundMuSecondH1", PxPartAroundMuSecondH1, &b_PxPartAroundMuSecondH1);
   fChain->SetBranchAddress("PyPartAroundMuSecondH1", PyPartAroundMuSecondH1, &b_PyPartAroundMuSecondH1);
   fChain->SetBranchAddress("PzPartAroundMuSecondH1", PzPartAroundMuSecondH1, &b_PzPartAroundMuSecondH1);
   fChain->SetBranchAddress("PtPartAroundMuSecondH1", PtPartAroundMuSecondH1, &b_PtPartAroundMuSecondH1);
   fChain->SetBranchAddress("EtaPartAroundMuSecondH1", EtaPartAroundMuSecondH1, &b_EtaPartAroundMuSecondH1);
   fChain->SetBranchAddress("PhiPartAroundMuSecondH1", PhiPartAroundMuSecondH1, &b_PhiPartAroundMuSecondH1);
   Notify();
}

Bool_t HiggsAnalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HiggsAnalysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HiggsAnalysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef HiggsAnalysis_cxx
