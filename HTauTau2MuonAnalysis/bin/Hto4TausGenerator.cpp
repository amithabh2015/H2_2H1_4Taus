#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <string> 

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRFIOFile.h"
#include "TH1F.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/EventWeight.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/RunLumiReader.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauUtils.h"

#include "TLorentzVector.h"

#include "TRandom.h"
using namespace std;
const float MuMass = 0.105658367;
const float EleMass = 0.000511;
const float PiMass = 0.13957018;

float invMass(float px1, float py1, float pz1, float m1,
	      float px2, float py2, float pz2, float m2) {

  TLorentzVector l1, l2;
  l1.SetXYZM(px1,py1,pz1,m1);
  l2.SetXYZM(px2,py2,pz2,m2);

  TLorentzVector l = l1 + l2;

  float mass = float(l.M());

  return mass;

}

int binNumber(float x, int nbins, float * bins) {

  int binN = 0;

  for (int iB=0; iB<nbins; ++iB) {
    if (x>=bins[iB]&&x<bins[iB+1]) {
      binN = iB;
      break;
    }
  }

  return binN;

}

float effBin(float x, int nbins, float * bins, float * eff) {

  int bin = binNumber(x, nbins, bins);

  return eff[bin];

}

double cosRestFrame(TLorentzVector boost, TLorentzVector vect) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double QToEta(double Q) {
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;
}

double EtaToQ(double Eta) {
  double Q = 2.0*TMath::ATan(TMath::Exp(-Eta));
  if (Q<0.0) Q += TMath::Pi();
  return Q;
}

double PtoEta(double Px, double Py, double Pz) {

  double P = TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
  double cosQ = Pz/P;
  double Q = TMath::ACos(cosQ);
  double Eta = - TMath::Log(TMath::Tan(0.5*Q));  
  return Eta;

}

double dPhiFrom2P(double Px1, double Py1,
		  double Px2, double Py2) {


  double prod = Px1*Px2 + Py1*Py2;
  double mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  double mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  double cosDPhi = prod/(mod1*mod2);
  
  return TMath::ACos(cosDPhi);

}

double deltaEta(double Px1, double Py1, double Pz1,
		double Px2, double Py2, double Pz2) {

  double eta1 = PtoEta(Px1,Py1,Pz1);
  double eta2 = PtoEta(Px2,Py2,Pz2);

  double dEta = eta1 - eta2;

  return dEta;

}

double deltaR(double Eta1, double Phi1,
	      double Eta2, double Phi2) {

  double Px1 = TMath::Cos(Phi1);
  double Py1 = TMath::Sin(Phi1);

  double Px2 = TMath::Cos(Phi2);
  double Py2 = TMath::Sin(Phi2);

  double dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  double dEta = Eta1 - Eta2;

  double dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}

double PtEtaToP(double Pt, double Eta) {

  //  double Q = EtaToQ(Eta);

  //double P = Pt/TMath::Sin(Q);
  double P = Pt*TMath::CosH(Eta);

  return P;
}
double Px(double Pt, double Phi){

  double Px=Pt*TMath::Cos(Phi);
  return Px;
}
double Py(double Pt, double Phi){

  double Py=Pt*TMath::Sin(Phi);
  return Py;
}
double Pz(double Pt, double Eta){

  double Pz=Pt*TMath::SinH(Eta);
  return Pz;
}
double InvariantMass(double energy,double Px,double Py, double Pz){

  double M_2=energy*energy-Px*Px-Py*Py-Pz*Pz;
  double M=TMath::Sqrt(M_2);
  return M;


}
double EFromPandM0(double M0,double Pt,double Eta){

  double E_2=M0*M0+PtEtaToP(Pt,Eta)*PtEtaToP(Pt,Eta);
  double E =TMath::Sqrt(E_2);
  return E;


}

const float muonMass = 0.10565837;
const float pionMass = 0.1396;
const float elecMass = 0.000511;
const float kaonMass = 0.494;

int main(int argc, char * argv[]) {

  // first argument - text file (list of files)  
  // second argument - isolation "DirIso" or "InvIso"
 
  // ***********************************************
  // ***************  steering cards ***************
  // ***********************************************

  // boolean cards
  bool applyMuonId=true;
  bool applyPUReweighting=true;
  bool reweightWithPUI=true;
  bool applyTrigger=true;
//PU reweighting files 
  TString puWeightsFileData = "Data_Pileup_2012_ReReco-600bins.root";
  TString puWeightsFileMc = "MC_Summer12_PU_S10-600bins.root";
	
  bool applyMuonIdSF = true;

  bool randomizeMuons = true;

  // kinematic cuts
  float ptHardLepCut = 17.0;
  float ptSoftLepCut = 10.0;
  float MuEtaTriggerCut = 2.1;
  float MuEtaCut = 2.1;


  // cuts on track
  float ptHardTrkCut = 2.5;
  float ptSoftTrkCut = 2.5;

  // loose bkg control region cuts
  float ptTrkLower = 1.0;
  float ptTrkUpper = 2.5;
  int NsoftTracks = 1;

  float ptTrkCut = 1.;
  float etaTrkCut = 99999;
  float ptGenPartCut = 0.5;
  float d0Cut = 0.03;
  float dzCut = 0.1;

  // standard
  float dzCutLoose = 1.0;
  float d0CutLoose = 1.0;
  // loose
  //  float dzCutLoose = 0.5;
  //  float d0CutLoose = 0.5;
  //  Super loose
  //  float dzCutLoose = 0.4;
  //  float d0CutLoose = 0.2;
  float d0CutTracks= 0.02; 
  float dzCutTracks= 0.04; 
  // looser conditions 
  //  float d0CutTracks= 0.1; 
  //  float dzCutTracks= 0.2;
  
  float d0CutLowerTracks = 0.02;
  float dzCutLowerTracks = 0.04;
  float dRmumuCut = 2.0;

  float jetEtaCut = 4.5;
  float bJetEtaCut = 2.4;
  float btagCSVL = 0.244;
  float btagCSVM = 0.679;

  
  if (argc!=2) {
    std::cout << "Usage of the program :  [file_list]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    exit(-1);
  }

  using namespace std;
 
  int run;
  int event;
  int lumi;

  int processId;


  int nProcessed;
  int nPV;
  int nPUI;
  float nPUITruth;
  float probPV[100];
  float ndofPV[100];
  unsigned int nTrkPV[100];
  float chi2PV[100];
  float xPV[100];
  float yPV[100];
  float zPV[100];
  float sumPtPV[100];

  // generator info

  float H2pt;
  float H2px;
  float H2py;
  float H2pz;
  float H2phi;
  float H2eta;
  float H2mass;

  float HiggsPt;


  float PtTau1;
  float PxTau1;
  float PyTau1;
  float PzTau1;
  float PhiTau1;
  float EtaTau1;
  int Tau1FinalStates;
  int FinalStateIdTau1[10];

  float PtTau2;
  float PxTau2;
  float PyTau2;
  float PzTau2;
  float PhiTau2;
  float EtaTau2;
  int Tau2FinalStates;
  int FinalStateIdTau2[10];

  float PtTau3;
  float PxTau3;
  float PyTau3;
  float PzTau3;
  float PhiTau3;
  float EtaTau3;
  int Tau3FinalStates;
  int FinalStateIdTau3[10];

  float PtTau4;
  float PxTau4;
  float PyTau4;
  float PzTau4;
  float PhiTau4;
  float EtaTau4;
  int Tau4FinalStates;
  int FinalStateIdTau4[10];
 
  float PtFirstH1;
  float PxFirstH1;
  float PyFirstH1;
  float PzFirstH1;
  float PhiFirstH1;
  float EtaFirstH1;
  float MassFirstH1;
  
  float PtSecondH1;
  float PxSecondH1;
  float PySecondH1;
  float PzSecondH1;
  float PhiSecondH1;
  float EtaSecondH1;
  float MassSecondH1;


  int NPartFirstH1;
  int IdPartFirstH1[10000];
  int QPartFirstH1[10000];
  float PxPartFirstH1[10000];
  float PyPartFirstH1[10000];
  float PzPartFirstH1[10000];
  float PtPartFirstH1[10000];
  float EtaPartFirstH1[10000];
  float PhiPartFirstH1[10000];

  int NPartSecondH1;
  int IdPartSecondH1[10000];
  int QPartSecondH1[10000];
  float PxPartSecondH1[10000];
  float PyPartSecondH1[10000];
  float PzPartSecondH1[10000];
  float PtPartSecondH1[10000];
  float EtaPartSecondH1[10000];
  float PhiPartSecondH1[10000];

  int NPartAroundMuFirstH1;
  int IdPartAroundMuFirstH1[10];
  int QPartAroundMuFirstH1[10];
  float PxPartAroundMuFirstH1[10];
  float PyPartAroundMuFirstH1[10];
  float PzPartAroundMuFirstH1[10];
  float PtPartAroundMuFirstH1[10];
  float EtaPartAroundMuFirstH1[10];
  float PhiPartAroundMuFirstH1[10];

  int NPartAroundMuSecondH1;
  int IdPartAroundMuSecondH1[10];
  float QPartAroundMuSecondH1[10];
  float PxPartAroundMuSecondH1[10];
  float PyPartAroundMuSecondH1[10];
  float PzPartAroundMuSecondH1[10];
  float PtPartAroundMuSecondH1[10];
  float EtaPartAroundMuSecondH1[10];
  float PhiPartAroundMuSecondH1[10];

  int _bosonPDG;
  int _decayPDG;
  bool _bosonExist;	
  float _ZMass;
  float _ZPx;
  float _ZPy;
  float _ZPz;
  
  // pos muon --->
  int posMuonGenpdg;
  bool posMuonGenExist;
  float posMuonGenQ;
  float posMuonGenpx;
  float posMuonGenpy;
  float posMuonGenpz;
  float posMuonGenpt;
  float posMuonGeneta;
  float posMuonGenphi;

  bool PMmotherExist;
  int PMmotherpdg;
  float PMmothereta;
  float PMmotherphi;
  float PMmotherpt;
  float PMmotherpx;
  float PMmotherpy;
  float PMmotherpz;

  bool PMGrmotherExist;
  int PMGrmotherpdg;
  float PMGrmothereta;
  float PMGrmotherphi;
  float PMGrmotherpt;
  float PMGrmotherpx;
  float PMGrmotherpy;
  float PMGrmotherpz;

  // neg muon --->
  int negMuonGenpdg;
  bool negMuonGenExist;
  float negMuonGenQ;
  float negMuonGenpx;
  float negMuonGenpy;
  float negMuonGenpz;
  float negMuonGenpt;
  float negMuonGeneta;
  float negMuonGenphi;

  int NMmotherpdg;
  bool NMmotherExist;
  float NMmothereta;
  float NMmotherphi;
  float NMmotherpt;
  float NMmotherpx;
  float NMmotherpy;
  float NMmotherpz;
 
  bool NMGrmotherExist;
  int NMGrmotherpdg;
  float NMGrmothereta;
  float NMGrmotherphi;
  float NMGrmotherpt;
  float NMGrmotherpx;
  float NMGrmotherpy;
  float NMGrmotherpz;
 


  float weight;

  float XSection;

  std::ifstream fileList(argv[1]);
  
  int nFiles;

  fileList >> nFiles;

  //  std::string chainNameGen("generatorNMSSM/Gen4Taus");
  std::string chainNameGen("DiMuonAnalysis/Gen4Taus");

  TString TStrName(argv[1]);
  std::cout <<TStrName <<std::endl;  
  bool isData = false;
  if (TStrName.Contains("Data")) {
    isData = true;
   
    std::cout << "=============================" << std::endl;
    std::cout << "==== Running on data ========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;    
    processId = 0;
  }
  else {
   
    std::cout << "=============================" << std::endl;
    std::cout << "== Running on Monte Carlo ===" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl; 
  }

  bool isQCD = false;
  if (TStrName.Contains("QCD")) {
    isQCD = true;
    std::cout << "=============================" << std::endl;
    std::cout << "====== QCD Monte Carlo ======" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    processId = 5;
    XSection = 84679;
  }
  
  bool isM10to20 = false;
  bool isMadgraph = false;

  bool isTTJets = false;
  if (TStrName.Contains("TTJets")) {
    std::cout << "=============================" << std::endl;
    std::cout << "========== TTJets ===========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    isTTJets = true;
    processId = 6;
    XSection = 158;
  }

  bool isDYToTauTau = false;
  if (TStrName.Contains("DYToTauTau")||TStrName.Contains("ZToTauTau")) {
    isDYToTauTau = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->tautau ======" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    processId = 1;
    XSection = 1667;
    if (TStrName.Contains("M10to20")) {
      isM10to20 = true;
      processId = 2;
      XSection = 3893;
    }
    if (TStrName.Contains("Madgraph")||TStrName.Contains("madgraph")) {
      isMadgraph = true;
      XSection = 3048;
    }
  }

  bool isDYToMuMu = false;
  if (TStrName.Contains("DYToMuMu")||TStrName.Contains("ZToMuMu")) {
    isDYToMuMu = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->mumu ========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;    
    processId = 3;
    XSection = 1667;
    if (TStrName.Contains("M10to20")) {
      isM10to20 = true;
      processId = 4;
      XSection = 3893;
    }
    if (TStrName.Contains("Madgraph")||TStrName.Contains("madgraph")) {
      isMadgraph = true;
      XSection = 3048;
    }
  }

  bool isDYToLL = false;
  if (TStrName.Contains("DYToLL")||TStrName.Contains("ZToLL")) {
    isDYToLL = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->l+l- ========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;    
    processId = 3;
    XSection = 1667;
    if (TStrName.Contains("M10to20")) {
      isM10to20 = true;
      processId = 4;
      XSection = 3893;
    }
    if (TStrName.Contains("Madgraph")||TStrName.Contains("madgraph")) {
      isMadgraph = true;
      XSection = 3048;
    }
  }

  if (isMadgraph) {
    std::cout << "=============================" << std::endl;
    std::cout << "==== Madgraph Monte Carlo ===" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
  }

  bool isHto4taus = false;
  if (TStrName.Contains("NMSSM_H2ToH1H1")) {
    isHto4taus = true;
    std::cout << "===================================" << std::endl;
    std::cout << "==== H2->2H1->4taus Monte Carlo ===" << std::endl;
    std::cout << "===================================" << std::endl;
  }

  
  TString FullName(argv[1]) ;
  
  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");

  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);

  // *** additional gen info
  TH1F * nMuonsFirstGenH = new TH1F("nMuonsFirstGenH","",4,-0.5,3.5);
  TH1F * nMuonsSecondGenH = new TH1F("nMuonsSecondGenH","",4,-0.5,3.5);

  TH1F * nTrkFirstGenH = new TH1F("nTrkFirstGenH","",4,-0.5,3.5);
  TH1F * nTrkSecondGenH = new TH1F("nTrkSecondGenH","",4,-0.5,3.5);
  
  TH1F * nEleFirstGenH = new TH1F("nEleFirstGenH","",4,-0.5,3.5);
  TH1F * nEleSecondGenH = new TH1F("nEleSecondGenH","",4,-0.5,3.5);

  TH1F * ptMuHardGenH = new TH1F("ptMuHardGenH","",100,0.,100.);
  TH1F * ptMuSoftGenH = new TH1F("ptMuSoftGenH","",100,0.,100.);

  TH1F * etaMuHardGenH = new TH1F("etaMuHardGenH","",80,-4.,4.);
  TH1F * etaMuSoftGenH = new TH1F("etaMuSoftGenH","",80,-4.,4.);

  TH1F * ptTrkFirstGenH = new TH1F("ptTrkFirstGenH","",200,0.,100.);
  TH1F * ptTrkSecondGenH = new TH1F("ptTrkSecondGenH","",200,0.,100.);

  TH1F * etaTrkFirstGenH = new TH1F("etaTrkFirstGenH","",80,-4.,4.);
  TH1F * etaTrkSecondGenH = new TH1F("etaTrkSecondGenH","",80,-4.,4.);

  TH1F * dRMuonsGenH = new TH1F("dRMuonsGenH","",100.,0.,5.);
  TH1F * dRMuonTrkFirstGenH = new TH1F("dRMuonTrkFirstGenH","",100.,0.,5.);
  TH1F * dRMuonTrkSecondGenH = new TH1F("dRMuonTrkSecondGenH","",100.,0.,5.);

  TH1F * muTrkMassFirstH = new TH1F("muTrkMassFirstH","",20.,0.,10.);
  TH1F * muTrkMassFirstCutsH = new TH1F("muTrkMassFirstCutsH","",20.,0.,10.);

  TH1F * muTrkMassSecondH = new TH1F("muTrkMassSecondH","",20.,0.,10.);
  TH1F * muTrkMassSecondCutsH = new TH1F("muTrkMassSecondCutsH","",20.,0.,10.);

  TH1F * muTrkMassH = new TH1F("muTrkMassH","",20.,0.,10.);
  TH1F * muTrkMassCutsH = new TH1F("muTrkMassCutsH","",20.,0.,10.);

  TH1F * HiggsPtH = new TH1F("HiggsPtH","",100,0.,200.);


  TH1F * puiweightMC;
  TH1F * puiweightData;

  for (int iFiles=0; iFiles<nFiles;++iFiles) {

   std::string filen;
   fileList >> filen;
   std::cout << iFiles+1 << " out of " << nFiles << "  ===> filename : " << filen << std::endl;
   TFile * file_ = TFile::Open(TString(filen));
   
   TChain * _genTree = new TChain(TString(chainNameGen));
    
   TH1F * histoInputEvents = NULL;
   
   histoInputEvents = (TH1F*)file_->Get("eventCountPostTrigger/EventCount");

   if (histoInputEvents==NULL) continue;

   int NE = int(histoInputEvents->GetEntries());

   std::cout << "      number of input events = " << NE << std::endl;

   for (int iE=0;iE<NE;++iE)
     inputEventsH->Fill(0.);


   _genTree->Add(TString(filen));

   
   _genTree->SetBranchAddress("H2pt",&H2pt);
   _genTree->SetBranchAddress("H2px",&H2px);
   _genTree->SetBranchAddress("H2py",&H2py);
   _genTree->SetBranchAddress("H2pz",&H2pz);
   _genTree->SetBranchAddress("H2phi",&H2phi);
   _genTree->SetBranchAddress("H2eta",&H2eta);
   _genTree->SetBranchAddress("H2mass",&H2mass);
   
   _genTree->SetBranchAddress("PtTau1",&PtTau1);
   _genTree->SetBranchAddress("PxTau1",&PxTau1);
   _genTree->SetBranchAddress("PyTau1",&PyTau1);
   _genTree->SetBranchAddress("PzTau1",&PzTau1);
   _genTree->SetBranchAddress("PhiTau1",&PhiTau1);
   _genTree->SetBranchAddress("EtaTau1",&EtaTau1);
   _genTree->SetBranchAddress("Tau1FinalStates",&Tau1FinalStates);
   _genTree->SetBranchAddress("FinalStateIdTau1",FinalStateIdTau1);
   
   _genTree->SetBranchAddress("PtTau2",&PtTau2);
   _genTree->SetBranchAddress("PxTau2",&PxTau2);
   _genTree->SetBranchAddress("PyTau2",&PyTau2);
   _genTree->SetBranchAddress("PzTau2",&PzTau2);
   _genTree->SetBranchAddress("PhiTau2",&PhiTau2);
   _genTree->SetBranchAddress("EtaTau2",&EtaTau2);
   _genTree->SetBranchAddress("Tau2FinalStates",&Tau2FinalStates);
   _genTree->SetBranchAddress("FinalStateIdTau2",FinalStateIdTau2);
   
   _genTree->SetBranchAddress("PtTau3",&PtTau3);
   _genTree->SetBranchAddress("PxTau3",&PxTau3);
   _genTree->SetBranchAddress("PyTau3",&PyTau3);
   _genTree->SetBranchAddress("PzTau3",&PzTau3);
   _genTree->SetBranchAddress("PhiTau3",&PhiTau3);
   _genTree->SetBranchAddress("EtaTau3",&EtaTau3);
   _genTree->SetBranchAddress("Tau3FinalStates",&Tau3FinalStates);
   _genTree->SetBranchAddress("FinalStateIdTau3",FinalStateIdTau3);
   
   _genTree->SetBranchAddress("PtTau4",&PtTau4);
   _genTree->SetBranchAddress("PxTau4",&PxTau4);
   _genTree->SetBranchAddress("PyTau4",&PyTau4);
   _genTree->SetBranchAddress("PzTau4",&PzTau4);
   _genTree->SetBranchAddress("PhiTau4",&PhiTau4);
   _genTree->SetBranchAddress("EtaTau4",&EtaTau4);
   _genTree->SetBranchAddress("Tau4FinalStates",&Tau4FinalStates);
   _genTree->SetBranchAddress("FinalStateIdTau4",FinalStateIdTau4);
   
   _genTree->SetBranchAddress("PtFirstH1",&PtFirstH1);
   _genTree->SetBranchAddress("PxFirstH1",&PxFirstH1);
   _genTree->SetBranchAddress("PyFirstH1",&PyFirstH1);
   _genTree->SetBranchAddress("PzFirstH1",&PzFirstH1);
   _genTree->SetBranchAddress("PhiFirstH1",&PhiFirstH1);
   _genTree->SetBranchAddress("EtaFirstH1",&EtaFirstH1);
   _genTree->SetBranchAddress("MassFirstH1",&MassFirstH1);

   _genTree->SetBranchAddress("PtSecondH1",&PtSecondH1);
   _genTree->SetBranchAddress("PxSecondH1",&PxSecondH1);
   _genTree->SetBranchAddress("PySecondH1",&PySecondH1);
   _genTree->SetBranchAddress("PzSecondH1",&PzSecondH1);
   _genTree->SetBranchAddress("PhiSecondH1",&PhiSecondH1);
   _genTree->SetBranchAddress("EtaSecondH1",&EtaSecondH1);
   _genTree->SetBranchAddress("MassSecondH1",&MassSecondH1);
   
   _genTree->SetBranchAddress("NPartFirstH1",&NPartFirstH1);
   _genTree->SetBranchAddress("IdPartFirstH1",IdPartFirstH1);
   _genTree->SetBranchAddress("QPartFirstH1", QPartFirstH1);
   _genTree->SetBranchAddress("PxPartFirstH1",PxPartFirstH1);
   _genTree->SetBranchAddress("PyPartFirstH1",PyPartFirstH1);
   _genTree->SetBranchAddress("PzPartFirstH1",PzPartFirstH1);
   _genTree->SetBranchAddress("PtPartFirstH1",PtPartFirstH1);
   _genTree->SetBranchAddress("EtaPartFirstH1",EtaPartFirstH1);
   _genTree->SetBranchAddress("PhiPartFirstH1",PhiPartFirstH1);
   
   _genTree->SetBranchAddress("NPartAroundMuFirstH1",&NPartAroundMuFirstH1);
   _genTree->SetBranchAddress("IdPartAroundMuFirstH1",IdPartAroundMuFirstH1);
   _genTree->SetBranchAddress("QPartAroundMuFirstH1", QPartAroundMuFirstH1);
   _genTree->SetBranchAddress("PxPartAroundMuFirstH1",PxPartAroundMuFirstH1);
   _genTree->SetBranchAddress("PyPartAroundMuFirstH1",PyPartAroundMuFirstH1);
   _genTree->SetBranchAddress("PzPartAroundMuFirstH1",PzPartAroundMuFirstH1);
   _genTree->SetBranchAddress("PtPartAroundMuFirstH1",PtPartAroundMuFirstH1);
   _genTree->SetBranchAddress("EtaPartAroundMuFirstH1",EtaPartAroundMuFirstH1);
   _genTree->SetBranchAddress("PhiPartAroundMuFirstH1",PhiPartAroundMuFirstH1);
   
   _genTree->SetBranchAddress("NPartSecondH1",&NPartSecondH1);
   _genTree->SetBranchAddress("IdPartSecondH1",IdPartSecondH1);
   _genTree->SetBranchAddress("QPartSecondH1", QPartSecondH1);
   _genTree->SetBranchAddress("PxPartSecondH1",PxPartSecondH1);
   _genTree->SetBranchAddress("PyPartSecondH1",PyPartSecondH1);
   _genTree->SetBranchAddress("PzPartSecondH1",PzPartSecondH1);
   _genTree->SetBranchAddress("PtPartSecondH1",PtPartSecondH1);
   _genTree->SetBranchAddress("EtaPartSecondH1",EtaPartSecondH1);
   _genTree->SetBranchAddress("PhiPartSecondH1",PhiPartSecondH1);
   
   _genTree->SetBranchAddress("NPartAroundMuSecondH1",&NPartAroundMuSecondH1);
   _genTree->SetBranchAddress("IdPartAroundMuSecondH1",IdPartAroundMuSecondH1);
   _genTree->SetBranchAddress("QPartAroundMuSecondH1", QPartAroundMuSecondH1);
   _genTree->SetBranchAddress("PxPartAroundMuSecondH1",PxPartAroundMuSecondH1);
   _genTree->SetBranchAddress("PyPartAroundMuSecondH1",PyPartAroundMuSecondH1);
   _genTree->SetBranchAddress("PzPartAroundMuSecondH1",PzPartAroundMuSecondH1);
   _genTree->SetBranchAddress("PtPartAroundMuSecondH1",PtPartAroundMuSecondH1);
   _genTree->SetBranchAddress("EtaPartAroundMuSecondH1",EtaPartAroundMuSecondH1);
   _genTree->SetBranchAddress("PhiPartAroundMuSecondH1",PhiPartAroundMuSecondH1);
   
   
   int numberOfCandidates = _genTree->GetEntries();
  
   std::cout << "      number of selected candidates = " << numberOfCandidates << std::endl;
   
   for (int iCand=0; iCand<numberOfCandidates; iCand++) { 
    
     _genTree->GetEntry(iCand);


     HiggsPtH->Fill(H2pt);

     float ptMuFirstGen = 0;
     float etaMuFirstGen = 0;
     float phiMuFirstGen = 0;     

     float ptMuSecondGen = 0;
     float etaMuSecondGen = 0;
     float phiMuSecondGen = 0;

     float ptTrkFirstGen = 0;
     float phiTrkFirstGen = 0;
     float etaTrkFirstGen = 0;
     
     float ptTrkSecondGen = 0;
     float phiTrkSecondGen = 0;
     float etaTrkSecondGen = 0;

     int nMuFirstGen = 0;
     int nMuSecondGen = 0;
     
     int nTrkFirstGen = 0;
     int nTrkSecondGen = 0;
     
     int nEleFirstGen = 0;
     int nEleSecondGen = 0;	

     int nChargedFirst = 0;
     float qTotalFirst = 0;
     float pxChargedFirst[10];
     float pyChargedFirst[10];
     float pzChargedFirst[10];
     float ptChargedFirst[10];
     float etaChargedFirst[10];
     float massChargedFirst[10];

     int nChargedSecond = 0;
     float qTotalSecond = 0;
     float pxChargedSecond[10];
     float pyChargedSecond[10];
     float pzChargedSecond[10];
     float ptChargedSecond[10];
     float etaChargedSecond[10];
     float massChargedSecond[10];

     for (int iP=0; iP<NPartFirstH1; ++iP) {
       

       bool isPosMuon = IdPartFirstH1[iP]==-13;
       bool isNegMuon = IdPartFirstH1[iP]==13;
       bool isMuon = isPosMuon || isNegMuon;
       
       bool isPosElec = IdPartFirstH1[iP]==-11;
       bool isNegElec = IdPartFirstH1[iP]==11;
       bool isElec = isPosElec || isNegElec;
       
       bool isPosPion = IdPartFirstH1[iP]==211;
       bool isNegPion = IdPartFirstH1[iP]==-211;
       bool isPion = isPosPion || isNegPion;
       
       bool isPosKaon = IdPartFirstH1[iP]==321;
       bool isNegKaon = IdPartFirstH1[iP]==-321;
       bool isKaon = isPosKaon || isNegKaon;

       float mass = muonMass;
       
       if (isElec)
	 mass = elecMass;

       if (isPion)
	 mass = pionMass;

       if (isKaon)
	 mass = kaonMass;

       float etaAbs = TMath::Abs(EtaPartFirstH1[iP]);

       if (isMuon) {
	 if (etaAbs<2.1) {
	   nMuFirstGen++;
	   if (PtPartFirstH1[iP]>ptMuFirstGen) {
	     ptMuFirstGen = PtPartFirstH1[iP];
	     etaMuFirstGen = EtaPartFirstH1[iP];
	     phiMuFirstGen = PhiPartFirstH1[iP];
	   }
	 }

       }
       if (isPion||isKaon) {
	 nTrkFirstGen++;
	 ptTrkFirstGen = PtPartFirstH1[iP];
	 etaTrkFirstGen = EtaPartFirstH1[iP];
	 phiTrkFirstGen = PhiPartFirstH1[iP];
       }
       if (isElec) {
	 nEleFirstGen++;
       }

       if ((QPartFirstH1[iP]<-0.5||QPartFirstH1[iP]>0.5)&&etaAbs<2.4) {

	 qTotalFirst += QPartFirstH1[iP];

	 pxChargedFirst[nChargedFirst] = PxPartFirstH1[iP];
	 pyChargedFirst[nChargedFirst] = PyPartFirstH1[iP];
	 pzChargedFirst[nChargedFirst] = PzPartFirstH1[iP];
	 massChargedFirst[nChargedFirst] = mass;
	 
	 ptChargedFirst[nChargedFirst] = PtPartFirstH1[iP];
	 etaChargedFirst[nChargedFirst] = EtaPartFirstH1[iP];

	 nChargedFirst++;

       }

       
     }
     
     for (int iP=0; iP<NPartSecondH1; ++iP) {
	  
       bool isPosMuon = IdPartSecondH1[iP]==-13;
       bool isNegMuon = IdPartSecondH1[iP]==13;
       bool isMuon = isPosMuon || isNegMuon;
       
       bool isPosElec = IdPartSecondH1[iP]==-11;
       bool isNegElec = IdPartSecondH1[iP]==11;
       bool isElec = isPosElec || isNegElec;
       
       bool isPosPion = IdPartSecondH1[iP]==211;
       bool isNegPion = IdPartSecondH1[iP]==-211;
       bool isPion = isPosPion || isNegPion;
       
       bool isPosKaon = IdPartSecondH1[iP]==321;
       bool isNegKaon = IdPartSecondH1[iP]==-321;
       bool isKaon = isPosKaon || isNegKaon;

       float mass = muonMass;

       if (isElec)
         mass = elecMass;

       if (isPion)
         mass = pionMass;

       if (isKaon)
         mass = kaonMass;

       float etaAbs = TMath::Abs(EtaPartSecondH1[iP]);
 
       if (isMuon) {
         if (etaAbs<2.1) {
           nMuSecondGen++;
           if (PtPartSecondH1[iP]>ptMuSecondGen) {
	     ptMuSecondGen = PtPartSecondH1[iP];
	     etaMuSecondGen = EtaPartSecondH1[iP];
	     phiMuSecondGen = PhiPartSecondH1[iP];
	   }
	 }
       }
       if (isPion||isKaon) {
	 nTrkSecondGen++;
	 ptTrkSecondGen = PtPartSecondH1[iP];
	 etaTrkSecondGen = EtaPartSecondH1[iP];
	 phiTrkSecondGen = PhiPartSecondH1[iP];
       }
       if (isElec) {
	 nEleSecondGen++;
       }
       
       if ((QPartSecondH1[iP]<-0.5||QPartSecondH1[iP]>0.5)&&etaAbs<2.4) {

         qTotalSecond += QPartSecondH1[iP];

         pxChargedSecond[nChargedSecond] = PxPartSecondH1[iP];
         pyChargedSecond[nChargedSecond] = PyPartSecondH1[iP];
         pzChargedSecond[nChargedSecond] = PzPartSecondH1[iP];
         massChargedSecond[nChargedSecond] = mass;

         ptChargedSecond[nChargedSecond] = PtPartSecondH1[iP];
         etaChargedSecond[nChargedSecond] = EtaPartSecondH1[iP];

         nChargedSecond++;

       }

     }


     
     float ptMuHardGen = ptMuFirstGen;
     float ptMuSoftGen = ptMuSecondGen;
     float etaMuHardGen = etaMuFirstGen;
     float etaMuSoftGen = etaMuSecondGen;
     
     if (ptMuSecondGen>ptMuFirstGen) {
       ptMuSoftGen = ptMuFirstGen;
       ptMuHardGen = ptMuSecondGen;
       etaMuSoftGen = etaMuFirstGen;
       etaMuHardGen = etaMuSecondGen;
     }

     

     if (nChargedFirst==2 && (qTotalFirst<1e-2 && qTotalFirst>-1e-2) && nMuFirstGen>0) {

       TLorentzVector firstPart; firstPart.SetXYZM(pxChargedFirst[0],
						   pyChargedFirst[0],
						   pzChargedFirst[0],
						   massChargedFirst[0]);

       TLorentzVector secondPart; secondPart.SetXYZM(pxChargedFirst[1],
						     pyChargedFirst[1],
						     pzChargedFirst[1],
						     massChargedFirst[1]);

       TLorentzVector muTrk = firstPart + secondPart;

       float mass = float(muTrk.M());

       muTrkMassFirstH->Fill(mass);
       muTrkMassH->Fill(mass);

       if (ptMuHardGen>17 && ptMuSoftGen>10 && ptChargedFirst[0]>2.5 && ptChargedFirst[1]>2.5) {
	 muTrkMassFirstCutsH->Fill(mass);
	 muTrkMassCutsH->Fill(mass);
       }

     }

     if (nChargedSecond==2 && (qTotalSecond<1e-2 && qTotalSecond>-1e-2) && nMuSecondGen>0) {

       TLorentzVector firstPart; firstPart.SetXYZM(pxChargedSecond[0],
                                                   pyChargedSecond[0],
                                                   pzChargedSecond[0],
                                                   massChargedSecond[0]);

       TLorentzVector secondPart; secondPart.SetXYZM(pxChargedSecond[1],
                                                     pyChargedSecond[1],
                                                     pzChargedSecond[1],
                                                     massChargedSecond[1]);

       TLorentzVector muTrk = firstPart + secondPart;

       float mass = float(muTrk.M());

       muTrkMassSecondH->Fill(mass);
       muTrkMassH->Fill(mass);

       if (ptMuHardGen>17 && ptMuSoftGen>10 && ptChargedSecond[0]>2.5 && ptChargedSecond[1]>2.5) {
         muTrkMassSecondCutsH->Fill(mass);
         muTrkMassCutsH->Fill(mass);
       }

     }


     float DeltaRmuons =  float(deltaR(etaMuFirstGen,phiMuFirstGen,
				       etaMuSecondGen,phiMuSecondGen));
     
     float DeltaRMuTrkFirst = float(deltaR(etaMuFirstGen,phiMuFirstGen,
					   etaTrkFirstGen,phiTrkFirstGen));
     
     float DeltaRMuTrkSecond = float(deltaR(etaMuSecondGen,phiMuSecondGen,
					    etaTrkSecondGen,phiTrkSecondGen));
     
     nMuonsFirstGenH->Fill(float(nMuFirstGen));
     nMuonsSecondGenH->Fill(float(nMuSecondGen));
     
     nEleFirstGenH->Fill(float(nEleFirstGen));
     nEleSecondGenH->Fill(float(nEleSecondGen));
     
     nTrkFirstGenH->Fill(float(nTrkFirstGen));
     nTrkSecondGenH->Fill(float(nTrkSecondGen));
     
     etaMuHardGenH->Fill(etaMuHardGen);
     etaMuSoftGenH->Fill(etaMuSoftGen);
     
     float absEtaMuHardGen = TMath::Abs(etaMuHardGen);
     float absEtaMuSoftGen = TMath::Abs(etaMuSoftGen);

     float absEtaTrkFirstGen = TMath::Abs(etaTrkFirstGen);
     float absEtaTrkSecondGen = TMath::Abs(etaTrkSecondGen);

     if (nMuFirstGen==1 && nMuSecondGen==1) {

       if (absEtaMuHardGen<2.1) 
	 ptMuHardGenH->Fill(ptMuHardGen);
       
       if (absEtaMuSoftGen<2.1) 
	 ptMuSoftGenH->Fill(ptMuSoftGen);
       
       if (absEtaMuHardGen<2.1 && absEtaMuSoftGen<2.1 && ptMuHardGen>17. && ptMuSoftGen>10.) { // muon kinematics
	 
	 dRMuonsGenH->Fill(DeltaRmuons);
       
	 if (nTrkFirstGen==1) {
	   etaTrkFirstGenH->Fill(etaTrkFirstGen);
	   if (absEtaTrkFirstGen<2.4)
	     ptTrkFirstGenH->Fill(ptTrkFirstGen);
	   if (absEtaTrkFirstGen<2.4 && ptTrkFirstGen>2.5) {
	     dRMuonTrkFirstGenH->Fill(DeltaRMuTrkFirst);
	   }
	 }
	 
	 if (nTrkSecondGen==1) {
	   etaTrkSecondGenH->Fill(etaTrkSecondGen);
	   if (absEtaTrkSecondGen<2.4)
	     ptTrkSecondGenH->Fill(ptTrkSecondGen);
	   if (absEtaTrkSecondGen<2.4 && ptTrkSecondGen>2.5) {
	     dRMuonTrkSecondGenH->Fill(DeltaRMuTrkSecond);
	   }
	 }
	 
       }
       
     }

   }
  }
  
  nProcessed = int(inputEventsH->GetEntries());
  std::cout <<nProcessed<< std::endl;
  
  file->cd("");
  file->Write();
  file->Close();
  delete file;
 
}

 

