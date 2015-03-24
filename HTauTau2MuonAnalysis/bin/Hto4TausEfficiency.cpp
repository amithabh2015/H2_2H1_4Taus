#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

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
#include "H2to2H1to4Taus/TauTauSkimming/interface/EventWeight.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/RunLumiReader.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauUtils.h"

#include "TLorentzVector.h"

#include "TRandom.h"

const float MuMass = 0.105658367;


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

int main(int argc, char * argv[]) {

  // first argument - text file (list of files)  
  // second argument - isolation "DirIso" or "InvIso"
 
  // ***********************************************
  // ***************  steering cards ***************
  // ***********************************************

  // boolean cards
  bool applyPUReweighting=true;
  bool reweightWithPUI=true;
  bool applyTrigger=true;
  // trigger
  // 0 : HLT_Mu24
  // 1 : HLT_Mu40
  // 2 : HLT_IsoMu24
  int triggerType = 0;

  //PU reweighting files
 
  TString puWeightsFileData = "Data_Pileup_2012_ReReco-600bins.root";
  TString puWeightsFileMc = "MC_Summer12_PU_S10-600bins.root";

  bool randomizeMuons = true;

  // kinematic cuts
  float ptTagLepCut   = 24.0;
  float ptProbeLepCut = 10.0;
  float etaTagLepCut = 2.1;
  float etaProbeLepCut = 2.4;

  // cuts on track
  //float ptHardTrkCut = 2.5;
  float ptSoftTrkCut = 2.5;

  // loose bkg control region cuts
  //float ptTrkLower = 1.0;
  //float ptTrkUpper = 1.5;

  float ptTrkCut = 1.0;
  float etaTrkCut = 2.4;
  //float ptGenPartCut = 0.5;
  float d0Cut = 0.03;
  float dzCut = 0.1;
  float dzCutLoose = 1.;
  float d0CutLoose = 1.;
  float d0CutTracks= 0.02; 
  float dzCutTracks= 0.04; 

  float relIsoCut = 0.15;

  /*float jetEtaCut = 4.5;
  float bJetEtaCut = 2.4;
  float btagCSVL = 0.244;
  float btagCSVM = 0.679;*/

  
  std::string firstArg(argv[1]);
  TString FirstArg(argv[1]);
  if (argc!=4) {
    std::cout << "Usage of the program : Hto4TausEfficiency [file_list] [apply_trigger] [trigger]" << std::endl;
    std::cout << "file_list      : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    std::cout << "apply_trigger  : True or False" << std::endl;
    std::cout << "trigger        : HLT_Mu24 , HLT_Mu40, HLT_IsoMu24" << std::endl;
    exit(-1);
  }
  else if (FirstArg=="--help"||FirstArg=="-help") {
    std::cout << "Usage of the program : Hto4TausEfficiency [file_list] [apply_trigger] [trigger]" << std::endl;
    std::cout << "file_list      : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    std::cout << "apply_trigger  : True or False" << std::endl;
    std::cout << "trigger        : HLT_Mu24 , HLT_Mu40, HLT_IsoMu24" << std::endl;
  }

  if (argc!=4)
    exit(-1);

  TString ApplyTrigger(argv[2]);
  if (ApplyTrigger=="False"||ApplyTrigger=="false") {
    applyTrigger = false;
  }
  else if (ApplyTrigger=="True"||ApplyTrigger=="true") {
    applyTrigger = true;
  }
  else {
    std::cout << "Specify for apply_trigger option False or True" << std::endl;
    exit(-1);
  }
  TString Trigger(argv[3]);
  if (Trigger=="HLT_Mu24") {
    triggerType = 0;
  }
  else if (Trigger=="HLT_Mu40") {
    triggerType = 1;
  }
  else if (Trigger=="HLT_IsoMu24") {
    triggerType = 2;
  }
  else {
    std::cout << "Specified trigger options : HLT_Mu24,  HLT_Mu40,  HLT_IsoMu24" << std::endl;
    exit(-1);
  }


  using namespace std;
 
  int run;
  int event;
  int lumi;

  //int processId;


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
 
  //Tracks from PV with highest pt sum
  int NPVtracks;
  //int ALLNPVtracks;
  float PtTrack[10000];
  float PxTrack[10000];
  float PyTrack[10000];
  float PzTrack[10000];
  //float PhiTrack;
  float EtaTrack[10000];
  int ChargeTrack[10000];
  double Chi2Track[10000];
  
  // tracks around Pos muon
  /*int PVtracksDRpos;
  float ptPVTracksPos[10000];
  float pxPVTracksPos[10000];
  float pyPVTracksPos[10000];
  float pzPVTracksPos[10000];
  int qPVTracksPos[10000];*/
  
  // tracks around Neg muon
  /*int PVtracksDRneg;
  float ptPVTracksNeg[10000];
  float pxPVTracksNeg[10000];
  float pyPVTracksNeg[10000];
  float pzPVTracksNeg[10000];
  int qPVTracksNeg[10000];*/

  // pfNoPu around Pos muon
 /* int NpfNoPuPos;
  float PtPfNoPuPos[10000];
  float PxPfNoPuPos[10000];
  float PyPfNoPuPos[10000];
  float PzPfNoPuPos[10000];
  float EnPfNoPuPos[10000];
  float IPPfNoPuPos[10000];
  float IPZPfNoPuPos[10000];
  int QPfNoPuPos[10000];*/
  
  // pfNoPu around Neg muon
  /*int NpfNoPuNeg;
  float PtPfNoPuNeg[10000];
  float PxPfNoPuNeg[10000];
  float PyPfNoPuNeg[10000];
  float PzPfNoPuNeg[10000];
  float EnPfNoPuNeg[10000];
  float IPPfNoPuNeg[10000];
  float IPZPfNoPuNeg[10000];
  int QPfNoPuNeg[10000];*/

  //  pfNoPu around Pos muon (Ip cuts)
  /*int NpfNoPuIpPos;
  int NpfNoPuIpPosAll;
  float PtPfNoPuIpPos[10000];
  float PxPfNoPuIpPos[10000];
  float PyPfNoPuIpPos[10000];
  float PzPfNoPuIpPos[10000];
  float EnPfNoPuIpPos[10000];
  float IPPfNoPuIpPos[10000];
  float IPZPfNoPuIpPos[10000];
  int QPfNoPuIpPos[10000];
  int IdNoPuIpPos[10000];*/

  // pfNoPu around Neg muon (Ip cuts)
  /*int NpfNoPuIpNeg;
  int NpfNoPuIpNegAll;
  float PtPfNoPuIpNeg[10000];
  float PxPfNoPuIpNeg[10000];
  float PyPfNoPuIpNeg[10000];
  float PzPfNoPuIpNeg[10000];
  float EnPfNoPuIpNeg[10000];
  float IPPfNoPuIpNeg[10000];
  float IPZPfNoPuIpNeg[10000];
  int QPfNoPuIpNeg[10000];
  int IdNoPuIpNeg[10000];*/

  //PF objects from the No PU collection
  int NpfNoPU;
  //int AllNpfNoPU;
  //int AllNpfNoPU_ip;
  float PxpfNoPU[10000];
  float PypfNoPU[10000];
  float PzpfNoPU[10000];
  float PtpfNoPU[10000];
  float EnpfNoPU[10000];
  float EtapfNoPU[10000];
  float ChargepfNoPU[10000];
  float ippfNoPU[10000];
  float ipSigpfNoPU[10000];
  float ipZpfNoPU[10000];
  float ipZSigpfNoPU[10000];
  float ip3DpfNoPU[10000];
  float ip3DSigpfNoPU[10000];
  int   idpfNoPU[10000];
  //Outer cone
  /*int Ntracks_OutConePos;
  int Ntracks_OutConeNeg;
  int Ntracks_OutConePosIP;
  int Ntracks_OutConeNegIP;
  float SumPtPosCone;
  float SumPtNegCone;
  float SumPtPosConeIP;
  float SumPtNegConeIP;*/
  // Dilepton kinematics

  float met;
  float metPx;
  float metPy;
  float metPz;
  float metEn;
  float metCovXX;
  float metCovXY;
  float metCovYX;
  float metCovYY;

  float metMVAPx;
  float metMVAPy;
  float metMVAPz;
  float metMVAEn;
  float metMVACovXX;
  float metMVACovXY;
  float metMVACovYX;
  float metMVACovYY;
  
  /*float InvMassH2;
  float InvMassTrackPlusPosMuon;
  float InvMassTrackPlusNegMuon;*/

  float _diLepEta;
  float DiLepPt;
  float DiLepPhi;
  float _diLepDPhi;
  float _diLepDEta;
  float _diLepDR;
 

  float _posLepMETDPhi;
  float _negLepMETDPhi;
 

  float _negLepPt;
  float _posLepPt;
  float _negLepEta;
  float _posLepEta;
  float _negLepPhi;
  float _posLepPhi;
  float _negLepQ;
  float _posLepQ;

  bool NegLep_IsoMu24;
  bool PosLep_IsoMu24;
  bool NegLep_Mu24;
  bool PosLep_Mu24;
  //bool NegLep_Mu30;
  //bool PosLep_Mu30;
  bool NegLep_Mu40;
  bool PosLep_Mu40;

  bool NegLep_Mu17_Mu8_Leg17;
  bool PosLep_Mu17_Mu8_Leg17;
  bool NegLep_Mu17_Mu8_Leg8;
  bool PosLep_Mu17_Mu8_Leg8;
  bool NegLep_Mu17_Mu8_DzFilter;
  bool PosLep_Mu17_Mu8_DzFilter;

  bool NegLep_Mu17_TkMu8_Leg17;
  bool PosLep_Mu17_TkMu8_Leg17;
  bool NegLep_Mu17_TkMu8_Leg8;
  bool PosLep_Mu17_TkMu8_Leg8;
  bool NegLep_Mu17_TkMu8_DzFilter;
  bool PosLep_Mu17_TkMu8_DzFilter;


  bool HLT_Mu17_Mu8;
  bool HLT_Mu17_TkMu8;
  bool HLT_Mu24;
  bool HLT_Mu40;
  bool HLT_IsoMu24;

  int NTracksPos;
  /*int NTracksPos_NoMuons;
  int NTracksPos_NoMuons_ip;
  int NTracksPos_ip;*/

  float PxTrackPos[10000];
  float PyTrackPos[10000];
  float PzTrackPos[10000];
  float PtTrackPos[10000];
  float EnergyTrackPos[10000];
  float QTrackPos[10000];
  int IdTrackPos[10000];

  int NTracksNeg;
  /*int NTracksNeg_NoMuons;
  int NTracksNeg_NoMuons_ip;
  int NTracksNeg_ip;*/
  float PxTrackNeg[10000];
  float PyTrackNeg[10000];
  float PzTrackNeg[10000];
  float PtTrackNeg[10000];
  float EnergyTrackNeg[10000];
  float QTrackNeg[10000];
  int IdTrackNeg[10000];

  // Jets
  int _nJets;

  float jetPt[50];
  float jetEta[50];
  float jetPhi[50];

  int chargedMultiplicity[50];
  bool jetIDLoose[50];
  bool jetPUIDLoose[50];
  float combSVBJetTag[50];
 


  // Collinear Approx
  bool _validDiTau;
  float _diTauMass;

  float _posLepETauRest;
  float _negLepETauRest;
  float _posLepEdiTauRest;
  float _negLepEdiTauRest;
  float _posLepTauERatio;
  float _negLepTauERatio;

  float _cosNegDiTau;

  float _diTauPx;
  float _diTauPy;
  float _diTauPz;
  float _diTauPt;
  float _diTauEta;
  float _diTauE;
 
  // Lepton isolation
  float IsoNeg;
  float IsoPos;
  float IsoTrkNeg;
  float IsoTrkPos;
  float IsoECalNeg;
  float IsoECalPos;
  float IsoHCalNeg;
  float IsoHCalPos;
  float IsoPFsNeg;
  float IsoPFsPos;

  /*float IsoChargedHadPFsNeg;
  float IsoChargedHadPFsPos;
  float IsoNeutralHadPFsNeg;
  float IsoNeutralHadPFsPos;
  float IsoPhotonsPFsNeg;
  float IsoPhotonsPFsPos;*/
  

  bool isGlobalMuPos;
  bool isGlobalMuNeg;

  bool isTrackerMuPos;
  bool isTrackerMuNeg;

  bool isPFMuPos;
  bool isPFMuNeg;

  int nPixelHitsPos;
  int nPixelHitsNeg;

  int nTrackerHitsPos;
  int nTrackerHitsNeg;

  int nMuonHitsPos;
  int nMuonHitsNeg;

  int nMuonStationsPos;
  int nMuonStationsNeg;

  float Chi2PosMu;
  float Chi2NegMu;

  float NdofPosMu;
  float NdofNegMu;

  // impact parameter significance
  float IpSig2DPos[100];
  float IpSig2DNeg[100];
  float IpSig3DPos[100];
  float IpSig3DNeg[100];
  float IpSigZPos[100];
  float IpSigZNeg[100];  
  
  // impact parameters
  float Ip2DPos[100];
  float Ip2DNeg[100];
  float Ip3DPos[100];
  float Ip3DNeg[100];
  float IpZPos[100];
  float IpZNeg[100];  
  //for TMVA
 /* float IPxyPos=0;
  float IPxyNeg=0;
  float IPzPos=0;
  float IPzNeg=0;*/
  // distance between muon tracks
  float TwoMuonDist3D;
  float TwoMuonDist3DE;
  float TwoMuonDist2D;
  float TwoMuonDist2DE;
  float TwoMuonDistRPhi3D;
  float TwoMuonDistRPhi3DE;
  float TwoMuonDistRPhi2D;
  float TwoMuonDistRPhi2DE;


  float DiLeptonMass;


  float weight;

  //float XSection;

  std::ifstream fileList(argv[1]);
  
  std::string rootFileName(argv[1]);
  //bool PFNoPUIP=true;

  std::string chainName("DiMuonAnalysis/HTauTau2Muons");
  std::string chainNameGen("DiMuonAnalysis/HTauTau2Muons");

  TString TStrName(rootFileName);
  std::cout <<TStrName <<std::endl;  
  bool isData = false;
  if (TStrName.Contains("Data")) {
    isData = true;
   
    std::cout << "=============================" << std::endl;
    std::cout << "==== Running on data ========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;    
    //processId = 0;
  }
  else {
  
    std::cout << "=============================" << std::endl;
    std::cout << "== Running on Monte Carlo ===" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl; 
  }

  //bool isQCD = false;
  if (TStrName.Contains("QCD")) {
    //isQCD = true;
    std::cout << "=============================" << std::endl;
    std::cout << "====== QCD Monte Carlo ======" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    //processId = 5;
    //XSection = 84679;
  }
  
  //bool isM10to20 = false;
  bool isMadgraph = false;

  //bool isTTJets = false;
  if (TStrName.Contains("TTJets")) {
    std::cout << "=============================" << std::endl;
    std::cout << "========== TTJets =====-=====" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    //isTTJets = true;
    //processId = 6;
    //XSection = 158;
  }

  //bool isDYToTauTau = false;
  if (TStrName.Contains("DYToTauTau")||TStrName.Contains("ZToTauTau")) {
    //isDYToTauTau = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->tautau ======" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;
    //processId = 1;
    //XSection = 1667;
    if (TStrName.Contains("M10to20")) {
      //isM10to20 = true;
      //processId = 2;
      //XSection = 3893;
    }
    if (TStrName.Contains("Madgraph")||TStrName.Contains("madgraph")) {
      isMadgraph = true;
      //XSection = 3048;
    }
  }

  //bool isDYToMuMu = false;
  if (TStrName.Contains("DYToMuMu")||TStrName.Contains("ZToMuMu")) {
    //isDYToMuMu = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->mumu ========" << std::endl;
    std::cout << "=============================" << std::endl;
    std::cout << std::endl;    
    //processId = 3;
    //XSection = 1667;
    if (TStrName.Contains("M10to20")) {
      //isM10to20 = true;
      //processId = 4;
      //XSection = 3893;
    }
    if (TStrName.Contains("Madgraph")||TStrName.Contains("madgraph")) {
      isMadgraph = true;
      //XSection = 3048;
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

  TString FullName = TStrName+"_Trigger"+ApplyTrigger+"_"+Trigger;

  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");
  EventWeight * eventWeight = new EventWeight("/afs/desy.de/user/r/rasp/public/Analysis",0);
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  TH1F * counterTriggerH = new TH1F("counterTriggerH","",1,-0.5,0.5);
  TH1F * counterKineCutsH = new TH1F("counterKineCutsH","",1,-0.5,0.5);

  /////*******************
  ///// Efficiency Studies
  /////
  ///// Z, J/Psi Upsilon decay
  /////*******************
  
  // muon id
  TH1F * ZInvMassPosTag_PtEta[5][3][2]; 
  TH1F * ZInvMassNegTag_PtEta[5][3][2];
  
  TH1F * JPsiInvMassPosTag_PtEta[5][3][2];
  TH1F * JPsiInvMassNegTag_PtEta[5][3][2];

  TH1F * UpsilonInvMassPosTag_PtEta[5][3][2];
  TH1F * UpsilonInvMassNegTag_PtEta[5][3][2];

  // HLT_Mu17_Mu8 Leg 17
  TH1F * ZInvMassPosTag_PtEta_Leg17[5][3][2];
  TH1F * ZInvMassNegTag_PtEta_Leg17[5][3][2];

  TH1F * JPsiInvMassPosTag_PtEta_Leg17[5][3][2];
  TH1F * JPsiInvMassNegTag_PtEta_Leg17[5][3][2];

  TH1F * UpsilonInvMassPosTag_PtEta_Leg17[5][3][2];
  TH1F * UpsilonInvMassNegTag_PtEta_Leg17[5][3][2];

  // HLT_Mu17_Mu8 Leg 8
  TH1F * ZInvMassPosTag_PtEta_Leg8[5][3][2];
  TH1F * ZInvMassNegTag_PtEta_Leg8[5][3][2];

  TH1F * JPsiInvMassPosTag_PtEta_Leg8[5][3][2];
  TH1F * JPsiInvMassNegTag_PtEta_Leg8[5][3][2];

  TH1F * UpsilonInvMassPosTag_PtEta_Leg8[5][3][2];
  TH1F * UpsilonInvMassNegTag_PtEta_Leg8[5][3][2];


  // HLT_Mu17_Mu8 dZ filter
  TH1F * ZInvMassPosTag_PtEta_dZ[2];
  TH1F * ZInvMassNegTag_PtEta_dZ[2];

  TH1F * JPsiInvMassPosTag_PtEta_dZ[2];
  TH1F * JPsiInvMassNegTag_PtEta_dZ[2];

  TH1F * UpsilonInvMassPosTag_PtEta_dZ[2];
  TH1F * UpsilonInvMassNegTag_PtEta_dZ[2];



  int nPtBins = 5;
  int nEtaBins = 3;
  float ptBins[6] = {10, 15, 20, 25, 30, 10000};
  float etaBins[4] = {-0.1, 0.8, 1.6, 2.4};

  int nPtBinsTrigger = 5;
  int nEtaBinsTrigger = 3;
  float ptBinsTrigger[6] = {10, 14, 17, 22, 30, 10000};
  float etaBinsTrigger[4] = {-0.1, 0.8, 1.6, 2.1};


  TString ptBinNames[5] = {"10To15",
			   "15To20",
			   "20To25",
			   "25To30",
			   "30ToInf"};

  TString ptBinNamesTrigger[5] = {"10To14",
				  "14To17",
				  "17To22",
				  "22To30",
				  "30ToInf"};

  TString etaBinNames[3] = {"0To0p8",
			    "0p8To1p6",
			    "1p6To2p4"};

  TString etaBinNamesTrigger[3] = {"0To0p8",
				   "0p8To1p6",
				   "1p6To2p1"};


  TString statusNames[2] = {"passed",
			    "failed"};

  for (int iPt=0; iPt<nPtBins; ++iPt) {
    for (int iEta=0; iEta<nEtaBins; ++iEta) {
      for (int iStatus=0; iStatus<2; ++iStatus) {
	// muon id studies

	ZInvMassPosTag_PtEta[iPt][iEta][iStatus] = new TH1F("ZInvMassPosTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
	ZInvMassNegTag_PtEta[iPt][iEta][iStatus] = new TH1F("ZInvMassNegTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
	JPsiInvMassPosTag_PtEta[iPt][iEta][iStatus] = new TH1F("JPsiInvMassPosTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
	JPsiInvMassNegTag_PtEta[iPt][iEta][iStatus] = new TH1F("JPsiInvMassNegTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
	UpsilonInvMassPosTag_PtEta[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassPosTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);
        UpsilonInvMassNegTag_PtEta[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassNegTag_Pt"+ptBinNames[iPt]+"_Eta"+etaBinNames[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);

	// trigger studies
	// leg 17
	ZInvMassPosTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("ZInvMassPosTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
        ZInvMassNegTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("ZInvMassNegTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
        JPsiInvMassPosTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("JPsiInvMassPosTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
        JPsiInvMassNegTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("JPsiInvMassNegTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
        UpsilonInvMassPosTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassPosTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);
        UpsilonInvMassNegTag_PtEta_Leg17[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassNegTag_Leg17_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);
	// leg 8
	ZInvMassPosTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("ZInvMassPosTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
        ZInvMassNegTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("ZInvMassNegTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",300,60.,120.);
        JPsiInvMassPosTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("JPsiInvMassPosTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
        JPsiInvMassNegTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("JPsiInvMassNegTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",160,2.5,4.1);
        UpsilonInvMassPosTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassPosTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);
        UpsilonInvMassNegTag_PtEta_Leg8[iPt][iEta][iStatus] = new TH1F("UpsilonInvMassNegTag_Leg8_Pt"+ptBinNamesTrigger[iPt]+"_Eta"+etaBinNamesTrigger[iEta]+"_"+statusNames[iStatus],"",200,8.0,12.0);
      }
    }
  }
  
  // dz filter
  for (int iStatus=0; iStatus<2; ++iStatus) {

    ZInvMassPosTag_PtEta_dZ[iStatus] = new TH1F("ZInvMassPosTag_dZ_"+statusNames[iStatus],"",300,60.,120.);
    ZInvMassNegTag_PtEta_dZ[iStatus] = new TH1F("ZInvMassNegTag_dZ_"+statusNames[iStatus],"",300,60.,120.);

    JPsiInvMassPosTag_PtEta_dZ[iStatus] = new TH1F("JPsiInvMassPosTag_dZ_"+statusNames[iStatus],"",160,2.5,4.1);
    JPsiInvMassNegTag_PtEta_dZ[iStatus] = new TH1F("JPsiInvMassNegTag_dZ_"+statusNames[iStatus],"",160,2.5,4.1);

    UpsilonInvMassPosTag_PtEta_dZ[iStatus] = new TH1F("UpsilonInvMassPosTag_dZ_"+statusNames[iStatus],"",200,8.0,12.0);
    UpsilonInvMassNegTag_PtEta_dZ[iStatus] = new TH1F("UpsilonInvMassNegTag_dZ_"+statusNames[iStatus],"",200,8.0,12.0);

  }


  TH1F * ZInvMassPosTag_passed = new TH1F("ZInvMassPosTag_passed","",300,60.,120.);
  TH1F * ZInvMassNegTag_passed = new TH1F("ZInvMassNegTag_passed","",300,60.,120.);

  TH1F * JPsiInvMassPosTag_passed = new TH1F("JPsiInvMassPosTag_passed","",160,2.5,4.1);
  TH1F * JPsiInvMassNegTag_passed = new TH1F("JPsiInvMassNegTag_passed","",160,2.5,4.1);

  TH1F * UpsilonInvMassPosTag_passed = new TH1F("UpsilonInvMassPosTag_passed","",200,8.,12.);
  TH1F * UpsilonInvMassNegTag_passed = new TH1F("UpsilonInvMassNegTag_passed","",200,8.,12.);

  TH1F * ZInvMassPosTag_failed = new TH1F("ZInvMassPosTag_failed","",300,60.,120.);
  TH1F * ZInvMassNegTag_failed = new TH1F("ZInvMassNegTag_failed","",300,60.,120.);

  TH1F * JPsiInvMassPosTag_failed = new TH1F("JPsiInvMassPosTag_failed","",160,2.5,4.1);
  TH1F * JPsiInvMassNegTag_failed = new TH1F("JPsiInvMassNegTag_failed","",160,2.5,4.1);

  TH1F * UpsilonInvMassPosTag_failed = new TH1F("UpsilonInvMassPosTag_failed","",200,8.,12.);
  TH1F * UpsilonInvMassNegTag_failed = new TH1F("UpsilonInvMassNegTag_failed","",200,8.,12.);



  //TH1F * puiweightMC;
  //TH1F * puiweightData;
  
  std::vector<std::string> jsonFiles;
  jsonFiles.push_back("/afs/desy.de/user/r/rasp/public/Analysis/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");
 
  RunLumiSelector runLumiSelector;
  runLumiSelector = RunLumiSelector(jsonFiles);


  //bool PUInitialised=eventWeight->InitPUWeightsS10Truth(puWeightsFileData, puWeightsFileMc);


  int nFiles;
  fileList >> nFiles;

  for (int iFiles=0; iFiles<nFiles;++iFiles) {

   std::string filen;
   fileList >> filen;

   std::cout << iFiles+1 << " out of " << nFiles << "  ===> filename : " << filen << std::endl;
   TFile * file_ = TFile::Open(TString(filen));
   
   TChain * _tree = NULL;
   _tree = new TChain(TString(chainName));
  
   if (_tree==NULL) continue;
    
   TH1F * histoInputEvents = NULL;
   
   histoInputEvents = (TH1F*)file_->Get("eventCountPostTrigger/EventCount");

   if (histoInputEvents==NULL) continue;

   int NE = int(histoInputEvents->GetEntries());

   std::cout << "      number of input events = " << NE << std::endl;

   for (int iE=0;iE<NE;++iE)
     inputEventsH->Fill(0.);

   _tree->Add(TString(filen));
    

   _tree->SetBranchAddress("Run",&run);
   _tree->SetBranchAddress("Event",&event);
   _tree->SetBranchAddress("Lumi",&lumi);
   
   
   _tree->SetBranchAddress("nPV",&nPV);
   
   _tree->SetBranchAddress("nPUI",&nPUI);
   _tree->SetBranchAddress("nPUTruth",&nPUITruth);
   
   _tree->SetBranchAddress("chi2PV",chi2PV);
   _tree->SetBranchAddress("ndofPV",ndofPV);
   _tree->SetBranchAddress("probPV",probPV);
   _tree->SetBranchAddress("xPV",xPV);
   _tree->SetBranchAddress("yPV",yPV);
   _tree->SetBranchAddress("zPV",zPV);
   _tree->SetBranchAddress("nTrkPV",nTrkPV);
   _tree->SetBranchAddress("sumPtPV",sumPtPV);
   
   _tree->SetBranchAddress("isGlobalMuPos",&isGlobalMuPos);
   _tree->SetBranchAddress("isGlobalMuNeg",&isGlobalMuNeg);
   
   _tree->SetBranchAddress("isTrackerMuPos",&isTrackerMuPos);
   _tree->SetBranchAddress("isTrackerMuNeg",&isTrackerMuNeg);
   
   _tree->SetBranchAddress("isPFMuPos",&isPFMuPos);
   _tree->SetBranchAddress("isPFMuNeg",&isPFMuNeg);
   
   _tree->SetBranchAddress("nPixelHitsPos",&nPixelHitsPos);
   _tree->SetBranchAddress("nPixelHitsNeg",&nPixelHitsNeg);
   
   _tree->SetBranchAddress("nTrackerHitsPos",&nTrackerHitsPos);
   _tree->SetBranchAddress("nTrackerHitsNeg",&nTrackerHitsNeg);
   
   _tree->SetBranchAddress("nMuonHitsPos",&nMuonHitsPos);
   _tree->SetBranchAddress("nMuonHitsNeg",&nMuonHitsNeg);
   
   _tree->SetBranchAddress("nMuonStationsPos",&nMuonStationsPos);
   _tree->SetBranchAddress("nMuonStationsNeg",&nMuonStationsNeg);
   
   _tree->SetBranchAddress("Chi2PosMu",&Chi2PosMu);
   _tree->SetBranchAddress("Chi2NegMu",&Chi2NegMu);
   
   _tree->SetBranchAddress("NdofPosMu",&NdofPosMu);
   _tree->SetBranchAddress("NdofNegMu",&NdofNegMu);
   
    
   TChain * _genTree = NULL;
   if (isHto4taus){
     TChain * _genTree = NULL;
     _genTree = new TChain(TString(chainNameGen));

     if (_genTree==NULL) continue;
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

     _tree->SetBranchAddress("NPartFirstH1",&NPartFirstH1);
     _tree->SetBranchAddress("IdPartFirstH1",IdPartFirstH1);
     _tree->SetBranchAddress("QPartFirstH1", QPartFirstH1);
     _tree->SetBranchAddress("PxPartFirstH1",PxPartFirstH1);
     _tree->SetBranchAddress("PyPartFirstH1",PyPartFirstH1);
     _tree->SetBranchAddress("PzPartFirstH1",PzPartFirstH1);
     _tree->SetBranchAddress("PtPartFirstH1",PtPartFirstH1);
     _tree->SetBranchAddress("EtaPartFirstH1",EtaPartFirstH1);
     _tree->SetBranchAddress("PhiPartFirstH1",PhiPartFirstH1);

     _tree->SetBranchAddress("NPartSecondH1",&NPartSecondH1);
     _tree->SetBranchAddress("IdPartSecondH1",IdPartSecondH1);
     _tree->SetBranchAddress("QPartSecondH1", QPartSecondH1);
     _tree->SetBranchAddress("PxPartSecondH1",PxPartSecondH1);
     _tree->SetBranchAddress("PyPartSecondH1",PyPartSecondH1);
     _tree->SetBranchAddress("PzPartSecondH1",PzPartSecondH1);
     _tree->SetBranchAddress("PtPartSecondH1",PtPartSecondH1);
     _tree->SetBranchAddress("EtaPartSecondH1",EtaPartSecondH1);
     _tree->SetBranchAddress("PhiPartSecondH1",PhiPartSecondH1);

     // positive generated muon
     _tree->SetBranchAddress("posMuonGenExist",&posMuonGenExist);
     _tree->SetBranchAddress("posMuonGenpdg",&posMuonGenpdg);
     _tree->SetBranchAddress("posMuonGenQ",&posMuonGenQ);
     _tree->SetBranchAddress("posMuonGenpx",&posMuonGenpx);
     _tree->SetBranchAddress("posMuonGenpy",&posMuonGenpy);
     _tree->SetBranchAddress("posMuonGenpz",&posMuonGenpz);
     _tree->SetBranchAddress("posMuonGenpt",&posMuonGenpt);
     _tree->SetBranchAddress("posMuonGenphi",&posMuonGenphi);
     _tree->SetBranchAddress("posMuonGeneta",&posMuonGeneta);

     _tree->SetBranchAddress("PMmotherExist",&PMmotherExist);    
     _tree->SetBranchAddress("PMmotherpdg",&PMmotherpdg);
     _tree->SetBranchAddress("PMmotherphi",&PMmotherphi);
     _tree->SetBranchAddress("PMmothereta",&PMmothereta);
     _tree->SetBranchAddress("PMmotherpx",&PMmotherpx);
     _tree->SetBranchAddress("PMmotherpy",&PMmotherpy);
     _tree->SetBranchAddress("PMmotherpz",&PMmotherpz);
     _tree->SetBranchAddress("PMmotherpt",&PMmotherpt);
     
     _tree->SetBranchAddress("PMGrmotherExist",&PMGrmotherExist);
     _tree->SetBranchAddress("PMGrmotherpdg",&PMGrmotherpdg);
     _tree->SetBranchAddress("PMGrmotherphi",&PMGrmotherphi);
     _tree->SetBranchAddress("PMGrmothereta",&PMGrmothereta);
     _tree->SetBranchAddress("PMGrmotherpx",&PMGrmotherpx);
     _tree->SetBranchAddress("PMGrmotherpy",&PMGrmotherpy);
     _tree->SetBranchAddress("PMGrmotherpz",&PMGrmotherpz);
     _tree->SetBranchAddress("PMGrmotherpt",&PMGrmotherpt);
    

     // negative generated muon
     _tree->SetBranchAddress("negMuonGenExist",&negMuonGenExist);
     _tree->SetBranchAddress("negMuonGenpdg",&negMuonGenpdg);
     _tree->SetBranchAddress("negMuonGenQ",&negMuonGenQ);
     _tree->SetBranchAddress("negMuonGenpx",&negMuonGenpx);
     _tree->SetBranchAddress("negMuonGenpy",&negMuonGenpy);
     _tree->SetBranchAddress("negMuonGenpz",&negMuonGenpz);
     _tree->SetBranchAddress("negMuonGenpt",&negMuonGenpt);
     _tree->SetBranchAddress("negMuonGenphi",&negMuonGenphi);
     _tree->SetBranchAddress("negMuonGeneta",&negMuonGeneta);

     _tree->SetBranchAddress("NMmotherExist",&NMmotherExist);  
     _tree->SetBranchAddress("NMmotherpdg",&NMmotherpdg);
     _tree->SetBranchAddress("NMmotherphi",&NMmotherphi);
     _tree->SetBranchAddress("NMmothereta",&NMmothereta);
     _tree->SetBranchAddress("NMmotherpx",&NMmotherpx);
     _tree->SetBranchAddress("NMmotherpy",&NMmotherpy);
     _tree->SetBranchAddress("NMmotherpz",&NMmotherpz);
     _tree->SetBranchAddress("NMmotherpt",&NMmotherpt);
    
     _tree->SetBranchAddress("NMGrmotherExist",&NMGrmotherExist);
     _tree->SetBranchAddress("NMGrmotherpdg",&NMGrmotherpdg);
     _tree->SetBranchAddress("NMGrmotherphi",&NMGrmotherphi);
     _tree->SetBranchAddress("NMGrmothereta",&NMGrmothereta);
     _tree->SetBranchAddress("NMGrmotherpx",&NMGrmotherpx);
     _tree->SetBranchAddress("NMGrmotherpy",&NMGrmotherpy);
     _tree->SetBranchAddress("NMGrmotherpz",&NMGrmotherpz);
     _tree->SetBranchAddress("NMGrmotherpt",&NMGrmotherpt);
     
     _tree->SetBranchAddress("Boson",&_bosonPDG);
     _tree->SetBranchAddress("Decay",&_decayPDG);
     
     _tree->SetBranchAddress("ZExist",&_bosonExist); 
     _tree->SetBranchAddress("ZMass",&_ZMass);
     _tree->SetBranchAddress("ZPx",&_ZPx);
     _tree->SetBranchAddress("ZPy",&_ZPy);
     _tree->SetBranchAddress("ZPz",&_ZPz);
   }
   

   // triggers =====>
   _tree->SetBranchAddress("HLT_Mu17_Mu8",&HLT_Mu17_Mu8);
   _tree->SetBranchAddress("HLT_Mu17_TkMu8",&HLT_Mu17_TkMu8);
   _tree->SetBranchAddress("HLT_Mu24",&HLT_Mu24);
   _tree->SetBranchAddress("HLT_Mu40",&HLT_Mu40);
   _tree->SetBranchAddress("HLT_IsoMu24",&HLT_IsoMu24);

   // Trigger legs
   _tree->SetBranchAddress("HLT_Mu24_TrigMatch_Leg24_PosMuon",&PosLep_Mu24);
   _tree->SetBranchAddress("HLT_Mu24_TrigMatch_Leg24_NegMuon",&NegLep_Mu24);
   _tree->SetBranchAddress("HLT_IsoMu24_TrigMatch_Leg24_PosMuon",&PosLep_IsoMu24);
   _tree->SetBranchAddress("HLT_IsoMu24_TrigMatch_Leg24_NegMuon",&NegLep_IsoMu24);
   _tree->SetBranchAddress("HLT_Mu40_TrigMatch_Leg40_PosMuon",&PosLep_Mu40);
   _tree->SetBranchAddress("HLT_Mu40_TrigMatch_Leg40_NegMuon",&NegLep_Mu40);
   
   // DoubleMu triggers

   // HLT_Mu17_Mu8
   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_Leg17_PosMuon",&PosLep_Mu17_Mu8_Leg17);
   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_Leg17_NegMuon",&NegLep_Mu17_Mu8_Leg17);

   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_Leg8_PosMuon",&PosLep_Mu17_Mu8_Leg8);
   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_Leg8_NegMuon",&NegLep_Mu17_Mu8_Leg8);

   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_DzFilter_PosMuon",&PosLep_Mu17_Mu8_DzFilter);
   _tree->SetBranchAddress("HLT_Mu17_Mu8_TrigMatch_DzFilter_NegMuon",&NegLep_Mu17_Mu8_DzFilter);


   // HLT_Mu17_TkMu8
   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_Leg17_PosMuon",&PosLep_Mu17_TkMu8_Leg17);
   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_Leg17_NegMuon",&NegLep_Mu17_TkMu8_Leg17);

   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_Leg8_PosMuon",&PosLep_Mu17_TkMu8_Leg8);
   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_Leg8_NegMuon",&NegLep_Mu17_TkMu8_Leg8);

   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_DzFilter_PosMuon",&PosLep_Mu17_TkMu8_DzFilter);
   _tree->SetBranchAddress("HLT_Mu17_TkMu8_TrigMatch_DzFilter_NegMuon",&NegLep_Mu17_TkMu8_DzFilter);

   


   // collinear approximation =====>
   _tree->SetBranchAddress("ValidDiTau",&_validDiTau);
   _tree->SetBranchAddress("DiTauMass",&_diTauMass);
   _tree->SetBranchAddress("PosLepETauRF",&_posLepETauRest);
   _tree->SetBranchAddress("NegLepETauRF",&_negLepETauRest);
   _tree->SetBranchAddress("PosLepEDiTauRF",&_posLepEdiTauRest);
   _tree->SetBranchAddress("NegLepEDiTauRF",&_negLepEdiTauRest);
   _tree->SetBranchAddress("PosLepTauERatio",&_posLepTauERatio);
   _tree->SetBranchAddress("NegLepTauERatio",&_negLepTauERatio);
   _tree->SetBranchAddress("CosNegDiTau",&_cosNegDiTau );

   // di-tau kinematics from collinear approximation
   _tree->SetBranchAddress("DiTauPx",&_diTauPx);
   _tree->SetBranchAddress("DiTauPy",&_diTauPy);
   _tree->SetBranchAddress("DiTauPz",&_diTauPz);
   _tree->SetBranchAddress("DiTauEta",&_diTauEta);
   _tree->SetBranchAddress("DiTauPt",&_diTauPt);
   _tree->SetBranchAddress("DiTauE",&_diTauE);

   // dilepton kinematics
   _tree->SetBranchAddress("DiLepDPhi",&_diLepDPhi);
   _tree->SetBranchAddress("DiLepDEta",&_diLepDEta);
   _tree->SetBranchAddress("DiLepDR",&_diLepDR);
   _tree->SetBranchAddress("PosMETDPhi",&_posLepMETDPhi);
   _tree->SetBranchAddress("NegMETDPhi",&_negLepMETDPhi);
   _tree->SetBranchAddress("DiLepEta",&_diLepEta);

   _tree->SetBranchAddress("NegLepPt",&_negLepPt);
   _tree->SetBranchAddress("PosLepPt",&_posLepPt);
   _tree->SetBranchAddress("NegLepEta",&_negLepEta);
   _tree->SetBranchAddress("PosLepEta",&_posLepEta);
   _tree->SetBranchAddress("NegLepPhi",&_negLepPhi);
   _tree->SetBranchAddress("PosLepPhi",&_posLepPhi);
   _tree->SetBranchAddress("NegLepQ",&_negLepQ);
   _tree->SetBranchAddress("PosLepQ",&_posLepQ);

   _tree->SetBranchAddress("DiLepPt",&DiLepPt);
   _tree->SetBranchAddress("DiLepPhi",&DiLepPhi);
   _tree->SetBranchAddress("DiLeptonMass",&DiLeptonMass);

   
   // Tracks around Pos muon
   _tree->SetBranchAddress("NTracksPos",&NTracksPos);
   _tree->SetBranchAddress("PxTrackPos",PxTrackPos);
   _tree->SetBranchAddress("PyTrackPos",PyTrackPos);
   _tree->SetBranchAddress("PzTrackPos",PzTrackPos);
   _tree->SetBranchAddress("PtTrackPos",PtTrackPos);
   _tree->SetBranchAddress("EnergyTrackPos",EnergyTrackPos);
   _tree->SetBranchAddress("QTrackPos", QTrackPos);
   _tree->SetBranchAddress("IdTrackPos",IdTrackPos);

   // Tracks around Neg muon
   _tree->SetBranchAddress("NTracksNeg",&NTracksNeg);
   _tree->SetBranchAddress("PxTrackNeg",PxTrackNeg);
   _tree->SetBranchAddress("PyTrackNeg",PyTrackNeg);
   _tree->SetBranchAddress("PzTrackNeg",PzTrackNeg);
   _tree->SetBranchAddress("PtTrackNeg",PtTrackNeg);
   _tree->SetBranchAddress("EnergyTrackNeg",EnergyTrackNeg);
   _tree->SetBranchAddress("QTrackNeg",QTrackNeg);
   _tree->SetBranchAddress("IdTrackNeg",IdTrackNeg);

   //Tracks from PV with highest pt sum
   _tree->SetBranchAddress("NPVtracks",&NPVtracks);
   _tree->SetBranchAddress("PtTrack",PtTrack);
   _tree->SetBranchAddress("PxTrack",PxTrack);
   _tree->SetBranchAddress("PyTrack",PyTrack);
   _tree->SetBranchAddress("PzTrack",PzTrack);
   _tree->SetBranchAddress("EtaTrack",EtaTrack);
   _tree->SetBranchAddress("ChargeTrack",ChargeTrack);
   _tree->SetBranchAddress("Chi2Track",Chi2Track);


   //PF objects from the No PU collection
   _tree->SetBranchAddress("NpfNoPU",&NpfNoPU);
   _tree->SetBranchAddress("PxpfNoPU",PxpfNoPU);
   _tree->SetBranchAddress("PypfNoPU",PypfNoPU);
   _tree->SetBranchAddress("PzpfNoPU",PzpfNoPU);
   _tree->SetBranchAddress("PtpfNoPU",PtpfNoPU);
   _tree->SetBranchAddress("EnpfNoPU",EnpfNoPU);
   _tree->SetBranchAddress("EtapfNoPU",EtapfNoPU);
   _tree->SetBranchAddress("ChargepfNoPU",ChargepfNoPU);
   _tree->SetBranchAddress("IppfNoPU",ippfNoPU);
   _tree->SetBranchAddress("IpSigpfNoPU",ipSigpfNoPU);
   _tree->SetBranchAddress("IpZpfNoPU",ipZpfNoPU);
   _tree->SetBranchAddress("IpZSigpfNoPU",ipZSigpfNoPU);
   _tree->SetBranchAddress("Ip3DpfNoPU",ip3DpfNoPU);
   _tree->SetBranchAddress("Ip3DSigpfNoPU",ip3DSigpfNoPU);
   _tree->SetBranchAddress("IdpfNoPU",idpfNoPU);


   // MET variables
   _tree->SetBranchAddress("MET",&met);
   _tree->SetBranchAddress("METPx",&metPx);
   _tree->SetBranchAddress("METPy",&metPy);
   _tree->SetBranchAddress("METPz",&metPz);
   _tree->SetBranchAddress("METEn",&metEn);
   _tree->SetBranchAddress("METCovXX",&metCovXX);
   _tree->SetBranchAddress("METCovXY",&metCovXY);
   _tree->SetBranchAddress("METCovYX",&metCovYX);
   _tree->SetBranchAddress("METCovYY",&metCovYY);
  
   // MVA MET variables
   _tree->SetBranchAddress("METMVAPx",&metMVAPx);
   _tree->SetBranchAddress("METMVAPy",&metMVAPy);
   _tree->SetBranchAddress("METMVAPz",&metMVAPz);
   _tree->SetBranchAddress("METMVAEn",&metMVAEn);
   _tree->SetBranchAddress("METMVACovXX",&metMVACovXX);
   _tree->SetBranchAddress("METMVACovXY",&metMVACovXY);
   _tree->SetBranchAddress("METMVACovYX",&metMVACovYX);
   _tree->SetBranchAddress("METMVACovYY",&metMVACovYY);
  
   // isolation variables
   _tree->SetBranchAddress("IsoNeg",&IsoNeg);
   _tree->SetBranchAddress("IsoPos",&IsoPos);
   _tree->SetBranchAddress("IsoTrkNeg",&IsoTrkNeg);
   _tree->SetBranchAddress("IsoTrkPos",&IsoTrkPos);
   _tree->SetBranchAddress("IsoECalNeg",&IsoECalNeg);
   _tree->SetBranchAddress("IsoECalPos",&IsoECalPos);
   _tree->SetBranchAddress("IsoHCalNeg",&IsoHCalNeg);
   _tree->SetBranchAddress("IsoHCalPos",&IsoHCalPos);
   _tree->SetBranchAddress("IsoPFsNeg",&IsoPFsNeg);
   _tree->SetBranchAddress("IsoPFsPos",&IsoPFsPos);

   // distance between muons 
   _tree->SetBranchAddress("TwoMuonDist3D",&TwoMuonDist3D);
   _tree->SetBranchAddress("TwoMuonDist3DE",&TwoMuonDist3DE);
   _tree->SetBranchAddress("TwoMuonDist2D",&TwoMuonDist2D);
   _tree->SetBranchAddress("TwoMuonDist2DE",&TwoMuonDist2DE);
   _tree->SetBranchAddress("TwoMuonDistRPhi3D",&TwoMuonDistRPhi3D);
   _tree->SetBranchAddress("TwoMuonDistRPhi3DE",&TwoMuonDistRPhi3DE);
   _tree->SetBranchAddress("TwoMuonDistRPhi2D",&TwoMuonDistRPhi2D);
   _tree->SetBranchAddress("TwoMuonDistRPhi2DE",&TwoMuonDistRPhi2DE);

   // impact parameter significance
   _tree->SetBranchAddress("IpSig2DPos",IpSig2DPos);
   _tree->SetBranchAddress("IpSig2DNeg",IpSig2DNeg);
   _tree->SetBranchAddress("IpSig3DPos",IpSig3DPos);
   _tree->SetBranchAddress("IpSig3DNeg",IpSig3DNeg);
   _tree->SetBranchAddress("IpSigZPos",IpSigZPos);
   _tree->SetBranchAddress("IpSigZNeg",IpSigZNeg);  
  
   // impact parmeaters
   _tree->SetBranchAddress("Ip2DPos",Ip2DPos);
   _tree->SetBranchAddress("Ip2DNeg",Ip2DNeg);
   _tree->SetBranchAddress("Ip3DPos",Ip3DPos);
   _tree->SetBranchAddress("Ip3DNeg",Ip3DNeg);
   _tree->SetBranchAddress("IpZPos",IpZPos);
   _tree->SetBranchAddress("IpZNeg",IpZNeg);  

   _tree->SetBranchAddress("noJets",&_nJets);
   _tree->SetBranchAddress("jetIDLoose",jetIDLoose);
   _tree->SetBranchAddress("jetPt",jetPt);
   _tree->SetBranchAddress("jetEta",jetEta);
   _tree->SetBranchAddress("jetPhi",jetPhi);
   _tree->SetBranchAddress("combSVBJetTag",combSVBJetTag);
   _tree->SetBranchAddress("puJetFullLoose", jetPUIDLoose);
   _tree->SetBranchAddress("chargedMultiplicity",chargedMultiplicity);

    
   int numberOfCandidates = _tree->GetEntries();
  
   std::cout << "      number of selected candidates = " << numberOfCandidates << std::endl;
   
   for (int iCand=0; iCand<numberOfCandidates; iCand++) { 
    
      _tree->GetEntry(iCand);

      if (isData && !runLumiSelector.accept(run, lumi))
	{
	  continue;
	}

      //----------------------
      // weight initialization 
      //---------------------- 
      weight = 1;

      // Vertex selection
      int iPVtx=0;
      float sumPtMax = -1;
      
      for (int iPV=0; iPV<nPV;++iPV) {
	float dVtx = TMath::Sqrt(xPV[iPV]*xPV[iPV]+yPV[iPV]*yPV[iPV]);
	if (ndofPV[iPV]>4.0 && probPV[iPV]>0.01 && dVtx < 2.0 && TMath::Abs(zPV[iPV])<24.0) {
	  if (sumPtPV[iPV]>sumPtMax) {
	    sumPtMax = sumPtPV[iPV];
	    iPVtx = iPV;
	  }
	  
	}
      }

      // cut on Z vertex;
      if (TMath::Abs(zPV[iPVtx])>24.0) continue;

      if (randomizeMuons) {

	float posLepPhiPlusPi2 = _posLepPhi + 0.5*TMath::Pi();
	float unitX = _posLepPt * TMath::Cos( posLepPhiPlusPi2 );
	float unitY = _posLepPt * TMath::Sin( posLepPhiPlusPi2 );
	float pxNeg = _negLepPt * TMath::Cos( _negLepPhi );
	float pyNeg = _negLepPt * TMath::Sin( _negLepPhi ); 
	float prod  = pxNeg*unitX + pyNeg*unitY;

	if (prod<0.) { // swap muons

	  bool tmpBool = isGlobalMuPos;
	  isGlobalMuPos = isGlobalMuNeg; 
	  isGlobalMuNeg = tmpBool; 

	  tmpBool = isTrackerMuPos;
	  isTrackerMuPos = isTrackerMuNeg;
	  isTrackerMuNeg = tmpBool;

	  tmpBool = isPFMuPos;
	  isPFMuPos = isPFMuNeg;
	  isPFMuNeg = tmpBool;

	  int tmpInt = nPixelHitsPos;
	  nPixelHitsPos = nPixelHitsNeg;
	  nPixelHitsNeg = tmpInt;

	  tmpInt = nTrackerHitsPos;
	  nTrackerHitsPos = nTrackerHitsNeg;
	  nTrackerHitsNeg = tmpInt;
	  
	  tmpInt = nMuonHitsPos;
	  nMuonHitsPos = nMuonHitsNeg;
	  nMuonHitsNeg = tmpInt;

	  tmpInt = nMuonStationsPos;
	  nMuonStationsPos = nMuonStationsNeg;
	  nMuonStationsNeg = tmpInt;

	  float tmpFloat = Chi2PosMu;
	  Chi2PosMu = Chi2NegMu;
	  Chi2NegMu = tmpFloat;

	  tmpFloat = NdofPosMu;
	  NdofPosMu = NdofNegMu;
	  NdofNegMu = tmpFloat;


	  tmpFloat = _posLepPt;
	  _posLepPt = _negLepPt;
	  _negLepPt = tmpFloat;

	  tmpFloat = _posLepEta;
	  _posLepEta = _negLepEta;
	  _negLepEta = tmpFloat;

	  tmpFloat = _posLepPhi;
	  _posLepPhi = _negLepPhi;
	  _negLepPhi = tmpFloat;

	  tmpFloat = _posLepQ;
	  _posLepQ = _negLepQ;
	  _negLepQ = tmpFloat;

	  tmpFloat = IpSig2DPos[iPVtx];
	  IpSig2DPos[iPVtx] = IpSig2DNeg[iPVtx];
	  IpSig2DNeg[iPVtx] = tmpFloat;

	  tmpFloat = IpSig3DPos[iPVtx];
	  IpSig3DPos[iPVtx] = IpSig3DNeg[iPVtx];
	  IpSig3DNeg[iPVtx] = tmpFloat;
	  
	  tmpFloat = IpSigZPos[iPVtx];
	  IpSigZPos[iPVtx] = IpSigZNeg[iPVtx];
	  IpSigZNeg[iPVtx] = tmpFloat;
	  

	  tmpFloat = Ip2DPos[iPVtx];
	  Ip2DPos[iPVtx] = Ip2DNeg[iPVtx];
	  Ip2DNeg[iPVtx] = tmpFloat;

	  tmpFloat = Ip3DPos[iPVtx];
	  Ip3DPos[iPVtx] = Ip3DNeg[iPVtx];
	  Ip3DNeg[iPVtx] = tmpFloat;
	  
	  tmpFloat = IpZPos[iPVtx];
	  IpZPos[iPVtx] = IpZNeg[iPVtx];
	  IpZNeg[iPVtx] = tmpFloat;
	  
	}

      }

      TLorentzVector MuonPos4;
      MuonPos4.SetPtEtaPhiM(_posLepPt,_posLepEta,_posLepPhi,muonMass);
      TLorentzVector MuonNeg4;
      MuonNeg4.SetPtEtaPhiM(_negLepPt,_negLepEta,_negLepPhi,muonMass);

      if (applyPUReweighting && !isData)
	{
	  float ratioF = 1;
	  ratioF = eventWeight->PUWeightS10Truth(float(reweightWithPUI ? nPUITruth : nPV));
	  //  std::cout << "apply PU weight: "<< ratioF <<std::endl;
	  weight *= ratioF;
	}
	
     
      //Trigger selection
      bool firedTrigger = HLT_Mu24;
      bool posLepTriggerMatch = PosLep_Mu24;
      bool negLepTriggerMatch = NegLep_Mu24;
      if (triggerType==1) {
	firedTrigger = HLT_Mu40;
	posLepTriggerMatch = PosLep_Mu40;
        negLepTriggerMatch = NegLep_Mu40;
      }
      if (triggerType==2) {
	firedTrigger = HLT_IsoMu24;
        posLepTriggerMatch = PosLep_IsoMu24;
        negLepTriggerMatch = NegLep_IsoMu24;
      }

      if (applyTrigger && !firedTrigger) continue;

      counterTriggerH->Fill(1.);

      float absPosMuEta = TMath::Abs(_posLepEta);
      float absNegMuEta = TMath::Abs(_negLepEta);

      if (_posLepPt<ptProbeLepCut) continue;
      if (_negLepPt<ptProbeLepCut) continue;

      if (absPosMuEta>etaProbeLepCut) continue;
      if (absNegMuEta>etaProbeLepCut) continue;

      counterKineCutsH->Fill(1.);

      float sumPtPos = 0;
      float sumPtNeg = 0;

      int nAroundPosAll = 0;
      int nAroundPosSig = 0;

      int nAroundNegAll = 0;
      int nAroundNegSig = 0;


      for (int itrack=0; itrack<NpfNoPU; itrack++){
        if (PtpfNoPU[itrack]>ptTrkCut && (ChargepfNoPU[itrack]>0.5 ||ChargepfNoPU[itrack]<-0.5) && TMath::Abs(EtapfNoPU[itrack])<etaTrkCut) {
	  if ( TMath::Abs(ipZpfNoPU[itrack])<dzCutLoose && TMath::Abs(ippfNoPU[itrack])<d0CutLoose) {
	    TLorentzVector pfCand;
	    pfCand.SetXYZM(PxpfNoPU[itrack],PypfNoPU[itrack],PzpfNoPU[itrack],pionMass);
	    float EtaPfCand = pfCand.Eta();
	    float PhiPfCand = pfCand.Phi();
	    float dRpos = float(deltaR(_posLepEta,_posLepPhi,EtaPfCand,PhiPfCand));
	    TLorentzVector diffPos = MuonPos4 - pfCand;
	    if (dRpos<0.5 && diffPos.P()>0.1) {
	      sumPtPos += PtpfNoPU[itrack];
	      nAroundPosAll++;
	      if (PtpfNoPU[itrack]>ptSoftTrkCut && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks && TMath::Abs(ippfNoPU[itrack])<d0CutTracks)
		nAroundPosSig++;
	    }
	    float dRneg = float(deltaR(_negLepEta,_negLepPhi,EtaPfCand,PhiPfCand));
            TLorentzVector diffNeg = MuonNeg4 - pfCand;
	    if (dRneg<0.5 && diffNeg.P()>0.1) {
              sumPtNeg += PtpfNoPU[itrack];
	      nAroundNegAll++;
	      if (PtpfNoPU[itrack]>ptSoftTrkCut && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks && TMath::Abs(ippfNoPU[itrack])<d0CutTracks)
                nAroundNegSig++;
	      
	    }
	  }
	}
      }

      float relIsoPos = sumPtPos / _posLepPt;
      float relIsoNeg = sumPtNeg / _negLepPt;


      //######### Beginning of cuts ###########################################

      // Cosmics rejection
      //      float deltaPt = 2*fabs(_posLepPt-_negLepPt)/(_posLepPt+_negLepPt);
      //      float deltaEta = fabs(_posLepEta+_negLepEta);
      //      if ( (deltaPt<0.25) && (deltaEta<0.002) && (_diLepDPhi>3.14)) continue;
      
      ////////////////

      TLorentzVector MuonPos4PlusMuonNeg4;
      MuonPos4PlusMuonNeg4=MuonPos4+MuonNeg4;
      float InvMassDiMuon = float(MuonPos4PlusMuonNeg4.M());

      ///////////////

      /////////// Efficiency Studies ////////////

      int iEtaPos = eventWeight->binNumber(absPosMuEta,nEtaBins,etaBins);
      int iPtPos  = eventWeight->binNumber(_posLepPt,nPtBins,ptBins);

      int iEtaNeg = eventWeight->binNumber(absNegMuEta,nEtaBins,etaBins);
      int iPtNeg  = eventWeight->binNumber(_negLepPt,nPtBins,ptBins);

      int iEtaPosTrigger = eventWeight->binNumber(absPosMuEta,nEtaBinsTrigger,etaBinsTrigger);
      int iPtPosTrigger  = eventWeight->binNumber(_posLepPt,nPtBinsTrigger,ptBinsTrigger);
      
      int iEtaNegTrigger = eventWeight->binNumber(absNegMuEta,nEtaBinsTrigger,etaBinsTrigger);
      int iPtNegTrigger  = eventWeight->binNumber(_negLepPt,nPtBinsTrigger,ptBinsTrigger);
      
//       std::cout << "neg lep pt = " << _negLepPt << "  neg lep eta = " << absNegMuEta << std::endl;
//       std::cout << "iPt = " << iPtNegTrigger << "  iEta = " << iEtaNegTrigger << std::endl;

//       std::cout << "HLT_Mu17_Mu8 = " << HLT_Mu17_Mu8 << "  HLT_Mu17_TkMu8 = " << HLT_Mu17_TkMu8 <<std::endl;

//       std::cout << "Neg : HLT_Mu17_Mu8_Leg17 = " << NegLep_Mu17_Mu8_Leg17 << 
// 	"  HLT_Mu17_TkMu8_Leg17 = " << NegLep_Mu17_TkMu8_Leg17 << std::endl;

//       std::cout << "Neg : HLT_Mu17_Mu8_Leg8 = " << NegLep_Mu17_Mu8_Leg8 << 
// 	"  HLT_Mu17_TkMu8_Leg8 = " << NegLep_Mu17_TkMu8_Leg8 << std::endl;

//       std::cout << "Pos : HLT_Mu17_Mu8_Leg17 = " << PosLep_Mu17_Mu8_Leg17 << 
// 	"  HLT_Mu17_TkMu8_Leg17 = " << PosLep_Mu17_TkMu8_Leg17 << std::endl;

//       std::cout << "Pos : HLT_Mu17_Mu8_Leg8 = " << PosLep_Mu17_Mu8_Leg8 << 
// 	"  HLT_Mu17_TkMu8_Leg8 = " << PosLep_Mu17_TkMu8_Leg8 << std::endl;
      

//       std::cout << "Pos : HLT_Mu17_Mu8_dZ = " << PosLep_Mu17_Mu8_DzFilter  << 
// 	"  HLT_Mu17_TkMu8_dZ = " << PosLep_Mu17_TkMu8_DzFilter  << std::endl;

//       std::cout << "Neg : HLT_Mu17_Mu8_dZ = " << NegLep_Mu17_Mu8_DzFilter  << 
// 	"  HLT_Mu17_TkMu8_dZ = " << NegLep_Mu17_TkMu8_DzFilter  << std::endl;

//       std::cout << std::endl;

      // and now apply muon Id --->
      
      bool PosMuId = true;
      bool NegMuId = true;
	
      if (!isGlobalMuPos) PosMuId = false;
      if (!isGlobalMuNeg) NegMuId = false;  // do not apply muon id on neg mu which is probe leg
      
      if(!isPFMuPos) PosMuId = false;
      if(!isPFMuNeg) NegMuId = false;   // do not apply muon id on neg mu which is probe leg
      
      float chi2ndfPos = Chi2PosMu/NdofPosMu; 
      float chi2ndfNeg = Chi2NegMu/NdofNegMu; 
	
      if (chi2ndfPos>10) PosMuId = false;
      if (chi2ndfNeg>10) NegMuId = false;
      
      if (nMuonHitsPos<1) PosMuId = false;
      if (nMuonHitsNeg<1) NegMuId = false;

      if (nMuonStationsPos<2) PosMuId = false;
      if (nMuonStationsNeg<2) NegMuId = false;	

      //	if (Ip2DPos[iPVtx]>0.02||Ip2DPos[iPVtx]<-0.02) continue;
      //	if (Ip2DNeg[iPVtx]>0.02||Ip2DNeg[iPVtx]<-0.02) continue;
      // 	if (IpZPos[iPVtx]>0.2||IpZPos[iPVtx]<-0.2) continue; 
      // 	if (IpZNeg[iPVtx]>0.2||IpZNeg[iPVtx]<-0.2) continue; 
      
      if (TMath::Abs(Ip2DPos[iPVtx])>d0Cut) PosMuId = false;
      if (TMath::Abs(Ip2DNeg[iPVtx])>d0Cut) NegMuId = false;
      
      if (TMath::Abs(IpZPos[iPVtx])>dzCut) PosMuId = false;
      if (TMath::Abs(IpZNeg[iPVtx])>dzCut) NegMuId = false;
      
      
      if (nPixelHitsPos<1) PosMuId = false;
      if (nPixelHitsNeg<1) NegMuId = false;
	
      // 	if (nTrackerHitsPos<11)  continue; 
      // 	if (nTrackerHitsNeg<11)  continue; 

      if (nTrackerHitsPos<5) PosMuId = false;
      if (nTrackerHitsNeg<5) NegMuId = false;

      if (nAroundPosAll>0) PosMuId = false;
      if (nAroundNegAll>0) NegMuId = false;

      bool dZfilter = (PosLep_Mu17_Mu8_DzFilter && NegLep_Mu17_Mu8_DzFilter) || 
	(PosLep_Mu17_TkMu8_DzFilter && NegLep_Mu17_TkMu8_DzFilter);


      // pos muon tag
      if ( PosMuId && _posLepPt>ptTagLepCut && absPosMuEta<etaTagLepCut && posLepTriggerMatch) {
	if (NegMuId) {
	  if (relIsoPos<relIsoCut) {

	    ZInvMassPosTag_PtEta[iPtNeg][iEtaNeg][0]->Fill(InvMassDiMuon,weight);
	    ZInvMassPosTag_passed->Fill(InvMassDiMuon,weight);

	    if (NegLep_Mu17_Mu8_Leg17||NegLep_Mu17_TkMu8_Leg17)
	      ZInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
	    else 
	      ZInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);

	    if (NegLep_Mu17_Mu8_Leg8||NegLep_Mu17_TkMu8_Leg8)
              ZInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
            else
              ZInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);

	    if (dZfilter) 
	      ZInvMassPosTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
	    else
	      ZInvMassPosTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);

	  }
	  JPsiInvMassPosTag_PtEta[iPtNeg][iEtaNeg][0]->Fill(InvMassDiMuon,weight);
	  UpsilonInvMassPosTag_PtEta[iPtNeg][iEtaNeg][0]->Fill(InvMassDiMuon,weight);
          JPsiInvMassPosTag_passed->Fill(InvMassDiMuon,weight);
          UpsilonInvMassPosTag_passed->Fill(InvMassDiMuon,weight);

	  if (NegLep_Mu17_Mu8_Leg17||NegLep_Mu17_TkMu8_Leg17) {
	    JPsiInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
	    UpsilonInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
	  }
	  else {
	    JPsiInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassPosTag_PtEta_Leg17[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);
	  }

	  if (NegLep_Mu17_Mu8_Leg8||NegLep_Mu17_TkMu8_Leg8) {
            JPsiInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][0]->Fill(InvMassDiMuon,weight);
          }
          else {
            JPsiInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassPosTag_PtEta_Leg8[iPtNegTrigger][iEtaNegTrigger][1]->Fill(InvMassDiMuon,weight);
          }

	  if (dZfilter) {
	    JPsiInvMassPosTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassPosTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
	  }
	  else {
	    JPsiInvMassPosTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassPosTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);
	  }


	}
	else {
	  if (relIsoPos<relIsoCut) {
	    ZInvMassPosTag_PtEta[iPtNeg][iEtaNeg][1]->Fill(InvMassDiMuon,weight);
	    ZInvMassPosTag_failed->Fill(InvMassDiMuon,weight);
	  }
          JPsiInvMassPosTag_PtEta[iPtNeg][iEtaNeg][1]->Fill(InvMassDiMuon,weight);
          UpsilonInvMassPosTag_PtEta[iPtNeg][iEtaNeg][1]->Fill(InvMassDiMuon,weight);
          JPsiInvMassPosTag_failed->Fill(InvMassDiMuon,weight);
          UpsilonInvMassPosTag_failed->Fill(InvMassDiMuon,weight);
	}
      }

      // neg muon tag 
      if ( NegMuId && _negLepPt>ptTagLepCut && absNegMuEta<etaTagLepCut && negLepTriggerMatch) {
        if (PosMuId) {
	  if (relIsoNeg<relIsoCut) {

	    ZInvMassNegTag_PtEta[iPtPos][iEtaPos][0]->Fill(InvMassDiMuon,weight);
	    ZInvMassNegTag_passed->Fill(InvMassDiMuon,weight);


            if (PosLep_Mu17_Mu8_Leg17||PosLep_Mu17_TkMu8_Leg17)
              ZInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
            else
              ZInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);

            if (PosLep_Mu17_Mu8_Leg8||PosLep_Mu17_TkMu8_Leg8)
              ZInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
            else
              ZInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);

            if (dZfilter)
              ZInvMassNegTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
            else
              ZInvMassNegTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);


	  }

	  JPsiInvMassNegTag_PtEta[iPtPos][iEtaPos][0]->Fill(InvMassDiMuon,weight);
          UpsilonInvMassNegTag_PtEta[iPtPos][iEtaPos][0]->Fill(InvMassDiMuon,weight);
          JPsiInvMassNegTag_passed->Fill(InvMassDiMuon,weight);
          UpsilonInvMassNegTag_passed->Fill(InvMassDiMuon,weight);

	  if (PosLep_Mu17_Mu8_Leg17||PosLep_Mu17_TkMu8_Leg17) {
            JPsiInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
          }
          else {
            JPsiInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_Leg17[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);
          }

          if (PosLep_Mu17_Mu8_Leg8||PosLep_Mu17_TkMu8_Leg8) {
            JPsiInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][0]->Fill(InvMassDiMuon,weight);
          }
          else {
            JPsiInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_Leg8[iPtPosTrigger][iEtaPosTrigger][1]->Fill(InvMassDiMuon,weight);
          }

	  if (dZfilter) {
            JPsiInvMassNegTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_dZ[0]->Fill(InvMassDiMuon,weight);
	  }
	  else {
	    JPsiInvMassNegTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);
            UpsilonInvMassNegTag_PtEta_dZ[1]->Fill(InvMassDiMuon,weight);
	  }

        }
        else {
	  if (relIsoNeg<relIsoCut) {
	    ZInvMassNegTag_PtEta[iPtPos][iEtaPos][1]->Fill(InvMassDiMuon,weight);
	    ZInvMassNegTag_failed->Fill(InvMassDiMuon,weight);
	  }
	  JPsiInvMassNegTag_PtEta[iPtPos][iEtaPos][1]->Fill(InvMassDiMuon,weight);
          UpsilonInvMassNegTag_PtEta[iPtPos][iEtaPos][1]->Fill(InvMassDiMuon,weight);
          JPsiInvMassNegTag_failed->Fill(InvMassDiMuon,weight);
          UpsilonInvMassNegTag_failed->Fill(InvMassDiMuon,weight);
        } 
      }


   }
    
    
   delete _tree;
   if(!isData) delete _genTree;
   file_->Close();
  }
  
  nProcessed = int(inputEventsH->GetEntries());
  std::cout <<nProcessed<< std::endl;
  
  file->cd("");
  
  file->Write();
 
  file->Close();
  
  delete file;
 
}

