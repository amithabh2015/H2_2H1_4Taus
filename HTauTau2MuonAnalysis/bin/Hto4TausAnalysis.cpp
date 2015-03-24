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
const float pzeroMass = 0.134976;
const float kzeroMass = 0.497614;
const float kaonMass = 0.4937;
const float tauMass  = 1.7768;

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

  int randomizeMuons = 2;

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

  float isoCone = 0.4;

  float ptTrkCut = 1.;
  float ptTrkCutIso = 0.0;

  float etaTrkCut = 99999;
  float ptGenPartCut = 0.5;
  // loose cuts
  float d0Cut = 0.03;
  float dzCut = 0.1;
  // tight cuts
  //  float d0Cut = 0.02;
  //  float dzCut = 0.04;

  // standard
  float dzCutLoose = 1.0;
  float d0CutLoose = 1.0;
  // loose
  //  float dzCutLoose = 0.5;
  //  float d0CutLoose = 0.5;
  //  Super loose
  //  float dzCutLoose = 0.4;
  //  float d0CutLoose = 0.2;
  // tight cuts
  float d0CutTracks= 0.02; 
  float dzCutTracks= 0.04; 
  // looser conditions 
  //  float d0CutTracks= 0.03; 
  //  float dzCutTracks= 0.1;
  
  float d0CutLowerTracks = 0.02;
  float dzCutLowerTracks = 0.04;
  float dRmumuCut = 2.0;

  float dRcone = 0.5;

  float jetEtaCut = 4.5;
  float bJetEtaCut = 2.4;
  float btagCSVL = 0.244;
  float btagCSVM = 0.679;

  float muMomScale = 0.0;
  //float trkMomScale = 0.0;

  
  if (argc<10) {
    std::cout << "Usage of the program : Hto4TausAnalysis [file_list] [dR(mu-mu)] [Ip scale up] [Ip scale down] [NSoftTracks] [applyPtReweigh] [ptScale] [applySF]" << std::endl;
    std::cout << "file_list : file list of RooT files to be processed. To run on Data the string has to include the string \"Data\"." << std::endl;
    std::cout << "Type of cut        : dR or dPhi" << std::endl;
    std::cout << "dR/dPhi(mu-mu) cut : 0.5, 1.0 1.5, 2.0, 2.5 " << std::endl;
    std::cout << "Ip scale up        : 1, 2, 3, 4, 5" << std::endl;
    std::cout << "Ip scale down      : -1, 0.2, 0.5, 1" << std::endl;
    std::cout << "N Soft Tracks      : 1, 2, 3, 4" << std::endl;
    std::cout << "Apply Pt Reweighting : True, true, false, False, 0 or 1" << std::endl;
    std::cout << "Scale for the Pt Reweighing: nominal, up or down" << std::endl;
    std::cout << "Apply Muon Id and Iso SF : True, true, False, false" << std::endl;
    exit(1);
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
  int qPVTracksPos[10000];
  
  // tracks around Neg muon
  int PVtracksDRneg;
  float ptPVTracksNeg[10000];
  float pxPVTracksNeg[10000];
  float pyPVTracksNeg[10000];
  float pzPVTracksNeg[10000];
  int qPVTracksNeg[10000];*/

  // pfNoPu around Pos muon
  int NpfNoPuPos;
  //float PtPfNoPuPos[10000];
  /*float PxPfNoPuPos[10000];
  float PyPfNoPuPos[10000];
  float PzPfNoPuPos[10000];
  float EnPfNoPuPos[10000];
  float IPPfNoPuPos[10000];
  float IPZPfNoPuPos[10000];
  int QPfNoPuPos[10000];*/
  
  // pfNoPu around Neg muon
  int NpfNoPuNeg;
  /*float PtPfNoPuNeg[10000];
  float PxPfNoPuNeg[10000];
  float PyPfNoPuNeg[10000];
  float PzPfNoPuNeg[10000];
  float EnPfNoPuNeg[10000];
  float IPPfNoPuNeg[10000];
  float IPZPfNoPuNeg[10000];
  int QPfNoPuNeg[10000];*/

  //  pfNoPu around Pos muon (Ip cuts)
  int NpfNoPuIpPos;
  int NpfNoPuIpPosSoft;
  int NpfNoPuIpPosAll;
  float PtPfNoPuIpPos[10000];
  float PxPfNoPuIpPos[10000];
  float PyPfNoPuIpPos[10000];
  float PzPfNoPuIpPos[10000];
  //float EnPfNoPuIpPos[10000];
  float IPPfNoPuIpPos[10000];
  float IPZPfNoPuIpPos[10000];
  int QPfNoPuIpPos[10000];
  int IdNoPuIpPos[10000];

  // pfNoPu around Neg muon (Ip cuts)
  int NpfNoPuIpNeg;
  int NpfNoPuIpNegSoft;
  int NpfNoPuIpNegAll;
  float PtPfNoPuIpNeg[10000];
  float PxPfNoPuIpNeg[10000];
  float PyPfNoPuIpNeg[10000];
  float PzPfNoPuIpNeg[10000];
  //float EnPfNoPuIpNeg[10000];
  float IPPfNoPuIpNeg[10000];
  float IPZPfNoPuIpNeg[10000];
  int QPfNoPuIpNeg[10000];
  int IdNoPuIpNeg[10000];

  //PF objects from the No PU collection
  int NpfNoPU;
  int AllNpfNoPU;
  int AllNpfNoPU_ip;
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
  int Ntracks_OutConePos;
  int Ntracks_OutConeNeg;
  int Ntracks_OutConePosIP;
  int Ntracks_OutConeNegIP;
  float SumPtPosCone;
  float SumPtNegCone;
  float SumPtPosConeIP;
  float SumPtNegConeIP;
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
  
  float InvMassH2;
  float InvMassTrackPlusPosMuon;
  float InvMassTrackPlusNegMuon;

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
  bool NegLep_Mu20;
  bool PosLep_Mu20;
  bool NegLep_Mu30;
  bool PosLep_Mu30;
  bool NegLep_Mu40;
  bool PosLep_Mu40;
  bool NegLep_Mu17_Mu8;
  bool PosLep_Mu17_Mu8;
  bool NegLep_Mu17_TkMu8;
  bool PosLep_Mu17_TkMu8;
  
  bool HLT_Mu17_Mu8;
  bool HLT_Mu17_TkMu8;

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
  float jetPx[50];
  float jetPy[50];
  float jetPz[50];
  float jetEn[50];
  float jetMass[50];

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
  /*float IPxyPos=0;
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

  float weight;

  //float XSection;

  std::ifstream fileList(argv[1]);
  
  int nFiles;

  fileList >> nFiles;

  std::string rootFileName(argv[1]);
  TString cutTypeStr(argv[2]);
  std::string dRcut(argv[3]);
  std::string RelaxedIp(argv[4]);
  std::string IpDown(argv[5]);  
  std::string NSoftTracks(argv[6]);
  std::string PtReweighting(argv[7]);
  std::string PtReweightingScale(argv[8]);
  std::string ApplyMuonSF(argv[9]);

  bool applyPTReweighting =true;
  bool isPtFileValid=false;
  string scalein="nominal";

  bool cutType = true;
  if (cutTypeStr=="dR")
    cutType = true;
  else if (cutTypeStr=="dPhi")
    cutType = false;
  else {
    std::cout << "cut type should be dR or dPhi" << std::endl;
    exit(-1);
  }


  TH1F  *PythiaH2PtWeights=new TH1F("PythiaH2PtWeights", "PythiaH2PtWeights",30, 0, 300);;
  TFile *filein;
  TString fileName("/afs/desy.de/user/r/rasp/public/Analysis/HPt_Spectrum.root");
  char scale_[100];
  float bin=5.;
	
  if (PtReweighting=="True" || PtReweighting=="true" || PtReweighting=="1") 
    {
      if (PtReweightingScale=="nominal" || PtReweightingScale=="Nominal" || PtReweightingScale=="up" || PtReweightingScale=="Up" || PtReweightingScale ==	"down" || PtReweightingScale=="Down") 
	{ 
	  std::transform(PtReweightingScale.begin(), PtReweightingScale.end(), PtReweightingScale.begin(), ::tolower);
	  scalein=PtReweightingScale;
		        
	  cout<<" The scale was initiated properly.... I will use the "<<scalein<<" scale "<<endl;}
      
      else {cout<<" The scale was not initiated properly.... I will use the "<<scalein<<" scale "<<endl;}
	
      //if (applyPTReweighting){
      //TFile *filein = TFile::Read(fileName);
      filein = TFile::Open(fileName,"read");
      if (filein->IsZombie()) {
	cout << "Error opening Reweighting H2 Pt spectrum file....Will use weights =1..." << endl;
	//exit(-1);
      }
      else {
	cout << "Succesfully opened Reweighting H2 Pt spectrum file...." << endl;
	sprintf(scale_,"scale_%s_step%g/Weights",scalein.c_str() ,bin);
	
	PythiaH2PtWeights = (TH1F*) filein->Get(scale_);
	
	//PythiaH2PtWeights = (TH1F*) filein->Get("scale_nominal_step5/Weights");
	isPtFileValid=true;
      }
      //}
    }
  if (PtReweighting=="False" || PtReweighting=="false" || PtReweighting=="0") applyPTReweighting=false;

  TFile * muonEffFile = new TFile("/afs/desy.de/user/r/rasp/public/Analysis/muEff.root"); 
  TH1F * muonIdSF = (TH1F*)muonEffFile->Get("muEff");

  if (ApplyMuonSF=="True"||ApplyMuonSF=="true")
    applyMuonIdSF = true;
  else if (ApplyMuonSF=="False"||ApplyMuonSF=="false")
    applyMuonIdSF = false;
  else {
    std::cout << "applySF should be True or False..."  << std::endl;
    exit(1);
  }


  if (NSoftTracks=="1")
    NsoftTracks = 1;
  else if (NSoftTracks=="2")
    NsoftTracks = 2;
  else if (NSoftTracks=="3")
    NsoftTracks = 3;
  else if (NSoftTracks=="4")
    NsoftTracks = 4;
  else {
    std::cout << "Undefined number of soft tracks : " << NSoftTracks << std::endl;
    exit(-1);
  }

  float relaxedIp = 1;
  if (RelaxedIp=="1")
    relaxedIp = 1;
  else if (RelaxedIp=="2")
    relaxedIp = 2;
  else if (RelaxedIp=="3")
    relaxedIp = 3;
  else if (RelaxedIp=="4")
    relaxedIp = 4;
  else if (RelaxedIp=="5")
    relaxedIp = 5;
  else {
    std::cout << "Undefined scale of track ip (upper) : " << RelaxedIp << std::endl;
    exit(-1);
  }
  d0CutTracks *= relaxedIp;
  dzCutTracks *= relaxedIp;

  float ipDown = -1;
  if (IpDown=="-1")
    ipDown = -1;
  else if (IpDown=="0.1")
    ipDown = 0.1;
  else if (IpDown=="0.2")
    ipDown = 0.2;
  else if (IpDown=="0.5")
    ipDown = 0.5;
  else if (IpDown=="1")
    ipDown = 1;
  else {
    std::cout << "Undefined scale of track ip (lower) : " << IpDown << std::endl;
    exit(-1);
  }
  d0CutLowerTracks *= ipDown;
  dzCutLowerTracks *= ipDown;

  bool ipAndLogic = false;

  if (ipDown>0&&ipDown<0.4)
    ipAndLogic = true;

  if (dRcut=="0"||dRcut=="0.0")
    dRmumuCut = -0.01;
  else if (dRcut=="0.5")
    dRmumuCut = 0.5;
  else if (dRcut=="1"||dRcut=="1.0")
    dRmumuCut = 1.0;
  else if (dRcut=="1.5")
    dRmumuCut = 1.5;
  else if (dRcut=="2"||dRcut=="2.0")
    dRmumuCut = 2.0;
  else if (dRcut=="2.5")
    dRmumuCut = 2.5;
  else {
    std::cout << "Undefined cut on dR(mu-mu) : " << dRcut << std::endl;
    exit(-1);
  }

  std::string chainName("DiMuonAnalysisSS/HTauTau2Muons");
  std::string chainNameGen("DiMuonAnalysisSS/HTauTau2Muons");

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
    std::cout << "========== TTJets ===========" << std::endl;
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

  //bool isDYToLL = false;
  if (TStrName.Contains("DYToLL")||TStrName.Contains("ZToLL")) {
    //isDYToLL = true;
    
    std::cout << "=============================" << std::endl;
    std::cout << "===== Z/gamma*->l+l- ========" << std::endl;
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

  
  TString FullName = TStrName +"_"+cutTypeStr+dRcut+"_IpUp"+RelaxedIp+"_IpDown"+IpDown+"_NSoft"+NSoftTracks;

  if  (PtReweighting=="True" || PtReweighting=="true" || PtReweighting=="1") FullName += "_Scale"+PtReweightingScale;
  else {FullName += "_NoPtRwgt";}
      
  
  FullName = FullName + "_MuIdScale" + ApplyMuonSF;

  TFile * file = new TFile(FullName+TString(".root"),"recreate");

  file->cd("");
  //  EventWeight * eventWeight = new EventWeight("/afs/desy.de/user/a/alkaloge/xxl-af-cms/Higgs/CMSSW_5_3_9/src/H2to2H1to4Taus/HTauTau2MuonAnalysis/RooT",0);
  EventWeight * eventWeight = new EventWeight("/afs/desy.de/user/r/rasp/public/Analysis",0); 
  TH1F * inputEventsH = new TH1F("inputEventsH","",1,-0.5,0.5);
  
  TH1F *counterInverseH=new TH1F("counterInverseH","",1,0.,2.);
  TH1F *counterInverse2H=new TH1F("counterInverse2H","",1,0.,2.);
  TH1F *counterInverse3H=new TH1F("counterInverse3H","",1,0.,2.);
  TH1F *counterInverse4H=new TH1F("counterInverse4H","",1,0.,2.);

  TH1F *counterSwapQH=new TH1F("counterSwapQH","",1,0.,2.);

  TH1F *counterDirectH=new TH1F("counterDirectH","",1,0.,2.);
  TH1F *counterNegCutsH=new TH1F("counterNegCutsH","",1,0.,2.);
  TH1F *counterPosCutsH=new TH1F("counterPosCutsH","",1,0.,2.);

  TH1F *counterPreIsoFinalH= new TH1F("counterPreIsoFinalH","",1,0.,2.);
  TH1F *counterPostIsoFinalH= new TH1F("counterPostIsoFinalH","",1,0.,2.);
  TH1F *counterFinalH=new TH1F("counterFinalH","",1,0.,2.);
  TH1F *counterSemiFinalH=new TH1F("counterSemiFinalH","",1,0.,2.);
  TH1F *countermuonDRH=new TH1F("countermuonDRH","",1,0.,2.);
  TH1F *counterPtH=new TH1F("counterPtH","",1,0.,2.);
  TH1F *counterEtaH=new TH1F("counterEtaH","",1,0.,2.);
  TH1F *counterMuonIDH=new TH1F("counterMuonIDH","",1,0.,2.);
  TH1F *counterTriggerH=new TH1F("counterTriggerH","",1,0.,2.);

  int nPtTauBins = 5;
  float ptTauBins[6] = {1,2.5,5,10,20,100};

  int nPtFullTauBins = 6;
  float ptFullTauBins[7] = {5,10,15,20,30,50,100};
  
  // === efficiency and fakes (track)

  TH1F * counterPtPosMuTrkH    = new TH1F("counterPtPosMuTrkH","",nPtTauBins,ptTauBins);
  TH1F * counterPtPosMuTrkIdH  = new TH1F("counterPtPosMuTrkIdH","",nPtTauBins,ptTauBins);
  TH1F * counterPtPosMuTrkIsoH = new TH1F("counterPtPosMuTrkIsoH","",nPtTauBins,ptTauBins);
  TH1F * counterPtPosMuTrkNonLepH = new TH1F("counterPtPosMuTrkNonLepH","",nPtTauBins,ptTauBins);

  TH1F * counterPtNegMuTrkH    = new TH1F("counterPtNegMuTrkH","",nPtTauBins,ptTauBins);
  TH1F * counterPtNegMuTrkIdH  = new TH1F("counterPtNegMuTrkIdH","",nPtTauBins,ptTauBins);
  TH1F * counterPtNegMuTrkIsoH = new TH1F("counterPtNegMuTrkIsoH","",nPtTauBins,ptTauBins);
  TH1F * counterPtNegMuTrkNonLepH = new TH1F("counterPtNegMuTrkNonLepH","",nPtTauBins,ptTauBins);

  TH1F * counterPtTrkH    = new TH1F("counterPtTrkH","",nPtTauBins,ptTauBins);
  TH1F * counterPtTrkIdH  = new TH1F("counterPtTrkIdH","",nPtTauBins,ptTauBins);
  TH1F * counterPtTrkIsoH = new TH1F("counterPtTrkIsoH","",nPtTauBins,ptTauBins);
  TH1F * counterPtTrkNonLepH = new TH1F("counterPtTrkNonLepH","",nPtTauBins,ptTauBins);

  // === efficiency and fakes (tau)

  TH1F * counterPtPosMuTauH    = new TH1F("counterPtPosMuTauH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtPosMuTauIdH  = new TH1F("counterPtPosMuTauIdH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtPosMuTauIsoH = new TH1F("counterPtPosMuTauIsoH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtPosMuTauNonLepH = new TH1F("counterPtPosMuTauNonLepH","",nPtFullTauBins,ptFullTauBins);

  TH1F * counterPtNegMuTauH    = new TH1F("counterPtNegMuTauH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtNegMuTauIdH  = new TH1F("counterPtNegMuTauIdH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtNegMuTauIsoH = new TH1F("counterPtNegMuTauIsoH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtNegMuTauNonLepH = new TH1F("counterPtNegMuTauNonLepH","",nPtFullTauBins,ptFullTauBins);

  TH1F * counterPtTauH    = new TH1F("counterPtTauH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtTauIdH  = new TH1F("counterPtTauIdH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtTauIsoH = new TH1F("counterPtTauIsoH","",nPtFullTauBins,ptFullTauBins);
  TH1F * counterPtTauNonLepH = new TH1F("counterPtTauNonLepH","",nPtFullTauBins,ptFullTauBins);

  // === efficiency and fakes (track and tau)

  TH2F * counterPtPosMuTrkTauH = new TH2F("counterPtPosMuTrkTauH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins); 
  TH2F * counterPtPosMuTrkTauIdH = new TH2F("counterPtPosMuTrkTauIdH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtPosMuTrkTauIsoH = new TH2F("counterPtPosMuTrkTauIsoH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtPosMuTrkTauNonLepH = new TH2F("counterPtPosMuTrkTauNonLepH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);

  TH2F * counterPtNegMuTrkTauH = new TH2F("counterPtNegMuTrkTauH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtNegMuTrkTauIdH = new TH2F("counterPtNegMuTrkTauIdH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtNegMuTrkTauIsoH = new TH2F("counterPtNegMuTrkTauIsoH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtNegMuTrkTauNonLepH = new TH2F("counterPtNegMuTrkTauNonLepH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);

  TH2F * counterPtTrkTauH = new TH2F("counterPtTrkTauH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtTrkTauIdH = new TH2F("counterPtTrkTauIdH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtTrkTauIsoH = new TH2F("counterPtTrkTauIsoH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);
  TH2F * counterPtTrkTauNonLepH = new TH2F("counterPtTrkTauNonLepH","",nPtTauBins,ptTauBins,nPtFullTauBins,ptFullTauBins);

  TH1F * tauPtHardDimuonH = new TH1F("tauPtHardDimuonH","",20,0,100);
  TH1F * tauPtSoftDimuonH = new TH1F("tauPtSoftDimuonH","",20,0,100);

  TH1F * tauPtHardDimuonIdH = new TH1F("tauPtHardDimuonIdH","",20,0,100);
  TH1F * tauPtSoftDimuonIdH = new TH1F("tauPtSoftDimuonIdH","",20,0,100);

  TH1F * tauPtHardDimuonIsoH = new TH1F("tauPtHardDimuonIsoH","",20,0,100);
  TH1F * tauPtSoftDimuonIsoH = new TH1F("tauPtSoftDimuonIsoH","",20,0,100);


  TH1F * trkPtHardDimuonH = new TH1F("trkPtHardDimuonH","",100,0,50);
  TH1F * trkPtSoftDimuonH = new TH1F("trkPtSoftDimuonH","",100,0,50);

  TH1F * trkPtHardDimuonIdH = new TH1F("trkPtHardDimuonIdH","",100,0,50);
  TH1F * trkPtSoftDimuonIdH = new TH1F("trkPtSoftDimuonIdH","",100,0,50);

  TH1F * trkPtHardDimuonIsoH = new TH1F("trkPtHardDimuonIsoH","",100,0,50);
  TH1F * trkPtSoftDimuonIsoH = new TH1F("trkPtSoftDimuonIsoH","",100,0,50);




  // --------

  TH1F * counterPreIsolationH = new TH1F("counterPreIsolationH","",1,0.,2.);
  TH1F * counterPostIsolationH = new TH1F("counterPostIsolationH","",1,0.,2.);
  
  TH1F * hardMuIsoTracksH = new TH1F("hardMuIsoTracksH","",11,-0.5,10.5);
  TH1F * softMuIsoTracksH = new TH1F("softMuIsoTracksH","",11,-0.5,10.5);
  
  TH1F * posMuIsoTracksH = new TH1F("posMuIsoTracksH","",11,-0.5,10.5);
  TH1F * negMuIsoTracksH = new TH1F("negMuIsoTracksH","",11,-0.5,10.5);

  

  //
  // Muon impact parameters after muon selection ->
  //
  TH1F * IpSig2DPosH=new TH1F("IpSig2DPosH","",50,-5.,5.);
  TH1F * IpSig2DNegH=new TH1F("IpSig2DNegH","",50,-5.,5.);
  TH2F * IpSig2D_2DH=new TH2F("IpSig2D_2DH","",50,-5.,5.,50,-5.,5.);
  TH1F * IpSig3DPosH=new TH1F("IpSig3DPosH","",50,-5.,5.);
  TH1F * IpSig3DNegH=new TH1F("IpSig3DNegH","",50,-5.,5.);
  TH1F * IpSigZPosH=new TH1F("IpSigZPosH","",50,-5.,5.);
  TH1F * IpSigZNegH=new TH1F("IpSigZNegH","",50,-5.,5.);

  TH1F * Ip2DPosH=new TH1F("Ip2DPosH","",100,-0.5,0.5);
  TH1F * Ip2DNegH=new TH1F("Ip2DNegH","",100,-0.5,0.5);
  TH1F * Ip2DBothH=new TH1F("Ip2DBothH","",100,-0.5,0.5);
  TH1F * IpZPosH=new TH1F("IpZPosH","",200,-1,1);
  TH1F * IpZNegH=new TH1F("IpZNegH","",200,-1,1);
  TH1F * IpZBothH=new TH1F("IpZBothH","",200,-1.,1.);

  //
  // Signal tracks --->
  // 
  // transverse impact parameters
  TH1F * Ip2DPosTrackH = new TH1F("Ip2DPosTrackH","",100,-0.5,0.5);
  TH1F * Ip2DNegTrackH = new TH1F("Ip2DNegTrackH","",100,-0.5,0.5);
  TH1F * Ip2DBothTrackH = new TH1F("Ip2DBothTrackH","",100,-0.5,0.5);
  
  TH1F * Ip2DPosTrackSideBandH = new TH1F("Ip2DPosTrackSideBandH","",100,-0.5,0.5);
  TH1F * Ip2DNegTrackSideBandH = new TH1F("Ip2DNegTrackSideBandH","",100,-0.5,0.5);
  TH1F * Ip2DBothTrackSideBandH = new TH1F("Ip2DBothTrackSideBandH","",100,-0.5,0.5);
  
  // longitudinal impact parameters
  TH1F * IpZPosTrackH = new TH1F("IpZPosTrackH","",200,-1.,1.);
  TH1F * IpZNegTrackH = new TH1F("IpZNegTrackH","",200,-1.,1.);
  TH1F * IpZBothTrackH = new TH1F("IpZBothTrackH","",200,-1.,1.);
  
  TH1F * IpZPosTrackSideBandH = new TH1F("IpZPosTrackSideBandH","",200,-1.,1.);
  TH1F * IpZNegTrackSideBandH = new TH1F("IpZNegTrackSideBandH","",200,-1.,1.);
  TH1F * IpZBothTrackSideBandH = new TH1F("IpZBothTrackSideBandH","",200,-1.,1.);

  // track pt 
  TH1F * PtPosTrackH = new TH1F("PtPosTrackH","",500,0.,50.);
  TH1F * PtNegTrackH = new TH1F("PtNegTrackH","",500,0.,50.);
  TH1F * PtBothTrackH = new TH1F("PtBothTrackH","",500,0.,50.);

  // track pt with Ip cuts 
  TH1F * PtPosTrackIpCutsH = new TH1F("PtPosTrackIpCutsH","",500,0.,50.);
  TH1F * PtNegTrackIpCutsH = new TH1F("PtNegTrackIpCutsH","",500,0.,50.);
  TH1F * PtBothTrackIpCutsH = new TH1F("PtBothTrackIpCutsH","",500,0.,50.);
  
  // track pt (side band)
  TH1F * PtPosTrackSideBandH = new TH1F("PtPosTrackSideBandH","",500,0.,50.);
  TH1F * PtNegTrackSideBandH = new TH1F("PtNegTrackSideBandH","",500,0.,50.);
  TH1F * PtBothTrackSideBandH = new TH1F("PtBothTrackSideBandH","",500,0.,50.);
  
  
  // track pt with Ip cuts (side band)
  TH1F * PtPosTrackIpCutsSideBandH = new TH1F("PtPosTrackIpCutsSideBandH","",500,0.,50.);
  TH1F * PtNegTrackIpCutsSideBandH = new TH1F("PtNegTrackIpCutsSideBandH","",500,0.,50.);
  TH1F * PtBothTrackIpCutsSideBandH = new TH1F("PtBothTrackIpCutsSideBandH","",500,0.,50.);

  // ************ Signal region kinematics
  TH1F * PosMuonPtH = new TH1F("PosMuonPtH","",100,0.,50.);
  TH1F * NegMuonPtH = new TH1F("NegMuonPtH","",100,0.,50.);

  TH1F * PosMuonEtaH = new TH1F("PosMuonEtaH","",50,-2.5,2.5);
  TH1F * NegMuonEtaH = new TH1F("NegMuonEtaH","",50,-2.5,2.5);

  TH1F * PosTrackPtH = new TH1F("PosTrackPtH","",100,0,50.0);
  TH1F * NegTrackPtH = new TH1F("NegTrackPtH","",100,0,50.0);
  
  TH1F * PosTrackEtaH = new TH1F("PosTrackEtaH","",50,-2.5,2.5);
  TH1F * NegTrackEtaH = new TH1F("NegTrackEtaH","",50,-2.5,2.5);

  TH1F * dRPosMuonTrackH = new TH1F("dRPosMuonTrackH","",40,0,2);
  TH1F * dRNegMuonTrackH = new TH1F("dRNegMuonTrackH","",40,0,2);

  TH1F * dRSoftMuonTrackH = new TH1F("dRSoftMuonTrackH","",40,0,2);
  TH1F * dRHardMuonTrackH = new TH1F("dRHardMuonTrackH","",40,0,2);



  // ********** Side band with non-isolated mu-trk pair *********

  TH1F * PosmuonPt_2QH = new TH1F("PosmuonPt_2QH","",100,0,50.0);
  TH1F * NegmuonPt_2QH = new TH1F("NegmuonPt_2QH","",100,0,50.0);
  //TH1F * PosmuonPt_2H = new TH1F("PosmuonPt_2H","",100,0,50.0);
  //TH1F * NegmuonPt_2H = new TH1F("NegmuonPt_2H","",100,0,50.0);
  TH1F * PosmuonPt_3QH = new TH1F("PosmuonPt_3QH","",100,0,50.0);
  TH1F * NegmuonPt_3QH = new TH1F("NegmuonPt_3QH","",100,0,50.0);
  //TH1F * PosmuonPt_3H = new TH1F("PosmuonPt_3H","",100,0,50.0);
  //TH1F * NegmuonPt_3H = new TH1F("NegmuonPt_3H","",100,0,50.0);
  TH1F * PosmuonPt_4QH = new TH1F("PosmuonPt_4QH","",100,0,50.0);
  TH1F * NegmuonPt_4QH = new TH1F("NegmuonPt_4QH","",100,0,50.0);
  //TH1F * PosmuonPt_4H = new TH1F("PosmuonPt_4H","",100,0,50.0);
  //TH1F * NegmuonPt_4H = new TH1F("NegmuonPt_4H","",100,0,50.0);

  TH1F * PosmuonPt_InvQH = new TH1F("PosmuonPt_InvQH","",100,0,50.0);
  TH1F * NegmuonPt_InvQH = new TH1F("NegmuonPt_InvQH","",100,0,50.0);
  //TH1F * PosmuonPt_InvH = new TH1F("PosmuonPt_InvH","",100,0,50.0);
  //TH1F * NegmuonPt_InvH = new TH1F("NegmuonPt_InvH","",100,0,50.0);

  TH1F * HardmuonPt_InvQH = new TH1F("HardmuonPt_InvQH","",100,0,50.0);
  TH1F * SoftmuonPt_InvQH = new TH1F("SoftmuonPt_InvQH","",100,0,50.0);

  //TH2F * HardSoftmuonPt_InvQH = new TH2F("HardSoftmuonPt_InvQH","",100,0,50.0,100,0,50.0); 

  TH1F * PosmuonEta_2QH = new TH1F("PosmuonEta_2QH","",50,-2.5,2.5);
  TH1F * PosmuonEta_3QH = new TH1F("PosmuonEta_3QH","",50,-2.5,2.5);
  TH1F * PosmuonEta_4QH = new TH1F("PosmuonEta_4QH","",50,-2.5,2.5);
  TH1F * PosmuonEta_InvQH = new TH1F("PosmuonEta_InvQH","",50,-2.5,2.5);

  TH1F * NegmuonEta_2QH = new TH1F("NegmuonEta_2QH","",50,-2.5,2.5);
  TH1F * NegmuonEta_3QH = new TH1F("NegmuonEta_3QH","",50,-2.5,2.5);
  TH1F * NegmuonEta_4QH = new TH1F("NegmuonEta_4QH","",50,-2.5,2.5);
  TH1F * NegmuonEta_InvQH = new TH1F("NegmuonEta_InvQH","",50,-2.5,2.5);

  TH1F * HardmuonEta_InvQH = new TH1F("HardmuonEta_InvQH","",50,-2.5,2.5);
  TH1F * SoftmuonEta_InvQH = new TH1F("SoftmuonEta_InvQH","",50,-2.5,2.5);

  TH1F * counterPos_InvQH = new TH1F("counterPos_InvQH","",1,0.,2.);
  TH1F * counterNeg_InvQH = new TH1F("counterNeg_InvQH","",1,0.,2.);

  TH1F * counterPos_2QH = new TH1F("counterPos_2QH","",1,0.,2.);
  TH1F * counterNeg_2QH = new TH1F("counterNeg_2QH","",1,0.,2.);

  TH1F * counterPos_3QH = new TH1F("counterPos_3QH","",1,0.,2.);
  TH1F * counterNeg_3QH = new TH1F("counterNeg_3QH","",1,0.,2.);

  TH1F * counterPos_4QH = new TH1F("counterPos_4QH","",1,0.,2.);
  TH1F * counterNeg_4QH = new TH1F("counterNeg_4QH","",1,0.,2.);

  TH1F * counterPos_SwapQH = new TH1F("counterPos_SwapQH","",1,0.,2.);
  TH1F * counterNeg_SwapQH = new TH1F("counterNeg_SwapQH","",1,0.,2.);


  //TH1F * PostrackPtH= new TH1F("PostrackPtH","",100,0,50.0);
  //TH1F * NegtrackPtH= new TH1F("NegtrackPtH","",100,0,50.0);

  //TH1F * PostrackPt_2H= new TH1F("PostrackPt_2H","",100,0,50.0);
  //TH1F * NegtrackPt_2H= new TH1F("NegtrackPt_2H","",100,0,50.0);
  TH1F * PostrackPt_2QH= new TH1F("PostrackPt_2QH","",100,0,50.0);
  TH1F * NegtrackPt_2QH= new TH1F("NegtrackPt_2QH","",100,0,50.0);
  //TH1F * PostrackPt_3H= new TH1F("PostrackPt_3H","",100,0,50.0);
  //TH1F * NegtrackPt_3H= new TH1F("NegtrackPt_3H","",100,0,50.0);
  TH1F * PostrackPt_3QH= new TH1F("PostrackPt_3QH","",100,0,50.0);
  TH1F * NegtrackPt_3QH= new TH1F("NegtrackPt_3QH","",100,0,50.0);
  //TH1F * PostrackPt_4H= new TH1F("PostrackPt_4H","",100,0,50.0);
  //TH1F * NegtrackPt_4H= new TH1F("NegtrackPt_4H","",100,0,50.0);
  TH1F * PostrackPt_4QH= new TH1F("PostrackPt_4QH","",100,0,50.0);
  TH1F * NegtrackPt_4QH= new TH1F("NegtrackPt_4QH","",100,0,50.0);
  //TH1F * PostrackPt_InvH= new TH1F("PostrackPt_InvH","",100,0,50.0);
  //TH1F * NegtrackPt_InvH= new TH1F("NegtrackPt_InvH","",100,0,50.0);
  TH1F * PostrackPt_InvQH= new TH1F("PostrackPt_InvQH","",100,0,50.0);
  TH1F * NegtrackPt_InvQH= new TH1F("NegtrackPt_InvQH","",100,0,50.0);

  TH1F * PostrackEta_2QH = new TH1F("PostrackEta_2QH","",50,-2.5,2.5);
  TH1F * PostrackEta_3QH = new TH1F("PostrackEta_3QH","",50,-2.5,2.5);
  TH1F * PostrackEta_4QH = new TH1F("PostrackEta_4QH","",50,-2.5,2.5);
  TH1F * PostrackEta_InvQH = new TH1F("PostrackEta_InvQH","",50,-2.5,2.5);

  TH1F * NegtrackEta_2QH = new TH1F("NegtrackEta_2QH","",50,-2.5,2.5);
  TH1F * NegtrackEta_3QH = new TH1F("NegtrackEta_3QH","",50,-2.5,2.5);
  TH1F * NegtrackEta_4QH = new TH1F("NegtrackEta_4QH","",50,-2.5,2.5);
  TH1F * NegtrackEta_InvQH = new TH1F("NegtrackEta_InvQH","",50,-2.5,2.5);

  //
  //TH1F * HardtrackPt_InvQH= new TH1F("HardtrackPt_InvQH","",100,0,50.0);
  //TH1F * SofttrackPt_InvQH= new TH1F("SofttrackPt_InvQH","",100,0,50.0);

  //TH2F * HardSofttrackPt_InvQH= new TH2F("HardSofttrackPt_InvQH","",100,0,50.0,100,0,50.);

  //TH1F * HardtrackEta_InvQH = new TH1F("HardtrackEta_InvQH","",50,-2.5,2.5);
  //TH1F * SofttrackEta_InvQH = new TH1F("SofttrackEta_InvQH","",50,-2.5,2.5);

  //
  TH1F * HardtrackMuPt_InvQH= new TH1F("HardtrackMuPt_InvQH","",100,0,50.0);
  TH1F * SofttrackMuPt_InvQH= new TH1F("SofttrackMuPt_InvQH","",100,0,50.0);

  //TH2F * HardSofttrackMuPt_InvQH= new TH2F("HardSofttrackMuPt_InvQH","",100,0,50.0,100,0,50.);

  TH1F * HardtrackMuEta_InvQH = new TH1F("HardtrackMuEta_InvQH","",50,-2.5,2.5);
  TH1F * SofttrackMuEta_InvQH = new TH1F("SofttrackMuEta_InvQH","",50,-2.5,2.5);


  TH1F * dRPosmuonTrk_2QH = new TH1F("dRPosmuonTrk_2QH","",40,0,2);
  TH1F * dRPosmuonTrk_3QH = new TH1F("dRPosmuonTrk_3QH","",40,0,2);
  TH1F * dRPosmuonTrk_4QH = new TH1F("dRPosmuonTrk_4QH","",40,0,2);
  TH1F * dRPosmuonTrk_InvQH = new TH1F("dRPosmuonTrk_InvQH","",40,0,2);

  TH1F * dRNegmuonTrk_2QH = new TH1F("dRNegmuonTrk_2QH","",40,0,2);
  TH1F * dRNegmuonTrk_3QH = new TH1F("dRNegmuonTrk_3QH","",40,0,2);
  TH1F * dRNegmuonTrk_4QH = new TH1F("dRNegmuonTrk_4QH","",40,0,2);
  TH1F * dRNegmuonTrk_InvQH = new TH1F("dRNegmuonTrk_InvQH","",40,0,2);

  TH1F * dRHardmuonTrk_InvQH = new TH1F("dRHardmuonTrk_InvQH","",40,0,2);
  TH1F * dRSoftmuonTrk_InvQH = new TH1F("dRSoftmuonTrk_InvQH","",40,0,2);


  TH1F * ippfNoPUH = new TH1F("ippfNoPUH","",100,-0.5,0.5);
  TH1F * ipZpfNoPUH = new TH1F("ipZpfNoPUH","",200,-1,1);
  
  //TH1F * sumIpAllH = new TH1F("sumIpAllH","",100,0,2);
  //TH1F * sumIpZAllH = new TH1F("sumIpZAllH","",100,0,4);


  /*TH1F * diLepPtH = new TH1F("diLepPtH","",100,0,500.0);
  TH1F * dPhiPosLepMetH = new TH1F("dPhiPosLepMetH","",50,0.,TMath::Pi());
  TH1F * dPhiNegLepMetH = new TH1F("dPhiNegLepMetH","",50,0.,TMath::Pi());
  TH1F * diLepEtaH = new TH1F("diLepEtaH","",180,-9.,9.);
  TH1F * dimuonPtH = new TH1F("dimuonPtH","",200,0.,400.);*/
  
 
  TH1F * nPVH = new TH1F("nPVH","",41,-0.5,40.5);

  //****************PV tracks
  /*TH1F * NAllPVtracksH=new TH1F("NAllPVtracksH","",100,-0.5,99.5);

  TH1F * NTracksPos_PVH = new TH1F("NTracksPos_PVH","",10,-0.5,9.5);
  TH1F * NTracksNeg_PVH = new TH1F("NTracksNeg_PVH","",10,-0.5,9.5);
  TH1F * NTracksPos_PV_NoMuonsH = new TH1F("NTracksPos_PV_NoMuonsH","",10,-0.5,9.5);
  TH1F * NTracksNeg_PV_NoMuonsH = new TH1F("NTracksNeg_PV_NoMuonsH","",10,-0.5,9.5);*/

  //****************PF
  TH1F * NpfNoPuPosH = new TH1F("NpfNoPuPosH","",10,-0.5,9.5);
  TH1F * NpfNoPuNegH = new TH1F("NpfNoPuNegH","",10,-0.5,9.5);
  TH1F * NpfNoPuIpPosH = new TH1F("NpfNoPuIpPosH","",10,-0.5,9.5);
  TH1F * NpfNoPuIpNegH = new TH1F("NpfNoPuIpNegH","",10,-0.5,9.5);

/*  TH1F * NTracksPos_NoMuonsH = new TH1F("NTracksPos_NoMuonsH","",10,-0.5,9.5);
  TH1F * NTracksNeg_NoMuonsH = new TH1F("NTracksNeg_NoMuonsH","",10,-0.5,9.5);
  TH1F * NTracksPos_NoMuons_ipH = new TH1F("NTracksPos_NoMuons_ipH","",10,-0.5,9.5);
  TH1F * NTracksNeg_NoMuons_ipH = new TH1F("NTracksNeg_NoMuons_ipH","",10,-0.5,9.5);
  TH1F * NTracksPos_ipH = new TH1F("NTracksPos_ipH","",10,-0.5,9.5);
  TH1F * NTracksNeg_ipH = new TH1F("NTracksNeg_ipH","",10,-0.5,9.5);*/

  TH1F * AllNpfNoPUH=new TH1F("AllNpfNoPUH","",100,-0.5,99.5);
  TH1F * AllNpfNoPU_ipH=new TH1F("AllNpfNoPU_ipH","",100,-0.5,99.5);


  //*******************************************
  // 1D mass distributions in the signal region
  //*******************************************
  TH1F * InvMassTrackPosH= new TH1F("InvMassTrackPosH","",100,0,2.0);
  TH1F * InvMassTrackNegH= new TH1F("InvMassTrackNegH","",100,0,2.0);

  //*******************************************
  //******** General control plots ************
  //*******************************************

  // signal region
  TH2F * InvMassTrackPlusMuon2D_H= new TH2F("InvMassTrackPlusMuon2D_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_UpH= new TH2F("InvMassTrackPlusMuon2D_UpH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_DownH= new TH2F("InvMassTrackPlusMuon2D_DownH","",100,0.,10.,100,0.,10.);
  
  TH1F * metH = new TH1F("metH","",100,0.,500.);
  TH1F * InvMassTracksPlusMuonsH = new TH1F("InvMassTracksPlusMuonsH","",200,0.,200.);

  TH1F * HardMuonPtH = new TH1F("HardMuonPtH","",100,0,50.0);
  TH1F * SoftMuonPtH = new TH1F("SoftMuonPtH","",100,0,50.0);
  TH1F * HardMuonEtaH = new TH1F("HardMuonEtaH","",50,-2.5,2.5);
  TH1F * SoftMuonEtaH = new TH1F("SoftMuonEtaH","",50,-2.5,2.5);

  TH1F * SumOfInvMass2SystemsH=new TH1F("SumOfInvMass2SystemsH","",50,0.,25.);
  TH1F * InvMassH2H= new TH1F("InvMassH2H","",200,0.,200.);

  TH1F * PtTrackHardH=new TH1F("PtTrackHardH","",100,0.,50.);
  TH1F * PtTrackSoftH=new TH1F("PtTrackSoftH","",100,0.,50.);
  TH1F * EtaTrackHardH=new TH1F("EtaTrackHardH","",50,-2.5,2.5);
  TH1F * EtaTrackSoftH=new TH1F("EtaTrackSoftH","",50,-2.5,2.5);
 
  TH1F * PtTrackMuHardH=new TH1F("PtTrackMuHardH","",100,0.,50.);
  TH1F * PtTrackMuSoftH=new TH1F("PtTrackMuSoftH","",100,0.,50.);
  TH1F * EtaTrackMuHardH=new TH1F("EtaTrackMuHardH","",50,-2.5,2.5);
  TH1F * EtaTrackMuSoftH=new TH1F("EtaTrackMuSoftH","",50,-2.5,2.5);



  // control region (at least one muon has second soft track)
  TH2F * InvMassTrackPlusMuon2D_ControlH = new TH2F("InvMassTrackPlusMuon2D_ControlH","",100,0.,10.,100,0.,10.);
  TH1F * met_ControlH = new TH1F("met_ControlH","",100,0.,500.);
  TH1F * InvMassTracksPlusMuons_ControlH = new TH1F("InvMassTracksPlusMuons_ControlH","",200,0.,200.);
  TH1F * HardMuonPt_ControlH = new TH1F("HardMuonPt_ControlH","",100,0.,200.);
  TH1F * SoftMuonPt_ControlH = new TH1F("SoftMuonPt_ControlH","",100,0.,200.);
  TH1F * SumOfInvMass2Systems_ControlH=new TH1F("SumOfInvMass2Systems_ControlH","",50,0.,25.);
  TH1F * InvMassH2_ControlH= new TH1F("InvMassH2_ControlH","",200,0.,200.);
  TH1F * PtTrackHard_ControlH=new TH1F("PtTrackHard_ControlH","",100,0.,100.);
  TH1F * PtTrackSoft_ControlH=new TH1F("PtTrackSoft_ControlH","",100,0.,100.);
  TH1F * EtaTrackHard_ControlH=new TH1F("EtaTrackHard_ControlH","",50,-2.5,2.5);
  TH1F * EtaTrackSoft_ControlH=new TH1F("EtaTrackSoft_ControlH","",50,-2.5,2.5);
  TH1F * dRMuMu_ControlH=new TH1F("dRMuMu_ControlH","",50,0,5.);
  TH1F * dEtaMuMu_ControlH=new TH1F("dEtaMuMu_ControlH","",50,0,5.);
  TH1F * dPhiMuMu_ControlH=new TH1F("dPhiMuMu_ControlH","",50,0.,TMath::Pi()); 

  TH1F * InvMassTrackPlusPosMuon1D_ControlH = new TH1F("InvMassTrackPlusPosMuon1D_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNegMuon1D_ControlH = new TH1F("InvMassTrackPlusNegMuon1D_ControlH","",100,0.,10.); 
  TH1F * InvMassTrackPlusMuon1D_ControlH = new TH1F("InvMassTrackPlusMuon1D_ControlH","",100,0.,10.); 
  

  // no Pt cuts on tracks
  TH1F * PtTrackHardSemiFinalH=new TH1F("PtTrackHardSemiFinalH","",100,0.,100.);
  TH1F * PtTrackSoftSemiFinalH=new TH1F("PtTrackSoftSemiFinalH","",100,0.,100.);


  // Muon selection
  TH1F * HardMuonPt_DimuonsH = new TH1F("HardMuonPt_DimuonsH","",40,0.,200.);
  TH1F * SoftMuonPt_DimuonsH = new TH1F("SoftMuonPt_DimuonsH","",40,0.,200.);
  TH1F * HardMuonEta_DimuonsH = new TH1F("HardMuonEta_DimuonsH","",100,-2.5,2.5);
  TH1F * SoftMuonEta_DimuonsH = new TH1F("SoftMuonEta_DimuonsH","",100,-2.5,2.5);
  TH1F * dRMuMu_DimuonsH=new TH1F("dRMuMu_DimuonsH","",50,0,5.);
  TH1F * dPhiMuMu_DimuonsH=new TH1F("dPhiMuMu_DimuonsH","",50,0,TMath::Pi());
  TH2F * dRdPhiMuMu_DimuonsH=new TH2F("dRdPhiMuMu_DimuonsH","",100,0,5.,100,0.,TMath::Pi());
  TH2F * HardMuonPtSoftMuonPt_DimuonsH = new TH2F("HardMuonPtSoftMuonPt_DimuonsH","",40,0,200,40,0,200);
  
  TH1F * diMuonMass_DimuonsH=new TH1F("diMuonMass_DimuonsH","",150,0,300);

  // dR between muon and nearest track;
  TH1F * dRMuPosTrkMinH = new TH1F("dRMuPosTrkMinH","",300,0.,3.);
  TH1F * dRMuNegTrkMinH = new TH1F("dRMuNegTrkMinH","",300,0.,3.);
  TH1F * dRMuTrkMinH = new TH1F("dRMuTrkMinH","",300,0.,3.);
  
  // track isolations
  TH1F * trackIsoPosH = new TH1F("trackIsoPosH","",50,0,2);
  TH1F * trackIsoNegH = new TH1F("trackIsoNegH","",50,0,2);
  
  TH1F * trackIsoHardH = new TH1F("trackIsoHardH","",50,0,2);
  TH1F * trackIsoSoftH = new TH1F("trackIsoSoftH","",50,0,2);

  TH1F * counterBeforeIsoH = new TH1F("counterBeforeIsoH","",1,-0.5,0.5);
  TH1F * counterAfterIsoH = new TH1F("counterAfterIsoH","",1,-0.5,0.5);


  //*************Outer Cone **************
  TH1F * Ntracks_OutConePosH=new TH1F("Ntracks_OutConePosH","",20,-0.5,19.5);
  TH1F * Ntracks_OutConeNegH=new TH1F("Ntracks_OutConeNegH","",20,-0.5,19.5);
  TH1F * Ntracks_OutConePosIPH=new TH1F("Ntracks_OutConePosIPH","",20,-0.5,19.5);
  TH1F * Ntracks_OutConeNegIPH=new TH1F("Ntracks_OutConeNegIPH","",20,-0.5,19.5);
  TH1F * SumPtPosConeH= new TH1F("SumPtPosConeH","",100,0,10.0);
  TH1F * SumPtNegConeH= new TH1F("SumPtNegConeH","",100,0,10.0);
  TH1F * SumPtPosConeIPH= new TH1F("SumPtPosConeIPH","",100,0,10.0);
  TH1F * SumPtNegConeIPH= new TH1F("SumPtNegConeIPH","",100,0,10.0);
  //***************************************************

  
 
  // *********************
  // 1st muon+track system
  // *********************
 
  // 1D muon-track mass distribution in signal region
  TH1F * InvMassTrackPlusPosH= new TH1F("InvMassTrackPlusPosH","",100,0.,10.);

  TH1F * InvMassTrackPlusPos_MuH= new TH1F("InvMassTrackPlusPos_MuH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_EleH= new TH1F("InvMassTrackPlusPos_EleH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_HadH= new TH1F("InvMassTrackPlusPos_HadH","",100,0.,10.);

  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it) 
  // track of the probed muon-track system do not pass dxy and dZ impact parameter requirements
  TH1F * InvMassTrackPlusPos_2H= new TH1F("InvMassTrackPlusPos_2H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_3H= new TH1F("InvMassTrackPlusPos_3H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_4H= new TH1F("InvMassTrackPlusPos_4H","",100,0.,10.);
  //TH1F * InvMassTrackPlusPos_InvNegH= new TH1F("InvMassTrackPlusPos_InvNegH","",100,0.,10.);

  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it)
  // track of the probed muon-track system does pass dxy and dZ impact parameter requirements
  TH1F * InvMassTrackPlusPos_2QH= new TH1F("InvMassTrackPlusPos_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_3QH= new TH1F("InvMassTrackPlusPos_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_4QH= new TH1F("InvMassTrackPlusPos_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_InvNegQH= new TH1F("InvMassTrackPlusPos_InvNegQH","",100,0.,10.);
  // tag muon has one track of the same sign around it
  TH1F * InvMassTrackPlusPos_SwapQH= new TH1F("InvMassTrackPlusPos_SwapQH","",100,0.,10.);

  // mass distributions for different track Id's
  // Mu
  TH1F * InvMassTrackPlusPos_Mu_2QH= new TH1F("InvMassTrackPlusPos_Mu_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_3QH= new TH1F("InvMassTrackPlusPos_Mu_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_4QH= new TH1F("InvMassTrackPlusPos_Mu_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_InvNegQH= new TH1F("InvMassTrackPlusPos_Mu_InvNegQH","",100,0.,10.);
  // Ele
  TH1F * InvMassTrackPlusPos_Ele_2QH= new TH1F("InvMassTrackPlusPos_Ele_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_3QH= new TH1F("InvMassTrackPlusPos_Ele_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_4QH= new TH1F("InvMassTrackPlusPos_Ele_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_InvNegQH= new TH1F("InvMassTrackPlusPos_Ele_InvNegQH","",100,0.,10.);  
  // Had
  TH1F * InvMassTrackPlusPos_Had_2QH= new TH1F("InvMassTrackPlusPos_Had_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_3QH= new TH1F("InvMassTrackPlusPos_Had_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_4QH= new TH1F("InvMassTrackPlusPos_Had_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_InvNegQH= new TH1F("InvMassTrackPlusPos_Had_InvNegQH","",100,0.,10.);  

  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it)

  // muon of the probed muon-track system has second soft track around it
  TH1F * InvMassTrackPlusPos_2Q_ControlH= new TH1F("InvMassTrackPlusPos_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_3Q_ControlH= new TH1F("InvMassTrackPlusPos_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_4Q_ControlH= new TH1F("InvMassTrackPlusPos_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_InvNegQ_ControlH= new TH1F("InvMassTrackPlusPos_InvNegQ_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_ControlH = new TH1F("InvMassTrackPlusPos_ControlH","",100,0.,10.);

  // mass distributions for different track Id's
  // Mu
  TH1F * InvMassTrackPlusPos_Mu_2Q_ControlH= new TH1F("InvMassTrackPlusPos_Mu_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_3Q_ControlH= new TH1F("InvMassTrackPlusPos_Mu_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_4Q_ControlH= new TH1F("InvMassTrackPlusPos_Mu_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Mu_InvNegQ_ControlH= new TH1F("InvMassTrackPlusPos_Mu_InvNegQ_ControlH","",100,0.,10.);

  // Ele
  TH1F * InvMassTrackPlusPos_Ele_2Q_ControlH= new TH1F("InvMassTrackPlusPos_Ele_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_3Q_ControlH= new TH1F("InvMassTrackPlusPos_Ele_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_4Q_ControlH= new TH1F("InvMassTrackPlusPos_Ele_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Ele_InvNegQ_ControlH= new TH1F("InvMassTrackPlusPos_Ele_InvNegQ_ControlH","",100,0.,10.);  
  // Had
  TH1F * InvMassTrackPlusPos_Had_2Q_ControlH= new TH1F("InvMassTrackPlusPos_Had_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_3Q_ControlH= new TH1F("InvMassTrackPlusPos_Had_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_4Q_ControlH= new TH1F("InvMassTrackPlusPos_Had_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Had_InvNegQ_ControlH= new TH1F("InvMassTrackPlusPos_Had_InvNegQ_ControlH","",100,0.,10.);  

  // invariant mass of the 1st muon-track pair in bins of invariant mass of the 2nd muon-track pair
  // all SS di-muon events
  TH1F * InvMassTrackPlusPos_Neg1GeVH= new TH1F("InvMassTrackPlusPos_Neg1GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg2GeVH= new TH1F("InvMassTrackPlusPos_Neg2GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg3GeVH= new TH1F("InvMassTrackPlusPos_Neg3GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg4GeVH= new TH1F("InvMassTrackPlusPos_Neg4GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_NegGT4GeVH = new TH1F("InvMassTrackPlusPos_NegGT4GeVH","",100,0.,10.);

  // invariant mass of the 1st muon-track pair in bins of invariant mass of the 2nd muon-track pair
  // 2nd muon has second soft track around it
  TH1F * InvMassTrackPlusPos_Neg1GeV_ControlPosH= new TH1F("InvMassTrackPlusPos_Neg1GeV_ControlPosH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg2GeV_ControlPosH= new TH1F("InvMassTrackPlusPos_Neg2GeV_ControlPosH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg3GeV_ControlPosH= new TH1F("InvMassTrackPlusPos_Neg3GeV_ControlPosH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg4GeV_ControlPosH= new TH1F("InvMassTrackPlusPos_Neg4GeV_ControlPosH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_NegGT4GeV_ControlPosH = new TH1F("InvMassTrackPlusPos_NegGT4GeV_ControlPosH","",100,0.,10.);

  // invariant mass of the 1st muon-track pair in bins of invariant mass of the 2nd muon-track pair
  // selected signal sample
  TH1F * InvMassTrackPlusPos_Neg1GeV_SignalH= new TH1F("InvMassTrackPlusPos_Neg1GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg2GeV_SignalH= new TH1F("InvMassTrackPlusPos_Neg2GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg3GeV_SignalH= new TH1F("InvMassTrackPlusPos_Neg3GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_Neg4GeV_SignalH= new TH1F("InvMassTrackPlusPos_Neg4GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_NegGT4GeV_SignalH = new TH1F("InvMassTrackPlusPos_NegGT4GeV_SignalH","",100,0.,10.);
  
  // *********************
  // 2nd muon+track system
  // *********************

  // 1D muon-track mass distribution in signal region
  TH1F * InvMassTrackPlusNegH= new TH1F("InvMassTrackPlusNegH","",100,0.,10.);

  TH1F * InvMassTrackPlusNeg_MuH= new TH1F("InvMassTrackPlusNeg_MuH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_EleH= new TH1F("InvMassTrackPlusNeg_EleH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_HadH= new TH1F("InvMassTrackPlusNeg_HadH","",100,0.,10.);

  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it) 
  // track of the probed muon-track system do not pass dxy and dZ impact parameter requirements
  TH1F * InvMassTrackPlusNeg_2H= new TH1F("InvMassTrackPlusNeg_2H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_3H= new TH1F("InvMassTrackPlusNeg_3H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_4H= new TH1F("InvMassTrackPlusNeg_4H","",100,0.,10.);
  //TH1F * InvMassTrackPlusNeg_InvPosH= new TH1F("InvMassTrackPlusNeg_InvPosH","",100,0.,10.);

  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it)
  // track of the probed muon-track system does pass dxy and dZ impact parameter requirements
  TH1F * InvMassTrackPlusNeg_2QH= new TH1F("InvMassTrackPlusNeg_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_3QH= new TH1F("InvMassTrackPlusNeg_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_4QH= new TH1F("InvMassTrackPlusNeg_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_InvPosQH= new TH1F("InvMassTrackPlusNeg_InvPosQH","",100,0.,10.);

  // tag muon has one track of the same sign around it
  TH1F * InvMassTrackPlusNeg_SwapQH= new TH1F("InvMassTrackPlusNeg_SwapQH","",100,0.,10.);

  // mass distributions for different track Id's
  // Mu
  TH1F * InvMassTrackPlusNeg_Mu_2QH= new TH1F("InvMassTrackPlusNeg_Mu_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_3QH= new TH1F("InvMassTrackPlusNeg_Mu_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_4QH= new TH1F("InvMassTrackPlusNeg_Mu_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_InvPosQH= new TH1F("InvMassTrackPlusNeg_Mu_InvPosQH","",100,0.,10.);
  // Ele
  TH1F * InvMassTrackPlusNeg_Ele_2QH= new TH1F("InvMassTrackPlusNeg_Ele_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_3QH= new TH1F("InvMassTrackPlusNeg_Ele_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_4QH= new TH1F("InvMassTrackPlusNeg_Ele_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_InvPosQH= new TH1F("InvMassTrackPlusNeg_Ele_InvPosQH","",100,0.,10.);  
  // Had
  TH1F * InvMassTrackPlusNeg_Had_2QH= new TH1F("InvMassTrackPlusNeg_Had_2QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_3QH= new TH1F("InvMassTrackPlusNeg_Had_3QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_4QH= new TH1F("InvMassTrackPlusNeg_Had_4QH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_InvPosQH= new TH1F("InvMassTrackPlusNeg_Had_InvPosQH","",100,0.,10.);  


  // 1D mass distribution of probed muon-track system (tag muon has 2, 3 or 4 tracks around it)
  // muon of the probed muon-track system has second soft track around it
  TH1F * InvMassTrackPlusNeg_2Q_ControlH= new TH1F("InvMassTrackPlusNeg_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_3Q_ControlH= new TH1F("InvMassTrackPlusNeg_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_4Q_ControlH= new TH1F("InvMassTrackPlusNeg_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_InvPosQ_ControlH= new TH1F("InvMassTrackPlusNeg_InvPosQ_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_ControlH = new TH1F("InvMassTrackPlusNeg_ControlH","",100,0.,10.);

  // mass distributions for different track Id's
  // Mu
  TH1F * InvMassTrackPlusNeg_Mu_2Q_ControlH= new TH1F("InvMassTrackPlusNeg_Mu_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_3Q_ControlH= new TH1F("InvMassTrackPlusNeg_Mu_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_4Q_ControlH= new TH1F("InvMassTrackPlusNeg_Mu_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Mu_InvPosQ_ControlH= new TH1F("InvMassTrackPlusNeg_Mu_InvPosQ_ControlH","",100,0.,10.);
  // Ele
  TH1F * InvMassTrackPlusNeg_Ele_2Q_ControlH= new TH1F("InvMassTrackPlusNeg_Ele_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_3Q_ControlH= new TH1F("InvMassTrackPlusNeg_Ele_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_4Q_ControlH= new TH1F("InvMassTrackPlusNeg_Ele_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Ele_InvPosQ_ControlH= new TH1F("InvMassTrackPlusNeg_Ele_InvPosQ_ControlH","",100,0.,10.);  
  // Had
  TH1F * InvMassTrackPlusNeg_Had_2Q_ControlH= new TH1F("InvMassTrackPlusNeg_Had_2Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_3Q_ControlH= new TH1F("InvMassTrackPlusNeg_Had_3Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_4Q_ControlH= new TH1F("InvMassTrackPlusNeg_Had_4Q_ControlH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Had_InvPosQ_ControlH= new TH1F("InvMassTrackPlusNeg_Had_InvPosQ_ControlH","",100,0.,10.);  
  

  // invariant mass of the 2nd muon-track pair in bins of invariant mass of the 1st muon-track pair
  // all SS di-muon events
  TH1F * InvMassTrackPlusNeg_Pos1GeVH= new TH1F("InvMassTrackPlusNeg_Pos1GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos2GeVH= new TH1F("InvMassTrackPlusNeg_Pos2GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos3GeVH= new TH1F("InvMassTrackPlusNeg_Pos3GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos4GeVH= new TH1F("InvMassTrackPlusNeg_Pos4GeVH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_PosGT4GeVH = new TH1F("InvMassTrackPlusNeg_PosGT4GeVH","",100,0.,10.);

  // invariant mass of the 2nd muon-track pair in bins of invariant mass of the 1st muon-track pair
  // 1st muon has second soft track around it
  TH1F * InvMassTrackPlusNeg_Pos1GeV_ControlNegH= new TH1F("InvMassTrackPlusNeg_Pos1GeV_ControlNegH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos2GeV_ControlNegH= new TH1F("InvMassTrackPlusNeg_Pos2GeV_ControlNegH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos3GeV_ControlNegH= new TH1F("InvMassTrackPlusNeg_Pos3GeV_ControlNegH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos4GeV_ControlNegH= new TH1F("InvMassTrackPlusNeg_Pos4GeV_ControlNegH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_PosGT4GeV_ControlNegH = new TH1F("InvMassTrackPlusNeg_PosGT4GeV_ControlNegH","",100,0.,10.);

  // invariant mass of the 2nd muon-track pair in bins of invariant mass of the 1st muon-track pair
  // selected signal sample
  TH1F * InvMassTrackPlusNeg_Pos1GeV_SignalH= new TH1F("InvMassTrackPlusNeg_Pos1GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos2GeV_SignalH= new TH1F("InvMassTrackPlusNeg_Pos2GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos3GeV_SignalH= new TH1F("InvMassTrackPlusNeg_Pos3GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_Pos4GeV_SignalH= new TH1F("InvMassTrackPlusNeg_Pos4GeV_SignalH","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_PosGT4GeV_SignalH = new TH1F("InvMassTrackPlusNeg_PosGT4GeV_SignalH","",100,0.,10.);

  // 

  TH1F * InvMassTrackPlusMuon_QCDH=new TH1F("InvMassTrackPlusMuon_QCDH","",100,0.,10.);
  TH1F * InvMassTrackPlusMuon_OSOSH=new TH1F("InvMassTrackPlusMuon_OSOSH","",100,0.,10.);

  //Bjet veto control plots
  TH1F * InvMassTrackPlusPos_bothBJets20H = new TH1F("InvMassTrackPlusPos_bothBJets20H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_bothBJets20H = new TH1F("InvMassTrackPlusNeg_bothBJets20H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_posBJet20H= new TH1F("InvMassTrackPlusPos_posBJet20H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_posBJet20H= new TH1F("InvMassTrackPlusNeg_posBJet20H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_negBJet20H= new TH1F("InvMassTrackPlusPos_negBJet20H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_negBJet20H= new TH1F("InvMassTrackPlusNeg_negBJet20H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_NoBJets20H= new TH1F("InvMassTrackPlusPos_NoBJets20H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_NoBJets20H= new TH1F("InvMassTrackPlusNeg_NoBJets20H","",100,0.,10.);

  TH1F * InvMassTrackPlusPos_bothBJets15H = new TH1F("InvMassTrackPlusPos_bothBJets15H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_bothBJets15H = new TH1F("InvMassTrackPlusNeg_bothBJets15H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_posBJet15H= new TH1F("InvMassTrackPlusPos_posBJet15H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_posBJet15H= new TH1F("InvMassTrackPlusNeg_posBJet15H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_negBJet15H= new TH1F("InvMassTrackPlusPos_negBJet15H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_negBJet15H= new TH1F("InvMassTrackPlusNeg_negBJet15H","",100,0.,10.);
  TH1F * InvMassTrackPlusPos_NoBJets15H= new TH1F("InvMassTrackPlusPos_NoBJets15H","",100,0.,10.);
  TH1F * InvMassTrackPlusNeg_NoBJets15H= new TH1F("InvMassTrackPlusNeg_NoBJets15H","",100,0.,10.);


//*******************Background "side-bands" using b jets(START)*********************
   TH1F * InvMassTrackPlusPos_NegIso_bothBJets20CSVL_2QH= new TH1F("InvMassTrackPlusPos_NegIso_bothBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_NegIso_bothBJets20CSVL_3QH= new TH1F("InvMassTrackPlusPos_NegIso_bothBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_NegIso_bothBJets20CSVL_InvPosQH= new TH1F("InvMassTrackPlusPos_NegIso_bothBJets20CSVL_InvPosQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_InvPosQH= new TH1F("InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_InvPosQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_NegIso_posBJet20CSVL_2QH= new TH1F("InvMassTrackPlusPos_NegIso_posBJet20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_posBJet20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_NegIso_posBJet20CSVL_2QH","",100,0.,10.);  
   TH1F * InvMassTrackPlusPos_NegIso_posBJet20CSVL_3QH= new TH1F("InvMassTrackPlusPos_NegIso_posBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_posBJet20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_NegIso_posBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_NegIso_posBJet20CSVL_InvPosQH= new TH1F("InvMassTrackPlusPos_NegIso_posBJet20CSVL_InvPosQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_posBJet20CSVL_InvPosQH= new TH1F("InvMassTrackPlusNeg_NegIso_posBJet20CSVL_InvPosQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_NegIso_negBJet20CSVL_2QH= new TH1F("InvMassTrackPlusPos_NegIso_negBJet20CSVL_2QH","",100,0.,10.); 
   TH1F * InvMassTrackPlusNeg_NegIso_negBJet20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_NegIso_negBJet20CSVL_2QH","",100,0.,10.); 
   TH1F * InvMassTrackPlusPos_NegIso_negBJet20CSVL_3QH= new TH1F("InvMassTrackPlusPos_NegIso_negBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_negBJet20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_NegIso_negBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_NegIso_negBJet20CSVL_InvPosQH= new TH1F("InvMassTrackPlusPos_NegIso_negBJet20CSVL_InvPosQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_negBJet20CSVL_InvPosQH= new TH1F("InvMassTrackPlusNeg_NegIso_negBJet20CSVL_InvPosQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_NegIso_NoBJets20CSVL_2QH= new TH1F("InvMassTrackPlusPos_NegIso_NoBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_2QH","",100,0.,10.);  
   TH1F * InvMassTrackPlusPos_NegIso_NoBJets20CSVL_3QH= new TH1F("InvMassTrackPlusPos_NegIso_NoBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_NegIso_NoBJets20CSVL_InvPosQH= new TH1F("InvMassTrackPlusPos_NegIso_NoBJets20CSVL_InvPosQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_InvPosQH= new TH1F("InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_InvPosQH","",100,0.,10.);
   ///
   TH1F * InvMassTrackPlusPos_PosIso_bothBJets20CSVL_2QH= new TH1F("InvMassTrackPlusPos_PosIso_bothBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_PosIso_bothBJets20CSVL_3QH= new TH1F("InvMassTrackPlusPos_PosIso_bothBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_PosIso_bothBJets20CSVL_InvNegQH= new TH1F("InvMassTrackPlusPos_PosIso_bothBJets20CSVL_InvNegQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_InvNegQH= new TH1F("InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_InvNegQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_PosIso_posBJet20CSVL_2QH= new TH1F("InvMassTrackPlusPos_PosIso_posBJet20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_posBJet20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_PosIso_posBJet20CSVL_2QH","",100,0.,10.);  
   TH1F * InvMassTrackPlusPos_PosIso_posBJet20CSVL_3QH= new TH1F("InvMassTrackPlusPos_PosIso_posBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_posBJet20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_PosIso_posBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_PosIso_posBJet20CSVL_InvNegQH= new TH1F("InvMassTrackPlusPos_PosIso_posBJet20CSVL_InvNegQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_posBJet20CSVL_InvNegQH= new TH1F("InvMassTrackPlusNeg_PosIso_posBJet20CSVL_InvNegQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_PosIso_negBJet20CSVL_2QH= new TH1F("InvMassTrackPlusPos_PosIso_negBJet20CSVL_2QH","",100,0.,10.); 
   TH1F * InvMassTrackPlusNeg_PosIso_negBJet20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_PosIso_negBJet20CSVL_2QH","",100,0.,10.); 
   TH1F * InvMassTrackPlusPos_PosIso_negBJet20CSVL_3QH= new TH1F("InvMassTrackPlusPos_PosIso_negBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_negBJet20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_PosIso_negBJet20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_PosIso_negBJet20CSVL_InvNegQH= new TH1F("InvMassTrackPlusPos_PosIso_negBJet20CSVL_InvNegQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_negBJet20CSVL_InvNegQH= new TH1F("InvMassTrackPlusNeg_PosIso_negBJet20CSVL_InvNegQH","",100,0.,10.);
 
   TH1F * InvMassTrackPlusPos_PosIso_NoBJets20CSVL_2QH= new TH1F("InvMassTrackPlusPos_PosIso_NoBJets20CSVL_2QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_2QH= new TH1F("InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_2QH","",100,0.,10.);  
   TH1F * InvMassTrackPlusPos_PosIso_NoBJets20CSVL_3QH= new TH1F("InvMassTrackPlusPos_PosIso_NoBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_3QH= new TH1F("InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_3QH","",100,0.,10.);
   TH1F * InvMassTrackPlusPos_PosIso_NoBJets20CSVL_InvNegQH= new TH1F("InvMassTrackPlusPos_PosIso_NoBJets20CSVL_InvNegQH","",100,0.,10.);
   TH1F * InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_InvNegQH= new TH1F("InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_InvNegQH","",100,0.,10.);

   TH1F *  InvMassTrackPlusPos_BJet20CSVL_2QH=new TH1F("InvMassTrackPlusPos_BJet20CSVL_2QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusPos_BJet20CSVL_3QH=new TH1F("InvMassTrackPlusPos_BJet20CSVL_3QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusPos_BJet20CSVL_InvNegQH=new TH1F("InvMassTrackPlusPos_BJet20CSVL_InvNegQH","",100,0.,10.);

   TH1F *  InvMassTrackPlusPos_NotBJet20CSVL_2QH=new TH1F("InvMassTrackPlusPos_NotBJet20CSVL_2QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusPos_NotBJet20CSVL_3QH=new TH1F("InvMassTrackPlusPos_NotBJet20CSVL_3QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusPos_NotBJet20CSVL_InvNegQH=new TH1F("InvMassTrackPlusPos_NotBJet20CSVL_InvNegQH","",100,0.,10.);

   TH1F *  InvMassTrackPlusNeg_BJet20CSVL_2QH=new TH1F("InvMassTrackPlusNeg_BJet20CSVL_2QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusNeg_BJet20CSVL_3QH=new TH1F("InvMassTrackPlusNeg_BJet20CSVL_3QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusNeg_BJet20CSVL_InvPosQH=new TH1F("InvMassTrackPlusNeg_BJet20CSVL_InvPosQH","",100,0.,10.);

   TH1F *  InvMassTrackPlusNeg_NotBJet20CSVL_2QH=new TH1F("InvMassTrackPlusNeg_NotBJet20CSVL_2QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusNeg_NotBJet20CSVL_3QH=new TH1F("InvMassTrackPlusNeg_NotBJet20CSVL_3QH","",100,0.,10.);
   TH1F *  InvMassTrackPlusNeg_NotBJet20CSVL_InvPosQH=new TH1F("InvMassTrackPlusNeg_NotBJet20CSVL_InvPosQH","",100,0.,10.);

  // 2D distributions mass of the 1st muon vs 2nd muon
  // signal region

  TH2F * InvMassTrackPlusMuon2D_0CSVLH= new TH2F("InvMassTrackPlusMuon2D_0CSVLH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_1CSVLH= new TH2F("InvMassTrackPlusMuon2D_1CSVLH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_2CSVLH= new TH2F("InvMassTrackPlusMuon2D_2CSVLH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_0CSVMH= new TH2F("InvMassTrackPlusMuon2D_0CSVMH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_1CSVMH= new TH2F("InvMassTrackPlusMuon2D_1CSVMH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_2CSVMH= new TH2F("InvMassTrackPlusMuon2D_2CSVMH","",100,0.,10.,100,0.,10.);


  TH2F * InvMass2DTrackPlusProbe_2H= new TH2F("InvMass2DTrackPlusProbe_2H","",100,0.,10.,100,0.,10.);
  TH2F * InvMass2DTrackPlusProbe_3H= new TH2F("InvMass2DTrackPlusProbe_3H","",100,0.,10.,100,0.,10.);
  TH2F * InvMass2DTrackPlusProbe_4H= new TH2F("InvMass2DTrackPlusProbe_4H","",100,0.,10.,100,0.,10.);
  TH2F * InvMass2DTrackPlusMuon_QCDH= new TH2F("InvMass2DTrackPlusMuon_QCDH","",100,0.,10.,100,0.,10.);
  TH2F * InvMass2DTrackPlusMuon_OSOSH= new TH2F("InvMass2DTrackPlusMuon_OSOSH","",100,0.,10.,100,0.,10.);

  // control region (one of the muons has second close-by soft track)
  TH2F * InvMassTrackPlusMuon2D_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_ControlPosH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_bothBJets20_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_bothBJets20_ControlPosH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_bothBJets20_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_bothBJets20_ControlNegH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_posBJet20_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_posBJet20_ControlPosH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_posBJet20_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_posBJet20_ControlNegH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_negBJet20_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_negBJet20_ControlPosH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_negBJet20_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_negBJet20_ControlNegH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_NoBJets20_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_NoBJets20_ControlPosH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_NoBJets20_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_NoBJets20_ControlNegH","",100,0.,10.,100,0.,10.);

  ///  
  TH2F * InvMassTrackPlusMuon2D_MuMu_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_MuMu_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuMu_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_MuMu_ControlPosH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_MuEle_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_MuEle_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuEle_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_MuEle_ControlPosH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_MuHad_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_MuHad_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuHad_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_MuHad_ControlPosH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_EleEle_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_EleEle_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_EleEle_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_EleEle_ControlPosH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_EleHad_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_EleHad_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_EleHad_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_EleHad_ControlPosH","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_HadHad_ControlNegH = new TH2F("InvMassTrackPlusMuon2D_HadHad_ControlNegH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_HadHad_ControlPosH = new TH2F("InvMassTrackPlusMuon2D_HadHad_ControlPosH","",100,0.,10.,100,0.,10.);

  // control region (both muons have second close-by soft track)
  TH2F * InvMassTrackPlusMuon2D_ControlBothH = new TH2F("InvMassTrackPlusMuon2D_ControlBothH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_bothBJets20_ControlH = new TH2F("InvMassTrackPlusMuon2D_bothBJets20_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_posBJet20_ControlH = new TH2F("InvMassTrackPlusMuon2D_posBJet20_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_negBJet20_ControlH = new TH2F("InvMassTrackPlusMuon2D_negBJet20_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_NoBJets20_ControlH = new TH2F("InvMassTrackPlusMuon2D_NoBJets20_ControlH","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_MuMu_ControlH = new TH2F("InvMassTrackPlusMuon2D_MuMu_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuEle_ControlH = new TH2F("InvMassTrackPlusMuon2D_MuEle_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuHad_ControlH = new TH2F("InvMassTrackPlusMuon2D_MuHad_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_EleEle_ControlH = new TH2F("InvMassTrackPlusMuon2D_EleEle_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_EleHad_ControlH = new TH2F("InvMassTrackPlusMuon2D_EleHad_ControlH","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_HadHad_ControlH = new TH2F("InvMassTrackPlusMuon2D_HadHad_ControlH","",100,0.,10.,100,0.,10.);

  //2D bjet matching
  TH2F * InvMassTrackPlusMuon2D_bothBJets20_H= new TH2F("InvMassTrackPlusMuon2D_bothBJets20_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_posBJet20_H= new TH2F("InvMassTrackPlusMuon2D_posBJet20_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_negBJet20_H= new TH2F("InvMassTrackPlusMuon2D_negBJet20_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_NoBJets20_H= new TH2F("InvMassTrackPlusMuon2D_NoBJets20_H","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_bothBJets15_H= new TH2F("InvMassTrackPlusMuon2D_bothBJets15_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_posBJet15_H= new TH2F("InvMassTrackPlusMuon2D_posBJet15_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_negBJet15_H= new TH2F("InvMassTrackPlusMuon2D_negBJet15_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_NoBJets15_H= new TH2F("InvMassTrackPlusMuon2D_NoBJets15_H","",100,0.,10.,100,0.,10.);

  /// Id Track categories
  TH2F * InvMassTrackPlusMuon2D_MuMuH = new TH2F("InvMassTrackPlusMuon2D_MuMu_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuEleH = new TH2F("InvMassTrackPlusMuon2D_MuEle_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_MuHadH = new TH2F("InvMassTrackPlusMuon2D_MuHad_H","",100,0.,10.,100,0.,10.);

  TH2F * InvMassTrackPlusMuon2D_EleEleH = new TH2F("InvMassTrackPlusMuon2D_EleEle_H","",100,0.,10.,100,0.,10.);
  TH2F * InvMassTrackPlusMuon2D_EleHadH = new TH2F("InvMassTrackPlusMuon2D_EleHad_H","",100,0.,10.,100,0.,10.);
  
  TH2F * InvMassTrackPlusMuon2D_HadHadH = new TH2F("InvMassTrackPlusMuon2D_HadHad_H","",100,0.,10.,100,0.,10.);

  /// 

  TH1F * InvMassTrackPlusMuon_QCD_UnfH= new TH1F("InvMassTrackPlusMuon_QCD_UnfH","",100,0.,100.);
  TH1F * InvMassTrackPlusMuon2D_UnfH= new TH1F("InvMassTrackPlusMuon2D_UnfH","",100,0.,100.);
  TH1F * InvMassTrackPlusMuon_OSOS_UnfH=new TH1F("InvMassTrackPlusMuon_OSOS_UnfH","",100,0.,100.);
 




  
  //TH1F * PtPFNoPU1H=new TH1F("PtPFNoPU1H","",200,0.,400.);
  //TH1F * PtPFNoPU2H=new TH1F("PtPFNoPU2H","",200,0.,400.);
  
  // Generator information --->
  TH1F * GenMuonsInTrkColH = new TH1F("GenMuonsInTrkColH","",2,-0.5,1.5);
  TH1F * GenMuonsInPfNoPuColH = new TH1F("GenMuonsInPfNoPuColH","",2,-0.5,1.5);

  TH1F * GenMuonsMatchedRecoMuonsH = new TH1F("GenMuonsMatchedRecoMuonsH","",2,-0.5,1.5);
  TH1F * IpGenMuonsMatchedRecoMuonsH = new TH1F("IpGenMuonsMatchedRecoMuonsH","",100,-0.5,0.5);
  TH1F * IpZGenMuonsMatchedRecoMuonsH = new TH1F("IpZGenMuonsMatchedRecoMuonsH","",200,-1.0,1.0);
  
  TH1F * IpGenMuonsH = new TH1F("IpGenMuonsH","",100,-0.5,0.5);
  TH1F * IpZGenMuonsH = new TH1F("IpZGenMuonsH","",200,-1.,1.);
  
  TH1F * GenTracksInTrkColH = new TH1F("GenTracksInTrkColH","",2,-0.5,1.5);
  TH1F * GenTracksInPfNoPuColH = new TH1F("GenTracksInPfNoPuColH","",2,-0.5,1.5);

  TH1F * IpGenTracksH = new TH1F("IpGenTracksH","",100,-0.5,0.5);
  TH1F * IpZGenTracksH = new TH1F("IpZGenTracksH","",200,-1.,1.);

  TH1F * PtGenMuonsEtaCutH  = new TH1F("PtGenMuonsEtaCutH","",100,0.,50.);
  TH1F * PtGenTracksEtaCutH = new TH1F("PtGenTracksEtaCutH","",100,0.,50.);

  TH1F * dPGenMuonsTrkH = new TH1F("dPGenMuonsTrkH","",100,0.,0.1);
  TH1F * dPGenTracksTrkH = new TH1F("dPGenTracksTrkH","",100,0.,0.1);

  TH1F * dRGenMuonsTrkH = new TH1F("dRGenMuonsTrkH","",100,0.,0.1);
  TH1F * dRGenTracksTrkH = new TH1F("dRGenTracksTrkH","",100,0.,0.1);


  TH1F * dPGenMuonsPfNoPuH = new TH1F("dPGenMuonsPfNoPuH","",100,0.,0.1);
  TH1F * dPGenTracksPfNoPuH = new TH1F("dPGenTracksPfNoPuH","",100,0.,0.1);

  TH1F * dRGenMuonsPfNoPuH = new TH1F("dRGenMuonsPfNoPuH","",100,0.,0.1);
  TH1F * dRGenTracksPfNoPuH = new TH1F("dRGenTracksPfNoPuH","",100,0.,0.1);

  TH1F * muontracks_effH= new TH1F("muontrack_effH","",100,0.,1.);
  TH1F * recotracks_effH= new TH1F("recotracks_effH","",100,0.,1.);
  
  TH1F * muonPF_effH= new TH1F("muonPF_effH","",100,0.,1.);
  TH1F * trackPF_effH= new TH1F("recoPF_effH","",100,0.,1.);

  TH1F * dPGenMuonsMatchedRecoMuonsH = new TH1F("dPGenMuonsMatchedRecoMuonsH","",100,0.,0.1);
  //TH1F * dRGenMuonsMatchedRecoMuonsH = new TH1F("dRGenMuonsMatchedRecoMuonsH","",100,0.,0.1);

  int nPtBins = 5;
  int nEtaBins = 3;

  float ptBins[6]  = {0.5, 1.0, 2.0, 5.0, 10.0, 1000};
  float etaBins[4] = {0.0, 0.8, 1.5, 2.5};

  // Id of tracks 
  TH1F * idTrackPosH = new TH1F("idTrackPosH","",3,-0.5,2.5);
  TH1F * idTrackNegH = new TH1F("idTrackNegH","",3,-0.5,2.5);

  TH1F * idTrackPos_DoubleTagH = new TH1F("idTrackPos_DoubleTagH","",3,-0.5,2.5);
  TH1F * idTrackNeg_DoubleTagH = new TH1F("idTrackNeg_DoubleTagH","",3,-0.5,2.5);

  TH1F * idTrackPos_SingleTagH = new TH1F("idTrackPos_SingleTagH","",3,-0.5,2.5);
  TH1F * idTrackNeg_SingleTagH = new TH1F("idTrackNeg_SingleTagH","",3,-0.5,2.5);

  TH1F * idTrackPos_NoTagH = new TH1F("idTrackPos_NoTagH","",3,-0.5,2.5);
  TH1F * idTrackNeg_NoTagH = new TH1F("idTrackNeg_NoTagH","",3,-0.5,2.5);
  
  // 2D distributions of track Id in the signal region
  TH2F * idTrack2DH = new TH2F("idTrack2DH","",3,-0.5,2.5,3,-0.5,2.5);

  // 2D distributions of track Id in the background side-band region (one of muons has additional soft track)
  TH2F * idTrack_ControlPos_2DH = new TH2F("idTrack_ControlPos_2DH","",3,-0.5,2.5,3,-0.5,2.5);
  TH2F * idTrack_ControlNeg_2DH = new TH2F("idTrack_ControlNeg_2DH","",3,-0.5,2.5,3,-0.5,2.5);
  TH2F * idTrack_Control_2DH    = new TH2F("idTrack_Control_2DH",   "",3,-0.5,2.5,3,-0.5,2.5);
  
  // 1D distributions of track Id in the background side-band region (one of muons has additional soft track) 
  TH1F * idTrackPos_ControlPosH = new TH1F("idTrackPos_ControlPosH","",3,-0.5,2.5);
  TH1F * idTrackPos_ControlNegH = new TH1F("idTrackPos_ControlNegH","",3,-0.5,2.5);
  TH1F * idTrackPos_ControlH    = new TH1F("idTrackPos_ControlH",   "",3,-0.5,2.5);
  
  TH1F * idTrackNeg_ControlPosH = new TH1F("idTrackNeg_ControlPosH","",3,-0.5,2.5);
  TH1F * idTrackNeg_ControlNegH = new TH1F("idTrackNeg_ControlNegH","",3,-0.5,2.5);
  TH1F * idTrackNeg_ControlH    = new TH1F("idTrackNeg_ControlH",   "",3,-0.5,2.5);

  // jet multiplicities (signal region)
  TH1F * nJets30H = new TH1F("nJets30H","",10,-0.5,9.5);
  TH1F * nJets30EtaH = new TH1F("nJets30EtaH","",10,-0.5,9.5);
  TH1F * nJets20CSVLH = new TH1F("nJets20CSVLH","",10,-0.5,9.5);
  TH1F * nJets20CSVMH = new TH1F("nJets20CSVMH","",10,-0.5,9.5);
  TH1F * nJets15CSVLH = new TH1F("nJets15CSVLH","",10,-0.5,9.5);
  TH1F * nJets15CSVMH = new TH1F("nJets15CSVMH","",10,-0.5,9.5);

  // jet multiplicities (control region)
  TH1F * nJets30_ControlH = new TH1F("nJets30_ControlH","",10,-0.5,9.5);
  TH1F * nJets30Eta_ControlH = new TH1F("nJets30Eta_ControlH","",10,-0.5,9.5);
  TH1F * nJets20CSVL_ControlH = new TH1F("nJets20CSVL_ControlH","",10,-0.5,9.5);
  TH1F * nJets20CSVM_ControlH = new TH1F("nJets20CSVM_ControlH","",10,-0.5,9.5);
  TH1F * nJets15CSVL_ControlH = new TH1F("nJets15CSVL_ControlH","",10,-0.5,9.5);
  TH1F * nJets15CSVM_ControlH = new TH1F("nJets15CSVM_ControlH","",10,-0.5,9.5);



  TH1F * nJets20CSVL_DoubleTagH = new TH1F("nJets20CSVL_DoubleTagH","",10,-0.5,9.5);
  //  TH1F * nJets20CSVM_DoubleTagH = new TH1F("nJets20CSVM_DoubleTagH","",10,-0.5,9.5);
  
  TH1F * nJets20CSVL_SingleTagH = new TH1F("nJets20CSVL_SingleTagH","",10,-0.5,9.5);
  //  TH1F * nJets20CSVM_SingleTagH = new TH1F("nJets20CSVM_SingleTagH","",10,-0.5,9.5);

  TH1F * nJets20CSVL_NoTagH = new TH1F("nJets20CSVL_NoTagH","",10,-0.5,9.5);
  //  TH1F * nJets20CSVM_NoTagH = new TH1F("nJets20CSVM_NoTagH","",10,-0.5,9.5);

  TH1F * HiggsPtFinalH = new TH1F("HiggsPtFinalH","",100,0,200.);


  //float EfficiencyGenTracksInTrkColBinsPt[5];
  //float EfficiencyGenTracksInTrkColBinsEta[3];

 

  TH1F * GenTracksInTrkColPtBinsH[5];
  TH1F * GenTracksInTrkColEtaBinsH[3];
  TH1F * GenTracksInTrkColPtEtaBinsH[5][3];

  TH1F * GenTracksInPfNoPuColPtBinsH[5];
  TH1F * GenTracksInPfNoPuColEtaBinsH[3];
  TH1F * GenTracksInPfNoPuColPtEtaBinsH[5][3];



  TString ptBinStr[5] = {"0p5To1",
			 "1To2",
			 "2To5",
			 "5To10",
			 "10ToInf"};

  TString etaBinStr[3] = {"0To0p8",
			  "0p8To1p5",
			  "1p5To2p5"};


  for (int iPt=0; iPt<5; ++iPt) {
    TString histName = "GenTracksInTrkCol_Pt" + ptBinStr[iPt] + "H";
    GenTracksInTrkColPtBinsH[iPt] = new TH1F(histName,"",2,-0.5,1.5);
    histName = "GenTracksInPfNoPuCol_Pt" + ptBinStr[iPt] + "H";
    GenTracksInPfNoPuColPtBinsH[iPt] = new TH1F(histName,"",2,-0.5,1.5);
    for (int iEta=0; iEta<3; ++iEta) {
      histName = "GenTracksInTrkCol_Pt" + ptBinStr[iPt] + "_Eta" + etaBinStr[iEta] + "H";
      GenTracksInTrkColPtEtaBinsH[iPt][iEta] = new TH1F(histName,"",2,-0.5,1.5);
      histName = "GenTracksInPfNoPuCol_Pt" + ptBinStr[iPt] + "_Eta" + etaBinStr[iEta] + "H";
      GenTracksInPfNoPuColPtEtaBinsH[iPt][iEta] = new TH1F(histName,"",2,-0.5,1.5);
    }
  }

  for (int iEta=0; iEta<3; ++iEta) {
    TString histName = "GenTracksInTrkCol_Eta" + etaBinStr[iEta] + "H";
    GenTracksInTrkColEtaBinsH[iEta] = new TH1F(histName,"",2,-0.5,1.5);
    histName = "GenTracksInPfNoPuCol_Eta" + etaBinStr[iEta] + "H";
    GenTracksInPfNoPuColEtaBinsH[iEta] = new TH1F(histName,"",2,-0.5,1.5);
  }

  TH1F * massJPsiEleEleH = new TH1F("massJPsiEleEleH","",500,2.6,3.6);
  TH1F * massJPsiEleMuH  = new TH1F("massJPsiEleMuH","",500,2.6,3.6);
  TH1F * massJPsiEleHadH = new TH1F("massJPsiEleHadH","",500,2.6,3.6);

  TH1F * massJPsiMuEleH = new TH1F("massJPsiMuEleH","",500,2.6,3.6);
  TH1F * massJPsiMuMuH  = new TH1F("massJPsiMuMuH","",500,2.6,3.6);
  TH1F * massJPsiMuHadH = new TH1F("massJPsiMuHadH","",500,2.6,3.6);

  // ************* 

  TH1F * massYpsilonEleEleH = new TH1F("massYpsilonEleEleH","",2000,8.,12.);
  TH1F * massYpsilonEleMuH  = new TH1F("massYpsilonEleMuH","",2000,8.,12.);
  TH1F * massYpsilonEleHadH = new TH1F("massYpsilonEleHadH","",2000,8.,12.);

  TH1F * massYpsilonMuEleH = new TH1F("massYpsilonMuEleH","",2000,8.,12.);
  TH1F * massYpsilonMuMuH  = new TH1F("massYpsilonMuMuH","",2000,8.,12.);
  TH1F * massYpsilonMuHadH = new TH1F("massYpsilonMuHadH","",2000,8.,12.);

  // ****

  TH1F * massK0HadEleH = new TH1F("massK0HadEleH","",500,0.45,0.55);
  TH1F * massK0HadMuH  = new TH1F("massK0HadMuH","",500,0.45,0.55);
  TH1F * massK0HadHadH = new TH1F("massK0HadHadH","",500,0.45,0.55);



  TH1F * particleIdH = new TH1F("particleIdH","",5,0,5);

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

  //TH1F * puiweightMC;
  //TH1F * puiweightData;
  //TMVA Play
//   TFile * outputFileLikelihood=NULL;
//   TTree * likeTree=NULL;
//   TTree * likeTreeBkg=NULL;

//   TString likefilename=rootFileName+"_like.root";
//   outputFileLikelihood = new TFile(likefilename,"recreate");
//   outputFileLikelihood->cd("");
//   likeTree = new TTree("likeTree","likeTree");
//   likeTree->Branch("PtTrackPos",&PtTrackPos,"PtTrackPos/F");
//   likeTree->Branch("PtTrackNeg",&PtTrackNeg,"PtTrackNeg/F");
//   likeTree->Branch("met",&met,"met/F");
//   likeTree->Branch("InvMassH2",&InvMassH2,"InvMassH2/F");
//   likeTree->Branch("InvMassTrackPlusPosMuon",&InvMassTrackPlusPosMuon,"InvMassTrackPlusPosMuon/F");
//   likeTree->Branch("InvMassTrackPlusNegMuon",&InvMassTrackPlusNegMuon,"InvMassTrackPlusNegMuon/F");
//   likeTree->Branch("IPxyPos",&IPxyPos,"IPxyPos/F");
//   likeTree->Branch("IPxyNeg",&IPxyNeg,"IPxyNeg/F");
//   likeTree->Branch("IPzPos",&IPzPos,"IPzPos/F");
//   likeTree->Branch("IPzNeg",&IPzNeg,"IPzNeg/F");
//   likeTree->Branch("posLepPt",&_posLepPt,"posLepPt/F");
//   likeTree->Branch("negLepPt",&_negLepPt,"negLepPt/F");
//   likeTree->Branch("SumPtPosConeIP",&SumPtPosConeIP,"SumPtPosConeIP/F");
//   likeTree->Branch("SumPtNegConeIP",&SumPtNegConeIP,"SumPtNegConeIP/F");
  
//   likeTreeBkg = new TTree("likeTreeBkg","likeTreeBkg");
//   likeTreeBkg->Branch("PtTrackPos",&PtTrackPos,"PtTrackPos/F");
//   likeTreeBkg->Branch("PtTrackNeg",&PtTrackNeg,"PtTrackNeg/F");
//   likeTreeBkg->Branch("met",&met,"met/F");
//   likeTreeBkg->Branch("InvMassH2",&InvMassH2,"InvMassH2/F");
//   likeTreeBkg->Branch("InvMassTrackPlusPosMuon",&InvMassTrackPlusPosMuon,"InvMassTrackPlusPosMuon/F");
//   likeTreeBkg->Branch("InvMassTrackPlusNegMuon",&InvMassTrackPlusNegMuon,"InvMassTrackPlusNegMuon/F");
//   likeTreeBkg->Branch("IPxyPos",&IPxyPos,"IPxyPos/F");
//   likeTreeBkg->Branch("IPxyNeg",&IPxyNeg,"IPxyNeg/F");
//   likeTreeBkg->Branch("IPzPos",&IPzPos,"IPzPos/F");
//   likeTreeBkg->Branch("IPzNeg",&IPzNeg,"IPzNeg/F");
//   likeTreeBkg->Branch("posLepPt",&_posLepPt,"posLepPt/F");
//   likeTreeBkg->Branch("negLepPt",&_negLepPt,"negLepPt/F");
//   likeTreeBkg->Branch("SumPtPosConeIP",&SumPtPosConeIP,"SumPtPosConeIP/F");
//   likeTreeBkg->Branch("SumPtNegConeIP",&SumPtNegConeIP,"SumPtNegConeIP/F");
 
  std::vector<std::string> jsonFiles;
  jsonFiles.push_back("/afs/desy.de/user/r/rasp/public/Analysis/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");
 
  RunLumiSelector runLumiSelector;
  runLumiSelector = RunLumiSelector(jsonFiles);


  bool puInitialized = eventWeight->InitPUWeightsS10Truth(puWeightsFileData, puWeightsFileMc);
 
  if (!puInitialized) {
    std::cout << "PU weighting files are not found..." << std::endl;
    exit(-1);
  }

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

     _tree->SetBranchAddress("PtFirstH1",&PtFirstH1);
     _tree->SetBranchAddress("PxFirstH1",&PxFirstH1);
     _tree->SetBranchAddress("PyFirstH1",&PyFirstH1);
     _tree->SetBranchAddress("PzFirstH1",&PzFirstH1);
     _tree->SetBranchAddress("PhiFirstH1",&PhiFirstH1);
     _tree->SetBranchAddress("EtaFirstH1",&EtaFirstH1);
     _tree->SetBranchAddress("MassFirstH1",&MassFirstH1);

     _tree->SetBranchAddress("PtSecondH1",&PtSecondH1);
     _tree->SetBranchAddress("PxSecondH1",&PxSecondH1);
     _tree->SetBranchAddress("PySecondH1",&PySecondH1);
     _tree->SetBranchAddress("PzSecondH1",&PzSecondH1);
     _tree->SetBranchAddress("PhiSecondH1",&PhiSecondH1);
     _tree->SetBranchAddress("EtaSecondH1",&EtaSecondH1);
     _tree->SetBranchAddress("MassSecondH1",&MassSecondH1);


   }

   if (!isData) {
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

     if (isHto4taus)
       _tree->SetBranchAddress("H2pt",&HiggsPt);

   }
   

   // triggers =====>
   _tree->SetBranchAddress("HLT_Mu17_Mu8",&HLT_Mu17_Mu8);
   _tree->SetBranchAddress("HLT_Mu17_TkMu8",&HLT_Mu17_TkMu8);

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

   _tree->SetBranchAddress("NegLep_IsoMu24",&NegLep_IsoMu24);
   _tree->SetBranchAddress("PosLep_IsoMu24",&PosLep_IsoMu24);
   _tree->SetBranchAddress("NegLep_Mu20",&NegLep_Mu20);
   _tree->SetBranchAddress("PosLep_Mu20",&PosLep_Mu20);
   _tree->SetBranchAddress("NegLep_Mu30",&NegLep_Mu30);
   _tree->SetBranchAddress("PosLep_Mu30",&PosLep_Mu30);
   _tree->SetBranchAddress("NegLep_Mu40",&NegLep_Mu40);
   _tree->SetBranchAddress("PosLep_Mu40",&PosLep_Mu40);
   _tree->SetBranchAddress("NegLep_Mu17_Mu8",&NegLep_Mu17_Mu8);
   _tree->SetBranchAddress("PosLep_Mu17_Mu8",&PosLep_Mu17_Mu8);
   _tree->SetBranchAddress("NegLep_Mu17_TkMu8",&NegLep_Mu17_TkMu8);
   _tree->SetBranchAddress("PosLep_Mu17_TkMu8",&PosLep_Mu17_TkMu8);

   _tree->SetBranchAddress("DiLepPt",&DiLepPt);
   _tree->SetBranchAddress("DiLepPhi",&DiLepPhi);
  
   
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
   _tree->SetBranchAddress("jetPx",jetPx);
   _tree->SetBranchAddress("jetPy",jetPy);
   _tree->SetBranchAddress("jetPz",jetPz);
   _tree->SetBranchAddress("jetEn",jetEn);
   _tree->SetBranchAddress("jetMass",jetMass);
   _tree->SetBranchAddress("combSVBJetTag",combSVBJetTag);
   _tree->SetBranchAddress("puJetFullLoose", jetPUIDLoose);
   _tree->SetBranchAddress("chargedMultiplicity",chargedMultiplicity);

    
   int numberOfCandidates = _tree->GetEntries();
  
   std::cout << "      number of selected candidates = " << numberOfCandidates << std::endl;
   
   for (int iCand=0; iCand<numberOfCandidates; iCand++) { 
    
      _tree->GetEntry(iCand);

      if (!isData) {
	if (muMomScale<-0.001||muMomScale>0.001) {
	  //      	std::cout << std::endl;
	  //      	std::cout << "Before : Pos = " << _posLepPt << "  Neg = " << _negLepPt << std::endl;
	  //     	std::cout << "Muon momentum rescaled" << std::endl;
	  _negLepPt *= (1+muMomScale);
	  _posLepPt *= (1+muMomScale);
	  //      	std::cout << "After : Pos = " << _posLepPt << "  Neg = " << _negLepPt << std::endl;
	  //      	std::cout << std::endl;
	  for (int itrack=0; itrack<NpfNoPU; itrack++){
	    PxpfNoPU[itrack] *= (1+muMomScale);
	    PypfNoPU[itrack] *= (1+muMomScale);
	    PzpfNoPU[itrack] *= (1+muMomScale);                                                                                       
	    PtpfNoPU[itrack] *= (1+muMomScale);                                                                                        
	  }                                                                                                                             
	}

	//      if (trkMomScale<-0.001||trkMomScale>0.001) {
	//      	std::cout << "Track momentum rescaled" << std::endl;
	//      	for (int itrack=0; itrack<NpfNoPU; itrack++){
	//      	  PxpfNoPU[itrack] *= (1+trkMomScale);
	//      	  PypfNoPU[itrack] *= (1+trkMomScale);
	//      	  PzpfNoPU[itrack] *= (1+trkMomScale);
	//	  PtpfNoPU[itrack] *= (1+trkMomScale);
	//      	}
      }
      
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


      if (randomizeMuons>0) {

	float posLepPhiPlusPi2 = _posLepPhi + 0.5*TMath::Pi();
	float unitX = _posLepPt * TMath::Cos( posLepPhiPlusPi2 );
	float unitY = _posLepPt * TMath::Sin( posLepPhiPlusPi2 );
	float pxNeg = _negLepPt * TMath::Cos( _negLepPhi );
	float pyNeg = _negLepPt * TMath::Sin( _negLepPhi ); 
	float prod  = pxNeg*unitX + pyNeg*unitY;

	bool randomize = prod<0.0;

	if (randomizeMuons==2) {

	  if (_posLepPt<_negLepPt && (event%2==0))
	    randomize = true;
	  if (_posLepPt>_negLepPt && (event%2!=0))
	    randomize = true;

	}

	if (randomize) { // swap muons

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

      TLorentzVector DiMuon4 = MuonPos4 + MuonNeg4;
      float dimuonMass = float(DiMuon4.M());

      int muonPFyes=0;
      int muonPFno=0;
      int muontrackyes=0;
      int muontrackno=0;

      int trackPFyes=0;
      int trackPFno=0;
      int recotrackyes = 0;
      int recotrackno = 0;

      // generator studies

      TLorentzVector genMuonFirstH1;
      TLorentzVector genTrackFirstH1;
      TLorentzVector genVisTauFirstH1; genVisTauFirstH1.SetXYZT(0,0,0,0);

      bool firstMuonNotFound = true;
      bool firstTrackNotFound = true;

      TLorentzVector genMuonSecondH1;
      TLorentzVector genTrackSecondH1;
      TLorentzVector genVisTauSecondH1; genVisTauSecondH1.SetXYZT(0,0,0,0);

      bool secondMuonNotFound = true;
      bool secondTrackNotFound = true;

     // bool is1prongFirst = true;
      int nChargedFirst = 0;

      //bool is1prongSecond = true;
      int nChargedSecond  = 0;

      if (isHto4taus) {

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

	//	std::cout << std::endl;
	//	std::cout << "NPartFirstH1 = " << NPartFirstH1 << std::endl;

	float ptTrackMax = 0;

	for (int iP=0; iP<NPartFirstH1; ++iP) {
	  
	  //	  std::cout << iP << " pdgId = " << IdPartFirstH1[iP] << std::endl;
	  
	  if (QPartFirstH1[iP]!=0)
	    nChargedFirst++;

	  int partId = TMath::Abs(IdPartFirstH1[iP]);

	  if ( (partId!=11) && (partId!=13) && (partId!=15) &&
	       (partId!=12) && (partId!=14) && (partId!=16)) {

	    float massPart = pionMass;
	    if (partId==321)
	      massPart = kaonMass;
	    if (partId==130)
	      massPart = kzeroMass;
	    if (partId==310)
	      massPart = kzeroMass;
	    if (partId==311)
	      massPart = kzeroMass;
	    if (partId==111)
	      massPart = pzeroMass;
	    if (partId==22)
	      massPart = 0;

	    TLorentzVector part4mom; part4mom.SetXYZM(PxPartFirstH1[iP],
						      PyPartFirstH1[iP],
						      PzPartFirstH1[iP],
						      massPart);

	    genVisTauFirstH1 += part4mom;
	    
	  }

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

	  if (isMuon) {
	    nMuFirstGen++;
	    ptMuFirstGen = PtPartFirstH1[iP];
	    etaMuFirstGen = EtaPartFirstH1[iP];
	    phiMuFirstGen = PhiPartFirstH1[iP];
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

          if (isMuon && firstMuonNotFound) {
            genMuonFirstH1.SetXYZM(PxPartFirstH1[iP],
                                   PyPartFirstH1[iP],
                                   PzPartFirstH1[iP],
                                   muonMass);
            firstMuonNotFound = false;
          }

          if (isPion && PtPartFirstH1[iP]>ptTrackMax) {
	    genTrackFirstH1.SetXYZM(PxPartFirstH1[iP],
				    PyPartFirstH1[iP],
				    PzPartFirstH1[iP],
				    pionMass);
	    ptTrackMax = PtPartFirstH1[iP];
            firstTrackNotFound = false;
          }

          if (isKaon && PtPartFirstH1[iP]>ptTrackMax) {
            genTrackFirstH1.SetXYZM(PxPartFirstH1[iP],
				    PyPartFirstH1[iP],
				    PzPartFirstH1[iP],
				    kaonMass);
	    ptTrackMax = PtPartFirstH1[iP];
            firstTrackNotFound = false;
          }
	  
	}

	//	std::cout << "nChargedFirst = " << nChargedFirst << std::endl;
	//	std::cout << std::endl;
	//	std::cout << "NPartSecondH1 = " << NPartSecondH1 << std::endl;

	ptTrackMax = 0;

	for (int iP=0; iP<NPartSecondH1; ++iP) {
	  
	  //	  std::cout << iP << " pdgId = " << IdPartSecondH1[iP] << std::endl;

          int partId = TMath::Abs(IdPartSecondH1[iP]);

	  if (QPartSecondH1[iP]!=0)
	    nChargedSecond++;

          if ( (partId!=11) && (partId!=13) && (partId!=15) &&
               (partId!=12) && (partId!=14) && (partId!=16)) {

            float massPart = pionMass;
            if (partId==321)
              massPart = kaonMass;
            if (partId==130)
              massPart = kzeroMass;
            if (partId==310)
              massPart = kzeroMass;
            if (partId==311)
              massPart = kzeroMass;
            if (partId==111)
              massPart = pzeroMass;
	    if (partId==22)
	      massPart = 0;

            TLorentzVector part4mom; part4mom.SetXYZM(PxPartSecondH1[iP],
                                                      PyPartSecondH1[iP],
                                                      PzPartSecondH1[iP],
                                                      massPart);

            genVisTauSecondH1 += part4mom;

          }

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

	  if (isMuon) {
	    nMuSecondGen++;
	    ptMuSecondGen = PtPartSecondH1[iP];
	    etaMuSecondGen = EtaPartSecondH1[iP];
	    phiMuSecondGen = PhiPartSecondH1[iP];
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

          if (isMuon && secondMuonNotFound) {
            genMuonSecondH1.SetXYZM(PxPartSecondH1[iP],
				    PyPartSecondH1[iP],
				    PzPartSecondH1[iP],
				    muonMass);
            secondMuonNotFound = false;
          }

          if (isPion && PtPartSecondH1[iP]>ptTrackMax) {
            genTrackSecondH1.SetXYZM(PxPartSecondH1[iP],
				     PyPartSecondH1[iP],
				     PzPartSecondH1[iP],
				     pionMass);
	    ptTrackMax = PxPartSecondH1[iP];
            secondTrackNotFound = false;
          }

          if (isKaon && PtPartSecondH1[iP]>ptTrackMax) {
            genTrackSecondH1.SetXYZM(PxPartSecondH1[iP],
				     PyPartSecondH1[iP],
				     PzPartSecondH1[iP],
				     kaonMass);
	    ptTrackMax = PxPartSecondH1[iP];
            secondTrackNotFound = false;
          }

	}

	//	std::cout << std::endl;
	//	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	//	std::cout << std::endl;

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
 
	for (int iP=0; iP<NPartFirstH1; ++iP) {

	  float idPart = 0.5;

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

	  if (isMuon)
	    idPart = 1.5;
	  if (isElec)
	    idPart = 2.5;
	  if (isPion)
	    idPart = 3.5;
	  if (isKaon)
	    idPart = 4.5;

	  //	  std::cout << iP << " Id = " << IdPartFirstH1[iP] 
	  //		    << " ; Q = " << QPartFirstH1[iP] 
	  //		    << " ; Pt = " << PtPartFirstH1[iP]
	  //		    << " ; Eta = " << EtaPartFirstH1[iP] << std::endl;

	  if (EtaPartFirstH1[iP]<2.4 && EtaPartFirstH1[iP]>-2.4 && QPartFirstH1[iP]!=0 ) {
	    float absEtaPart = TMath::Abs(EtaPartFirstH1[iP]);
	    if (isMuon) 
	      PtGenMuonsEtaCutH->Fill(PtPartFirstH1[iP]);
	    else
	      PtGenTracksEtaCutH->Fill(PtPartFirstH1[iP]); 

	    if (PtPartFirstH1[iP]>ptGenPartCut) {

	      particleIdH->Fill(idPart);

	      int ptBin = binNumber(PtPartFirstH1[iP],nPtBins,ptBins);
	      int etaBin = binNumber(absEtaPart,nEtaBins,etaBins);

	      TLorentzVector genPart;
	      if (isMuon) 
		genPart.SetXYZM(PxPartFirstH1[iP],PyPartFirstH1[iP],PzPartFirstH1[iP],muonMass);
	      else
		genPart.SetXYZM(PxPartFirstH1[iP],PyPartFirstH1[iP],PzPartFirstH1[iP],pionMass);
	      
	      // ***************************
	      // running over Trk collection
	      // ***************************
	      float Diff = 1;
	      float DR = 1;
	      for (int iT=0; iT<NPVtracks; ++iT) {
		if (PtTrack[iT]>ptTrkCut) {
		  TLorentzVector track;
		  track.SetXYZM(PxTrack[iT],PyTrack[iT],PzTrack[iT],pionMass);
		  float phiTrk = track.Phi();
		  float etaTrk = track.Eta();
		  TLorentzVector diff = genPart - track;
		  float relDiff = diff.P()/genPart.P();
		  if (relDiff<Diff) {
		    Diff = relDiff;
		    DR = float(deltaR(EtaPartFirstH1[iP],PhiPartFirstH1[iP],
				      etaTrk,phiTrk));
		  }
		}
	      }
	      
	      if (isMuon) {
		if (Diff<0.05) {
		  GenMuonsInTrkColH->Fill(1.);
		  muontrackyes++;
		}
		else {
		  GenMuonsInTrkColH->Fill(0.);
		  muontrackno++;
		}
		dPGenMuonsTrkH->Fill(Diff);
		dRGenMuonsTrkH->Fill(DR);
	      }
	      else if (isPion || isKaon) {
		if (Diff<0.05){
		  GenTracksInTrkColH->Fill(1.);
		  recotrackyes++;
		  GenTracksInTrkColPtBinsH[ptBin]->Fill(1.);
		  GenTracksInTrkColEtaBinsH[etaBin]->Fill(1.);
		  GenTracksInTrkColPtEtaBinsH[ptBin][etaBin]->Fill(1.);
		}
		else{
		  GenTracksInTrkColH->Fill(0.);
                  GenTracksInTrkColPtBinsH[ptBin]->Fill(0.);
                  GenTracksInTrkColEtaBinsH[etaBin]->Fill(0.);
                  GenTracksInTrkColPtEtaBinsH[ptBin][etaBin]->Fill(0.);
		  recotrackno++;
		}
		dPGenTracksTrkH->Fill(Diff);
		dRGenTracksTrkH->Fill(DR);
	      }
	      // **********************************
	      // end of running over Trk collection
	      // **********************************
	      

	      // **********************************
	      // running over PfNoPileUp collection
	      // **********************************
	      Diff = 1;
	      DR   = 1; 
	      int iMatchedPfNoPu = 0;
	      for (int iT=0; iT<NpfNoPU; ++iT) {
		if (PtpfNoPU[iT]>ptTrkCut && (ChargepfNoPU[iT]<-0.5||ChargepfNoPU[iT]>0.5)) {
		  TLorentzVector track;
		  track.SetXYZM(PxpfNoPU[iT],PypfNoPU[iT],PzpfNoPU[iT],pionMass);
		  float phiTrk = track.Phi();
		  float etaTrk = track.Eta();
		  TLorentzVector diff = genPart - track;
		  float relDiff = diff.P()/genPart.P();
		  if (relDiff<Diff) {
		    Diff = relDiff;
		    DR = float(deltaR(EtaPartFirstH1[iP],PhiPartFirstH1[iP],
				      etaTrk,phiTrk));
		    iMatchedPfNoPu = iT;
		  }
		}
	      }
	      
	      float IP = ippfNoPU[iMatchedPfNoPu];
	      float IPZ = ipZpfNoPU[iMatchedPfNoPu];
	      if (IP>0.5) IP = 0.499;
	      if (IP<-0.5) IP = -0.499;
	      if (IPZ>1) IPZ = 0.999;
	      if (IPZ<-1) IPZ = -0.999;
	      

	      if (isMuon) {
		if (Diff<0.05) {
		  GenMuonsInPfNoPuColH->Fill(1.);
		  IpGenMuonsH->Fill(IP);
		  IpZGenMuonsH->Fill(IPZ); 
		  muonPFyes++;
		}
		else {
		  GenMuonsInPfNoPuColH->Fill(0.);
		  muonPFno++;
		}
		dPGenMuonsPfNoPuH->Fill(Diff);
		dRGenMuonsPfNoPuH->Fill(DR);
	      }
	      else if ( isPion || isKaon ) {
		if (Diff<0.05) {
		  GenTracksInPfNoPuColH->Fill(1.);
                  GenTracksInPfNoPuColPtBinsH[ptBin]->Fill(1.);
                  GenTracksInPfNoPuColEtaBinsH[etaBin]->Fill(1.);
                  GenTracksInPfNoPuColPtEtaBinsH[ptBin][etaBin]->Fill(1.);
		  IpGenTracksH->Fill(IP);
		  IpZGenTracksH->Fill(IPZ);
		  trackPFyes++;
		}
		else{
		  GenTracksInPfNoPuColH->Fill(0.);
		  GenTracksInPfNoPuColPtBinsH[ptBin]->Fill(0.);
                  GenTracksInPfNoPuColEtaBinsH[etaBin]->Fill(0.);
                  GenTracksInPfNoPuColPtEtaBinsH[ptBin][etaBin]->Fill(0.);
		  trackPFno++;
		}
		dPGenTracksPfNoPuH->Fill(Diff);
		dRGenTracksPfNoPuH->Fill(DR);
	      }

	      // *****************************************
	      // end of running over PfNoPileUp collection
	      // *****************************************
	      

	      if (isMuon) {
		Diff = 1;
	      
		TLorentzVector diffPos = MuonPos4 - genPart;
		TLorentzVector diffNeg = MuonNeg4 - genPart;
		
		float DiffPos = diffPos.P()/genPart.P();
		float DiffNeg = diffNeg.P()/genPart.P();
		
		if (DiffPos<DiffNeg) {
		  if (DiffPos<0.05) {
		    GenMuonsMatchedRecoMuonsH->Fill(1.);
		    float muIP = Ip2DPos[iPVtx];
		    float muIPZ = IpZPos[iPVtx];
		    if (muIP>0.5) muIP = 0.499;
		    if (muIP<-0.5) muIP = -0.499;
		    if (muIPZ>1) muIPZ = 0.999;
		    if (muIPZ<-1) muIPZ = -0.999;
		    IpGenMuonsMatchedRecoMuonsH->Fill(muIP);
		    IpZGenMuonsMatchedRecoMuonsH->Fill(muIPZ);
		  }
		  else 
		    GenMuonsMatchedRecoMuonsH->Fill(0.);
		  dPGenMuonsMatchedRecoMuonsH->Fill(DiffPos);
		}
		else {
		  if (DiffNeg<0.05) {
		    GenMuonsMatchedRecoMuonsH->Fill(1.);
		    float muIP = Ip2DNeg[iPVtx];
		    float muIPZ = IpZNeg[iPVtx];
		    if (muIP>0.5) muIP = 0.499;
		    if (muIP<-0.5) muIP = -0.499;
		    if (muIPZ>1) muIPZ = 0.999;
		    if (muIPZ<-1) muIPZ = -0.999;
		    IpGenMuonsMatchedRecoMuonsH->Fill(muIP);
		    IpZGenMuonsMatchedRecoMuonsH->Fill(muIPZ);
		  }
		  else
		    GenMuonsMatchedRecoMuonsH->Fill(0.);
		  dPGenMuonsMatchedRecoMuonsH->Fill(DiffNeg);
		}
	      }
	    }
	  }
	}
      }

      // cut on Z vertex;
      if (TMath::Abs(zPV[iPVtx])>24.0) continue;

      double muontracks_eff=(double)muontrackyes/(double)(muontrackyes+muontrackno);
      double recotracks_eff=(double)recotrackyes/(double)(recotrackyes+recotrackno);
      muontracks_effH->Fill(muontracks_eff);
      recotracks_effH->Fill(recotracks_eff);
      double muonPF_eff=(double)muonPFyes/(double)(muonPFyes+muonPFno);
      double trackPF_eff=(double)trackPFyes/(double)(trackPFyes+trackPFno);
      muonPF_effH->Fill(muonPF_eff);
      trackPF_effH->Fill(trackPF_eff);
	
	
      if (applyPUReweighting && !isData) // applying PU reweighting
	{
	  float ratioF = 1;
	  
	  ratioF = eventWeight->PUWeightS10Truth(float(reweightWithPUI ? nPUITruth : nPV));
	  
	  //std::cout << "apply PU weight: "<< ratioF <<std::endl;
	  weight *= ratioF;
	  
	}
		
      if (!isData && applyMuonIdSF) {

	float posMuPt = _posLepPt;
	if (posMuPt>1000)
	  posMuPt = 999;

	float negMuPt = _negLepPt;
        if (negMuPt>1000)
          negMuPt = 999;

	float posMuEta = TMath::Abs(_posLepEta);
	if (posMuEta>2.4)
	  posMuEta = 2.3;

        float negMuEta = TMath::Abs(_negLepEta);
        if (negMuEta>2.4)
          negMuEta = 2.3;

	float wPos = muonIdSF->GetBinContent(muonIdSF->FindBin(posMuPt,posMuEta));
	float wNeg = muonIdSF->GetBinContent(muonIdSF->FindBin(negMuPt,negMuEta));

	//	std::cout << "posMuPt = " << _posLepPt << " ; posMuEta = " << _posLepEta << " ; sf = " << wPos << std::endl;
	//	std::cout << "negMuPt = " << _negLepPt << " ; negMuEta = " << _negLepEta << " ; sf = " << wNeg << std::endl;
	//	std::cout << std::endl;

	weight *= wPos;
	weight *= wNeg;

      }


      if (applyPTReweighting && isPtFileValid && !isData)
	{
	  //cout<< " Spectrum  "<<HiggsPt<<"  "<<isPtFileValid<<"  "<<endl;
	  
	  int lastbin = PythiaH2PtWeights->FindLastBinAbove(0);
	  float range_ = lastbin * PythiaH2PtWeights->GetBinWidth(lastbin);
	  float wgt = 1;
	  if (HiggsPt > range_) wgt=1;
	  else {
	    
	    int FBin = PythiaH2PtWeights->FindBin(HiggsPt);
	    wgt = PythiaH2PtWeights->GetBinContent(FBin);
	    //cout<<  "   Pt reweighting  "<<wgt <<" PU weight "<<weight<<"  "<<HiggsPt<<"  "<<FBin<<"  "<<range_<<"  "<<scale_<<endl;
	    weight *=wgt;
	  }
	}
      
      //Trigger selection

      bool oneOfTwoFired = HLT_Mu17_Mu8 || HLT_Mu17_TkMu8;

      if (applyTrigger && !oneOfTwoFired) continue;
      counterTriggerH->Fill(1.);

      //######### Beginning of cuts ###########################################

      // Cosmics rejection
      //      float deltaPt = 2*fabs(_posLepPt-_negLepPt)/(_posLepPt+_negLepPt);
      //      float deltaEta = fabs(_posLepEta+_negLepEta);
      //      if ( (deltaPt<0.25) && (deltaEta<0.002) && (_diLepDPhi>3.14)) continue;
      

      // dR between two muons
      float muonDR=(float)deltaR((double)_negLepEta,(double)_negLepPhi,(double)_posLepEta,(double)_posLepPhi);
      
	
      float ptMax = _posLepPt;
      float ptMin = _negLepPt;

      bool posHard = true;
      
      float etaMax = _posLepEta;
      float etaMin = _negLepEta;
      
      if (ptMax<ptMin) {
	posHard = false;
	ptMax = _negLepPt;
	ptMin = _posLepPt;
	etaMax = _negLepEta;
	etaMin = _posLepEta;
      }
	
      float absEtaMax = TMath::Abs(etaMax);
      float absEtaMin = TMath::Abs(etaMin);

      // and now apply muon Id --->

      if (applyMuonId) {
	
	if (!isGlobalMuPos) continue;
	if (!isGlobalMuNeg) continue;
	
	if(!isPFMuPos)continue;
	if(!isPFMuNeg)continue;
	
	float chi2ndfPos = Chi2PosMu/NdofPosMu; 
	float chi2ndfNeg = Chi2NegMu/NdofNegMu; 
	
	if (chi2ndfPos>10) continue;
	if (chi2ndfNeg>10) continue;
	
	if (nMuonHitsPos<1) continue;
	if (nMuonHitsNeg<1) continue;

	if (nMuonStationsPos<2) continue;
	if (nMuonStationsNeg<2) continue;	

	if (nPixelHitsPos<1) continue;
	if (nPixelHitsNeg<1) continue;
	
	// 	if (nTrackerHitsPos<11)  continue; 
	// 	if (nTrackerHitsNeg<11)  continue; 

	if (nTrackerHitsPos<5) continue;
	if (nTrackerHitsNeg<5) continue;

	//	if (Ip2DPos[iPVtx]>0.02||Ip2DPos[iPVtx]<-0.02) continue;
	//	if (Ip2DNeg[iPVtx]>0.02||Ip2DNeg[iPVtx]<-0.02) continue;
	// 	if (IpZPos[iPVtx]>0.2||IpZPos[iPVtx]<-0.2) continue; 
	// 	if (IpZNeg[iPVtx]>0.2||IpZNeg[iPVtx]<-0.2) continue; 

	bool selectMuons = 
	  (muonDR > dRmumuCut) && 
	  (absEtaMax<MuEtaTriggerCut) && 
	  (absEtaMin<MuEtaCut) &&
	  (ptMax>ptHardLepCut) && 
	  (ptMin>ptSoftLepCut) ;
	if (selectMuons) {

	  float ipSig2DPos = IpSig2DPos[iPVtx];
	  float ipSig2DNeg = IpSig2DNeg[iPVtx];

	  float ipSig3DPos = IpSig3DPos[iPVtx];
	  float ipSig3DNeg = IpSig3DNeg[iPVtx];

	  float ipSigZPos = IpSigZPos[iPVtx];
	  float ipSigZNeg = IpSigZNeg[iPVtx];

	  float ip2DPos = Ip2DPos[iPVtx];
	  float ip2DNeg = Ip2DNeg[iPVtx];

	  float ipZPos = IpZPos[iPVtx];
	  float ipZNeg = IpZNeg[iPVtx];
	  
	  if(ipSig2DPos<-5.0)ipSig2DPos=-4.99;
	  if(ipSig2DPos>5.0)ipSig2DPos=4.99;
	  if(ipSig2DNeg<-5.0)ipSig2DNeg=-4.99;
	  if(ipSig2DNeg>5.0)ipSig2DNeg=4.99;
	  IpSig2DPosH->Fill(ipSig2DPos,weight);
	  IpSig2DNegH->Fill(ipSig2DNeg,weight);
	  IpSig2D_2DH->Fill(ipSig2DNeg,ipSig2DPos,weight);
	  
	  if(ipSig3DPos<-5.0)ipSig3DPos=-4.99;
	  if(ipSig3DPos>5.0)ipSig3DPos=4.99;
	  if(ipSig3DNeg<-5.0)ipSig3DNeg=-4.99;
	  if(ipSig3DNeg>5.0)ipSig3DNeg=4.99;
	  IpSig3DPosH->Fill(ipSig3DPos,weight);
	  IpSig3DNegH->Fill(ipSig3DNeg,weight);
	  
	  if(ipSigZPos<-5.0)ipSigZPos=-4.99;
	  if(ipSigZPos>5.0)ipSigZPos=4.99;
	  if(ipSigZNeg<-5.0)ipSigZNeg=-4.99;
	  if(ipSigZNeg>5.0)ipSigZNeg=4.99;
	  IpSigZPosH->Fill(ipSigZPos,weight);
	  IpSigZNegH->Fill(ipSigZNeg,weight);

	  if (ip2DPos<-0.5) ip2DPos = -0.499;
	  if (ip2DNeg>0.5) ip2DNeg = 0.499;
	  Ip2DPosH->Fill(ip2DPos,weight);
	  Ip2DNegH->Fill(ip2DNeg,weight);
	  Ip2DBothH->Fill(ip2DNeg,weight);
	  Ip2DBothH->Fill(ip2DPos,weight);

	  if (ipZPos<-1.) ipZPos = -0.999;
	  if (ipZNeg>1.) ipZNeg = 0.999;
	  IpZPosH->Fill(ipZPos,weight);
	  IpZNegH->Fill(ipZNeg,weight);
	  IpZBothH->Fill(ipZNeg,weight);
	  IpZBothH->Fill(ipZPos,weight);


	}

	if (TMath::Abs(Ip2DPos[iPVtx])>d0Cut) continue;
	if (TMath::Abs(Ip2DNeg[iPVtx])>d0Cut) continue;

	if (TMath::Abs(IpZPos[iPVtx])>dzCut) continue;
	if (TMath::Abs(IpZNeg[iPVtx])>dzCut) continue;
	
      }

      HTauTauUtils utils;

      float muonDEta = TMath::Abs(_posLepEta-_negLepEta);
      float muonDPhi = utils.dPhiFrom2Phi(_posLepPhi,_negLepPhi);

      counterMuonIDH->Fill(1.);

      if ( absEtaMax>MuEtaTriggerCut ) continue;
      if ( absEtaMin>MuEtaCut ) continue;
      counterEtaH->Fill(1.);
      
      // cut on pt of muons
      if (ptMax<ptHardLepCut) continue;
      if (ptMin<ptSoftLepCut) continue;
      counterPtH->Fill(1.);
      
      dRdPhiMuMu_DimuonsH->Fill(muonDR,muonDPhi,weight);
      dRMuMu_DimuonsH->Fill(muonDR,weight);
      dPhiMuMu_DimuonsH->Fill(muonDPhi,weight);
      diMuonMass_DimuonsH->Fill(dimuonMass,weight);

      if (cutType) {
	if (muonDR<dRmumuCut) continue;
      }
      else {
	if (muonDPhi<dRmumuCut) continue;
      }
	  
      countermuonDRH->Fill(1.);


      HardMuonPt_DimuonsH->Fill(ptMax,weight);
      SoftMuonPt_DimuonsH->Fill(ptMin,weight);
      HardMuonEta_DimuonsH->Fill(etaMax,weight);
      SoftMuonEta_DimuonsH->Fill(etaMin,weight);
      HardMuonPtSoftMuonPt_DimuonsH->Fill(ptMax,ptMin,weight);
      


      //*****************************************No PU collection********************************* 
      
      AllNpfNoPU=0;
      AllNpfNoPU_ip=0;
      NpfNoPuPos = 0;
      NpfNoPuNeg = 0;
      NpfNoPuIpPos = 0;
      NpfNoPuIpNeg = 0;
      NpfNoPuIpPosAll  = 0;
      NpfNoPuIpNegAll  = 0;
      NpfNoPuIpPosSoft = 0;
      NpfNoPuIpNegSoft = 0;
      

      Ntracks_OutConePos= 0;
      Ntracks_OutConeNeg= 0;
      Ntracks_OutConePosIP= 0;
      Ntracks_OutConeNegIP= 0;
      SumPtPosCone = 0;
      SumPtNegCone = 0;
      SumPtPosConeIP = 0;
      SumPtNegConeIP = 0;


      float totPosIso = 0;
      float totNegIso = 0;

      float dRminPos = 1e+6;
      float dRminNeg = 1e+6;

      for (int itrack=0; itrack<NpfNoPU; itrack++){
	//	std::cout << PtpfNoPU[itrack] << std::endl;
	if (PtpfNoPU[itrack]>ptTrkCutIso && (ChargepfNoPU[itrack]>0.5 ||ChargepfNoPU[itrack]<-0.5) && TMath::Abs(EtapfNoPU[itrack])<etaTrkCut) {
	  
	  TLorentzVector pfCand;
	  pfCand.SetXYZM(PxpfNoPU[itrack],PypfNoPU[itrack],PzpfNoPU[itrack],pionMass);
	  float EtaPfCand = pfCand.Eta();
	  float PhiPfCand = pfCand.Phi();
	  float dRpos = float(deltaR(_posLepEta,_posLepPhi,EtaPfCand,PhiPfCand));
	  TLorentzVector diffPos = MuonPos4 - pfCand;
	  if (dRpos<isoCone && diffPos.P()>0.1) {
	    totPosIso += PtpfNoPU[itrack];
	  }
	  float dRneg = float(deltaR(_negLepEta,_negLepPhi,EtaPfCand,PhiPfCand));
	  TLorentzVector diffNeg = MuonNeg4 - pfCand;
	  if (dRneg<isoCone && diffNeg.P()>0.1) {
	    totNegIso += PtpfNoPU[itrack];
	  }



	}
	if (PtpfNoPU[itrack]>ptTrkCut && (ChargepfNoPU[itrack]>0.5 ||ChargepfNoPU[itrack]<-0.5) && TMath::Abs(EtapfNoPU[itrack])<etaTrkCut) { 
	  AllNpfNoPU++;
	  int idTrkI = TMath::Abs(idpfNoPU[itrack]);
	  float IP = ippfNoPU[itrack];
	  if (IP>0.5) IP = 0.499;
	  if (IP<-0.5) IP = -0.499;
	  ippfNoPUH->Fill(IP,weight);
	  float IPZ = ipZpfNoPU[itrack];
	  if (IPZ>1) IPZ = 0.999;
	  if (IPZ<-1) IPZ = -0.999;
	  ipZpfNoPUH->Fill(IPZ,weight);
	  float IPabs = TMath::Abs(ippfNoPU[itrack]);
	  float IPZabs = TMath::Abs(ipZpfNoPU[itrack]);
	  float px1 = PxpfNoPU[itrack];
	  float py1 = PypfNoPU[itrack];
	  float pz1 = PzpfNoPU[itrack];
	  if (IPabs<0.5 && IPZabs<0.5) {
	    for (int jtrack=0; jtrack<NpfNoPU; jtrack++){
	      if (itrack!=jtrack) {
		if (ChargepfNoPU[itrack]*ChargepfNoPU[jtrack]<-0.5) {
		  float IPabsX = TMath::Abs(ippfNoPU[jtrack]);
		  float IPZabsX = TMath::Abs(ipZpfNoPU[jtrack]);
		  if (IPabsX<0.5&&IPZabsX<0.5) {
		    int idTrkJ = TMath::Abs(idpfNoPU[jtrack]);
		    float px2 = PxpfNoPU[jtrack];
		    float py2 = PypfNoPU[jtrack];
		    float pz2 = PzpfNoPU[jtrack];

		    float massEleEle = invMass(px1,py1,pz1,EleMass,px2,py2,pz2,EleMass);
		    //float massEleMu  = invMass(px1,py1,pz1,EleMass,px2,py2,pz2,MuMass);
		    //float massEleHad = invMass(px1,py1,pz1,EleMass,px2,py2,pz2,PiMass);

		    //float massMuEle  = invMass(px1,py1,pz1,MuMass,px2,py2,pz2,EleMass);
		    float massMuMu   = invMass(px1,py1,pz1,MuMass,px2,py2,pz2,MuMass);
		    //float massMuHad  = invMass(px1,py1,pz1,MuMass,px2,py2,pz2,PiMass);

		    //float massHadEle = invMass(px1,py1,pz1,PiMass,px2,py2,pz2,EleMass);
		    //float massHadMu  = invMass(px1,py1,pz1,PiMass,px2,py2,pz2,MuMass);
		    float massHadHad = invMass(px1,py1,pz1,PiMass,px2,py2,pz2,PiMass);
		    
		    if (idTrkI==11) {
		      if (idTrkJ==11) {
			massJPsiEleEleH->Fill(massEleEle,weight);
			massYpsilonEleEleH->Fill(massEleEle,weight);
		      }
		      else if (idTrkJ==13) {
			massJPsiEleMuH->Fill(massEleEle,weight);
			massYpsilonEleMuH->Fill(massEleEle,weight);
		      }
		      else {
			massJPsiEleHadH->Fill(massEleEle,weight);
			massYpsilonEleHadH->Fill(massEleEle,weight);
		      }
		    }
		    else if (idTrkI==13) {
		      if (idTrkJ==11) {
			massJPsiMuEleH->Fill(massMuMu,weight);
			massYpsilonMuEleH->Fill(massMuMu,weight);
		      }
		      else if (idTrkJ==13) {
			massJPsiMuMuH->Fill(massMuMu,weight);
			massYpsilonMuMuH->Fill(massMuMu,weight);
		      }
		      else {
			massJPsiMuHadH->Fill(massMuMu,weight);
			massYpsilonMuHadH->Fill(massMuMu,weight);
		      }
		    }
		    else {
		      if (idTrkJ==11) {
			massK0HadEleH->Fill(massHadHad,weight);
		      }
		      else if (idTrkJ==13) {
			massK0HadMuH->Fill(massHadHad,weight);
		      }
		      else {
			massK0HadHadH->Fill(massHadHad,weight);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  TLorentzVector pfCand;
	  pfCand.SetXYZM(PxpfNoPU[itrack],PypfNoPU[itrack],PzpfNoPU[itrack],pionMass);
	  float EtaPfCand = pfCand.Eta();
	  float PhiPfCand = pfCand.Phi();
	  float dRpos = float(deltaR(_posLepEta,_posLepPhi,EtaPfCand,PhiPfCand));
	  TLorentzVector diffPos = MuonPos4 - pfCand;
	  if (dRpos<dRcone && diffPos.P()>0.1) {
	    //PxPfNoPuPos[NpfNoPuPos] = PxpfNoPU[itrack];
	    //PxPfNoPuPos[NpfNoPuPos] = PxpfNoPU[itrack];
	    /*PyPfNoPuPos[NpfNoPuPos] = PypfNoPU[itrack];
	    PzPfNoPuPos[NpfNoPuPos] = PzpfNoPU[itrack];
	    QPfNoPuPos[NpfNoPuPos] = ChargepfNoPU[itrack];
	    EnPfNoPuPos[NpfNoPuPos] = EnpfNoPU[itrack];
	    IPPfNoPuPos[NpfNoPuPos] = ippfNoPU[itrack];
	    IPZPfNoPuPos[NpfNoPuPos] = ipZpfNoPU[itrack];*/
	    NpfNoPuPos++;
	  }

	  //************* Second Cone for bjet rejection **************
	  if (dRpos>0.5 && dRpos<0.7 && diffPos.P()>0.1){
	    Ntracks_OutConePos++;
	    SumPtPosCone += PtpfNoPU[itrack];
	  }


	  float dRneg = float(deltaR(_negLepEta,_negLepPhi,EtaPfCand,PhiPfCand));
	  TLorentzVector diffNeg = MuonNeg4 - pfCand;
	  /*if (dRneg<dRcone&& diffNeg.P()>0.1) {
	    PtPfNoPuNeg[NpfNoPuNeg] = PtpfNoPU[itrack];
	    PxPfNoPuNeg[NpfNoPuNeg] = PxpfNoPU[itrack];
	    PyPfNoPuNeg[NpfNoPuNeg] = PypfNoPU[itrack];
	    PzPfNoPuNeg[NpfNoPuNeg] = PzpfNoPU[itrack];
	    QPfNoPuNeg[NpfNoPuNeg] = ChargepfNoPU[itrack];
	    EnPfNoPuNeg[NpfNoPuNeg] = EnpfNoPU[itrack];
	    IPPfNoPuNeg[NpfNoPuNeg] = ippfNoPU[itrack];
	    IPZPfNoPuNeg[NpfNoPuNeg] = ipZpfNoPU[itrack];
	    NpfNoPuNeg++;
	  }*/
	  //************* Second Cone for bjet rejection **************
	  if (dRneg>0.5 && dRneg<0.7 && diffNeg.P()>0.1){
	    Ntracks_OutConeNeg++;
	    SumPtNegCone += PtpfNoPU[itrack];
	  }

	  //************************* Tracks with cuts on the IP ***********************
	  if ( TMath::Abs(ipZpfNoPU[itrack])<dzCutLoose && TMath::Abs(ippfNoPU[itrack])<d0CutLoose)  { 
	    AllNpfNoPU_ip++;

	    if (diffPos.P()>0.1 && PtpfNoPU[itrack]>ptSoftTrkCut && (ChargepfNoPU[itrack]*_posLepQ<0) && 
		TMath::Abs(ippfNoPU[itrack])<d0CutTracks && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks) {
	      if (dRpos<dRminPos)
		dRminPos = dRpos;
	    }
	    if (diffNeg.P()>0.1 && PtpfNoPU[itrack]>ptSoftTrkCut && (ChargepfNoPU[itrack]*_negLepQ<0) && 
		TMath::Abs(ippfNoPU[itrack])<d0CutTracks && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks) {
	      if (dRneg<dRminNeg)
		dRminNeg = dRneg;
	    }
	    
	    if (dRpos<dRcone&& diffPos.P()>0.1) {
 	      if (PtpfNoPU[itrack]>ptSoftTrkCut && (ChargepfNoPU[itrack]*_posLepQ<0) && 
 		  TMath::Abs(ippfNoPU[itrack])<d0CutTracks && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks) {
//	      if (ChargepfNoPU[itrack]*_posLepQ<0) {
		PtPfNoPuIpPos[NpfNoPuIpPos] = PtpfNoPU[itrack];
		PxPfNoPuIpPos[NpfNoPuIpPos] = PxpfNoPU[itrack];
		PyPfNoPuIpPos[NpfNoPuIpPos] = PypfNoPU[itrack];
		PzPfNoPuIpPos[NpfNoPuIpPos] = PzpfNoPU[itrack];
		QPfNoPuIpPos[NpfNoPuIpPos] = ChargepfNoPU[itrack];
		//EnPfNoPuIpPos[NpfNoPuIpPos] = EnpfNoPU[itrack];
		IPPfNoPuIpPos[NpfNoPuIpPos] = ippfNoPU[itrack];
		IPZPfNoPuIpPos[NpfNoPuIpPos] = ipZpfNoPU[itrack];
		IdNoPuIpPos[NpfNoPuIpPos] = idpfNoPU[itrack];
		NpfNoPuIpPos++;
	      }
	      else {
		if ( PtpfNoPU[itrack]>ptTrkLower && PtpfNoPU[itrack]<ptTrkUpper )
		  NpfNoPuIpPosSoft++;
	      }
	      NpfNoPuIpPosAll++;
	    }
	    if (dRneg<dRcone && diffNeg.P()>0.1) {
 	      if (PtpfNoPU[itrack]>ptSoftTrkCut && (ChargepfNoPU[itrack]*_negLepQ<0) &&
 		  TMath::Abs(ippfNoPU[itrack])<d0CutTracks && TMath::Abs(ipZpfNoPU[itrack])<dzCutTracks) {
		//	      if (ChargepfNoPU[itrack]*_negLepQ<0) {
		PtPfNoPuIpNeg[NpfNoPuIpNeg] = PtpfNoPU[itrack];
		PxPfNoPuIpNeg[NpfNoPuIpNeg] = PxpfNoPU[itrack];
		PyPfNoPuIpNeg[NpfNoPuIpNeg] = PypfNoPU[itrack];
		PzPfNoPuIpNeg[NpfNoPuIpNeg] = PzpfNoPU[itrack];
		QPfNoPuIpNeg[NpfNoPuIpNeg] = ChargepfNoPU[itrack];
		//EnPfNoPuIpNeg[NpfNoPuIpNeg] = EnpfNoPU[itrack];
		IPPfNoPuIpNeg[NpfNoPuIpNeg] = ippfNoPU[itrack];
		IPZPfNoPuIpNeg[NpfNoPuIpNeg] = ipZpfNoPU[itrack];
		IdNoPuIpNeg[NpfNoPuIpNeg] = idpfNoPU[itrack];
		NpfNoPuIpNeg++;
	      }
	      else {
		if ( PtpfNoPU[itrack]>ptTrkLower && PtpfNoPU[itrack]<ptTrkUpper )
		  NpfNoPuIpNegSoft++;
	      }
	      NpfNoPuIpNegAll++;
	    }


	    //************* Second Cone for bjet rejection **************
	    if (dRpos>0.5 && dRpos<0.7 && diffPos.P()>0.1){
	      Ntracks_OutConePosIP++;
	      SumPtPosConeIP += PtpfNoPU[itrack];
	    }
	    
	    if (dRneg>0.5 && dRneg<0.7 && diffNeg.P()>0.1){
	      Ntracks_OutConeNegIP++;
	      SumPtNegConeIP += PtpfNoPU[itrack];
	    }

	  }
	}
      }

      TLorentzVector diffPosFirstMu  = MuonPos4 - genMuonFirstH1; 
      TLorentzVector diffPosSecondMu = MuonPos4 - genMuonSecondH1;

      TLorentzVector diffNegFirstMu  = MuonNeg4 - genMuonFirstH1;
      TLorentzVector diffNegSecondMu = MuonNeg4 - genMuonSecondH1;

      float resoPosFirst  = float(diffPosFirstMu.P())/float(MuonPos4.P());
      float resoPosSecond = float(diffPosSecondMu.P())/float(MuonPos4.P()); 

      float resoNegFirst  = float(diffNegFirstMu.P())/float(MuonNeg4.P());
      float resoNegSecond = float(diffNegSecondMu.P())/float(MuonNeg4.P());

      bool PosFirst = false;
      bool PosSecond = false;

      bool NegFirst = false;
      bool NegSecond = false;

      TLorentzVector posTrkGen4;
      TLorentzVector negTrkGen4;

      //      std::cout << std::endl;

      //      if (PMmotherExist) {
      //	std::cout << "Pos muon mother --->" << std::endl;
      //	std::cout << "PdgId :" << PMmotherpdg << std::endl;
      //	std::cout << "(Pt,Eta,Phi) = (" << PMmotherpt << "," << PMmothereta << "," << PMmotherphi << ")" << std::endl; 
      //      }
      //      else {
      //	std::cout << "mother of pos muon not found" << std::endl;
      //      }
	
      //      std::cout << std::endl;
      //      if (NMmotherExist) {
      //	std::cout << "Neg muon mother --->" << std::endl;
      //	std::cout << "PdgId :" << NMmotherpdg << std::endl;
      //	std::cout << "(Pt,Eta,Phi) = (" << NMmotherpt << "," << NMmothereta << "," << NMmotherphi << ")" << std::endl;
      //      }
      //      else {
      //	std::cout << "mother of pos muon not found" << std::endl;
      //      }

      //      std::cout << std::endl;
      //      std::cout << "+++++++++++++++++++++++++++++" << std::endl;

      TLorentzVector posTauGen4;
      TLorentzVector negTauGen4;

      if ( (resoPosSecond<resoPosFirst) && resoPosSecond<0.1) {
      	PosSecond = true;
      	posTrkGen4 = genTrackSecondH1;
	
      }

      if ( (resoPosFirst<resoPosSecond) && resoPosFirst<0.1) {
        PosFirst = true;
	posTrkGen4 = genTrackFirstH1;
      }

      if ( (resoNegSecond<resoNegFirst) && resoNegSecond<0.1) {
        NegSecond = true;
	negTrkGen4 = genTrackSecondH1;
      }

      if ( (resoNegFirst<resoNegSecond) && resoNegFirst<0.1) {
        NegFirst = true;
	negTrkGen4 = genTrackFirstH1;
      }

      // efficiency study --->
      if (isHto4taus ) {

	TLorentzVector firstH1; firstH1.SetXYZM(PxFirstH1,PyFirstH1,PzFirstH1,MassFirstH1);
	TLorentzVector secondH1; secondH1.SetXYZM(PxSecondH1,PySecondH1,PzSecondH1,MassSecondH1);

	TLorentzVector genPosMuTau4; genPosMuTau4.SetXYZM(PMmotherpx,PMmotherpy,PMmotherpz,tauMass);
	TLorentzVector genNegMuTau4; genNegMuTau4.SetXYZM(NMmotherpx,NMmotherpy,NMmotherpz,tauMass);

	if ( (PosFirst && !firstMuonNotFound && !firstTrackNotFound ) || (PosSecond && !secondMuonNotFound && !secondTrackNotFound ) ) {
	  float ptTrk = posTrkGen4.Pt();

	  if (ptTrk>99.9) ptTrk = 99.9; 

	  float ptTau = 0;
	  float etaTau = -10;
	  if (PosFirst) {
	    //	    ptTau  = (firstH1 - genPosMuTau4).Pt();
	    //	    etaTau = TMath::Abs((firstH1 - genPosMuTau4).Eta());
	    ptTau  = genVisTauFirstH1.Pt();
	    etaTau = genVisTauFirstH1.Eta();
	  }
	  if (PosSecond) {
	    //	    ptTau  = (secondH1 - genPosMuTau4).Pt(); 
	    //	    etaTau = TMath::Abs((secondH1 - genPosMuTau4).Eta());
	    ptTau  = genVisTauSecondH1.Pt();
            etaTau = genVisTauSecondH1.Eta();
	  }

	  if (ptTau>99.9) ptTau = 99.9;

	  if (etaTau<2.4) {
	  
	    counterPtPosMuTrkH->Fill(ptTrk,weight);
	    counterPtTrkH->Fill(ptTrk,weight);

	    counterPtPosMuTauH->Fill(ptTau,weight);
	    counterPtTauH->Fill(ptTau,weight); 

	    counterPtPosMuTrkTauH->Fill(ptTrk,ptTau,weight);
	    counterPtTrkTauH->Fill(ptTrk,ptTau,weight);

	    if (_posLepPt > _negLepPt) {
	      tauPtHardDimuonH->Fill(ptTau,weight);
	      trkPtHardDimuonH->Fill(ptTrk,weight);
	    }
	    else {
	      tauPtSoftDimuonH->Fill(ptTau,weight);
	      trkPtSoftDimuonH->Fill(ptTrk,weight);
	    }

	    bool foundTrack = false;

	    for (int iPart=0; iPart<NpfNoPuIpPos; iPart++) {
	      TLorentzVector part4; part4.SetXYZM(PxPfNoPuIpPos[iPart],
						  PyPfNoPuIpPos[iPart],
						  PzPfNoPuIpPos[iPart],
						  pionMass);
	      TLorentzVector diffPart = posTrkGen4 - part4;
	      float resoDR = deltaR(part4.Eta(),part4.Phi(),
				    posTrkGen4.Eta(),posTrkGen4.Phi());
	      if (resoDR<0.3)
		foundTrack = true;
	      
	    }

	    if (foundTrack) {
	      counterPtPosMuTrkIdH->Fill(ptTrk,weight);
	      counterPtTrkIdH->Fill(ptTrk,weight);
	      counterPtPosMuTauIdH->Fill(ptTau,weight);
	      counterPtTauIdH->Fill(ptTau,weight);
	      counterPtPosMuTrkTauIdH->Fill(ptTrk,ptTau,weight);
	      counterPtTrkTauIdH->Fill(ptTrk,ptTau,weight);
	      if (_posLepPt > _negLepPt) {
		tauPtHardDimuonIdH->Fill(ptTau,weight);
		trkPtHardDimuonIdH->Fill(ptTrk,weight);
	      }
	      else {
		tauPtSoftDimuonIdH->Fill(ptTau,weight);
		trkPtSoftDimuonIdH->Fill(ptTrk,weight);
	      }
	      if (NpfNoPuIpPosAll==1 && NpfNoPuIpPos==1) {
		counterPtPosMuTrkIsoH->Fill(ptTrk,weight);
		counterPtTrkIsoH->Fill(ptTrk,weight);
		counterPtPosMuTauIsoH->Fill(ptTau,weight);
		counterPtTauIsoH->Fill(ptTau,weight);
		counterPtPosMuTrkTauIsoH->Fill(ptTrk,ptTau,weight);
		counterPtTrkTauIsoH->Fill(ptTrk,ptTau,weight);
		if (_posLepPt > _negLepPt) {
		  tauPtHardDimuonIsoH->Fill(ptTau,weight);
		  trkPtHardDimuonIsoH->Fill(ptTrk,weight);
		}
		else {
		  tauPtSoftDimuonIsoH->Fill(ptTau,weight);
		  trkPtSoftDimuonIsoH->Fill(ptTrk,weight);
		}	    
		int idtrk = TMath::Abs(IdNoPuIpPos[0]);
		if (idtrk!=11 && idtrk!=13) {
		  counterPtPosMuTrkNonLepH->Fill(ptTrk,weight);
		  counterPtTrkNonLepH->Fill(ptTrk,weight);
		  counterPtPosMuTauNonLepH->Fill(ptTau,weight);
		  counterPtTauNonLepH->Fill(ptTau,weight);
		  counterPtPosMuTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
		  counterPtTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
		}
	      }
	    }
	  }
	}

        if ( (NegFirst && !firstMuonNotFound && !firstTrackNotFound ) || (NegSecond && !secondMuonNotFound && !secondTrackNotFound ) ) {
          float ptTrk = negTrkGen4.Pt();

	  if (ptTrk>99.9) ptTrk = 99.9;

          float ptTau = 0;
	  float etaTau = -10;
          if (NegFirst) {
	    //            ptTau  = (firstH1 - genNegMuTau4).Pt();
	    //	    etaTau = TMath::Abs((firstH1 - genNegMuTau4).Eta());
	    ptTau  = genVisTauFirstH1.Pt();
            etaTau = genVisTauFirstH1.Eta();
	  }
          if (NegSecond) {
	    //            ptTau  = (secondH1 - genNegMuTau4).Pt();
	    //	    etaTau = TMath::Abs((secondH1 - genNegMuTau4).Eta());
	    ptTau  = genVisTauSecondH1.Pt();
            etaTau = genVisTauSecondH1.Eta();
	  }

          if (ptTau>99.9) ptTau = 99.9;

	  if (etaTau<2.4) {

	    counterPtNegMuTrkH->Fill(ptTrk,weight);
	    counterPtTrkH->Fill(ptTrk,weight);

	    counterPtNegMuTauH->Fill(ptTau,weight);
	    counterPtTauH->Fill(ptTau,weight);

	    counterPtNegMuTrkTauH->Fill(ptTrk,ptTau,weight);
	    counterPtTrkTauH->Fill(ptTrk,ptTau,weight);

	    if (_negLepPt > _posLepPt) {
	      tauPtHardDimuonH->Fill(ptTau,weight);
	      trkPtHardDimuonH->Fill(ptTrk,weight);
	    }
	    else {
	      tauPtSoftDimuonH->Fill(ptTau,weight);
	      trkPtSoftDimuonH->Fill(ptTrk,weight);
	    }

	    bool foundTrack = false;

	    for (int iPart=0; iPart<NpfNoPuIpNeg; iPart++) {
	      TLorentzVector part4; part4.SetXYZM(PxPfNoPuIpNeg[iPart],
						  PyPfNoPuIpNeg[iPart],
						  PzPfNoPuIpNeg[iPart],
						  pionMass);
	      TLorentzVector diffPart = negTrkGen4 - part4;
	      float resoDR = deltaR(part4.Eta(),part4.Phi(),
				    negTrkGen4.Eta(),negTrkGen4.Phi());
	      if (resoDR<0.3)
		foundTrack = true;
	      
	    }

	    if (foundTrack) {
	      counterPtNegMuTrkIdH->Fill(ptTrk,weight);
	      counterPtTrkIdH->Fill(ptTrk,weight);
	      counterPtNegMuTauIdH->Fill(ptTau,weight);
	      counterPtTauIdH->Fill(ptTau,weight);
	      counterPtNegMuTrkTauIdH->Fill(ptTrk,ptTau,weight);
	      counterPtTrkTauIdH->Fill(ptTrk,ptTau,weight);
	      if (_negLepPt > _posLepPt) {
		tauPtHardDimuonIdH->Fill(ptTau,weight);
		trkPtHardDimuonIdH->Fill(ptTrk,weight);
	      }
	      else {
		tauPtSoftDimuonIdH->Fill(ptTau,weight);
		trkPtSoftDimuonIdH->Fill(ptTrk,weight);
	      }
	      if (NpfNoPuIpNegAll==1 && NpfNoPuIpNeg==1) {
		counterPtNegMuTrkIsoH->Fill(ptTrk,weight);
		counterPtTrkIsoH->Fill(ptTrk,weight);
		counterPtNegMuTauIsoH->Fill(ptTau,weight);
		counterPtTauIsoH->Fill(ptTau,weight);
		counterPtNegMuTrkTauIsoH->Fill(ptTrk,ptTau,weight);
		counterPtTrkTauIsoH->Fill(ptTrk,ptTau,weight);
		if (_negLepPt > _posLepPt) {
		  tauPtHardDimuonIsoH->Fill(ptTau,weight);
		  trkPtHardDimuonIsoH->Fill(ptTrk,weight);
		}
		else {
		  tauPtSoftDimuonIsoH->Fill(ptTau,weight);
		  trkPtSoftDimuonIsoH->Fill(ptTrk,weight);
		}
                int idtrk = TMath::Abs(IdNoPuIpNeg[0]);
                if (idtrk!=11 && idtrk!=13) {
                  counterPtNegMuTrkNonLepH->Fill(ptTrk,weight);
                  counterPtTrkNonLepH->Fill(ptTrk,weight);
                  counterPtNegMuTauNonLepH->Fill(ptTau,weight);
                  counterPtTauNonLepH->Fill(ptTau,weight);
                  counterPtNegMuTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
                  counterPtTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
                }
	      }
	    }
	  }
	}

      }

      if (isData) {

	TLorentzVector posMuJet4;
	TLorentzVector negMuJet4;
	
	float dRPosMuJetMin = 1000;
	float dRNegMuJetMin = 1000;

	for (int iJ=0; iJ<_nJets; ++iJ) {
	  
	  float dRPosMuJet = deltaR(_posLepEta,_posLepPhi,jetEta[iJ],jetPhi[iJ]);

	  if (dRPosMuJet<dRPosMuJetMin) {
	    dRPosMuJetMin = dRPosMuJet;
	    posMuJet4.SetXYZM(jetPx[iJ],jetPy[iJ],jetPz[iJ],jetEn[iJ]);
	  }

	  float dRNegMuJet = deltaR(_negLepEta,_negLepPhi,jetEta[iJ],jetPhi[iJ]);

          if (dRNegMuJet<dRNegMuJetMin) {
            dRNegMuJetMin = dRNegMuJet;
            negMuJet4.SetXYZM(jetPx[iJ],jetPy[iJ],jetPz[iJ],jetEn[iJ]); 
          }


	}
	  

	if (dRPosMuJetMin<0.5 && posMuJet4.Eta()>-2.4 && posMuJet4.Eta()<2.4) {

	  float ptTau = posMuJet4.Pt();

	  if (ptTau>99.9) ptTau = 99.9;
	  
	  counterPtPosMuTauH->Fill(ptTau,weight);
	  counterPtTauH->Fill(ptTau,weight);

	  float ptTrk = 1.0;
	  if (NpfNoPuIpPos>0) 
	    ptTrk = PtPfNoPuIpPos[0];

	  if (ptTrk>99.9) ptTrk = 99.9;
	    
	  counterPtPosMuTrkH->Fill(ptTrk,weight);
	  counterPtTrkH->Fill(ptTrk,weight);

	  counterPtPosMuTrkTauH->Fill(ptTrk,ptTau,weight);
          counterPtTrkTauH->Fill(ptTrk,ptTau,weight);


	  if (NpfNoPuIpPos>0) {
            counterPtPosMuTrkIdH->Fill(ptTrk,weight);
            counterPtTrkIdH->Fill(ptTrk,weight);
	    counterPtPosMuTauIdH->Fill(ptTau,weight);
	    counterPtTauIdH->Fill(ptTau,weight);
            counterPtPosMuTrkTauIdH->Fill(ptTrk,ptTau,weight);
            counterPtTrkTauIdH->Fill(ptTrk,ptTau,weight);
	    if (NpfNoPuIpPosAll==1 && NpfNoPuIpPos==1) {
	      counterPtPosMuTrkIsoH->Fill(ptTrk,weight);
              counterPtTrkIsoH->Fill(ptTrk,weight);
 	      counterPtPosMuTauIsoH->Fill(ptTau,weight);
	      counterPtTauIsoH->Fill(ptTau,weight);
	      counterPtPosMuTrkTauIsoH->Fill(ptTrk,ptTau,weight);
              counterPtTrkTauIsoH->Fill(ptTrk,ptTau,weight);
	      int idtrk = TMath::Abs(IdNoPuIpPos[0]);
	      if (idtrk!=11 && idtrk!=13) {
		counterPtPosMuTrkNonLepH->Fill(ptTrk,weight);
		counterPtTrkNonLepH->Fill(ptTrk,weight);
		counterPtPosMuTauNonLepH->Fill(ptTau,weight);
		counterPtTauNonLepH->Fill(ptTau,weight);
		counterPtPosMuTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
		counterPtTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
	      }
	    }
	  }
	}

	if (dRNegMuJetMin<0.5 && negMuJet4.Eta()>-2.4 && negMuJet4.Eta()<2.4 ) {

          float ptTau = negMuJet4.Pt();

          if (ptTau>99.9) ptTau = 99.9;

          counterPtNegMuTauH->Fill(ptTau,weight);
          counterPtTauH->Fill(ptTau,weight);

          float ptTrk = 1.0;
          if (NpfNoPuIpNeg>0)
            ptTrk = PtPfNoPuIpNeg[0];

          if (ptTrk>99.9) ptTrk = 99.9;

          counterPtNegMuTrkH->Fill(ptTrk,weight);
          counterPtTrkH->Fill(ptTrk,weight);

          counterPtNegMuTrkTauH->Fill(ptTrk,ptTau,weight);
          counterPtTrkTauH->Fill(ptTrk,ptTau,weight);

          if (NpfNoPuIpNeg>0) {
            counterPtNegMuTrkIdH->Fill(ptTrk,weight);
            counterPtTrkIdH->Fill(ptTrk,weight);
            counterPtNegMuTauIdH->Fill(ptTau,weight);
            counterPtTauIdH->Fill(ptTau,weight);
            counterPtNegMuTrkTauIdH->Fill(ptTrk,ptTau,weight);
            counterPtTrkTauIdH->Fill(ptTrk,ptTau,weight);
            if (NpfNoPuIpNegAll==1 && NpfNoPuIpNeg==1) {
              counterPtNegMuTrkIsoH->Fill(ptTrk,weight);
              counterPtTrkIsoH->Fill(ptTrk,weight);
              counterPtNegMuTauIsoH->Fill(ptTau,weight);
              counterPtTauIsoH->Fill(ptTau,weight);
              counterPtNegMuTrkTauIsoH->Fill(ptTrk,ptTau,weight);
              counterPtTrkTauIsoH->Fill(ptTrk,ptTau,weight);
	      int idtrk = TMath::Abs(IdNoPuIpPos[0]);
	      if (idtrk!=11 && idtrk!=13) {
		counterPtPosMuTrkNonLepH->Fill(ptTrk,weight);
		counterPtTrkNonLepH->Fill(ptTrk,weight);
		counterPtPosMuTauNonLepH->Fill(ptTau,weight);
		counterPtTauNonLepH->Fill(ptTau,weight);
		counterPtPosMuTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
		counterPtTrkTauNonLepH->Fill(ptTrk,ptTau,weight);
	      }
            }
          }
        }
      }
 

      dRMuPosTrkMinH->Fill(dRminPos,weight);
      dRMuNegTrkMinH->Fill(dRminNeg,weight);
      dRMuTrkMinH->Fill(dRminPos,weight);
      dRMuTrkMinH->Fill(dRminNeg,weight);
      
      AllNpfNoPUH->Fill(AllNpfNoPU,weight);
      AllNpfNoPU_ipH->Fill(AllNpfNoPU_ip,weight);

      NpfNoPuPosH->Fill(NpfNoPuPos,weight);
      NpfNoPuNegH->Fill(NpfNoPuNeg,weight);

      NpfNoPuIpPosH->Fill(NpfNoPuIpPos,weight);
      NpfNoPuIpNegH->Fill(NpfNoPuIpNeg,weight);

     //****************************************************************************************
      TLorentzVector MET4;
      MET4.SetPxPyPzE(metPx,metPy,metPz,metEn);

      TLorentzVector trackPos4;
      TLorentzVector trackNeg4;
      
      //      bool PosTrackCutsTMVA = false;
      //      bool NegTrackCutsTMVA = false;
      
      int IdTrackPosSelected = -1;
      int IdTrackNegSelected = -1;

      float relIsoPos = totPosIso/_posLepPt;
      float relIsoNeg = totNegIso/_negLepPt;
      

      trackIsoPosH->Fill(relIsoPos,weight);
      trackIsoNegH->Fill(relIsoNeg,weight);

      if (_posLepPt>_negLepPt) {
	trackIsoHardH->Fill(relIsoPos,weight);
	trackIsoSoftH->Fill(relIsoNeg,weight);
      }
      else {
	trackIsoHardH->Fill(relIsoNeg,weight);
	trackIsoSoftH->Fill(relIsoPos,weight);
      }

      counterBeforeIsoH->Fill(float(0.),weight);

      //      if (relIsoPos>0.4) continue;
      //      if (relIsoNeg>0.4) continue;

      counterAfterIsoH->Fill(float(0.),weight);
      

      // Id of pos track : IdNoPuIpPos
      if (IdNoPuIpPos[0]==11||IdNoPuIpPos[0]==-11)
	IdTrackPosSelected = 0;
      else if (IdNoPuIpPos[0]==13||IdNoPuIpPos[0]==-13)
	IdTrackPosSelected = 1;
      else
	IdTrackPosSelected = 2;
      // Id of neg track : IdNoPuIpNeg
      if (IdNoPuIpNeg[0]==11||IdNoPuIpNeg[0]==-11)
	IdTrackNegSelected = 0;
      else if (IdNoPuIpNeg[0]==13||IdNoPuIpNeg[0]==-13)
	IdTrackNegSelected = 1;
      else
	IdTrackNegSelected = 2;
      
      
      trackPos4.SetXYZM(PxPfNoPuIpPos[0],PyPfNoPuIpPos[0],PzPfNoPuIpPos[0],pionMass);
      trackNeg4.SetXYZM(PxPfNoPuIpNeg[0],PyPfNoPuIpNeg[0],PzPfNoPuIpNeg[0],pionMass);
      //Quality cuts on track
      bool PosTrackCuts = NpfNoPuIpPos==1;
      bool NegTrackCuts = NpfNoPuIpNeg==1;
      bool PosMuonTrackQne0 = (QPfNoPuIpPos[0]*_posLepQ>0. && trackPos4.Pt()>ptSoftTrkCut && TMath::Abs(IPPfNoPuIpPos[0])<d0CutTracks && TMath::Abs(IPZPfNoPuIpPos[0])<dzCutTracks && NpfNoPuIpPosAll==1);
      bool NegMuonTrackQne0 = (QPfNoPuIpNeg[0]*_negLepQ>0. && trackNeg4.Pt()>ptSoftTrkCut && TMath::Abs(IPPfNoPuIpNeg[0])<d0CutTracks && TMath::Abs(IPZPfNoPuIpNeg[0])<dzCutTracks && NpfNoPuIpNegAll==1);
      
      /*IPxyPos=IPPfNoPuIpPos[0];
      IPxyNeg=IPPfNoPuIpNeg[0];
      IPzPos=IPZPfNoPuIpPos[0];
      IPzNeg=IPZPfNoPuIpNeg[0];*/
      //

      int NtracksPos = NpfNoPuIpPosAll;
      int NtracksNeg = NpfNoPuIpNegAll;


      int NpfNoPuIpHardAll = NpfNoPuIpPosAll;
      int NpfNoPuIpSoftAll = NpfNoPuIpNegAll;

      if (_negLepPt > _posLepPt) {
	NpfNoPuIpHardAll = NpfNoPuIpNegAll;
	NpfNoPuIpSoftAll = NpfNoPuIpPosAll;
      }

      counterPreIsolationH->Fill(1.,weight);
      
      if (NpfNoPuIpHardAll==1 && NpfNoPuIpSoftAll==1)
	counterPostIsolationH->Fill(1.,weight);

      hardMuIsoTracksH->Fill(float(NpfNoPuIpHardAll),weight);
      softMuIsoTracksH->Fill(float(NpfNoPuIpSoftAll),weight);

      posMuIsoTracksH->Fill(float(NpfNoPuIpPosAll),weight);
      negMuIsoTracksH->Fill(float(NpfNoPuIpNegAll),weight);
      

      //****************************************************************************************

      double InvMassTrackPosMuon=trackPos4.M();
      double InvMassTrackNegMuon=trackNeg4.M();
	
      TLorentzVector TrackPlusMuonPos4;
      TrackPlusMuonPos4=MuonPos4+trackPos4;
      InvMassTrackPlusPosMuon=TrackPlusMuonPos4.M();

      TLorentzVector TrackPlusMuonNeg4;
      TrackPlusMuonNeg4=MuonNeg4+trackNeg4;
      InvMassTrackPlusNegMuon=TrackPlusMuonNeg4.M();
      
      TLorentzVector TwoTrackPlusMuons4;
      TwoTrackPlusMuons4=MuonPos4+trackPos4+MuonNeg4+trackNeg4;
      
      double InvMassTracksPlusMuons=TwoTrackPlusMuons4.M();

      TLorentzVector TracksMuonsMET;
      TracksMuonsMET=TwoTrackPlusMuons4+MET4;
      InvMassH2 = TracksMuonsMET.M();

      double SumOfInvMass2Systems=InvMassTrackPlusNegMuon+InvMassTrackPlusPosMuon; //also to be used for QCD background estimation

      //cuts on Pt of tracks
      float PtTrackPos = trackPos4.Pt();
      float PtTrackNeg = trackNeg4.Pt();
      float EtaTrackPos = trackPos4.Eta();
      float EtaTrackNeg = trackNeg4.Eta();

      float PtTrackHard=PtTrackPos;
      float PtTrackSoft=PtTrackNeg;
      float EtaTrackHard=EtaTrackPos;
      float EtaTrackSoft=EtaTrackNeg;
      
      if (PtTrackHard<PtTrackSoft){
	PtTrackHard=PtTrackNeg;
	PtTrackSoft=PtTrackPos;
	EtaTrackHard=EtaTrackNeg;
	EtaTrackSoft=EtaTrackPos;
      }

      bool PtTrackPosCut = PtTrackPos > ptHardTrkCut;
      bool PtTrackNegCut = PtTrackNeg > ptSoftTrkCut;

      if (event%2==0) {
	PtTrackPosCut = PtTrackPos > ptSoftTrkCut;
	PtTrackNegCut = PtTrackNeg > ptHardTrkCut;
      }

      // **************************
      // counting jets
      // and matching b-tagged 
      // jets to muons
      // **************************
      int nJets30 = 0;
      int nJets30Eta = 0;
      int nJets20CSVL  = 0;
      int nJets20CSVM  = 0;
      int nJets15CSVL  = 0;
      int nJets15CSVM  = 0;
      
      //bool posMuMatched20CSVM = false;
      bool posMuMatched20CSVL = false;
      //bool posMuMatched15CSVM = false;
      bool posMuMatched15CSVL = false;
      
      //bool negMuMatched20CSVM = false;
      bool negMuMatched20CSVL = false;
      //bool negMuMatched15CSVM = false;
      bool negMuMatched15CSVL = false;

      for (int iJ=0; iJ<_nJets; ++iJ) { 

        //      std::cout << "Jet  Pt = " << jetPt[iJ] << "   Eta = " << jetEta[iJ] << "   Phi = " << jetPhi[iJ] << std::endl;

        float absJetEta = TMath::Abs(jetEta[iJ]);
        bool jetId = jetIDLoose[iJ];
        bool jetPuId = jetPUIDLoose[iJ];

	if (chargedMultiplicity[iJ]<2) continue;
        if (!jetId) continue;
        if (!jetPuId) continue;
        if (absJetEta>jetEtaCut) continue;

	float dRPosLepJet = utils.deltaR(_posLepEta,_posLepPhi,jetEta[iJ],jetPhi[iJ]);
	float dRNegLepJet = utils.deltaR(_negLepEta,_negLepPhi,jetEta[iJ],jetPhi[iJ]);

        if (jetPt[iJ]>30. && dRPosLepJet>0.5 && dRNegLepJet>0.5) {
          nJets30++;
          if (absJetEta<bJetEtaCut)
            nJets30Eta++;
        }


        if (jetPt[iJ]>20. && absJetEta < bJetEtaCut) {
          if (combSVBJetTag[iJ]>btagCSVL) {
	    if (dRPosLepJet>0.5 && dRNegLepJet>0.5) 
	      nJets20CSVL++;
	    if (dRPosLepJet<0.4)
	      posMuMatched20CSVL = true;
	    if (dRNegLepJet<0.4)
	      negMuMatched20CSVL = true;
	  }
          if (combSVBJetTag[iJ]>btagCSVM) {
	    if (dRPosLepJet>0.5 && dRNegLepJet>0.5)
	      nJets20CSVM++;
	    /*if (dRPosLepJet<0.4)
	      posMuMatched20CSVM = true;
	    if (dRNegLepJet<0.4)
	      negMuMatched20CSVM = true;*/
	  }
        }

	if (jetPt[iJ]>15. && absJetEta < bJetEtaCut) {
          if (combSVBJetTag[iJ]>btagCSVL) {
	    if (dRPosLepJet>0.5 && dRNegLepJet>0.5)
	      nJets15CSVL++;
	    if (dRPosLepJet<0.4)
	      posMuMatched15CSVL = true;
	    if (dRNegLepJet<0.4)
	      negMuMatched15CSVL = true;
	  }
          if (combSVBJetTag[iJ]>btagCSVM) {
	    if (dRPosLepJet>0.5 && dRNegLepJet>0.5)
	      nJets15CSVM++;
	    /*if (dRPosLepJet<0.4)
	      posMuMatched15CSVM = true;
	    if (dRNegLepJet<0.4)
	      negMuMatched15CSVM = true;*/
	  }
        }


      }

      bool isPosHasSoftTrk = NpfNoPuIpPos==1 && ((NpfNoPuIpPosSoft + 1) == NpfNoPuIpPosAll) && NpfNoPuIpPosSoft>0 && NpfNoPuIpPosSoft<=NsoftTracks;
      bool isNegHasSoftTrk = NpfNoPuIpNeg==1 && ((NpfNoPuIpNegSoft + 1) == NpfNoPuIpNegAll) && NpfNoPuIpNegSoft>0 && NpfNoPuIpNegSoft<=NsoftTracks;

      bool bkgdControlNeg = (NpfNoPuIpNegAll==1 && NegTrackCuts) && isPosHasSoftTrk;
      bool bkgdControlPos = (NpfNoPuIpPosAll==1 && PosTrackCuts) && isNegHasSoftTrk; 

      bool invIpPos = TMath::Abs(IPPfNoPuIpPos[0])>d0CutLowerTracks && TMath::Abs(IPZPfNoPuIpPos[0])>dzCutLowerTracks;
      bool invIpNeg = TMath::Abs(IPPfNoPuIpNeg[0])>d0CutLowerTracks && TMath::Abs(IPZPfNoPuIpNeg[0])>dzCutLowerTracks;	

      
      // Signal tracks

      float ip2DTrackPos = IPPfNoPuIpPos[0];
      float ip2DTrackNeg = IPPfNoPuIpNeg[0];
      float ipZTrackNeg  = IPZPfNoPuIpPos[0];
      float ipZTrackPos  = IPZPfNoPuIpNeg[0];

      if (ip2DTrackPos>0.5) ip2DTrackPos = 0.499;
      if (ip2DTrackPos<-0.5) ip2DTrackPos = -0.499;
      if (ip2DTrackNeg>0.5) ip2DTrackNeg = 0.499;
      if (ip2DTrackNeg<-0.5) ip2DTrackNeg = -0.499;
      
      if (ipZTrackPos>1) ipZTrackPos = 0.999;
      if (ipZTrackPos<-1) ipZTrackPos = -0.999;
      if (ipZTrackNeg>1) ipZTrackNeg = 0.999;
      if (ipZTrackNeg<-1) ipZTrackNeg = -0.999;


      if (NpfNoPuIpPos==1) {
	Ip2DPosTrackH->Fill(ip2DTrackPos,weight);
	IpZPosTrackH->Fill(ipZTrackPos,weight);
	Ip2DBothTrackH->Fill(ip2DTrackPos,weight);
	IpZBothTrackH->Fill(ipZTrackPos,weight);
	PtPosTrackH->Fill(PtPfNoPuIpPos[0],weight);
	PtBothTrackH->Fill(PtPfNoPuIpPos[0],weight);
	if (TMath::Abs(ipZTrackPos)<dzCutTracks && TMath::Abs(ip2DTrackPos)<d0CutTracks) {
	  PtPosTrackIpCutsH->Fill(PtPfNoPuIpPos[0],weight);
	  PtBothTrackIpCutsH->Fill(PtPfNoPuIpPos[0],weight);
	}
	if (NtracksNeg<4) {
	  Ip2DPosTrackSideBandH->Fill(ip2DTrackPos,weight);
	  IpZPosTrackSideBandH->Fill(ipZTrackPos,weight);
	  Ip2DBothTrackSideBandH->Fill(ip2DTrackPos,weight);
	  IpZBothTrackSideBandH->Fill(ipZTrackPos,weight);
	  PtPosTrackSideBandH->Fill(PtPfNoPuIpPos[0],weight);
	  PtBothTrackSideBandH->Fill(PtPfNoPuIpPos[0],weight);
	  if (TMath::Abs(ipZTrackPos)<dzCutTracks && TMath::Abs(ip2DTrackPos)<d0CutTracks) {
	    PtPosTrackIpCutsSideBandH->Fill(PtPfNoPuIpPos[0],weight);
	    PtBothTrackIpCutsSideBandH->Fill(PtPfNoPuIpPos[0],weight);
	  }
	}
      }

      if (NpfNoPuIpNeg==1) {
	Ip2DNegTrackH->Fill(ip2DTrackNeg,weight);
	IpZNegTrackH->Fill(ipZTrackNeg,weight);
	Ip2DBothTrackH->Fill(ip2DTrackNeg,weight);
	IpZBothTrackH->Fill(ipZTrackNeg,weight);
	PtNegTrackH->Fill(PtPfNoPuIpNeg[0],weight);
	PtBothTrackH->Fill(PtPfNoPuIpNeg[0],weight);
	if (TMath::Abs(ipZTrackNeg)<dzCutTracks && TMath::Abs(ip2DTrackNeg)<d0CutTracks) {
	  PtNegTrackIpCutsH->Fill(PtPfNoPuIpNeg[0],weight);
	  PtBothTrackIpCutsH->Fill(PtPfNoPuIpNeg[0],weight);
	}
	if (NtracksPos<4) {
	  Ip2DNegTrackSideBandH->Fill(ip2DTrackNeg,weight);
	  IpZNegTrackSideBandH->Fill(ipZTrackNeg,weight);
	  Ip2DBothTrackSideBandH->Fill(ip2DTrackNeg,weight);
	  IpZBothTrackSideBandH->Fill(ipZTrackNeg,weight);
	  PtNegTrackSideBandH->Fill(PtPfNoPuIpNeg[0],weight);
	  PtBothTrackSideBandH->Fill(PtPfNoPuIpNeg[0],weight);
	  if (TMath::Abs(ipZTrackNeg)<dzCutTracks && TMath::Abs(ip2DTrackNeg)<d0CutTracks) {
	    PtNegTrackIpCutsSideBandH->Fill(PtPfNoPuIpNeg[0],weight);
	    PtBothTrackIpCutsSideBandH->Fill(PtPfNoPuIpNeg[0],weight);
	  }
	}

      }


     
      bool AdditionalInverse = invIpPos || invIpNeg ;
      if (ipAndLogic)
	AdditionalInverse = invIpPos && invIpNeg ;

      //background evaluation studies

      //Check if Invariant mass of one system is depended on the other

      // both muons are non-isolated
      // first muon -->
      if (InvMassTrackPlusPosMuon<=1.)
	InvMassTrackPlusNeg_Pos1GeVH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>1. && InvMassTrackPlusPosMuon<=2.)
	InvMassTrackPlusNeg_Pos2GeVH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>2. && InvMassTrackPlusPosMuon<=3.)
	InvMassTrackPlusNeg_Pos3GeVH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>3. && InvMassTrackPlusPosMuon<=4. )
	InvMassTrackPlusNeg_Pos4GeVH->Fill(InvMassTrackPlusNegMuon,weight);

      if (InvMassTrackPlusPosMuon>4)
	InvMassTrackPlusNeg_PosGT4GeVH->Fill(InvMassTrackPlusNegMuon,weight);
      // second muon -->
      if (InvMassTrackPlusNegMuon<=1.)
	InvMassTrackPlusPos_Neg1GeVH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>1. && InvMassTrackPlusNegMuon<=2.)
	InvMassTrackPlusPos_Neg2GeVH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>2. && InvMassTrackPlusNegMuon<=3.)
	InvMassTrackPlusPos_Neg3GeVH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>3. && InvMassTrackPlusNegMuon<=4.)
	InvMassTrackPlusPos_Neg4GeVH->Fill(InvMassTrackPlusPosMuon,weight);

      if (InvMassTrackPlusNegMuon>4)
	InvMassTrackPlusPos_NegGT4GeVH->Fill(InvMassTrackPlusPosMuon,weight);


      // one of the muons have additional soft track around it
      // first muon --->
      if (bkgdControlNeg && AdditionalInverse) {

	// id's
	idTrack_ControlNeg_2DH->Fill(float(IdTrackPosSelected),float(IdTrackNegSelected),weight);
	idTrackPos_ControlNegH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_ControlNegH->Fill(float(IdTrackNegSelected),weight);

	// general control plots
	InvMassTrackPlusMuon2D_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	InvMassTrackPlusMuon2D_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);	
	met_ControlH->Fill(met,weight);
	InvMassTracksPlusMuons_ControlH->Fill(InvMassTracksPlusMuons,weight);
	HardMuonPt_ControlH->Fill(ptMax,weight);
	SoftMuonPt_ControlH->Fill(ptMin,weight);
	SumOfInvMass2Systems_ControlH->Fill(SumOfInvMass2Systems,weight);
	InvMassH2_ControlH->Fill(InvMassH2,weight);
	PtTrackHard_ControlH->Fill(PtTrackHard,weight);
	PtTrackSoft_ControlH->Fill(PtTrackSoft,weight);
	EtaTrackHard_ControlH->Fill(EtaTrackHard,weight);
	EtaTrackSoft_ControlH->Fill(EtaTrackSoft,weight);
	nJets30_ControlH->Fill(nJets30,weight);
	nJets30Eta_ControlH->Fill(nJets30Eta,weight);
	nJets20CSVL_ControlH->Fill(nJets20CSVL,weight);
	nJets20CSVM_ControlH->Fill(nJets20CSVM,weight);
	nJets15CSVL_ControlH->Fill(nJets15CSVL,weight);
	nJets15CSVM_ControlH->Fill(nJets15CSVM,weight);

	InvMassTrackPlusPosMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNegMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);

	dRMuMu_ControlH->Fill(muonDR,weight);
	dEtaMuMu_ControlH->Fill(muonDEta,weight);
	dPhiMuMu_ControlH->Fill(muonDPhi,weight);



	if (InvMassTrackPlusPosMuon<=1.)
	  InvMassTrackPlusNeg_Pos1GeV_ControlNegH->Fill(InvMassTrackPlusNegMuon,weight);
      
	if (InvMassTrackPlusPosMuon>1. && InvMassTrackPlusPosMuon<=2.)
	  InvMassTrackPlusNeg_Pos2GeV_ControlNegH->Fill(InvMassTrackPlusNegMuon,weight);
      
	if (InvMassTrackPlusPosMuon>2. && InvMassTrackPlusPosMuon<=3.)
	  InvMassTrackPlusNeg_Pos3GeV_ControlNegH->Fill(InvMassTrackPlusNegMuon,weight);
      
	if (InvMassTrackPlusPosMuon>3. && InvMassTrackPlusPosMuon<=4.)
	  InvMassTrackPlusNeg_Pos4GeV_ControlNegH->Fill(InvMassTrackPlusNegMuon,weight);

        if (InvMassTrackPlusPosMuon>4.)
          InvMassTrackPlusNeg_PosGT4GeV_ControlNegH->Fill(InvMassTrackPlusNegMuon,weight);


	if (posMuMatched20CSVL && negMuMatched20CSVL) 
	  InvMassTrackPlusMuon2D_bothBJets20_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_posBJet20_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (!posMuMatched20CSVL && negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_negBJet20_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (!posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_NoBJets20_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	// Tracks id's
	if (IdTrackPosSelected==1&&IdTrackNegSelected==1) // MuMu
	  InvMassTrackPlusMuon2D_MuMu_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==0) // MuEle
	  InvMassTrackPlusMuon2D_MuEle_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==1) // MuEle (swapped)
	  InvMassTrackPlusMuon2D_MuEle_ControlNegH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==2) // MuHad
	  InvMassTrackPlusMuon2D_MuHad_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==1) // MuHad (swapped)
	  InvMassTrackPlusMuon2D_MuHad_ControlNegH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
	
	if (IdTrackPosSelected==0&&IdTrackNegSelected==0) // EleEle
	  InvMassTrackPlusMuon2D_EleEle_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==2) // EleHad
	  InvMassTrackPlusMuon2D_EleHad_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (IdTrackPosSelected==2&&IdTrackNegSelected==0) // EleHad (swapped)
	  InvMassTrackPlusMuon2D_EleHad_ControlNegH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==2) // HadHad 
	  InvMassTrackPlusMuon2D_HadHad_ControlNegH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	

      }
      // second muon --->
      if (bkgdControlPos && AdditionalInverse) {

	// track id's
	idTrack_ControlPos_2DH->Fill(float(IdTrackPosSelected),float(IdTrackNegSelected),weight);
	idTrackPos_ControlPosH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_ControlPosH->Fill(float(IdTrackNegSelected),weight);

	// general control plots
	InvMassTrackPlusMuon2D_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	InvMassTrackPlusMuon2D_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);	
	met_ControlH->Fill(met,weight);
	InvMassTracksPlusMuons_ControlH->Fill(InvMassTracksPlusMuons,weight);
	HardMuonPt_ControlH->Fill(ptMax,weight);
	SoftMuonPt_ControlH->Fill(ptMin,weight);
	SumOfInvMass2Systems_ControlH->Fill(SumOfInvMass2Systems,weight);
	InvMassH2_ControlH->Fill(InvMassH2,weight);
	PtTrackHard_ControlH->Fill(PtTrackHard,weight);
	PtTrackSoft_ControlH->Fill(PtTrackSoft,weight);
	EtaTrackHard_ControlH->Fill(EtaTrackHard,weight);
	EtaTrackSoft_ControlH->Fill(EtaTrackSoft,weight);
	nJets30_ControlH->Fill(nJets30,weight);
	nJets30Eta_ControlH->Fill(nJets30Eta,weight);
	nJets20CSVL_ControlH->Fill(nJets20CSVL,weight);
	nJets20CSVM_ControlH->Fill(nJets20CSVM,weight);
	nJets15CSVL_ControlH->Fill(nJets15CSVL,weight);
	nJets15CSVM_ControlH->Fill(nJets15CSVM,weight);
	
	InvMassTrackPlusPosMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNegMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);

        dRMuMu_ControlH->Fill(muonDR,weight);
        dEtaMuMu_ControlH->Fill(muonDEta,weight);
        dPhiMuMu_ControlH->Fill(muonDPhi,weight);

	if (InvMassTrackPlusNegMuon<=1.)
	  InvMassTrackPlusPos_Neg1GeV_ControlPosH->Fill(InvMassTrackPlusPosMuon,weight);
      
	if (InvMassTrackPlusNegMuon>1. && InvMassTrackPlusNegMuon<=2.)
	  InvMassTrackPlusPos_Neg2GeV_ControlPosH->Fill(InvMassTrackPlusPosMuon,weight);
      
	if (InvMassTrackPlusNegMuon>2. && InvMassTrackPlusNegMuon<=3.)
	  InvMassTrackPlusPos_Neg3GeV_ControlPosH->Fill(InvMassTrackPlusPosMuon,weight);
      
	if (InvMassTrackPlusNegMuon>3. && InvMassTrackPlusNegMuon<=4.)
	  InvMassTrackPlusPos_Neg4GeV_ControlPosH->Fill(InvMassTrackPlusPosMuon,weight);

        if (InvMassTrackPlusNegMuon>4.)
          InvMassTrackPlusPos_NegGT4GeV_ControlPosH->Fill(InvMassTrackPlusPosMuon,weight);


	if (posMuMatched20CSVL && negMuMatched20CSVL) 
	  InvMassTrackPlusMuon2D_bothBJets20_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_posBJet20_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (!posMuMatched20CSVL && negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_negBJet20_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (!posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_NoBJets20_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	// Tracks id's
	if (IdTrackPosSelected==1&&IdTrackNegSelected==1) // MuMu
	  InvMassTrackPlusMuon2D_MuMu_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==0) // MuEle
	  InvMassTrackPlusMuon2D_MuEle_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==1) // MuEle (swapped)
	  InvMassTrackPlusMuon2D_MuEle_ControlPosH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==2) // MuHad
	  InvMassTrackPlusMuon2D_MuHad_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==1) // MuHad (swapped)
	  InvMassTrackPlusMuon2D_MuHad_ControlPosH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
	
	if (IdTrackPosSelected==0&&IdTrackNegSelected==0) // EleEle
	  InvMassTrackPlusMuon2D_EleEle_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==2) // EleHad
	  InvMassTrackPlusMuon2D_EleHad_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (IdTrackPosSelected==2&&IdTrackNegSelected==0) // EleHad (swapped)
	  InvMassTrackPlusMuon2D_EleHad_ControlPosH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==2) // HadHad 
	  InvMassTrackPlusMuon2D_HadHad_ControlPosH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      }
      // each muon has additional soft track around it
      if ( isPosHasSoftTrk && isNegHasSoftTrk && AdditionalInverse) {

	idTrack_Control_2DH->Fill(float(IdTrackPosSelected),float(IdTrackNegSelected),weight);
	idTrackPos_ControlH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_ControlH->Fill(float(IdTrackNegSelected),weight);

	InvMassTrackPlusMuon2D_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);	
	InvMassTrackPlusMuon2D_ControlBothH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);	
	met_ControlH->Fill(met,weight);
	InvMassTracksPlusMuons_ControlH->Fill(InvMassTracksPlusMuons,weight);
	HardMuonPt_ControlH->Fill(ptMax,weight);
	SoftMuonPt_ControlH->Fill(ptMin,weight);
	SumOfInvMass2Systems_ControlH->Fill(SumOfInvMass2Systems,weight);
	InvMassH2_ControlH->Fill(InvMassH2,weight);
	PtTrackHard_ControlH->Fill(PtTrackHard,weight);
	PtTrackSoft_ControlH->Fill(PtTrackSoft,weight);
	EtaTrackHard_ControlH->Fill(EtaTrackHard,weight);
	EtaTrackSoft_ControlH->Fill(EtaTrackSoft,weight);
	nJets30_ControlH->Fill(nJets30,weight);
	nJets30Eta_ControlH->Fill(nJets30Eta,weight);
	nJets20CSVL_ControlH->Fill(nJets20CSVL,weight);
	nJets20CSVM_ControlH->Fill(nJets20CSVM,weight);
	nJets15CSVL_ControlH->Fill(nJets15CSVL,weight);
	nJets15CSVM_ControlH->Fill(nJets15CSVM,weight);

	InvMassTrackPlusPosMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNegMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusMuon1D_ControlH->Fill(InvMassTrackPlusNegMuon,weight);

        dRMuMu_ControlH->Fill(muonDR,weight);
        dEtaMuMu_ControlH->Fill(muonDEta,weight);
        dPhiMuMu_ControlH->Fill(muonDPhi,weight);

	if (posMuMatched20CSVL && negMuMatched20CSVL) 
	  InvMassTrackPlusMuon2D_bothBJets20_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_posBJet20_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (!posMuMatched20CSVL && negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_negBJet20_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (!posMuMatched20CSVL && !negMuMatched20CSVL)
	  InvMassTrackPlusMuon2D_NoBJets20_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);


	// Tracks id's
	if (IdTrackPosSelected==1&&IdTrackNegSelected==1) // MuMu
	  InvMassTrackPlusMuon2D_MuMu_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==0) // MuEle
	  InvMassTrackPlusMuon2D_MuEle_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==1) // MuEle (swapped)
	  InvMassTrackPlusMuon2D_MuEle_ControlH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==1&&IdTrackNegSelected==2) // MuHad
	  InvMassTrackPlusMuon2D_MuHad_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==1) // MuHad (swapped)
	  InvMassTrackPlusMuon2D_MuHad_ControlH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
	
	if (IdTrackPosSelected==0&&IdTrackNegSelected==0) // EleEle
	  InvMassTrackPlusMuon2D_EleEle_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

	if (IdTrackPosSelected==0&&IdTrackNegSelected==2) // EleHad
	  InvMassTrackPlusMuon2D_EleHad_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	
	if (IdTrackPosSelected==2&&IdTrackNegSelected==0) // EleHad (swapped)
	  InvMassTrackPlusMuon2D_EleHad_ControlH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);

	if (IdTrackPosSelected==2&&IdTrackNegSelected==2) // HadHad 
	  InvMassTrackPlusMuon2D_HadHad_ControlH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	

      }


      if (isPosHasSoftTrk && invIpPos) {
	if (NtracksNeg==2) {
	  InvMassTrackPlusPos_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksNeg==3) {
	  InvMassTrackPlusPos_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksNeg==4) {
	  InvMassTrackPlusPos_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksNeg>1 && NtracksNeg<4) {
	  InvMassTrackPlusPos_InvNegQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_InvNegQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_InvNegQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_InvNegQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	InvMassTrackPlusPos_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
      }


      if (isNegHasSoftTrk && invIpNeg) {
	if (NtracksPos==2) {
	  InvMassTrackPlusNeg_2Q_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_2Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksPos==3) {
	  InvMassTrackPlusNeg_3Q_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_3Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksPos==4) {
	  InvMassTrackPlusNeg_4Q_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_4Q_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksPos>1 && NtracksPos<4) {
	  InvMassTrackPlusNeg_InvPosQ_ControlH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_InvPosQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_InvPosQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_InvPosQ_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	InvMassTrackPlusNeg_ControlH->Fill(InvMassTrackPlusPosMuon,weight);
      }


      //Inverse and direct isolation. CUTS ON TRACK NUMBER ONLY
      if ( NtracksPos>1 && NtracksPos<4 && NtracksNeg==1)counterPosCutsH->Fill(1.,weight);
      if ( NtracksNeg>1 && NtracksNeg<4 && NtracksPos==1)counterNegCutsH->Fill(1.,weight);
      if ( NtracksPos==2 && NtracksNeg==2) counterInverse2H->Fill(1.,weight);
      if ( NtracksPos==3 && NtracksNeg==3) counterInverse3H->Fill(1.,weight);
      if ( NtracksPos==4 && NtracksNeg==4) counterInverse4H->Fill(1.,weight);
      if ( NtracksPos>1 && NtracksNeg>1 && NtracksPos<4 && NtracksNeg<4)counterInverseH->Fill(1.,weight);
      if ( PosMuonTrackQne0 && NegMuonTrackQne0) counterSwapQH->Fill(1.,weight); 
      if ( NtracksPos==1 && NtracksNeg==1) counterDirectH->Fill(1.,weight);

     

      //Invariant mass when one muon is allowed to have more than one tracks around it

      float drpos = deltaR(_posLepEta,_posLepPhi,trackPos4.Eta(),trackPos4.Phi());
      float drneg = deltaR(_negLepEta,_negLepPhi,trackNeg4.Eta(),trackNeg4.Phi());

      if (NtracksPos==1 && PosTrackCuts){
	if (NtracksNeg==2){
	  InvMassTrackPlusPos_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	  PosmuonPt_2QH->Fill(_posLepPt,weight);
	  PosmuonEta_2QH->Fill(_posLepEta,weight);
	  PostrackPt_2QH->Fill(trackPos4.Pt(),weight);
	  PostrackEta_2QH->Fill(trackPos4.Eta(),weight);
	  dRPosmuonTrk_2QH->Fill(drpos,weight);
	  counterPos_2QH->Fill(1.,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksNeg==3){
	  InvMassTrackPlusPos_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	  PosmuonPt_3QH->Fill(_posLepPt,weight);
          PosmuonEta_3QH->Fill(_posLepEta,weight);
	  PostrackPt_3QH->Fill(trackPos4.Pt(),weight);
          PostrackEta_3QH->Fill(trackPos4.Eta(),weight);
          dRPosmuonTrk_3QH->Fill(drpos,weight);
	  counterPos_3QH->Fill(1.,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	}
	if (NtracksNeg==4){
	  InvMassTrackPlusPos_4QH->Fill(InvMassTrackPlusPosMuon,weight);
	  PosmuonPt_4QH->Fill(_posLepPt,weight);
          PosmuonEta_4QH->Fill(_posLepEta,weight);
	  PostrackPt_4QH->Fill(trackPos4.Pt(),weight);
          PostrackEta_4QH->Fill(trackPos4.Eta(),weight);
          dRPosmuonTrk_4QH->Fill(drpos,weight);
	  counterPos_4QH->Fill(1.,weight);
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_4QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_4QH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_4QH->Fill(InvMassTrackPlusPosMuon,weight);	    
	}
	if (NtracksNeg>1 && NtracksNeg<4){
	  InvMassTrackPlusPos_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	  PosmuonPt_InvQH->Fill(_posLepPt,weight);
	  PosmuonEta_InvQH->Fill(_posLepEta,weight);
	  PostrackPt_InvQH->Fill(trackPos4.Pt(),weight);
	  PostrackEta_InvQH->Fill(trackPos4.Eta(),weight);
	  dRPosmuonTrk_InvQH->Fill(drpos,weight);
	  counterPos_InvQH->Fill(1.,weight);
	  if (posHard) {
	    HardmuonPt_InvQH->Fill(_posLepPt,weight);
	    HardmuonEta_InvQH->Fill(_posLepEta,weight);
	    HardtrackMuPt_InvQH->Fill(trackPos4.Pt(),weight);
	    HardtrackMuEta_InvQH->Fill(trackPos4.Eta(),weight);
	    dRHardmuonTrk_InvQH->Fill(drpos,weight);
	  }
	  else {
	    SoftmuonPt_InvQH->Fill(_posLepPt,weight);
	    SoftmuonEta_InvQH->Fill(_posLepEta,weight);
	    SofttrackMuPt_InvQH->Fill(trackPos4.Pt(),weight);
	    SofttrackMuEta_InvQH->Fill(trackPos4.Eta(),weight);
	    dRSoftmuonTrk_InvQH->Fill(drpos,weight);
	  }
	  if (IdTrackPosSelected==1)
	    InvMassTrackPlusPos_Mu_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==0)
	    InvMassTrackPlusPos_Ele_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	  if (IdTrackPosSelected==2)
	    InvMassTrackPlusPos_Had_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);	  
	}
	if (NegMuonTrackQne0) {
	  counterPos_SwapQH->Fill(1.,weight);
	  InvMassTrackPlusPos_SwapQH->Fill(InvMassTrackPlusPosMuon,weight);
	}
      }
     
      if(NtracksNeg==1 && NegTrackCuts){
	if (NtracksPos==2){
	  InvMassTrackPlusNeg_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	  NegmuonPt_2QH->Fill(_negLepPt,weight);
	  NegmuonEta_2QH->Fill(_negLepEta,weight);
	  NegtrackPt_2QH->Fill(trackNeg4.Pt(),weight);
	  NegtrackEta_2QH->Fill(trackNeg4.Eta(),weight);
	  dRNegmuonTrk_2QH->Fill(drneg,weight);
	  counterNeg_2QH->Fill(1.,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	}
	if (NtracksPos==3){
	  InvMassTrackPlusNeg_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	  NegmuonPt_3QH->Fill(_negLepPt,weight);
          NegmuonEta_3QH->Fill(_negLepEta,weight);
	  NegtrackPt_3QH->Fill(trackNeg4.Pt(),weight);
	  NegtrackEta_3QH->Fill(trackNeg4.Eta(),weight);
	  dRNegmuonTrk_3QH->Fill(drneg,weight);
	  counterNeg_3QH->Fill(1.,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	}
	if (NtracksPos==4){
	  InvMassTrackPlusNeg_4QH->Fill(InvMassTrackPlusNegMuon,weight);
	  NegmuonPt_4QH->Fill(_negLepPt,weight);
          NegmuonEta_4QH->Fill(_negLepEta,weight);
	  NegtrackPt_4QH->Fill(trackNeg4.Pt(),weight);
          NegtrackEta_4QH->Fill(trackNeg4.Eta(),weight);
	  dRNegmuonTrk_4QH->Fill(drneg,weight);
	  counterNeg_4QH->Fill(1.,weight);
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_4QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_4QH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_4QH->Fill(InvMassTrackPlusNegMuon,weight);
	}
	if (NtracksPos>1 && NtracksPos<4){
	  InvMassTrackPlusNeg_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	  NegmuonPt_InvQH->Fill(_negLepPt,weight);
          NegmuonEta_InvQH->Fill(_negLepEta,weight);
	  NegtrackPt_InvQH->Fill(trackNeg4.Pt(),weight);
          NegtrackEta_InvQH->Fill(trackNeg4.Eta(),weight);
          dRNegmuonTrk_InvQH->Fill(drneg,weight);
	  counterNeg_InvQH->Fill(1.,weight);
	  if (posHard) {
	    SoftmuonPt_InvQH->Fill(_negLepPt,weight);
	    SoftmuonEta_InvQH->Fill(_negLepEta,weight);
	    SofttrackMuPt_InvQH->Fill(trackNeg4.Pt(),weight);
	    SofttrackMuEta_InvQH->Fill(trackNeg4.Eta(),weight);
	    dRSoftmuonTrk_InvQH->Fill(drneg,weight);
	  }
	  else {
	    HardmuonPt_InvQH->Fill(_negLepPt,weight);
	    HardmuonEta_InvQH->Fill(_negLepEta,weight);
	    HardtrackMuPt_InvQH->Fill(trackNeg4.Pt(),weight);
	    HardtrackMuEta_InvQH->Fill(trackNeg4.Eta(),weight);
	    dRHardmuonTrk_InvQH->Fill(drneg,weight);
	  }
	  if (IdTrackNegSelected==1)
	    InvMassTrackPlusNeg_Mu_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==0)
	    InvMassTrackPlusNeg_Ele_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	  if (IdTrackNegSelected==2)
	    InvMassTrackPlusNeg_Had_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	}
	if (PosMuonTrackQne0) {
	  counterNeg_SwapQH->Fill(1.,weight);
	  InvMassTrackPlusNeg_SwapQH->Fill(InvMassTrackPlusNegMuon,weight);
	}
      }

      // *****BTag Categories********
      if(posMuMatched20CSVL){
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2)InvMassTrackPlusPos_BJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	    
	    if (NtracksNeg==3)InvMassTrackPlusPos_BJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	     
	    if (NtracksNeg>1 && NtracksNeg<4)InvMassTrackPlusPos_BJet20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	  }
	}
      }

      if(!posMuMatched20CSVL){
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2)InvMassTrackPlusPos_NotBJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	    
	    if (NtracksNeg==3)InvMassTrackPlusPos_NotBJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	    
	    if (NtracksNeg>1 && NtracksNeg<4)InvMassTrackPlusPos_NotBJet20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	  }
	}
      }

      if(negMuMatched20CSVL){
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2)InvMassTrackPlusNeg_BJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    
	    if (NtracksPos==3)InvMassTrackPlusNeg_BJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	     
	    if (NtracksPos>1 && NtracksPos<4)InvMassTrackPlusNeg_BJet20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	  }
	}
      }
      if(!negMuMatched20CSVL){
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2)InvMassTrackPlusNeg_NotBJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    
	    if (NtracksPos==3)InvMassTrackPlusNeg_NotBJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	     
	    if (NtracksPos>1 && NtracksPos<4)InvMassTrackPlusNeg_NotBJet20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	  }
	}
      }
      


//********************Background "side bands" using b jets*********************   
      if(posMuMatched20CSVL && negMuMatched20CSVL){
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2){
	      InvMassTrackPlusPos_PosIso_bothBJets20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg==3){
	      InvMassTrackPlusPos_PosIso_bothBJets20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg>1 && NtracksNeg<4){
	      InvMassTrackPlusPos_PosIso_bothBJets20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_bothBJets20CSVL_InvNegQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(posMuMatched20CSVL && !negMuMatched20CSVL){  
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2){
	      InvMassTrackPlusPos_PosIso_posBJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_posBJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg==3){
	      InvMassTrackPlusPos_PosIso_posBJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_posBJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg>1 && NtracksNeg<4){
	      InvMassTrackPlusPos_PosIso_posBJet20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_posBJet20CSVL_InvNegQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(!posMuMatched20CSVL && negMuMatched20CSVL){
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2){
	      InvMassTrackPlusPos_PosIso_negBJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_negBJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg==3){
	      InvMassTrackPlusPos_PosIso_negBJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_negBJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg>1 && NtracksNeg<4){
	      InvMassTrackPlusPos_PosIso_negBJet20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_negBJet20CSVL_InvNegQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(!posMuMatched20CSVL && !negMuMatched20CSVL){
	if (NtracksPos==1 && PosTrackCuts){
	  if (PtTrackPosCut){ 
	    if (NtracksNeg==2){
	      InvMassTrackPlusPos_PosIso_NoBJets20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg==3){
	      InvMassTrackPlusPos_PosIso_NoBJets20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksNeg>1 && NtracksNeg<4){
	      InvMassTrackPlusPos_PosIso_NoBJets20CSVL_InvNegQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_PosIso_NoBJets20CSVL_InvNegQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      ////
      if(posMuMatched20CSVL && negMuMatched20CSVL){
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2){
	      InvMassTrackPlusPos_NegIso_bothBJets20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos==3){
	      InvMassTrackPlusPos_NegIso_bothBJets20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos>1 && NtracksPos<4){
	      InvMassTrackPlusPos_NegIso_bothBJets20CSVL_InvPosQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_bothBJets20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(posMuMatched20CSVL && !negMuMatched20CSVL){  
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2){
	      InvMassTrackPlusPos_NegIso_posBJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_posBJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos==3){
	      InvMassTrackPlusPos_NegIso_posBJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_posBJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos>1 && NtracksPos<4){
	      InvMassTrackPlusPos_NegIso_posBJet20CSVL_InvPosQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_posBJet20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(!posMuMatched20CSVL && negMuMatched20CSVL){
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2){
	      InvMassTrackPlusPos_NegIso_negBJet20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_negBJet20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos==3){
	      InvMassTrackPlusPos_NegIso_negBJet20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_negBJet20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos>1 && NtracksPos<4){
	      InvMassTrackPlusPos_NegIso_negBJet20CSVL_InvPosQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_negBJet20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      }
      if(!posMuMatched20CSVL && !negMuMatched20CSVL){
	if (NtracksNeg==1 && NegTrackCuts){
	  if (PtTrackNegCut){ 
	    if (NtracksPos==2){
	      InvMassTrackPlusPos_NegIso_NoBJets20CSVL_2QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_2QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos==3){
	      InvMassTrackPlusPos_NegIso_NoBJets20CSVL_3QH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_3QH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	    if (NtracksPos>1 && NtracksPos<4){
	      InvMassTrackPlusPos_NegIso_NoBJets20CSVL_InvPosQH->Fill(InvMassTrackPlusPosMuon,weight);
	      InvMassTrackPlusNeg_NegIso_NoBJets20CSVL_InvPosQH->Fill(InvMassTrackPlusNegMuon,weight);
	    }
	  }
	}
      } 
      

      counterPreIsoFinalH->Fill(1.,weight);

      //********************  cuts on track number ******************************
      if (NtracksNeg!=1)continue;
      if (NtracksPos!=1)continue;

      counterPostIsoFinalH->Fill(1.,weight);

      //******************** cuts on track quality ****************************** 
      if (!PosTrackCuts) continue;
      if (!NegTrackCuts) continue;

      counterSemiFinalH->Fill(1.,weight);

      PtTrackHardSemiFinalH->Fill(PtTrackHard,weight);
      PtTrackSoftSemiFinalH->Fill(PtTrackSoft,weight);

      if (PtTrackHard<ptHardTrkCut)continue;
      if (PtTrackSoft<ptSoftTrkCut)continue;

      counterFinalH->Fill(1.,weight);
     
      HiggsPtFinalH->Fill(HiggsPt,weight);

      PtTrackHardH->Fill(PtTrackHard,weight);
      PtTrackSoftH->Fill(PtTrackSoft,weight);
      EtaTrackHardH->Fill(EtaTrackHard,weight);
      EtaTrackSoftH->Fill(EtaTrackSoft,weight);
      
      HardMuonPtH->Fill(ptMax,weight);
      SoftMuonPtH->Fill(ptMin,weight);
      HardMuonEtaH->Fill(etaMax,weight);
      SoftMuonEtaH->Fill(etaMin,weight);


      if (posHard) {
	PtTrackMuHardH->Fill(PtTrackPos,weight);
	PtTrackMuSoftH->Fill(PtTrackNeg,weight);
	EtaTrackMuHardH->Fill(EtaTrackPos,weight);
	EtaTrackMuSoftH->Fill(EtaTrackNeg,weight);
      }
      else {
	PtTrackMuHardH->Fill(PtTrackNeg,weight);
        PtTrackMuSoftH->Fill(PtTrackPos,weight);
        EtaTrackMuHardH->Fill(EtaTrackNeg,weight);
        EtaTrackMuSoftH->Fill(EtaTrackPos,weight);
      }

      PosMuonPtH->Fill(_posLepPt,weight);
      NegMuonPtH->Fill(_negLepPt,weight);

      PosMuonEtaH->Fill(_posLepEta,weight);
      NegMuonEtaH->Fill(_negLepEta,weight);

      PosTrackPtH->Fill(PtTrackPos,weight);
      NegTrackPtH->Fill(PtTrackNeg,weight);

      PosTrackEtaH->Fill(EtaTrackPos,weight);
      NegTrackEtaH->Fill(EtaTrackNeg,weight);

      dRPosMuonTrackH->Fill(drpos,weight); 
      dRNegMuonTrackH->Fill(drneg,weight);

      if (posHard) {
	dRHardMuonTrackH->Fill(drpos,weight); 
	dRSoftMuonTrackH->Fill(drneg,weight);
      }
      else {
	dRHardMuonTrackH->Fill(drneg,weight); 
	dRSoftMuonTrackH->Fill(drpos,weight);
      }

      nJets30H->Fill(float(nJets30),weight);
      nJets30EtaH->Fill(float(nJets30Eta),weight);
      nJets20CSVLH->Fill(float(nJets20CSVL),weight);
      nJets20CSVMH->Fill(float(nJets20CSVM),weight);
      nJets15CSVLH->Fill(float(nJets15CSVL),weight);
      nJets15CSVMH->Fill(float(nJets15CSVM),weight);

      idTrackPosH->Fill(float(IdTrackPosSelected),weight);
      idTrackNegH->Fill(float(IdTrackNegSelected),weight);

      idTrack2DH->Fill(float(IdTrackPosSelected),float(IdTrackNegSelected),weight);

      //Filling histograms
      Ntracks_OutConePosH->Fill(Ntracks_OutConePos,weight);
      Ntracks_OutConeNegH->Fill(Ntracks_OutConeNeg,weight);
      Ntracks_OutConePosIPH->Fill(Ntracks_OutConePosIP,weight);
      Ntracks_OutConeNegIPH->Fill(Ntracks_OutConeNegIP,weight);
      SumPtPosConeH->Fill(SumPtPosCone,weight);
      SumPtNegConeH->Fill(SumPtNegCone,weight);
      SumPtPosConeIPH->Fill(SumPtPosConeIP,weight);
      SumPtNegConeIPH->Fill(SumPtNegConeIP,weight);

     
      SumOfInvMass2SystemsH->Fill(SumOfInvMass2Systems,weight);

      metH->Fill(met,weight);
    
      nPVH->Fill(nPV);
      InvMassTrackPlusPosH->Fill(InvMassTrackPlusPosMuon,weight);
      InvMassTrackPlusNegH->Fill(InvMassTrackPlusNegMuon,weight);
      InvMassTracksPlusMuonsH->Fill(InvMassTracksPlusMuons,weight);
      InvMassTrackPosH->Fill(InvMassTrackPosMuon,weight);
      InvMassTrackNegH->Fill(InvMassTrackNegMuon,weight);


      if (IdTrackPosSelected==0)
	InvMassTrackPlusPos_EleH->Fill(InvMassTrackPlusPosMuon,weight);
      if (IdTrackPosSelected==1)
	InvMassTrackPlusPos_MuH->Fill(InvMassTrackPlusPosMuon,weight);
      if (IdTrackPosSelected==2)
	InvMassTrackPlusPos_HadH->Fill(InvMassTrackPlusPosMuon,weight);      


      if (IdTrackNegSelected==0)
	InvMassTrackPlusNeg_EleH->Fill(InvMassTrackPlusNegMuon,weight);
      if (IdTrackNegSelected==1)
	InvMassTrackPlusNeg_MuH->Fill(InvMassTrackPlusNegMuon,weight);
      if (IdTrackNegSelected==2)
	InvMassTrackPlusNeg_HadH->Fill(InvMassTrackPlusNegMuon,weight);      


      InvMassH2H->Fill(InvMassH2,weight);

      //bjet matching combinations 1D
      if(posMuMatched20CSVL && negMuMatched20CSVL){
	InvMassTrackPlusPos_bothBJets20H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_bothBJets20H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(posMuMatched20CSVL && !negMuMatched20CSVL){
	InvMassTrackPlusPos_posBJet20H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_posBJet20H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(!posMuMatched20CSVL && negMuMatched20CSVL){
	InvMassTrackPlusPos_negBJet20H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_negBJet20H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(!posMuMatched20CSVL && !negMuMatched20CSVL){
	InvMassTrackPlusPos_NoBJets20H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_NoBJets20H->Fill(InvMassTrackPlusNegMuon,weight);
      }

      if(posMuMatched15CSVL && negMuMatched15CSVL){
	InvMassTrackPlusPos_bothBJets15H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_bothBJets15H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(posMuMatched15CSVL && !negMuMatched15CSVL){
	InvMassTrackPlusPos_posBJet15H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_posBJet15H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(!posMuMatched15CSVL && negMuMatched15CSVL){
	InvMassTrackPlusPos_negBJet15H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_negBJet15H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      if(!posMuMatched15CSVL && !negMuMatched15CSVL){
	InvMassTrackPlusPos_NoBJets15H->Fill(InvMassTrackPlusPosMuon,weight);
	InvMassTrackPlusNeg_NoBJets15H->Fill(InvMassTrackPlusNegMuon,weight);
      }
      ///


      // mass(trk,mu-) in bins of mass(trk,mu+)
      
      if (InvMassTrackPlusPosMuon<=1.)
	InvMassTrackPlusNeg_Pos1GeV_SignalH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>1. && InvMassTrackPlusPosMuon<=2.)
	InvMassTrackPlusNeg_Pos2GeV_SignalH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>2. && InvMassTrackPlusPosMuon<=3.)
	InvMassTrackPlusNeg_Pos3GeV_SignalH->Fill(InvMassTrackPlusNegMuon,weight);
      
      if (InvMassTrackPlusPosMuon>3. && InvMassTrackPlusPosMuon<=4.)
	InvMassTrackPlusNeg_Pos4GeV_SignalH->Fill(InvMassTrackPlusNegMuon,weight);

      if (InvMassTrackPlusPosMuon>4.)
	InvMassTrackPlusNeg_PosGT4GeV_SignalH->Fill(InvMassTrackPlusNegMuon,weight);


      // mass(trk,mu+) in bins of mass(trk,mu-)
      if (InvMassTrackPlusNegMuon<=1.)
	InvMassTrackPlusPos_Neg1GeV_SignalH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>1. && InvMassTrackPlusNegMuon<=2.)
	InvMassTrackPlusPos_Neg2GeV_SignalH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>2. && InvMassTrackPlusNegMuon<=3.)
	InvMassTrackPlusPos_Neg3GeV_SignalH->Fill(InvMassTrackPlusPosMuon,weight);
      
      if (InvMassTrackPlusNegMuon>3. && InvMassTrackPlusNegMuon<=4.)
	InvMassTrackPlusPos_Neg4GeV_SignalH->Fill(InvMassTrackPlusPosMuon,weight);

      if (InvMassTrackPlusNegMuon>4.)
	InvMassTrackPlusPos_NegGT4GeV_SignalH->Fill(InvMassTrackPlusPosMuon,weight);

      
      InvMassTrackPlusMuon2D_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      float sigmaTrkEffPos = 0.04;
      float sigmaTrkEffNeg = 0.04;

      if (IdTrackPosSelected==0)
	sigmaTrkEffPos = 0.10;
      else if (IdTrackPosSelected==1)
	sigmaTrkEffPos = 0.06;

      if (IdTrackNegSelected==0) 
	sigmaTrkEffNeg = 0.10;
      else if (IdTrackNegSelected==1)
	sigmaTrkEffNeg = 0.06;
      
      float sigmaTrkEff = sigmaTrkEffPos + sigmaTrkEffNeg;
      
      InvMassTrackPlusMuon2D_UpH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight*(1+sigmaTrkEff));
      InvMassTrackPlusMuon2D_DownH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight*(1-sigmaTrkEff));


      //2D for different jet multiplicities
      if (nJets20CSVL==0)
	InvMassTrackPlusMuon2D_0CSVLH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      else if (nJets20CSVL==1)
	InvMassTrackPlusMuon2D_1CSVLH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      else
	InvMassTrackPlusMuon2D_2CSVLH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      if (nJets20CSVM==0)
        InvMassTrackPlusMuon2D_0CSVMH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      else if (nJets20CSVM==1)
        InvMassTrackPlusMuon2D_1CSVMH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      else
        InvMassTrackPlusMuon2D_2CSVMH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      //2D for combinations of systems matching or not bjets
      if(posMuMatched20CSVL && negMuMatched20CSVL) { 
	InvMassTrackPlusMuon2D_bothBJets20_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	nJets20CSVL_DoubleTagH->Fill(nJets20CSVL,weight);
	idTrackPos_DoubleTagH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_DoubleTagH->Fill(float(IdTrackNegSelected),weight);
      }
      if(posMuMatched20CSVL && !negMuMatched20CSVL) { 
	InvMassTrackPlusMuon2D_posBJet20_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	nJets20CSVL_SingleTagH->Fill(nJets20CSVL,weight);
	idTrackPos_SingleTagH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_SingleTagH->Fill(float(IdTrackNegSelected),weight);
      }
      if(!posMuMatched20CSVL && negMuMatched20CSVL) { 
	InvMassTrackPlusMuon2D_negBJet20_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	nJets20CSVL_SingleTagH->Fill(nJets20CSVL,weight);
	idTrackPos_SingleTagH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_SingleTagH->Fill(float(IdTrackNegSelected),weight);
	
      }
      if(!posMuMatched20CSVL && !negMuMatched20CSVL) { 
	InvMassTrackPlusMuon2D_NoBJets20_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
	nJets20CSVL_NoTagH->Fill(nJets20CSVL,weight);
	idTrackPos_NoTagH->Fill(float(IdTrackPosSelected),weight);
	idTrackNeg_NoTagH->Fill(float(IdTrackNegSelected),weight);
      }


      // Tracks id's
      if (IdTrackPosSelected==1&&IdTrackNegSelected==1) // MuMu
	InvMassTrackPlusMuon2D_MuMuH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      if (IdTrackPosSelected==1&&IdTrackNegSelected==0) // MuEle
	InvMassTrackPlusMuon2D_MuEleH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);

      if (IdTrackPosSelected==0&&IdTrackNegSelected==1) // MuEle (swapped)
	InvMassTrackPlusMuon2D_MuEleH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
      
      if (IdTrackPosSelected==1&&IdTrackNegSelected==2) // MuHad
	InvMassTrackPlusMuon2D_MuHadH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      
      if (IdTrackPosSelected==2&&IdTrackNegSelected==1) // MuHad (swapped)
	InvMassTrackPlusMuon2D_MuHadH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
	
      if (IdTrackPosSelected==0&&IdTrackNegSelected==0) // EleEle
	InvMassTrackPlusMuon2D_EleEleH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      
      if (IdTrackPosSelected==0&&IdTrackNegSelected==2) // EleHad
	InvMassTrackPlusMuon2D_EleHadH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      
      if (IdTrackPosSelected==2&&IdTrackNegSelected==0) // EleHad (swapped)
	InvMassTrackPlusMuon2D_EleHadH->Fill(InvMassTrackPlusNegMuon,InvMassTrackPlusPosMuon,weight);
      
      if (IdTrackPosSelected==2&&IdTrackNegSelected==2) // HadHad 
	InvMassTrackPlusMuon2D_HadHadH->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);



      if(posMuMatched15CSVL && negMuMatched15CSVL)InvMassTrackPlusMuon2D_bothBJets15_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      if(posMuMatched15CSVL && !negMuMatched15CSVL)InvMassTrackPlusMuon2D_posBJet15_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      if(!posMuMatched15CSVL && negMuMatched15CSVL)InvMassTrackPlusMuon2D_negBJet15_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);
      if(!posMuMatched15CSVL && !negMuMatched15CSVL)InvMassTrackPlusMuon2D_NoBJets15_H->Fill(InvMassTrackPlusPosMuon,InvMassTrackPlusNegMuon,weight);


   }
    
   if(iFiles==nFiles-1) std::cout <<"i think i am done"<< std::endl;
    
   delete _tree;
   if(!isData) delete _genTree;
   file_->Close();
  }
  
  for(int x=1; x<InvMassTrackPlusPos_2H->GetNbinsX()+1; x++){
    for(int y=1; y<InvMassTrackPlusNeg_2H->GetNbinsX()+1;y++){
      InvMass2DTrackPlusProbe_2H->SetBinContent(x,y,(InvMassTrackPlusPos_2H->GetBinContent(x))*(InvMassTrackPlusNeg_2H->GetBinContent(y)));
    }
    
  }
  for(int x=1; x<InvMassTrackPlusPos_3H->GetNbinsX()+1; x++){
    for(int y=1; y<InvMassTrackPlusNeg_3H->GetNbinsX()+1;y++){
      InvMass2DTrackPlusProbe_3H->SetBinContent(x,y,(InvMassTrackPlusPos_3H->GetBinContent(x))*(InvMassTrackPlusNeg_3H->GetBinContent(y)));
    }
    
  }
  for(int x=1; x<InvMassTrackPlusPos_4H->GetNbinsX()+1; x++){
    for(int y=1; y<InvMassTrackPlusNeg_4H->GetNbinsX()+1;y++){
      InvMass2DTrackPlusProbe_4H->SetBinContent(x,y,(InvMassTrackPlusPos_4H->GetBinContent(x))*(InvMassTrackPlusNeg_4H->GetBinContent(y)));
    }
    
  }
  

  for(int x=1; x<InvMassTrackPlusNeg_InvPosQH->GetNbinsX()+1; x++){
    InvMassTrackPlusMuon_QCDH->SetBinContent(x,(InvMassTrackPlusNeg_InvPosQH->GetBinContent(x))+(InvMassTrackPlusPos_InvNegQH->GetBinContent(x)));
  }
  
 for(int x=1; x<InvMassTrackPlusMuon_QCDH->GetNbinsX()+1; x++){
   for(int y=1; y<InvMassTrackPlusMuon_QCDH->GetNbinsX()+1;y++){
     InvMass2DTrackPlusMuon_QCDH->SetBinContent(x,y,(InvMassTrackPlusMuon_QCDH->GetBinContent(x))*(InvMassTrackPlusMuon_QCDH->GetBinContent(y)));
   }
 }

 //Unfolded histogram,"isolated muon+track" system constructed from QCD events
 int i=0;
 for(int x=1; x<InvMassTrackPlusMuon_QCDH->GetNbinsX()+1; x++){
   for(int y=1; y<InvMassTrackPlusMuon_QCDH->GetNbinsX()+1;y++){
     i++;
     InvMassTrackPlusMuon_QCD_UnfH->SetBinContent(i,InvMass2DTrackPlusMuon_QCDH->GetBinContent(x,y));
   }
 }
 i=0;
 for(int x=1; x<InvMassTrackPlusMuon2D_H->GetNbinsX()+1; x++){
   for(int y=1; y<InvMassTrackPlusMuon2D_H->GetNbinsX()+1;y++){
     i++;
     InvMassTrackPlusMuon2D_UnfH->SetBinContent(i,InvMassTrackPlusMuon2D_H->GetBinContent(x,y));
   }
 }


 std::cout <<"number of bins in unfolded "<<i<< std::endl;
 
 //OS vs OS contructed for QCD background
 for(int x=1; x<InvMassTrackPlusPos_SwapQH->GetNbinsX()+1; x++){
    InvMassTrackPlusMuon_OSOSH->SetBinContent(x,(InvMassTrackPlusPos_SwapQH->GetBinContent(x))+(InvMassTrackPlusNeg_SwapQH->GetBinContent(x)));
  }



 for(int x=1; x<InvMassTrackPlusMuon_OSOSH->GetNbinsX()+1; x++){
   for(int y=1; y<InvMassTrackPlusMuon_OSOSH->GetNbinsX()+1;y++){
     InvMass2DTrackPlusMuon_OSOSH->SetBinContent(x,y,(InvMassTrackPlusMuon_OSOSH->GetBinContent(x))*(InvMassTrackPlusMuon_OSOSH->GetBinContent(y)));
   }
   
 }
 //Unfolding
 i=0;
 for(int x=1; x<InvMass2DTrackPlusMuon_OSOSH->GetNbinsX()+1; x++){
   for(int y=1; y<InvMass2DTrackPlusMuon_OSOSH->GetNbinsX()+1;y++){
     i++;
     InvMassTrackPlusMuon_OSOS_UnfH->SetBinContent(i,InvMass2DTrackPlusMuon_OSOSH->GetBinContent(x,y));
   }
 }

//   outputFileLikelihood->cd("");
//   outputFileLikelihood->Write();
//   outputFileLikelihood->Close();
//   delete outputFileLikelihood;

// float xInputEvents = float(inputEventsH->GetEntries()); 
 nProcessed = int(inputEventsH->GetEntries());
 std::cout <<nProcessed<< std::endl;
 
 file->cd("");
 
 file->Write();
 
 file->Close();
 
 delete file;
 
}

 

