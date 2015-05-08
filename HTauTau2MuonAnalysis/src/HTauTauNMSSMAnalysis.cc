// -*- C++ -*-
// system include files
#include <memory>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


//#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauVertexSelector.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetMuons.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonSelector.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonTrackSelector.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonIsolator.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauKinematicsCalculator.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetPFParticles.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h>
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <string.h>
//#include <vector.h>
#include <TMath.h>

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/TrackReco/interface/Track.h>

#include <DataFormats/TrackReco/interface/HitPattern.h>

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauImpactParameter.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauUtils.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h" 
#include "DataFormats/PatCandidates/interface/Photon.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

#include "DataFormats/Math/interface/LorentzVector.h"

#include <algorithm>
#include <vector>
#include <utility>

#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGeneralInfo.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTriggerInfo.h"
//#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMassInfo.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauPileUpInfo.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauVertexInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

//#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"
//#include "RecoJets/JetProducers/interface/PileupJetIdentifier.h"
//#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

inline bool PFJetPairPtSorter(const std::pair<reco::PFJet, const reco::PFJet*>& pair1, const std::pair<reco::PFJet, const reco::PFJet*>& pair2) { return pair1.first.pt() > pair2.first.pt(); }


class HTauTauNMSSMAnalysis : public edm::EDAnalyzer {
   public:
      explicit HTauTauNMSSMAnalysis(const edm::ParameterSet&);
      ~HTauTauNMSSMAnalysis();


   private:
      void compX1X2byCollinearApprox(double& x1, double& x2, double pxLeg1, double pyLeg1, double pxLeg2, double pyLeg2, double pxMEt, double pyMEt)
      {
        double x1_numerator = pxLeg1*pyLeg2 - pxLeg2*pyLeg1;
        double x1_denominator = pyLeg2*(pxLeg1 + pxMEt) - pxLeg2*(pyLeg1 + pyMEt);
        x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
        //std::cout << "x1 = " << x1 << std::endl;

	double x2_numerator = x1_numerator;
	double x2_denominator = pxLeg1*(pyLeg2 + pyMEt) - pyLeg1*(pxLeg2 + pxMEt);
	x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
	//std::cout << "x2 = " << x2 << std::endl;
      };

     double getPhysX(double x, bool& isWithinPhysRange)
      { 
	double physX = x;

	isWithinPhysRange = true;
	
	if ( x < 0. ) {
	  physX = 0.;
	  isWithinPhysRange = false;
	}
	
	if ( x > 1. ) {
	  physX = 1.;
	  isWithinPhysRange = false;
	}
	
	return physX;
      };
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &setup);
      virtual void endJob() ;
      bool isDmeson(int pdg);
      bool isBmeson(int pdg);
      bool isJPsi(int pdg);
      bool isYpsilon(int pdg);
	// void fillMuonMVAVariables(
// 		const reco::Muon& posMuon,
// 		const reco::Muon& negMuon,
// 		const reco::Vertex& primaryVertex,
// 		const std::vector<reco::PFCandidate> pfCandidates,
// 		const edm::Event&, const edm::EventSetup&);
	bool isPFMuon(const pat::MuonRef& muon, const std::vector<reco::PFCandidate> pfCandidates);

  // ----------member data ---------------------------

// 	MuonMVAEstimator fMuonIsoMVA;
// 	MuonMVAEstimator fMuonIDMVA;

	//HTauTauGeneralInfo ntpGeneralInfo;
	HTauTauTriggerInfo ntpTriggerInfo;
  	//HTauTauMassInfo<pat::Muon, pat::Muon> ntpMassInfo;
	HTauTauPileUpInfo ntpPileUpInfo;
	HTauTauVertexInfo ntpVertexInfo;
  // ** sources
  edm::InputTag _vtxSrc;
  double _zVtxCut;
  double _dVtxCut;
  double _probVtxCut;

  edm::InputTag _muonSrc;
  //edm::InputTag _origMuonSrc;
  edm::InputTag _patMetSrc;
  //edm::InputTag _patMetMVASrc;
  edm::InputTag _jetSrc;
  //edm::InputTag _puJetIdSrc;
  //edm::InputTag _puJetMvaSrc;
  edm::InputTag _recoJetSrc;
  std::string _jetCorrService;
  edm::InputTag _photonSrc;
  edm::InputTag _genParticleSrc;
  //edm::InputTag _origGenParticleSrc;


  edm::InputTag _pfParticleSrc;
  edm::InputTag _pfNoPileUpSrc;
  edm::InputTag _pfPileUpSrc;


  // ** cuts and steering cards ----->

  double _ptMin;
  double _etaMax;
  double _ptMinHard;
  double _isoCut;
  double _jetEnergyMinCut;
  int    _nOfJetsCut;
  double _metCut;
  double _diLepDPhiCut;
  double _deltaR;

  bool _sameSignMuons;
  int _printOut;
  bool _doMC;

  int _year;
  std::string _period;

 //photons
  //  int nophotons;
  //  int numberOfPhotons;
  //  float photonPt[50];
  //  float photonEt[50];
  //  float photonEta[50];
  //  float photonPhi[50];
  //  float photonPx[50];
  //  float photonPy[50];
  //  float photonPz[50];
  //  float photonEn[50];
  
  // ** jets
  int _nJets;
  int noJets;
  int noRecoJets;

  bool jetIDLoose[50];
  bool jetIDMedium[50];
  bool jetIDTight[50];

  float jetPt[50];
  float jetEt[50];
  float jetEta[50];
  float jetPhi[50];
  float jetPx[50];
  float jetPy[50];
  float jetPz[50];
  float jetEn[50];
  float jetMass[50];
  float jetArea[50];
  float jetCorrFactor[50];

  int noConstituents[50];
  float neuEMEnergyFr[50];
  float neuHadronEnergyFr[50];
  int neuHadronMultiplicity[50];
  int neuMultiplicity[50];
  float photonEnergyFr[50];
  int photonMultiplicity[50];
  float hfEMEnergyFr[50];
  int hfEMMultiplicity[50];
  float hfHadronEnergyFr[50];
  int hfHadronMultiplicity[50];
  int noChargedConstituents[50];
  float chargedEnergyFr[50];
  int chargedMultiplicity[50];
  float electronEnergyFr[50];
  int electronMultiplicity[50];
  float muonEnergyFr[50];
  int muonMultiplicity[50];

  /*float puJetInputJetEta[50];
  float puJetInputJetPt[50];
  float puJetInputJetNCharged[50];
  float puJetInputJetNNeutrals[50];
  float puJetInputJetDZ[50];
  float puJetInputJetNParticles[50];
  float puJetInputJetDR2Mean[50];
  float puJetInputJetDRMean[50];
  float puJetInputJetFrac01[50];
  float puJetInputJetFrac02[50];
  float puJetInputJetFrac03[50];
  float puJetInputJetFrac04[50];
  float puJetInputJetFrac05[50];
  float puJetInputJetFrac06[50];
  float puJetInputJetFrac07[50];
  float puJetInputJetBeta[50];
  float puJetInputJetBetaStar[50];
  float puJetInputJetBetaClassic[50];
  float puJetInputJetBetaStarClassic[50];
  float puJetInputJetPtD[50];
  float puJetInputJetNvtx[50];

  float puJetFullMVA[50];
  bool puJetFullLoose[50];
  bool puJetFullMedium[50];
  bool puJetFullTight[50];*/

  bool genJetExist[50];
  float genJetPx[50];
  float genJetPy[50];
  float genJetPz[50];
  float genJetEn[50];
  int genJetCollisionId[50];

  bool genPartonExist[50];
  int  genPartonPdg[50];
  float genPartonPx[50];
  float genPartonPy[50];
  float genPartonPz[50];
  float genPartonEn[50];

  float jetBProbBJetTag[50];
  float jetProbBJetTag[50];
  float tcHPBJetTag[50];
  float tcHEBJetTag[50];
  float svHEBJetTag[50];
  float svHPBJetTag[50];
  float combSVBJetTag[50];
  float combSVMVABJetTag[50];

  bool recoJetIDLoose[50];
  bool recoJetIDMedium[50];
  bool recoJetIDTight[50];

  float recoJetPt[50];
  float recoJetEt[50];
  float recoJetEta[50];
  float recoJetPhi[50];
  float recoJetMass[50];
  float recoJetArea[50];

  // ***************************
  // ** Secondary vertices *****
  // ***************************

  int nSv;
  float fdSv[20];
  float fdSigSv[20];
  float massSv[20];
  int   nTrkSv[20];
  float pxTrkSv[20][10];
  float pyTrkSv[20][10];
  float pzTrkSv[20][10];
  int jetSv[20];

  // ***************************
  // ** PU related variables ***
  // ***************************


  int _nPtBinsPFSum;
  float _pfChargedSum[15];
  int   _nPfCharged[15]; 

  float _pfNeutralSum[15];
  int   _nPfNeutral[15];

  float _pfChargedPileUpSum[15];
  int   _nPfChargedPileUp[15];

  float _pfChargedNoPileUpSum[15];
  int   _nPfChargedNoPileUp[15];

  float _pfNeutralNoPileUpSum[15];
  int   _nPfNeutralNoPileUp[15];

  float _ptThresholds[15]; 

  float _isoPUChargedHadPFsNeg[15];
  float _isoPUChargedHadPFsPos[15];


  //tracks coming from the vertex with highest pt2 sum
  int nPVtracks;
  float pxTrack[1000];
  float pyTrack[1000];
  float pzTrack[1000];
  float ptTrack[1000];
  float phiTrack[1000];
  float etaTrack[1000];
  double chi2Track[1000];
  int chargeTrack[1000];

  int nPVPostracks;
  float ptPVTracksPos[1000];
  float phiPVTracksPos[1000];
  float etaPVTracksPos[1000];
  float pxPVTracksPos[1000];
  float pyPVTracksPos[1000];
  float pzPVTracksPos[1000];
  float  qPVTracksPos[1000];

  int nPVNegtracks;
  float ptPVTracksNeg[1000];
  float phiPVTracksNeg[1000];
  float etaPVTracksNeg[1000];
  float pxPVTracksNeg[1000];
  float pyPVTracksNeg[1000];
  float pzPVTracksNeg[1000];
  float  qPVTracksNeg[1000];

  //all NoPU candidates.
  int npfNoPU; 
  float pxpfNoPU[10000];
  float pypfNoPU[10000];
  float pzpfNoPU[10000];
  float ptpfNoPU[10000];
  float phipfNoPU[10000];
  float etapfNoPU[10000];
  float enpfNoPU[10000];
  float chargepfNoPU[10000];
  float ippfNoPU[10000];
  float ipSigpfNoPU[10000];
  float ipZpfNoPU[10000];
  float ipZSigpfNoPU[10000];
  float ip3DpfNoPU[10000];
  float ip3DSigpfNoPU[10000];
  int   idpfNoPU[10000];

  //NoPU PF around first muon
  int ntracksPos;
  float pxTrackPos[1000];
  float pyTrackPos[1000];
  float pzTrackPos[1000];
  float enTrackPos[1000];
  float ptTrackPos[1000];
  float etaTrackPos[1000];
  float phiTrackPos[1000];
  int idTrackPos[1000];
  float qTrackPos[1000];
  float ipTrackPos[1000];
  float ipSigTrackPos[1000];
  float ipZTrackPos[1000];
  float ipZSigTrackPos[1000];
  float ip3DTrackPos[1000];
  float ip3DSigTrackPos[1000];

  //PU PF around first muon
  int ntracksPileUpPos;
  float pxTrackPileUpPos[1000];
  float pyTrackPileUpPos[1000];
  float pzTrackPileUpPos[1000];
  float enTrackPileUpPos[1000];
  float ptTrackPileUpPos[1000];
  float etaTrackPileUpPos[1000];
  float phiTrackPileUpPos[1000];
  int idTrackPileUpPos[1000];
  float qTrackPileUpPos[1000];

  //NoPU PF around second muon
  int ntracksNeg;
  float pxTrackNeg[1000];
  float pyTrackNeg[1000];
  float pzTrackNeg[1000];
  float enTrackNeg[1000];
  float ptTrackNeg[1000];
  float etaTrackNeg[1000];
  float phiTrackNeg[1000];
  int idTrackNeg[1000];
  float qTrackNeg[1000];
  float ipTrackNeg[1000];
  float ipSigTrackNeg[1000];
  float ipZTrackNeg[1000];
  float ipZSigTrackNeg[1000];
  float ip3DTrackNeg[1000];
  float ip3DSigTrackNeg[1000];

  //PU PF around second muon
  int ntracksPileUpNeg;
  float pxTrackPileUpNeg[1000];
  float pyTrackPileUpNeg[1000];
  float pzTrackPileUpNeg[1000];
  float enTrackPileUpNeg[1000];
  float ptTrackPileUpNeg[1000];
  float etaTrackPileUpNeg[1000];
  float phiTrackPileUpNeg[1000];
  int idTrackPileUpNeg[1000];
  float qTrackPileUpNeg[1000];

  // ** histograms

  TH1F * h_cutFlow_;
  TH1F * h_numLep_;
  TH1F * h_ptLep_;
  TH1F * h_etaLep_;
  TH1F * h_isoLep_;
  TH1F * h_ptHard_;
  TH1F * h_jetEt_;
  TH1F * h_numJets_;
  TH1F * h_2LepDeltaPhi_;
  TH1F * h_MET_;

  TH1F * h_diLepMass_;
  TH1F * h_diLepMassEtaCut_;

  TH1F * h_cosNegDiLepton_;
  TH1F * h_mTDiLepton_;
  TH1F * h_mTDiLeptonMET_;

  TH1F * h_lepPosMETDPhi_;
  TH1F * h_lepNegMETDPhi_;

  TH1F * h_lepEDif_;
  TH1F * h_lepPtDif_;

  TH1F * h_diLepMETDPhi_;

  TH1F * h_diLeptonPL_;
  TH1F * h_diLeptonEta_;

  TH1F * h_ipLeptonPosTrk_;
  TH1F * h_ipLeptonNegTrk_;
  TH1F * h_ipSigLeptonPosTrk_;
  TH1F * h_ipSigLeptonNegTrk_;

  TH1F * h_diTauMass_;
  TH1F * h_diTaudiLepMass_;
  TH1F * h_lepETauRestFrame_;
  TH1F * h_lepEDifDiTauRestFrame_;
  TH1F * h_lepESumDiTauRestFrame_;
  TH1F * h_cosNegDiTau_;
  TH1F * h_lepTauERatio_;
  TH1F * h_diTauPL_;
  TH1F * h_diTauEta_;
  TH1F * h_xVtxStrange_;
  TH1F * h_yVtxStrange_;
  TH1F * h_xVtxNormal_;
  TH1F * h_yVtxNormal_;
  TH1F * h_nHitsStrange_;
  TH1F * h_nVtxHitsStrange_;
  TH1F * h_probTrkStrange_;
  TH1F * h_nHitsNormal_;
  TH1F * h_nVtxHitsNormal_;
  TH1F * h_probTrkNormal_;

  TH1F * h_pt_pfNeutral_;

  // ** input and passed events ---->

  int _event;
  int _selEvent;

  // ** trees --->

  TTree * _tree;

  TTree * _genTree; 

  float H2Pt;
  float H2Px;
  float H2Py;
  float H2Pz;
  float H2Phi;
  float H2Eta;
  float H2mass;

  float ptFirstH1;
  float pxFirstH1;
  float pyFirstH1;
  float pzFirstH1;
  float phiFirstH1;
  float etaFirstH1;
  float massFirstH1;

  float ptSecondH1; 
  float pxSecondH1;
  float pySecondH1;
  float pzSecondH1;
  float phiSecondH1;
  float etaSecondH1;
  float massSecondH1;

  float tau1pt;
  float tau1px;
  float tau1py;
  float tau1pz;
  float tau1eta;
  float tau1phi;
  int FSO1;
  int finalstateIdTau1[20];

  float tau2pt;
  float tau2px;
  float tau2py;
  float tau2pz;
  float tau2eta;
  float tau2phi;
  int FSO2;
  int finalstateIdTau2[20];

  float tau3pt;
  float tau3px;
  float tau3py;
  float tau3pz;
  float tau3eta;
  float tau3phi;
  int FSO3;
  int finalstateIdTau3[20];

  float tau4pt;
  float tau4px;
  float tau4py;
  float tau4pz;
  float tau4eta;
  float tau4phi;
  int FSO4;
  int finalstateIdTau4[20];

  int nPartFirstH1;
  int pdgIdPartFirstH1[20];
  int qPartFirstH1[20];
  float pxPartFirstH1[20];
  float pyPartFirstH1[20];
  float pzPartFirstH1[20];
  float ptPartFirstH1[20];
  float etaPartFirstH1[20];
  float phiPartFirstH1[20];

  int nPartAroundMuFirstH1;
  int pdgIdPartAroundMuFirstH1[20];
  int qPartAroundMuFirstH1[20];
  float pxPartAroundMuFirstH1[20];
  float pyPartAroundMuFirstH1[20];
  float pzPartAroundMuFirstH1[20];
  float ptPartAroundMuFirstH1[20];
  float etaPartAroundMuFirstH1[20];
  float phiPartAroundMuFirstH1[20];

  int nPartSecondH1;
  int pdgIdPartSecondH1[20];
  int qPartSecondH1[20];
  float pxPartSecondH1[20];
  float pyPartSecondH1[20];
  float pzPartSecondH1[20];
  float ptPartSecondH1[20];
  float etaPartSecondH1[20];
  float phiPartSecondH1[20];

  int nPartAroundMuSecondH1;
  int pdgIdPartAroundMuSecondH1[20];
  int qPartAroundMuSecondH1[20];
  float pxPartAroundMuSecondH1[20];
  float pyPartAroundMuSecondH1[20];
  float pzPartAroundMuSecondH1[20];
  float ptPartAroundMuSecondH1[20];
  float etaPartAroundMuSecondH1[20];
  float phiPartAroundMuSecondH1[20];

  // lepton isolation --->

  float _isoECalNeg;
  float _isoHCalNeg;
  float _isoTrkNeg;
  float _isoPFsNeg;
  float _isoChargedHadPFsNeg;
  float _isoNeutralHadPFsNeg;
  float _isoPhotonsPFsNeg;
  float _isoNeg; // relative Iso

  float _isoECalPos;
  float _isoHCalPos;
  float _isoTrkPos;
  float _isoPFsPos;
  float _isoChargedHadPFsPos;
  float _isoNeutralHadPFsPos;
  float _isoPhotonsPFsPos;
  float _isoPos; // relative Iso

  float _mvaIDMuPos;
  float _mvaIDMuNeg;
  
  float _mvaIsoMuPos;
  float _mvaIsoMuNeg;

  bool _isPosGlobalMu;
  bool _isNegGlobalMu;

  bool _isPosTrackerMu;
  bool _isNegTrackerMu;

  bool _isPosPFMu;
  bool _isNegPFMu;

  int _nPixelHitsPos;
  int _nPixelHitsNeg;

  int _nTrackerHitsPos;
  int _nTrackerHitsNeg;
  
  int _nMuonHitsPos;
  int _nMuonHitsNeg;

  int _nMuonStationsPos;
  int _nMuonStationsNeg;

  float _chi2PosMu;
  float _chi2NegMu;

  float _ndofPosMu;
  float _ndofNegMu;

 

  // lepton ip and ip sig

  float _ipSig2DPos[100];
  float _ipSig2DNeg[100];
  float _ip2DPos[100];
  float _ip2DNeg[100];


  float _ipSig3DPos[100];
  float _ipSig3DNeg[100];
  float _ip3DPos[100];
  float _ip3DNeg[100];

  float _ipSigZPos[100];
  float _ipSigZNeg[100];
  float _ipZPos[100];
  float _ipZNeg[100];

  float _twoMuonDist2D;
  float _twoMuonDist2DE;
  
  float _twoMuonDist3D;
  float _twoMuonDist3DE;

  float _twoMuonDistRPhi3D;
  float _twoMuonDistRPhi3DE;

  float _twoMuonDistRPhi2D;
  float _twoMuonDistRPhi2DE;

  // dilepton kinematics

  float _diLepMassVar;

  float _diLepDPhi;
  float _diLepDEta;
  float _diLepDR;
  float _met;

  float _posLepMETDPhi;
  float _negLepMETDPhi;

  float _diLepEta;
  float _negLepPt;
  float _posLepPt;
 float _negLepPx;
  float _posLepPx;
 float _negLepPy;
  float _posLepPy;
 float _negLepPz;
  float _posLepPz;
  float _negLepEta;
  float _posLepEta;
  float _negLepPhi;
  float _posLepPhi;
  float _posLepEn;
  float _negLepEn;
  float _negLepQ;
  float _posLepQ;
  bool _negLep_IsoMu24;
  bool _posLep_IsoMu24;
  bool _negLep_Mu20;
  bool _posLep_Mu20;
  bool _negLep_Mu30;
  bool _posLep_Mu30;
  bool _negLep_Mu40;
  bool _posLep_Mu40;
  bool _negLep_Mu17_Mu8;
  bool _posLep_Mu17_Mu8;
  bool _negLep_Mu17_TkMu8;
  bool _posLep_Mu17_TkMu8;



  float _diLepPt;
  float _diLepPhi;

  float _metPx;
  float _metPy;
  float _metPz;
  float _metEn;

  /*float _metCovXX;
  float _metCovXY;
  float _metCovYX;
  float _metCovYY;

  float _metMVAPx;
  float _metMVAPy;
  float _metMVAPz;
  float _metMVAEn;

  float _metMVACovXX;
  float _metMVACovXY;
  float _metMVACovYX;
  float _metMVACovYY;*/

  // original muons
  float _embedWeight;
  /*float _origDiLepMass;
  float _origNegLepPt;
  float _origPosLepPt;
  float _origNegLepEta;
  float _origPosLepEta;
  float _origNegLepPhi;
  float _origPosLepPhi;
  float _origNegLepChargedHadronIso;
  float _origPosLepChargedHadronIso;*/

  // ditau reconstrution --->

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

  // generator info

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

  int posMuonGenNoMothers;
  int PMmotherNoMothers;

  bool PMmotherExist;
  int PMmotherpdg;
  bool PMmomISmuon;
  float PMmothereta;
  float PMmotherphi;
  float PMmotherpt;
  float PMmotherpx;
  float PMmotherpy;
  float PMmotherpz;

  int posMuonGenNoGrmothers;
  int PMGrmotherNoMothers ;

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

  int negMuonGenNoMothers;
  int NMmotherNoMothers;

  bool NMmotherExist;
  int NMmotherpdg;
  bool NMmomISmuon;
  float NMmothereta;
  float NMmotherphi;
  float NMmotherpt;
  float NMmotherpx;
  float NMmotherpy;
  float NMmotherpz;

  int negMuonGenNoGrmothers;
  int NMGrmotherNoMothers ;

  bool NMGrmotherExist;
  int NMGrmotherpdg;
  float NMGrmothereta;
  float NMGrmotherphi;
  float NMGrmotherpt;
  float NMGrmotherpx;
  float NMGrmotherpy;
  float NMGrmotherpz;

  int _bosonPDG;
  int _decayPDG;
  int _origBosonPDG;
  int _origDecayPDG;

  bool _bosonExist;	
  float _ZMass;
  float _ZPx;
  float _ZPy;
  float _ZPz;
  int VBFid;
  bool isVBF;
  bool bbbarPhi;
  float _HMass;
  float _HPx;
  float _HPy;
  float _HPz;
 
  float benergy;
  float bPx;
  float bPy;
  float bPz;
  float bbarenergy;
  float bbarPx;
  float bbarPy;
  float bbarPz;

  float parton1pt;
  float parton1pdg;
  float parton1Px;
  float parton1Py;
  float parton1Pz;
  float parton1energy;
  float parton2pt;
  float parton2pdg;
  float parton2Px;
  float parton2Py;
  float parton2Pz;
  float parton2energy;

  bool _firstEvent;
  bool _firstTriggerEvent;

  float _posDecayMuonPx;
  float _posDecayMuonPy;
  float _posDecayMuonPz;
  int   _posDecayMuonPdg;

  float _negDecayMuonPx;
  float _negDecayMuonPy;
  float _negDecayMuonPz;
  int   _negDecayMuonPdg; 

  int _nFSRPos;
  int _nFSRNeg;

  int _FSRPosPdg[100];	
  float _FSRPosPx[100];
  float _FSRPosPy[100];
  float _FSRPosPz[100];

  int _FSRNegPdg[100];
  float _FSRNegPx[100];
  float _FSRNegPy[100];
  float _FSRNegPz[100];

	CLHEP::RandFlat * rndFlat;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HTauTauNMSSMAnalysis::HTauTauNMSSMAnalysis(const edm::ParameterSet& iConfig) : /*ntpGeneralInfo(iConfig),*/ ntpTriggerInfo(iConfig), ntpPileUpInfo(iConfig), ntpVertexInfo(iConfig)
{
   // initialization of cuts and steering cards ====>
  _vtxSrc         = iConfig.getParameter<edm::InputTag>  (  "PrimaryVertexSource" );
  _zVtxCut        = iConfig.getParameter<double>         (  "ZVertexCut" );
  _dVtxCut        = iConfig.getParameter<double>         (  "DVertexCut" );
  _probVtxCut     = iConfig.getParameter<double>         (  "ProbVertexCut" );

  _muonSrc        = iConfig.getParameter<edm::InputTag>  (  "MuonSource" );
  //_origMuonSrc    = iConfig.getParameter<edm::InputTag>  (  "OrigMuonSource" );
  _patMetSrc      = iConfig.getParameter<edm::InputTag>  (  "PatMETSource" );
  //_patMetMVASrc   = iConfig.getParameter<edm::InputTag>  (  "PatMETMVASource" );
  _jetSrc         = iConfig.getParameter<edm::InputTag>  (  "JetSource" );
  //_puJetIdSrc     = iConfig.getParameter<edm::InputTag>  (  "PuJetIdSource" );
  //_puJetMvaSrc    = iConfig.getParameter<edm::InputTag>  (  "PuJetMvaSource" );
  _recoJetSrc     = iConfig.getParameter<edm::InputTag>  (  "RecoJetSource" );
  _photonSrc      = iConfig.getParameter<edm::InputTag>  (  "PhotonSource" );
  _genParticleSrc = iConfig.getParameter<edm::InputTag>  (  "GenParticleSource" );
  //_origGenParticleSrc = iConfig.getParameter<edm::InputTag>  (  "OrigGenParticleSource" );

  _pfParticleSrc  = iConfig.getParameter<edm::InputTag>  (  "PFParticleSource" ); 
  _pfNoPileUpSrc  = iConfig.getParameter<edm::InputTag>  (  "pfNoPileUpSource" ); 
  _pfPileUpSrc  = iConfig.getParameter<edm::InputTag>    (  "pfPileUpSource" ); 

  _ptMin           = iConfig.getParameter<double>         (  "PtMinCut" );
  _etaMax          = iConfig.getParameter<double>         (  "EtaMaxCut" ); 
  _ptMinHard       = iConfig.getParameter<double>         (  "PtMinHardLeptonCut" );
  _isoCut          = iConfig.getParameter<double>         (  "IsolationCut" );
  _nOfJetsCut      = iConfig.getParameter<int>            (  "NumberOfJetsCut" );
  _jetEnergyMinCut = iConfig.getParameter<double>         (  "JetEnergyMinCut" );
  _metCut          = iConfig.getParameter<double>         (  "METCut" );
  _diLepDPhiCut    = iConfig.getParameter<double>         (  "DiLeptonDPhiCut"); 
  _deltaR          = iConfig.getParameter<double>         (  "IsoDeltaR" );

  _sameSignMuons   = iConfig.getParameter<bool>            (  "SameSignMuons" );
  _printOut        = iConfig.getParameter<int>            (  "PrintOut" );  
  _doMC            = iConfig.getParameter<bool>            (  "DoMC" );

  _year = iConfig.getParameter<int>("year");
  _period = iConfig.getParameter<std::string>("period");

  _nPtBinsPFSum = 15;
  
  float ptThresholds[15] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0};
 
  for (int iP=0; iP<_nPtBinsPFSum; ++iP) 
    _ptThresholds[iP] = ptThresholds[iP];

	std::vector<std::string> muonid_weightfiles;
	std::vector<std::string> muoniso_weightfiles;

	const char* CMSSW_BASE = getenv("CMSSW_BASE");
	if(!CMSSW_BASE) throw cms::Exception("Not in a CMSSW environment!");
	const std::string prefix = std::string(CMSSW_BASE) + "/src/";


// 	fMuonIDMVA.initialize("MuonID_BDTG", MuonMVAEstimator::kID, true, muonid_weightfiles);
// 	fMuonIsoMVA.initialize("MuonIso_BDTG_IsoRings", MuonMVAEstimator::kIsoRings, true, muoniso_weightfiles);
// 	fMuonIsoMVA.SetPrintMVADebug(kTRUE);

	// random number
	/*edm::Service<edm::RandomNumberGenerator> rng;
	if ( ! rng.isAvailable())
	{
		throw cms::Exception("Configuration")
			<< "EcalTBMCInfoProducer requires the RandomNumberGeneratorService\n"
			"which is not present in the configuration file.  You must add the service\n"
			"in the configuration file or remove the modules that require it.";
	}
	CLHEP::HepRandomEngine& engine = rng->getEngine();
	rndFlat = new CLHEP::RandFlat(engine);*/
}


HTauTauNMSSMAnalysis::~HTauTauNMSSMAnalysis()
{
}


//
// member functions
//

void HTauTauNMSSMAnalysis::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	//ntpGeneralInfo.beginRun(iRun, iSetup);
	ntpTriggerInfo.beginRun(iRun, iSetup);
	ntpPileUpInfo.beginRun(iRun, iSetup);
	ntpVertexInfo.beginRun(iRun, iSetup);
}


// ------------ method called to for each event  ------------
void
HTauTauNMSSMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;

  HTauTauUtils utils;
  //int numberoftaus=0;
  if (_doMC) {

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel(_genParticleSrc,genParticles);
    
    //	int nParticles = 50; 
    //	if (genParticles->size()<50)
    int nParticles = genParticles->size();
    bool foundH2=false;
    
    bool foundFirstHiggs = false;
    bool foundSecondHiggs = false;


    for (int iP=0; iP<nParticles;++iP) {
      GenParticleRef part(genParticles,iP);
      int pdg = part->pdgId();
      //int status = part->status();

      if (!foundH2) {
	if (pdg==35) {
	  foundH2= true;
	  //int nM = part->numberOfMothers();
	  /*for (int iM=0;iM<nM;++iM) {
	    const Candidate * mum = part->mother(iM);
	    //			std::cout <<"H2 mother number "<<iM<<" is  " <<mum->pdgId()<<std::endl;
	  }*/
	  //int nD = part->numberOfDaughters();
	  //	      std::cout<<"H2 daughters : "<<nD <<std::endl;
	  /*for (int iD=0;iD<nD;++iD) {
	    const Candidate * kid = part->daughter(iD);
	    //		std::cout <<"H2 daughter number "<<iD<<" is  " <<kid->pdgId()<<std::endl;
	    //		std::cout <<"H2 daughter status "<<iD<<" is  " <<kid->status()<<std::endl;
	    //		std::cout <<"H2 daughter pt "<<iD<<" is  " <<kid->pt()<<std::endl;
	  }*/
	  H2Pt = part->pt();
	  //	      std::cout <<"H2 pt store is  " <<part->pt()<<std::endl;
	  H2Px =  part->px();
	  H2Py=  part->py();
	  H2Pz=  part->pz();
	  H2Phi=  part->phi();
	  H2Eta=  part->eta();
	  H2mass=  part->mass();
	  
	}
      }
      
      if (foundFirstHiggs && (!foundSecondHiggs)) {
	if (pdg==25) {
	  foundSecondHiggs = true;
	  //std::cout <<"Second H1 info :  "<<std::endl;
	  //int nM = part->numberOfMothers();
	  /*for (int iM=0;iM<nM;++iM) {
	    const Candidate * mum = part->mother(iM);
	    //		std::cout <<"Second H1 mother number "<<iM<<" is  " <<mum->pdgId()<<std::endl;
	    //		std::cout <<"Second H1 mother status "<<iM<<" is  " <<mum->status()<<std::endl;
	    //		std::cout <<"Second H1 mother pt "<<iM<<" is  " <<mum->pt()<<std::endl;
	  }*/
	  int nD = part->numberOfDaughters();
	  for (int iD=0;iD<nD;++iD) {
	    const Candidate * kid = part->daughter(iD);
	    //		std::cout <<"Second H1 daughter number "<<iD<<" is  " <<kid->pdgId()<<std::endl;
	    //		std::cout <<"Second H1 daughter status "<<iD<<" is  " <<kid->status()<<std::endl;
	    //		std::cout <<"Second H1 daughter pt "<<iD<<" is  " <<kid->pt()<<std::endl;
	    if (kid->pdgId()==15){
	      //		  std::cout << "Hey WORLD!!!!!!!!!!!!!!I FOUND A TAU!!!!!!!!!!!!" <<std::endl;
	      //		  std::cout <<"with "<<  kid->numberOfDaughters()<<" kids"<<std::endl;
	      //		  std::cout <<"and status "<<  kid->status()<<std::endl;
	      //		  std::cout <<"and pt "<<  kid->pt()<<std::endl;
	      int ngd= kid->numberOfDaughters();
	      for (int igd=0;igd<ngd;++igd) {
		const Candidate * grandkid = kid->daughter(igd);		    
		//		  std::cout <<"of type "<<grandkid->pdgId()<<std::endl;
		//		    std::cout <<"and status"<<grandkid->status()<<std::endl;
		//		    std::cout <<"and pt"<<grandkid->pt()<<std::endl;
		int FSO=grandkid->numberOfDaughters();
		//		  std::cout <<"with "<<FSO<< "FSO"<<std::endl;
		
		if (grandkid->pdgId()==15){
		  for (int ifs=0;ifs<FSO;++ifs){
		    const Candidate * finalstate = grandkid->daughter(ifs);	
		    finalstateIdTau3[ifs] =finalstate->pdgId();
		    //		      std::cout <<"Final state object of type"<<finalstate->pdgId()<<std::endl;
		  }
		  FSO3=FSO;
		}
	      }
	      tau3pt=kid->pt();
	      tau3px=kid->px();
	      tau3py=kid->py();
	      tau3pz=kid->pz();
	      tau3eta=kid->eta();
	      tau3phi=kid->phi();
	      
	    }
	    if (kid->pdgId()==-15){
	      //std::cout << "Hey WORLD!!!!!!!!!!!!!!I FOUND AN ANTI- TAU!!!!!!!!!!!!" <<std::endl;
	      //		  std::cout <<"with "<<  kid->numberOfDaughters()<<" kids"<<std::endl;
	      //		  std::cout <<"and status "<<  kid->status()<<std::endl;
	      //		  std::cout <<"and pt "<<  kid->pt()<<std::endl;
	      int ngd= kid->numberOfDaughters();
	      for (int igd=0;igd<ngd;++igd) {
		const Candidate * grandkid = kid->daughter(igd);		    
		//		  std::cout <<"of type "<<grandkid->pdgId()<<std::endl;
		//		  std::cout <<"and status"<<grandkid->status()<<std::endl;
		//		  std::cout <<"and pt "<<  grandkid->pt()<<std::endl;
		int FSO=grandkid->numberOfDaughters();
		//		  std::cout <<"with "<<FSO<<"FSO"<<std::endl;
		
		if (grandkid->pdgId()==-15){
		  for (int ifs=0;ifs<FSO;++ifs){
		    const Candidate * finalstate = grandkid->daughter(ifs);	
		    finalstateIdTau4[ifs] =finalstate->pdgId();
		    //		      std::cout <<"Final state object of type"<<finalstate->pdgId()<<std::endl;
		  }
		  FSO4=FSO; 
		}
	      }
	      tau4pt=kid->pt();
	      tau4px=kid->px();
	      tau4py=kid->py();
	      tau4pz=kid->pz();
	      tau4eta=kid->eta();
	      tau4phi=kid->phi();
	    }
	  }
	  ptSecondH1= part->pt();
	  pxSecondH1 =  part->px();
	  pySecondH1 =  part->py();
	  pzSecondH1 =  part->pz();
	  phiSecondH1 =  part->phi();
	  etaSecondH1 =  part->eta();
	  massSecondH1 =  part->mass();
	}
      }
      
      if (!foundFirstHiggs) {
	if (pdg==25) {
	  //      std::cout <<"First H1 info :  "<<std::endl;
	  //int nM = part->numberOfMothers();
	  /*for (int iM=0;iM<nM;++iM) {
	    const Candidate * mum = part->mother(iM);
	    //		std::cout <<"First H1 mother number "<<iM<<" is  " <<mum->pdgId()<<std::endl;
	    //		std::cout <<"First H1 mother status "<<iM<<" is  " <<mum->status()<<std::endl;
	    //		std::cout <<"First H1 mother pt "<<iM<<" is  " <<mum->pt()<<std::endl;
	  }*/
	  int nD = part->numberOfDaughters();
	  for (int iD=0;iD<nD;++iD) {
	    const Candidate * kid = part->daughter(iD);
	    //		std::cout <<"First H1 daughter number "<<iD<<" is  " <<kid->pdgId()<<std::endl;
	    //		std::cout <<"First H1 daughter status "<<iD<<" is  " <<kid->status()<<std::endl;
	    //		std::cout <<"First H1 daughter pt "<<iD<<" is  " <<kid->pt()<<std::endl;
	    if (kid->pdgId()==15){
	      //		  std::cout << "Hey WORLD!!!!!!!!!!!!!!I FOUND A TAU!!!!!!!!!!!!" <<std::endl;
	      //		  std::cout <<"with "<<  kid->numberOfDaughters()<<" kids"<<std::endl;
	      //		  std::cout <<"and status "<<  kid->status()<<std::endl;
	      //		  std::cout <<"and pt "<<  kid->pt()<<std::endl;
	      int ngd= kid->numberOfDaughters();
	      for (int igd=0;igd<ngd;++igd) {
		const Candidate * grandkid = kid->daughter(igd);		    
		//		    std::cout <<"of type "<<grandkid->pdgId()<<std::endl;
		//		    std::cout <<"and status"<<grandkid->status()<<std::endl;
		//		    std::cout <<"and pt"<<grandkid->pt()<<std::endl;
		int FSO=grandkid->numberOfDaughters();
		//		    std::cout <<"with "<<FSO<<"FSO"<<std::endl;
		
		if (grandkid->pdgId()==15){
		  for (int ifs=0;ifs<FSO;++ifs){
		    const Candidate * finalstate = grandkid->daughter(ifs);	
		    finalstateIdTau1[ifs] =finalstate->pdgId();
		    //			std::cout <<"Final state object of type"<<finalstate->pdgId()<<std::endl;
		    
		  }
		  FSO1=FSO;
		}
		tau1pt=kid->pt();
		tau1px=kid->px();
		tau1py=kid->py();
		tau1pz=kid->pz();
		tau1eta=kid->eta();
		tau1phi=kid->phi();
		
	      }
	      if (kid->pdgId()==-15){
		//		  std::cout << "Hey WORLD!!!!!!!!!!!!!!I FOUND AN ANTI- TAU!!!!!!!!!!!!" <<std::endl;
		//		  std::cout <<"with "<<  kid->numberOfDaughters()<<" kids"<<std::endl;
		//		  std::cout <<"and status "<<  kid->status()<<std::endl;
		//		  std::cout <<"and pt "<<  kid->pt()<<std::endl;
		int ngd= kid->numberOfDaughters();
		for (int igd=0;igd<ngd;++igd) {
		  const Candidate * grandkid = kid->daughter(igd);		    
		  //		    std::cout <<"of type "<<grandkid->pdgId()<<std::endl;
		  //		    std::cout <<"and status"<<grandkid->status()<<std::endl;
		  //		    std::cout <<"and pt"<<grandkid->pt()<<std::endl;
		  int FSO=grandkid->numberOfDaughters();
		  //		    std::cout <<"with "<<FSO<<"FSO"<<std::endl;
		  
		  if (grandkid->pdgId()==-15){
		    for (int ifs=0;ifs<FSO;++ifs){
		      const Candidate* finalstate = grandkid->daughter(ifs);	
		      finalstateIdTau2[ifs] =finalstate->pdgId();
		      //			std::cout <<"Final state object of type"<<finalstate->pdgId()<<std::endl;
		    }
		    FSO2=FSO;
		  }
		}
		tau2pt=kid->pt();
		tau2px=kid->px();
		tau2py=kid->py();
		tau2pz=kid->pz();
		tau2eta=kid->eta();
		tau2phi=kid->phi();
	      }
	    }
	    
	    foundFirstHiggs = true;
	    ptFirstH1 =  part->pt();
	    pxFirstH1 =  part->px();
	    pyFirstH1 =  part->py();
	    pzFirstH1 =  part->pz();
	    phiFirstH1 =  part->phi();
	    etaFirstH1 =  part->eta();
	    massFirstH1 =  part->mass();
	  }
	}
	
	if (foundFirstHiggs && foundSecondHiggs) 
	  break;
	
      }
      
      int stableFirst[20];
      int stableSecond[20];
      int nStableFirst = 0;
      int nStableSecond = 0;
      nPartAroundMuFirstH1 = 0;
      nPartAroundMuSecondH1 = 0;
      bool foundFirstMuH1 = false;
      bool foundSecondMuH1 = false;
      
      for (int iP=0; iP<nParticles;++iP) {
	GenParticleRef part(genParticles,iP);
	int status = part->status();
	if (status != 1) continue;
	bool continueSearch = true;
	int nMothers = part->numberOfMothers();
	const Candidate * motherPart = NULL;
	if (nMothers!=1)
	  continueSearch = false;
	else 
	  motherPart = part->mother(0);
	while (continueSearch) {
	  int pdgMother = motherPart->pdgId();
	  float ptMother = motherPart->pt();
	  float ptDiffFirst = abs(ptMother-ptFirstH1);
	  float ptDiffSecond = abs(ptMother-ptSecondH1);
	  if (pdgMother==25 && (ptDiffFirst<1e-5 || ptDiffSecond<1e-5)) {
	    if (ptDiffFirst<1e-5) {
	      if (nStableFirst<19) {
		stableFirst[nStableFirst] = iP;
	      }
	      nStableFirst++;
	      continueSearch = false;
	    }
	    if (ptDiffSecond<1e-5) {
	      if (nStableSecond<19)
		stableSecond[nStableSecond] = iP;
	      nStableSecond++;
	      continueSearch = false;
	    }
	  }
	  else {
	    int nGrandMothers = motherPart->numberOfMothers();  
	    if (nGrandMothers!=1) 
	      continueSearch = false;
	    else {
	      motherPart = motherPart->mother(0);
	    }
	  }
	}
      }
      if (_printOut==-1) {
	std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout << std::endl;
	std::cout << "First H1 -----> " << std::endl;
      }
      if (nStableFirst>19)
	nStableFirst = 19;
      nPartFirstH1 = nStableFirst;
      for (int iP=0; iP<nStableFirst; iP++) {
	int iPlace = stableFirst[iP];
	GenParticleRef part(genParticles,iPlace);
	if (_printOut==-1)
	  printf("%2i -->  pdg = %7i   pt = %7.2f   eta = %5.2f   phi = %5.2f\n",
		 iP,part->pdgId(),part->pt(),part->eta(),part->phi());
	pdgIdPartFirstH1[iP] = part->pdgId();
	qPartFirstH1[iP]=0;
	if (part->charge()<-0.5)
	  qPartFirstH1[iP]=-1;
	else if (part->charge()>0.5)
	  qPartFirstH1[iP]=1;
	pxPartFirstH1[iP]=part->px();
	pyPartFirstH1[iP]=part->py();
	pzPartFirstH1[iP]=part->pz();
	ptPartFirstH1[iP]=part->pt();
	etaPartFirstH1[iP]=part->eta();
	phiPartFirstH1[iP]=part->phi();
	if ((!foundFirstMuH1) && (part->pdgId()==13||part->pdgId()==-13)) {
	  foundFirstMuH1 = true;
	  for (int jP=0; jP<nParticles; ++jP) {
	    if (nPartAroundMuFirstH1>19) break;
	    if (iPlace==jP) continue;
	    GenParticleRef partX(genParticles,jP);
	    if (partX->status()!=1) continue;
	    if (partX->pt()<1.0) continue;
	    int qPartX = 0;
	    if (partX->charge()<-0.5)
	      qPartX = -1;
	    else if (partX->charge()>0.5)
	      qPartX = 1;
	    if (qPartX == 0) continue;
	    float deltaRX = utils.deltaR(part->eta(),part->phi(),
					 partX->eta(),partX->phi());
	    if (deltaRX>0.5) continue;
	    pdgIdPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->pdgId();
	    qPartAroundMuFirstH1[nPartAroundMuFirstH1] = qPartX;
	    pxPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->px();
	    pyPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->py();
	    pzPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->pz();
	    ptPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->pt();
	    etaPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->eta();
	    phiPartAroundMuFirstH1[nPartAroundMuFirstH1] = partX->phi();
	    nPartAroundMuFirstH1++;
	  }
	}
      }
      
      if (_printOut==-1) {
	if (nPartAroundMuFirstH1>0) {
	  std::cout << "Particles around muon ---> " << std::endl;
	  for (int iP=0; iP<nPartAroundMuFirstH1; ++iP)
	    printf("%2i -->  pdg = %7i   pt = %7.2f   eta = %5.2f   phi = %5.2f\n",
		   iP,pdgIdPartAroundMuFirstH1[iP],
		   ptPartAroundMuFirstH1[iP],
		   etaPartAroundMuFirstH1[iP],
		   phiPartAroundMuFirstH1[iP]);
	}
      }
      
      if (_printOut==-1) {
	std::cout << std::endl;
	std::cout << "Second H1 ----> " << std::endl;
      }
      if (nStableSecond>19)
	nStableSecond = 19;
      nPartSecondH1 = nStableSecond;
      for (int iP=0; iP<nStableSecond; iP++) {
	int iPlace = stableSecond[iP];
	GenParticleRef part(genParticles,iPlace);
	if (_printOut==-1)
	  printf("%2i -->  pdg = %7i   pt = %7.2f   eta = %5.2f   phi = %5.2f\n",
		 iP,part->pdgId(),part->pt(),part->eta(),part->phi());
	pdgIdPartSecondH1[iP] = part->pdgId();
	qPartSecondH1[iP]=0;
	if (part->charge()<-0.5)
	  qPartSecondH1[iP]=-1;
	else if (part->charge()>0.5)
	  qPartSecondH1[iP]=1;
	pxPartSecondH1[iP]=part->px();
	pyPartSecondH1[iP]=part->py();
	pzPartSecondH1[iP]=part->pz();
	ptPartSecondH1[iP]=part->pt();
	etaPartSecondH1[iP]=part->eta();
	phiPartSecondH1[iP]=part->phi();
	if ((!foundSecondMuH1) && (part->pdgId()==13||part->pdgId()==-13)) {
	  foundSecondMuH1 = true;
	  for (int jP=0; jP<nParticles; ++jP) {
	    if (nPartAroundMuSecondH1>19) break;
	    if (iPlace==jP) continue;
	    GenParticleRef partX(genParticles,jP);
	    if (partX->status()!=1) continue;
	    if (partX->pt()<1.0) continue;
	    int qPartX = 0;
	    if (partX->charge()<-0.5)
	      qPartX = -1;
	    else if (partX->charge()>0.5)
	      qPartX = 1;
	    if (qPartX == 0) continue;
	    float deltaRX = utils.deltaR(part->eta(),part->phi(),
					 partX->eta(),partX->phi());
	    if (deltaRX>0.5) continue;
	    pdgIdPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->pdgId();
	    qPartAroundMuSecondH1[nPartAroundMuSecondH1] = qPartX;
	    pxPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->px();
	    pyPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->py();
	    pzPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->pz();
	    ptPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->pt();
	    etaPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->eta();
	    phiPartAroundMuSecondH1[nPartAroundMuSecondH1] = partX->phi();
	    nPartAroundMuSecondH1++;
	  }
	}
      }
      if (_printOut==-1) {
	if (nPartAroundMuSecondH1>0) {
	  std::cout << "Particles around muon ---> " << std::endl;
	  for (int iP=0; iP<nPartAroundMuSecondH1; ++iP)
	    printf("%2i -->  pdg = %7i   pt = %7.2f   eta = %5.2f   phi = %5.2f\n",
		   iP,pdgIdPartAroundMuSecondH1[iP],
		   ptPartAroundMuSecondH1[iP],
		   etaPartAroundMuSecondH1[iP],
		   phiPartAroundMuSecondH1[iP]);
	}
      }
      if (_printOut==-1) {
	std::cout << std::endl;
	std::cout << std::endl;
      }
    }
    _genTree->Fill();
  }

  _event++;
    
  _bosonPDG = 0;
  _decayPDG = 0;
  _origBosonPDG = 0;
  _origDecayPDG = 0;
    
  float step = 0.0;
  h_cutFlow_->Fill(step);
  
  // Trigger Information
  //ntpGeneralInfo.analyze(iEvent, iSetup);
  
  // Trigger Information
  ntpTriggerInfo.analyze(iEvent, iSetup);
  
  // Pile-Up Information
  ntpPileUpInfo.analyze(iEvent, iSetup);
  
  // Pile-Up Information
  ntpVertexInfo.analyze(iEvent, iSetup);
  reco::Vertex primaryVtx = ntpVertexInfo.getPrimaryVtx();
  std::vector<reco::Vertex> pVtxs = ntpVertexInfo.getVertices();
  
  if (ntpVertexInfo.primaryVertexFound())
    {
      step = 1.0;
      h_cutFlow_->Fill(step);
    }
  

  // get muons -->
  HTauTauGetMuons getMuon(_muonSrc);
  patMuonVector initialMuon;
  getMuon.getMuons(iEvent,initialMuon);
  int nInitialMuon = initialMuon.size();
  h_numLep_->Fill(float(nInitialMuon));

  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "Number of muons = " << nInitialMuon << std::endl;
  }

  if (nInitialMuon>1) {
    step = 2.0;
    h_cutFlow_->Fill(step);
    for (int iE=0;iE<nInitialMuon;++iE) {
      h_ptLep_->Fill(float(initialMuon[iE]->pt()));
      h_etaLep_->Fill(float(initialMuon[iE]->eta()));
      if (_printOut>0) 
	std::cout << "Muon " << iE << "  Q = " << initialMuon[iE]->charge() << " pt = " << initialMuon[iE]->pt()
		  << " eta = " << initialMuon[iE]->eta() << "  phi = " << initialMuon[iE]->phi() << std::endl;
    }
  }
  else {
    return;
  }

  /*// get orig muons (if it is set) -->
  if(!_origMuonSrc.label().empty())
  {
    edm::Handle<std::vector<reco::Muon> > muons;
    iEvent.getByLabel(_origMuonSrc, muons);

    assert(muons.isValid() && muons->size() >= 2);

    const reco::Muon* negOrigMuon = NULL;
    const reco::Muon* posOrigMuon = NULL;

    // Use hardest orig muons
    for(unsigned int i = 0; i < muons->size(); ++i)
    {
      if((*muons)[i].charge() > 0 && (!posOrigMuon || posOrigMuon->p4().pt() < (*muons)[i].p4().pt()))
        posOrigMuon = &(*muons)[i];
      if((*muons)[i].charge() < 0 && (!negOrigMuon || negOrigMuon->p4().pt() < (*muons)[i].p4().pt()))
        negOrigMuon = &(*muons)[i];
    }

    assert(negOrigMuon != NULL);
    assert(posOrigMuon != NULL);

    _origDiLepMass = (negOrigMuon->p4() + posOrigMuon->p4()).M();
    _origNegLepPt = negOrigMuon->p4().pt();
    _origPosLepPt = posOrigMuon->p4().pt();
    _origNegLepEta = negOrigMuon->p4().eta();
    _origPosLepEta = posOrigMuon->p4().eta();
    _origNegLepPhi = negOrigMuon->p4().phi();
    _origPosLepPhi = posOrigMuon->p4().phi();
    _origNegLepChargedHadronIso = negOrigMuon->pfIsolationR04().sumChargedHadronPt;
    _origPosLepChargedHadronIso = negOrigMuon->pfIsolationR04().sumChargedHadronPt;
  }*/

  // Get weight for embedded events
  edm::Handle<GenFilterInfo> hGenFilterInfo;
  iEvent.getByLabel(edm::InputTag("generator", "minVisPtFilter", "EmbeddedRECO"), hGenFilterInfo);
  if(hGenFilterInfo.isValid())
  {
    _embedWeight = hGenFilterInfo->filterEfficiency();
  }
  else
  {
    // Older samples have a simple double value
    edm::Handle<double> hUserWeight;
    iEvent.getByLabel(edm::InputTag("generator", "weight", "EmbeddedRECO"), hUserWeight);
    if(hUserWeight.isValid()) _embedWeight = *hUserWeight;
    else _embedWeight = 1.0f; // probably not an embedded sample
  }

  // Muon selector -->
  bool sameSign = false;
  
  if (_sameSignMuons) {
    sameSign = true;
  }
  
  HTauTauMuonSelector muonSelector(_ptMin, _etaMax, _ptMinHard, sameSign);
  patMuonVector selectedMuon;
  muonSelector.selectMuons(initialMuon,selectedMuon);
  int nSelMuon = selectedMuon.size();

  if (_printOut>0) {
    std::cout << "Number of selected muons = " << nSelMuon << std::endl;
  }

  if (nSelMuon>1) {
    step = 3.0;
    h_cutFlow_->Fill(step);
    for (int iE=0;iE<nSelMuon;++iE) {
      pat::MuonRef muonRef = selectedMuon[iE];
      double TrackIso = muonRef->ecalIso ();
      double EcalIso = muonRef->hcalIso();
      double HcalIso = muonRef->trackIso();
      float isolation = float((TrackIso+EcalIso+HcalIso)/(fmax(0.1,muonRef->pt())));
      h_isoLep_->Fill(isolation);
      if (_printOut>0) {
	std::cout << "Muon " << iE << " Isolation = " << isolation << std::endl;
      }
    }
  }
  else {
    return;
  }
   
  // Muon isolator -->
  HTauTauMuonIsolator muonIsolator(_isoCut);
  patMuonVector isolatedMuon;
  bool mode = true;
  muonIsolator.isolateMuonsRel(mode,selectedMuon,isolatedMuon);
  int nIsoMuon = isolatedMuon.size();
  if (_printOut>0) {
    std::cout << "Number of isolated muons = " << nIsoMuon << std::endl;
  }
  if (nIsoMuon>1) {
    step = 4.0;
    h_cutFlow_->Fill(step);
  }
  else {
    return;
  }
  
  // select pair -->
  patMuonVector selectedPair;
  double ptHard = muonSelector.selectPairMuons(isolatedMuon,selectedPair);
  h_ptHard_->Fill(float(ptHard));
  int nPair = selectedPair.size();
  if (_printOut>0) {
    std::cout << " Selected pair = " << nPair << " ptHard = " << ptHard <<  std::endl;
  }
  if (ptHard>_ptMinHard&&nPair==2) {
    step = 5.0;
    h_cutFlow_->Fill(step);
  }
  else {
    return;
  }
  
  pat::MuonRef negMuon = selectedPair.at(0);
  pat::MuonRef posMuon = selectedPair.at(1);

  Handle<reco::PFCandidateCollection> NoPileUp;
  Handle<reco::PFCandidateCollection> pflow;
  Handle<reco::PFCandidateCollection> pfPileUp;
 

  if (iEvent.getByLabel(_pfParticleSrc,pflow)) {

  }
  else {
    std::cout << "Collection " << _pfParticleSrc << "  not found " << std::endl; 
  }

  // 4-vectors;
  TLorentzVector NegLorentz(negMuon->px(),negMuon->py(),negMuon->pz(),negMuon->energy());
  TLorentzVector PosLorentz(posMuon->px(),posMuon->py(),posMuon->pz(),posMuon->energy());

  // run trigger matching

  ntpTriggerInfo.analyzeObject(posMuon->p4(), "PosMuon", iEvent, iSetup, 0.3);
  ntpTriggerInfo.analyzeObject(negMuon->p4(), "NegMuon", iEvent, iSetup, 0.3);

//   std::cout << std::endl;
//   std::cout << "NPart H1(1) = " << nPartFirstH1 << std::endl;
//   for (int iP=0; iP<nPartFirstH1; iP++) {
//     if (qPartFirstH1[iP]==-1 || qPartFirstH1[iP]==1) {
//       std::cout << "  q = " << qPartFirstH1[iP]  
// 		<< "  Id = " << pdgIdPartFirstH1[iP]
// 		<< "  pt = " << ptPartFirstH1[iP] 
// 		<< "  eta = " << etaPartFirstH1[iP] 
// 		<< "  phi = " << phiPartFirstH1[iP] << std::endl;
//     }
//   }
//   std::cout << std::endl;
//   std::cout << "NPart H1(2) = " << nPartSecondH1 << std::endl;
//   for (int iP=0; iP<nPartSecondH1; iP++) {
//     if (qPartSecondH1[iP]==-1 || qPartSecondH1[iP]==1) {
//       std::cout << "  q = " << qPartSecondH1[iP]  
// 		<< "  Id = " << pdgIdPartSecondH1[iP]
// 		<< "  pt = " << ptPartSecondH1[iP] 
// 		<< "  eta = " << etaPartSecondH1[iP] 
// 		<< "  phi = " << phiPartSecondH1[iP] << std::endl;
//     }
//   }

//   std::cout << std::endl;
//   std::cout << "Muons ---->" << std::endl;
//   std::cout << "Muon1 ->   pt = " << negMuon->pt() << "   q = " << negMuon->charge() << std::endl;
//   std::cout << "Muon2 ->   pt = " << posMuon->pt() << "   q = " << posMuon->charge() << std::endl;
//   std::cout << std::endl;

  // *************************************************************************
  // * Collection of Vertices : vertex with highest trk pt2 sum is chosen ====>
  // *************************************************************************
  Handle< std::vector<reco::Vertex> > vertexCol;
  reco::Vertex HighPTVtx;
  std::vector<reco::Vertex> vertices;
  vertices.clear();

  if (iEvent.getByLabel(_vtxSrc,vertexCol) && vertexCol.isValid())
    {
      int numberOfVertices = vertexCol->size();
      double ptmax = 0;
      
      for (int iVertex = 0; iVertex< std::min(100,numberOfVertices); iVertex++)
	{
	  reco::Vertex tmpVertex = vertexCol->at(iVertex);
	  vertices.push_back(tmpVertex);
	  double ptsum = 0;
	  for(reco::Vertex::trackRef_iterator iTrack  = tmpVertex.tracks_begin(); iTrack != tmpVertex.tracks_end();++iTrack) {
	    double pt=(*iTrack)->pt();
	    ptsum += (pt*pt);
	  }
		

	  if(ptsum > ptmax)
	    {
	      // primary vertex?
	      double chi2 = tmpVertex.chi2();
	      int ndf = int(tmpVertex.ndof());
	      double prob = TMath::Prob(chi2,ndf);
	      double xVtx = tmpVertex.x();
	      double yVtx = tmpVertex.y();
	      double d = sqrt(xVtx*xVtx+yVtx*yVtx);
	      if (fabs(tmpVertex.z())<_zVtxCut && d < _dVtxCut && prob > _probVtxCut)
		{
		  ptmax = ptsum;
		  HighPTVtx = tmpVertex;
		  
		}
	    }
	}
    }
  // **************************
  // No PileUp PFlow Collection
  // **************************
  npfNoPU = 0;
  ntracksPos = 0;
  ntracksNeg = 0;

  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);
  
  const TransientTrackBuilder * transientTrackBuilder = builder.product();
  HTauTauImpactParameter ipTool;

  std::cout << std::endl;

  if (iEvent.getByLabel(_pfNoPileUpSrc,NoPileUp))
    {
      int nallpfNoPU = NoPileUp->size();
      if (_printOut > 0 )
	std::cout << "   N(ALL NoPileUp) = " << nallpfNoPU << "\n";
      //      int nSelNoPUcand = 0;
      for (int iP = 0; iP < nallpfNoPU; ++iP)
	{
	  reco::PFCandidateRef pfCand(NoPileUp, iP);
	  float ptPfCand = pfCand->pt();
	  float chargePfCand = pfCand->charge();

	  float ipPfCand = 0;
	  float ipSigPfCand = 0;
	  float ipZPfCand = 0;
	  float ipZSigPfCand = 0;
	  float ip3DPfCand = 0;
	  float ip3DSigPfCand = 0;

	  if (chargePfCand<-0.5||chargePfCand>0.5) {
	    reco::TrackRef trkRef = pfCand->trackRef();
	    if (!trkRef.isNull()) {

	      //	      std::cout << "Particle " << iP << "  track quality = " << trkRef->quality(Track::highPurity) << std::endl;

	      TransientTrack transientTrk = transientTrackBuilder->build(trkRef);
	      IPMeasurement ipMeas = ipTool.impactParameter(transientTrk,HighPTVtx);
		  
	      Measurement1D trkIP   = ipMeas.ip;
	      Measurement1D trkIPZ  = ipMeas.ipZ;
	      Measurement1D trkIP3D = ipMeas.ip3D;
	      
	      ipPfCand = trkIP.value();
	      ipSigPfCand = trkIP.significance();
	      
	      ipZPfCand = trkIPZ.value();
	      ipZSigPfCand = trkIPZ.significance();
	      
	      ip3DPfCand = trkIP3D.value();
	      ip3DSigPfCand = trkIP3D.significance();
	    }
	  }

// 	  if (chargePfCand<-0.5||chargePfCand>0.5) {
// 	    if (ptPfCand>1.0) {
// 	      std::cout << "pfNoPU : " << nSelNoPUcand 
// 			<< " ---> pt = " << ptPfCand 
// 			<< "  eta = " << pfCand->eta()
// 			<< "  phi = " << pfCand->phi()
// 			<< "  charge = " <<  pfCand->charge() 
// 			<< "  ip = " << ipPfCand 
// 			<< "  ipZ = " << ipZPfCand 
// 			<< std::endl; 
// 	      nSelNoPUcand++;
// 	    }
// 	  }	  

	  if (ptPfCand>0.5 && npfNoPU<10000){
	   
	    pxpfNoPU[npfNoPU]  = pfCand->px();
	    pypfNoPU[npfNoPU]  = pfCand->py();
	    pzpfNoPU[npfNoPU]  = pfCand->pz();
	    ptpfNoPU[npfNoPU]  = pfCand->pt();
	    etapfNoPU[npfNoPU] = pfCand->eta();
	    enpfNoPU[npfNoPU]  = pfCand->energy();
	    chargepfNoPU[npfNoPU] = pfCand->charge();
	    ippfNoPU[npfNoPU] = ipPfCand;
	    ipSigpfNoPU[npfNoPU] = ipSigPfCand;
	    ipZpfNoPU[npfNoPU] = ipZPfCand;
	    ipZSigpfNoPU[npfNoPU] = ipZSigPfCand;
	    ip3DpfNoPU[npfNoPU] = ip3DPfCand;
	    ip3DSigpfNoPU[npfNoPU] = ip3DSigPfCand;
	    idpfNoPU[npfNoPU] = pfCand->pdgId();
	    npfNoPU ++;

	  }
	 
	  if (chargePfCand<-0.5||chargePfCand>0.5) {
	    if (ptPfCand>1.) {

	      TLorentzVector pfCandLorentz(pfCand->px(),pfCand->py(),pfCand->pz(),pfCand->energy());
	      
	      // count tracks around pos muon
	      float deltaRPos = utils.deltaR(posMuon->eta(),posMuon->phi(),
					     pfCand->eta(),pfCand->phi());
	      TLorentzVector diffPosLorentz = PosLorentz - pfCandLorentz; 
	      if (deltaRPos<0.5 && diffPosLorentz.P()>0.1 && ntracksPos<1000) {
		pxTrackPos[ntracksPos] = pfCand->px();
		pyTrackPos[ntracksPos] = pfCand->py();
		pzTrackPos[ntracksPos] = pfCand->pz();
		ptTrackPos[ntracksPos] = pfCand->pt();
		etaTrackPos[ntracksPos] = pfCand->eta();
		phiTrackPos[ntracksPos] = pfCand->phi();
		enTrackPos[ntracksPos] = pfCand->energy();
		qTrackPos[ntracksPos]  = pfCand->charge();
		idTrackPos[ntracksPos] = pfCand->pdgId();
		ipTrackPos[ntracksPos] = ipPfCand;
		ipSigTrackPos[ntracksPos] = ipSigPfCand;	
		ipZTrackPos[ntracksPos] = ipZPfCand;
		ipZSigTrackPos[ntracksPos] = ipZSigPfCand;	
		ip3DTrackPos[ntracksPos] = ip3DPfCand;
		ip3DSigTrackPos[ntracksPos] = ip3DSigPfCand;	
		ntracksPos++;
	      }

	      // count tracks around neg muon
	      float deltaRNeg = utils.deltaR(negMuon->eta(),negMuon->phi(),
					     pfCand->eta(),pfCand->phi());
	      TLorentzVector diffNegLorentz = NegLorentz - pfCandLorentz; 
	      if (deltaRNeg<0.5 && diffNegLorentz.P()>0.1 && ntracksNeg<1000) {
		pxTrackNeg[ntracksNeg] = pfCand->px();
		pyTrackNeg[ntracksNeg] = pfCand->py();
		pzTrackNeg[ntracksNeg] = pfCand->pz();
		ptTrackNeg[ntracksNeg] = pfCand->pt();
		etaTrackNeg[ntracksNeg] = pfCand->eta();
		phiTrackNeg[ntracksNeg] = pfCand->phi();
		enTrackNeg[ntracksNeg] = pfCand->energy();
		qTrackNeg[ntracksNeg]  = pfCand->charge();
		idTrackNeg[ntracksNeg] = pfCand->pdgId();
		ipTrackNeg[ntracksPos] = ipPfCand;
		ipSigTrackNeg[ntracksPos] = ipSigPfCand;	
		ipZTrackNeg[ntracksPos] = ipZPfCand;
		ipZSigTrackNeg[ntracksPos] = ipZSigPfCand;	
		ip3DTrackNeg[ntracksPos] = ip3DPfCand;
		ip3DSigTrackNeg[ntracksPos] = ip3DSigPfCand;
		ntracksNeg++;
	      }

	    }
	   
	  }

	}
      //      std::cout << "Selected pfNoPUcand = " << nSelNoPUcand << std::endl;
    }
  else
    {
      std::cout << "NoPileUp not found: " << "\n";
    }
 
  // ***********************
  // PileUp PFlow Collection
  // ***********************

  ntracksPileUpPos = 0;
  ntracksPileUpNeg = 0;

  if (iEvent.getByLabel(_pfPileUpSrc,pfPileUp))
    {
      int NpfPileUp = pfPileUp->size();
      if (_printOut > 0 )
	std::cout << "   N(pfPileUp) = " << NpfPileUp<< "\n";
      for (int iP = 0; iP < NpfPileUp; ++iP)
        {

	  reco::PFCandidateRef pfCand(pfPileUp, iP);
          float ptPfCand = pfCand->pt();
          float chargePfCand = pfCand->charge();

          if (chargePfCand<-0.5||chargePfCand>0.5) {
            if (ptPfCand>1.0) {

              TLorentzVector pfCandLorentz(pfCand->px(),pfCand->py(),pfCand->pz(),pfCand->energy());

              // count tracks around pos muon                                                                                                                                                                                               
              float deltaRPos = utils.deltaR(posMuon->eta(),posMuon->phi(),
                                             pfCand->eta(),pfCand->phi());
              TLorentzVector diffPosLorentz = PosLorentz - pfCandLorentz;
              if (deltaRPos<0.5 && diffPosLorentz.P()>0.1 && ntracksPileUpPos<1000) {
                pxTrackPileUpPos[ntracksPileUpPos] = pfCand->px();
                pyTrackPileUpPos[ntracksPileUpPos] = pfCand->py();
                pzTrackPileUpPos[ntracksPileUpPos] = pfCand->pz();
                ptTrackPileUpPos[ntracksPileUpPos] = pfCand->pt();
                etaTrackPileUpPos[ntracksPileUpPos] = pfCand->eta();
                phiTrackPileUpPos[ntracksPileUpPos] = pfCand->phi();
                enTrackPileUpPos[ntracksPileUpPos] = pfCand->energy();
                qTrackPileUpPos[ntracksPileUpPos]  = pfCand->charge();
                idTrackPileUpPos[ntracksPileUpPos] = pfCand->pdgId();
                ntracksPileUpPos++;
              }

              // count tracks around neg muon                                                                                                                                                                                               
              float deltaRNeg = utils.deltaR(negMuon->eta(),negMuon->phi(),
                                             pfCand->eta(),pfCand->phi());
              TLorentzVector diffNegLorentz = NegLorentz - pfCandLorentz;
              if (deltaRNeg<0.5 && diffNegLorentz.P()>0.1 && ntracksPileUpNeg<1000) {
                pxTrackPileUpNeg[ntracksPileUpNeg] = pfCand->px();
                pyTrackPileUpNeg[ntracksPileUpNeg] = pfCand->py();
                pzTrackPileUpNeg[ntracksPileUpNeg] = pfCand->pz();
                ptTrackPileUpNeg[ntracksPileUpNeg] = pfCand->pt();
                etaTrackPileUpNeg[ntracksPileUpNeg] = pfCand->eta();
                phiTrackPileUpNeg[ntracksPileUpNeg] = pfCand->phi();
                enTrackPileUpNeg[ntracksPileUpNeg] = pfCand->energy();
                qTrackPileUpNeg[ntracksPileUpNeg]  = pfCand->charge();
                idTrackPileUpNeg[ntracksPileUpNeg] = pfCand->pdgId();
                ntracksPileUpNeg++;
              }

            }

          }

        }
    }
  else
    {
      std::cout << "PileUp not found: " << _pfPileUpSrc << "\n";
    }

  if (_printOut==-1) {

    printf("1st Muon =>  pt=%7.2f  eta=%5.2f  phi=%5.2f\n",posMuon->pt(),posMuon->eta(),posMuon->phi());
    std::cout << "No PU particles ----> " << std::endl;
    for (int iTrk=0; iTrk<ntracksPos;++iTrk)
      printf("    pt=%7.2f  eta=%5.2f  phi=%5.2f  Q=%4.1f\n",ptTrackPos[iTrk],etaTrackPos[iTrk],phiTrackPos[iTrk],qTrackPos[iTrk]);
    std::cout << "PU particles ----> " << std::endl;
    for (int iTrk=0; iTrk<ntracksPileUpPos;++iTrk)
      printf("    pt=%7.2f  eta=%5.2f  phi=%5.2f  Q=%4.1f\n",
	     ptTrackPileUpPos[iTrk],
	     etaTrackPileUpPos[iTrk],
	     phiTrackPileUpPos[iTrk],
	     qTrackPileUpPos[iTrk]);
    
    std::cout << std::endl;
    printf("2nd Muon =>  pt=%7.2f  eta=%5.2f  phi=%5.2f\n",negMuon->pt(),negMuon->eta(),negMuon->phi());
    std::cout << "No PU particles ----> " << std::endl;
    for (int iTrk=0; iTrk<ntracksNeg;++iTrk)
      printf("    pt=%7.2f  eta=%5.2f  phi=%5.2f  Q=%4.1f\n",ptTrackNeg[iTrk],etaTrackNeg[iTrk],phiTrackNeg[iTrk],qTrackNeg[iTrk]);
    std::cout << "PU particles ----> " << std::endl;
    for (int iTrk=0; iTrk<ntracksPileUpNeg;++iTrk)
      printf("    pt=%7.2f  eta=%5.2f  phi=%5.2f  Q=%4.1f\n",
	     ptTrackPileUpNeg[iTrk],
	     etaTrackPileUpNeg[iTrk],
	     phiTrackPileUpNeg[iTrk],
	     qTrackPileUpNeg[iTrk]);
    
    
    std::cout << std::endl;

  }

  //  if (ntracksPos>3)
  //    return;

  //  if (ntracksPos<1)
  //    return;

  //  if (ntracksNeg>3)
  //    return;

  //  if (ntracksNeg<1)
  //    return;


  if (_printOut>0) {
    std::cout << "Pos Muon  Q = " << posMuon->charge() <<  "  (Px,Py,Pz) = (" 
	      << posMuon->px() << "," << posMuon->py() << "," << posMuon->pz() << ")" << std::endl;
    std::cout << "NTracks(Pos) = " << ntracksPos << std::endl;
    for (int iTrk=0; iTrk<ntracksPos; iTrk++) {
      std::cout << "  " << iTrk << " q = " 
		<< qTrackPos[iTrk] 
		<< "  id = " << idTrackPos[iTrk]
		<< "  (Px,Py,Pz) = (" 
		<< pxTrackPos[iTrk] << ","
		<< pyTrackPos[iTrk] << "," 
		<< pzTrackPos[iTrk] << ")" << std::endl;
      
    }
    
    std::cout << std::endl;
    std::cout << "Neg Muon  Q = " << negMuon->charge() <<  "  (Px,Py,Pz) = (" 
	      << negMuon->px() << "," << negMuon->py() << "," << posMuon->pz() << ")" << std::endl;

    std::cout << "NTracks(Neg) = " << ntracksNeg << std::endl;
    for (int iTrk=0; iTrk<ntracksNeg; iTrk++) {
      std::cout << "  " << iTrk << " q = " 
		<< qTrackNeg[iTrk] 
		<< "  id = " << idTrackNeg[iTrk]
		<< "  (Px,Py,Pz) = (" 
		<< pxTrackNeg[iTrk] << ","
		<< pyTrackNeg[iTrk] << "," 
		<< pzTrackNeg[iTrk] << ")" << std::endl;
      
    }
    std::cout << std::endl;
  }


  _posLep_Mu20     = false;
  _posLep_Mu30     = false;
  _posLep_Mu40     = false;
  _posLep_IsoMu24  = false;
  _posLep_Mu17_Mu8 = false;
  _posLep_Mu17_TkMu8 = false;
  

  _negLep_Mu20     = false;
  _negLep_Mu30     = false;
  _negLep_Mu40     = false;
  _negLep_IsoMu24  = false;
  _negLep_Mu17_Mu8 = false;
  _negLep_Mu17_TkMu8 = false; 

  // positive muon

  const pat::TriggerObjectStandAlone * matchedTrigObjPos = posMuon->triggerObjectMatch();

  if (matchedTrigObjPos==NULL) {
    // empty line
  }
  else {
    std::vector<std::string> pathNames = matchedTrigObjPos->pathNames();
    int nPaths = int(pathNames.size());
    if (_printOut>0)
      std::cout << "Matched triggers (Pos) : " << std::endl;
      
    for (int iP=0; iP<nPaths; ++iP) {
	  
      TString PathName(pathNames[iP]);
      
      if (_printOut>0)
	std::cout << "     " << PathName << std::endl;
      
      // setup triggers --->
      if (PathName.Contains("HLT_Mu20_v"))
	_posLep_Mu20 = true;

      if (PathName.Contains("HLT_Mu30_v"))
	_posLep_Mu30 = true;

      if (PathName.Contains("HLT_Mu40_v"))
	_posLep_Mu40 = true;

      if (PathName.Contains("HLT_IsoMu24_v"))
	_posLep_IsoMu24 = true;

      if (PathName.Contains("HLT_Mu17_Mu8_v"))
	_posLep_Mu17_Mu8 = true;

      if (PathName.Contains("HLT_Mu17_TkMu8_v"))
	_posLep_Mu17_TkMu8 = true;


    }
  }  

  if (_printOut>0)
    std::cout << std::endl;

  // negative muon --->

  const pat::TriggerObjectStandAlone * matchedTrigObjNeg = negMuon->triggerObjectMatch();

  if (matchedTrigObjNeg==NULL) {
    // empty line
  }
  else {
    std::vector<std::string> pathNames = matchedTrigObjNeg->pathNames();
    int nPaths = int(pathNames.size());
    if (_printOut>0)
      std::cout << "Matched triggers (Neg) : " << std::endl;
      
    for (int iP=0; iP<nPaths; ++iP) {
	  
      TString PathName(pathNames[iP]);
      
      if (_printOut>0)
	std::cout << "     " << PathName << std::endl;
      
      // setup triggers --->
      if (PathName.Contains("HLT_Mu20_v"))
	_negLep_Mu20 = true;

      if (PathName.Contains("HLT_Mu30_v"))
	_negLep_Mu30 = true;

      if (PathName.Contains("HLT_Mu40_v"))
	_negLep_Mu40 = true;

      if (PathName.Contains("HLT_IsoMu24_v"))
	_negLep_IsoMu24 = true;

      if (PathName.Contains("HLT_Mu17_Mu8_v"))
	_negLep_Mu17_Mu8 = true;

      if (PathName.Contains("HLT_Mu17_TkMu8_v"))
	_negLep_Mu17_TkMu8 = true;


    }
  }  

  if (_printOut>0)
    std::cout << std::endl;


  // Tracks from hard interaction vertex --->
  nPVtracks=0;
  //  int nSelVtxTrk = 0;
  //  std::cout << std::endl;
  for(reco::Vertex::trackRef_iterator iTrack  = HighPTVtx.tracks_begin(); iTrack != HighPTVtx.tracks_end();++iTrack) {
    if((*iTrack)->pt()>0.5 && nPVtracks<1000){

      ptTrack[nPVtracks]=(*iTrack)->pt();
      pxTrack[nPVtracks]=(*iTrack)->px();
      pyTrack[nPVtracks]=(*iTrack)->py();
      pzTrack[nPVtracks]=(*iTrack)->pz();
      etaTrack[nPVtracks]=(*iTrack)->eta();
      phiTrack[nPVtracks]=(*iTrack)->phi();
      chargeTrack[nPVtracks]=(*iTrack)->charge();
      chi2Track[nPVtracks]=(*iTrack)->chi2();
//       if (ptTrack[nPVtracks]>1.0) {
// 	if (chargeTrack[nPVtracks]<-0.5||chargeTrack[nPVtracks]>0.5) {
// 	  std::cout << "Trk : " << nSelVtxTrk 
// 		    << " --> pt = " << ptTrack[nPVtracks] 
// 		    << "  eta = " << etaTrack[nPVtracks]  
// 		    << "  phi = " << phiTrack[nPVtracks] 
// 		    << "   q = " <<  chargeTrack[nPVtracks] << std::endl;
// 	  nSelVtxTrk++;
// 	}
//       }
      nPVtracks++;
    }
  }

  // std::cout << "Number of selected vtx tracks = " << nSelVtxTrk << std::endl;
  // std::cout << std::endl;

 nPVPostracks=0;
 nPVNegtracks=0;
 float pionMass = 0.1396;
 for (int itrack=0; itrack<nPVtracks; itrack++)
   {

     if (chi2Track[itrack]!=0&&ptTrack[itrack]>1.0) {

       // count tracks around pos muon                                                                                                                                                                                               
       float deltaRPos = utils.deltaR(posMuon->eta(),posMuon->phi(),
				      etaTrack[itrack],phiTrack[itrack]);
   
       TLorentzVector trkLorentz;
       trkLorentz.SetXYZM(pxTrack[itrack],
			  pyTrack[itrack],
			  pzTrack[itrack],
			  pionMass);
          
       TLorentzVector diffPosLorentz = PosLorentz - trkLorentz;
       if ( deltaRPos<0.5 && diffPosLorentz.P()>0.1 && nPVPostracks<1000 ) {
	 ptPVTracksPos[nPVPostracks]=ptTrack[itrack];
	 phiPVTracksPos[nPVPostracks]=phiTrack[itrack];
	 etaPVTracksPos[nPVPostracks]=etaTrack[itrack];
	 pxPVTracksPos[nPVPostracks]=pxTrack[itrack];
	 pyPVTracksPos[nPVPostracks]=pyTrack[itrack];
	 pzPVTracksPos[nPVPostracks]=pzTrack[itrack];
	 qPVTracksPos[nPVPostracks]=chargeTrack[itrack];
	 nPVPostracks++;
       }

       // count tracks around neg muon                                                                                                                                                                                               
            
       float deltaRNeg = utils.deltaR(negMuon->eta(),negMuon->phi(),
				      etaTrack[itrack],phiTrack[itrack]);
       TLorentzVector diffNegLorentz = NegLorentz - trkLorentz;
       if ( deltaRNeg<0.5 && diffNegLorentz.P() > 0.1 && nPVNegtracks<1000 ) {
	 ptPVTracksNeg[nPVNegtracks]=ptTrack[itrack];
	 phiPVTracksNeg[nPVNegtracks]=phiTrack[itrack];
	 etaPVTracksNeg[nPVNegtracks]=etaTrack[itrack];
	 pxPVTracksNeg[nPVNegtracks]=pxTrack[itrack];
	 pyPVTracksNeg[nPVNegtracks]=pyTrack[itrack];
	 pzPVTracksNeg[nPVNegtracks]=pzTrack[itrack];
	 qPVTracksNeg[nPVNegtracks]=chargeTrack[itrack];
	 nPVNegtracks++; 
       }
       
     }
     
   }


  // ********************
  // accessing jets ====>
  // ********************
  edm::Handle <pat::JetCollection> jetHandle;
  iEvent.getByLabel( _jetSrc, jetHandle );
  int numberOfJets = jetHandle->size();
  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "Number of Jets = " << numberOfJets << std::endl;
  }

  _nJets = 0;
  for (int iJ=0;iJ<numberOfJets;++iJ) {
    pat::JetRef thisJet( jetHandle,iJ );
    double JetEt = thisJet->et();
    h_jetEt_->Fill(float(JetEt));
    if (_printOut>0)
      std::cout << "Jet " << iJ << " Et = " << JetEt << std::endl;
    if (JetEt>_jetEnergyMinCut)
      _nJets++;
  }

  if (_printOut>0)
    std::cout << "Number of selected jets = " << _nJets << std::endl;
  if (_nJets<=_nOfJetsCut) {
    step = 6;
    h_cutFlow_->Fill(step);
  }
  else {
    return;
  }
  
  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "Jet info ---> " << std::endl;
  }

  edm::ParameterSet paramsLoose;
  edm::ParameterSet paramsMedium;
  edm::ParameterSet paramsTight;

  paramsLoose.addParameter(std::string("version"),std::string("FIRSTDATA"));
  paramsLoose.addParameter(std::string("quality"),std::string("LOOSE"));

  paramsMedium.addParameter(std::string("version"),std::string("FIRSTDATA"));
  paramsMedium.addParameter(std::string("quality"),std::string("MEDIUM"));

  paramsTight.addParameter(std::string("version"),std::string("FIRSTDATA"));
  paramsTight.addParameter(std::string("quality"),std::string("TIGHT"));

  PFJetIDSelectionFunctor jetIDLooseFunctor( paramsLoose );
  PFJetIDSelectionFunctor jetIDMediumFunctor( paramsMedium );
  PFJetIDSelectionFunctor jetIDTightFunctor( paramsTight );

  pat::strbitset ret_loose = jetIDLooseFunctor.getBitTemplate();
  pat::strbitset ret_medium = jetIDMediumFunctor.getBitTemplate();
  pat::strbitset ret_tight = jetIDTightFunctor.getBitTemplate();

  /*Handle<ValueMap<StoredPileupJetIdentifier> > puJetId;
  iEvent.getByLabel(_puJetIdSrc, puJetId);

  Handle<ValueMap<float> > puJetMVAFull;
  iEvent.getByLabel(edm::InputTag(_puJetMvaSrc.label(), "fullDiscriminant"), puJetMVAFull);*/

  Handle<ValueMap<int> > puJetIdFull;
  //iEvent.getByLabel(edm::InputTag(_puJetMvaSrc.label(), "fullId"), puJetIdFull);

  noJets=0;
  nSv=0;
  for (int iJ=0;iJ<numberOfJets;++iJ) {    
    pat::JetRef jet (jetHandle,iJ);
    if (jet->et()>_jetEnergyMinCut) {
      if (noJets<50) {
	jetPt[noJets]=jet->pt();
	jetEt[noJets]=jet->et();
	jetEta[noJets]=jet->eta();
	jetPhi[noJets]=jet->phi();
	jetPx[noJets]=jet->px();
	jetPy[noJets]=jet->py();
	jetPz[noJets]=jet->pz();
	jetEn[noJets]=jet->energy();
	jetMass[noJets]=jet->mass();
	jetArea[noJets]=jet->jetArea();
	jetCorrFactor[noJets]= jet->jecFactor(0);

	//	std::cout << "Jet : " << iJ << std::endl;

	if (_printOut>1)
	{
		if (iJ==0)
		{
			std::cout << "Levels:\n";
			for (std::vector<std::string>::const_iterator itStr = jet->availableJECLevels().begin(); itStr != jet->availableJECLevels().end(); ++itStr)
				std::cout << "\t" << *itStr << "\n";
		}
		std::cout << "Jet:" << jet->pt() << "\n";
		for (size_t idxJEC = 0; idxJEC < jet->availableJECLevels().size(); ++idxJEC)
			std::cout << "\t" << idxJEC << " " << jet->jecFactor(idxJEC) << "\t" << jet->correctedP4(idxJEC).pt() << "\n";
	}

	ret_loose.set(false);
	jetIDLoose[noJets] = jetIDLooseFunctor(*jet, ret_loose);

	ret_medium.set(false);
	jetIDMedium[noJets] = jetIDMediumFunctor(*jet, ret_medium);

	ret_tight.set(false);
	jetIDTight[noJets] = jetIDTightFunctor(*jet, ret_tight);

	noConstituents[noJets]=int(jet->getPFConstituents().size());
	neuEMEnergyFr[noJets]=jet->neutralEmEnergyFraction();
	neuHadronEnergyFr[noJets]=jet->neutralHadronEnergyFraction();
	neuHadronMultiplicity[noJets]=jet->neutralHadronMultiplicity();
	neuMultiplicity[noJets]=jet->neutralMultiplicity();
	photonEnergyFr[noJets]=jet->photonEnergyFraction();
	photonMultiplicity[noJets]=jet->photonMultiplicity();
	hfEMEnergyFr[noJets]=jet->HFEMEnergyFraction();
	hfEMMultiplicity[noJets]=jet->HFEMMultiplicity();
	hfHadronEnergyFr[noJets]=jet->HFHadronEnergyFraction();
	hfHadronMultiplicity[noJets]=jet->HFHadronMultiplicity();
	noChargedConstituents[noJets]=int(jet->chargedMultiplicity());
	chargedEnergyFr[noJets]=jet->chargedEmEnergyFraction()+jet->chargedHadronEnergyFraction();
	chargedMultiplicity[noJets]=jet->chargedMultiplicity();
	electronEnergyFr[noJets]=jet->electronEnergy()/jet->energy();
	electronMultiplicity[noJets]=jet->electronMultiplicity();
	muonEnergyFr[noJets]=jet->muonEnergy()/jet->energy();
	muonMultiplicity[noJets]=jet->muonMultiplicity();

	/*puJetInputJetEta[noJets] = (*puJetId)[jet->originalObjectRef()].jetEta();
	puJetInputJetPt[noJets] = (*puJetId)[jet->originalObjectRef()].jetPt();
	puJetInputJetNCharged[noJets] = (*puJetId)[jet->originalObjectRef()].nCharged();
	puJetInputJetNNeutrals[noJets] = (*puJetId)[jet->originalObjectRef()].nNeutrals();
	puJetInputJetDZ[noJets] = (*puJetId)[jet->originalObjectRef()].dZ();
	puJetInputJetNParticles[noJets] = (*puJetId)[jet->originalObjectRef()].nParticles();
	puJetInputJetDR2Mean[noJets] = (*puJetId)[jet->originalObjectRef()].dR2Mean();
	puJetInputJetDRMean[noJets] = (*puJetId)[jet->originalObjectRef()].dRMean();
	puJetInputJetFrac01[noJets] = (*puJetId)[jet->originalObjectRef()].frac01();
	puJetInputJetFrac02[noJets] = (*puJetId)[jet->originalObjectRef()].frac02();
	puJetInputJetFrac03[noJets] = (*puJetId)[jet->originalObjectRef()].frac03();
	puJetInputJetFrac04[noJets] = (*puJetId)[jet->originalObjectRef()].frac04();
	puJetInputJetFrac05[noJets] = (*puJetId)[jet->originalObjectRef()].frac05();
	puJetInputJetFrac06[noJets] = (*puJetId)[jet->originalObjectRef()].frac06();
	puJetInputJetFrac07[noJets] = (*puJetId)[jet->originalObjectRef()].frac07();
	puJetInputJetBeta[noJets] = (*puJetId)[jet->originalObjectRef()].beta();
	puJetInputJetBetaStar[noJets] = (*puJetId)[jet->originalObjectRef()].betaStar();
	puJetInputJetBetaClassic[noJets] = (*puJetId)[jet->originalObjectRef()].betaClassic();
	puJetInputJetBetaStarClassic[noJets] = (*puJetId)[jet->originalObjectRef()].betaStarClassic();
	puJetInputJetPtD[noJets] = (*puJetId)[jet->originalObjectRef()].ptD();
	puJetInputJetNvtx[noJets] = (*puJetId)[jet->originalObjectRef()].nvtx();

	//puJetFullMVA[noJets] = (*puJetMVAFull)[jet->originalObjectRef()];
	puJetFullLoose[noJets] = PileupJetIdentifier::passJetId((*puJetIdFull)[jet->originalObjectRef()], PileupJetIdentifier::kLoose);
	puJetFullMedium[noJets] = PileupJetIdentifier::passJetId((*puJetIdFull)[jet->originalObjectRef()], PileupJetIdentifier::kMedium);
	puJetFullTight[noJets] = PileupJetIdentifier::passJetId((*puJetIdFull)[jet->originalObjectRef()], PileupJetIdentifier::kTight);*/


	// checking generator level info -->
	if (_doMC) {

	  const reco::GenJet * genJet = jet->genJet();
	  if (genJet == NULL) {
	    genJetExist[noJets] = false;	    
	    genJetPx[noJets] = 0.;
	    genJetPy[noJets] = 0.;
	    genJetPz[noJets] = 0.;
	    genJetEn[noJets] = 0.;
	    genJetCollisionId[noJets] = -1.;
	  }
	  else {
	    genJetExist[noJets] = true;	    
	    genJetPx[noJets] = genJet->px();
	    genJetPy[noJets] = genJet->py();
	    genJetPz[noJets] = genJet->pz();
	    genJetEn[noJets] = genJet->energy();

	    std::vector<const reco::GenParticle*> constituents = genJet->getGenConstituents();
	    genJetCollisionId[noJets] = constituents.at(0)->collisionId();
	  }
	  const reco::GenParticle * genParton = jet->genParton();
	  if (genParton == NULL) {
	    genPartonExist[noJets] = false;
	    genPartonPdg[noJets] = 0.;
	    genPartonPx[noJets] = 0.;
	    genPartonPy[noJets] = 0.;
	    genPartonPz[noJets] = 0.;
	    genPartonEn[noJets] = 0.;
	  }
	  else {
	    genPartonExist[noJets] = true;
	    genPartonPdg[noJets] = genParton->pdgId();
	    genPartonPx[noJets] = genParton->px();
	    genPartonPy[noJets] = genParton->py();
	    genPartonPz[noJets] = genParton->pz();
	    genPartonEn[noJets] = genParton->energy();
	  }

	}



	jetBProbBJetTag[noJets] = -10000.0;
	jetProbBJetTag[noJets] = -10000.0;
	tcHPBJetTag[noJets] = -10000.0;
	tcHEBJetTag[noJets] = -10000.0;
	svHEBJetTag[noJets] = -10000.0;
	svHPBJetTag[noJets] = -10000.0;
	combSVBJetTag[noJets] = -10000.0;
	combSVMVABJetTag[noJets] = -10000.0;
	const std::vector< std::pair< std::string, float > > pairDiscriVector = jet->getPairDiscri();
	int nDiscri = pairDiscriVector.size();
	if (_printOut>0) {
	  std::cout << " Jet" << noJets << " Et=" << jet->et() 
		    << " #const.=" << noConstituents[noJets]
		    << " #ch.=" << noChargedConstituents[noJets] 
		    << " NeutHadFr=" << neuHadronEnergyFr[noJets]
		    << " GammaFr=" << photonEnergyFr[noJets]
		    << " ChargedF=" << chargedEnergyFr[noJets]
		    << " ElecFr=" << electronEnergyFr[noJets] << std::endl;
	  std::cout << "        number of B-Jet Discriminants = " << nDiscri << std::endl;
	}
	for (int iD=0;iD<nDiscri;++iD) {
	  std::pair<std::string, float> pairDiscri = pairDiscriVector[iD];
	  TString nameOfDiscri(pairDiscri.first);
	  float tag = pairDiscri.second;

	  if (_printOut>0) 
	    std::cout << "     " << nameOfDiscri << " : " << tag << std::endl;    

	  if (nameOfDiscri==TString("jetBProbabilityBJetTags")) {
	    jetBProbBJetTag[noJets] = tag;
	    //	    std::cout << " jetBProbabilityBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("jetProbabilityBJetTags")) {
	    jetProbBJetTag[noJets] = tag;
	    //	    std::cout << " jetProbabilityBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("trackCountingHighPurBJetTags")) {
	    tcHPBJetTag[noJets] = tag;
	    //	    std::cout << " trackCountingHighPurBJetTags OK" << std::endl;
	  }
	  if (nameOfDiscri==TString("trackCountingHighEffBJetTags")) {
	    tcHEBJetTag[noJets] = tag;
	    //	    std::cout << " trackCountingHighEffBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("simpleSecondaryVertexHighEffBJetTags")) {
	    svHEBJetTag[noJets] = tag;
	    //	    std::cout << " simpleSecondaryVertexHighEffBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("simpleSecondaryVertexHighPurBJetTags")) {
	    svHPBJetTag[noJets] = tag;
	    //	    std::cout << " simpleSecondaryVertexHighPurBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("combinedSecondaryVertexBJetTags")) {
	    combSVBJetTag[noJets] = tag;
	    //	    std::cout << " combinedSecondaryVertexBJetTags OK " << std::endl;
	  }
	  if (nameOfDiscri==TString("combinedSecondaryVertexMVABJetTags")) {
	    combSVMVABJetTag[noJets] = tag;	  
	    //	    std::cout << " combinedSecondaryVertexBJetTags OK " << std::endl;
	  }
	}
	if (nSv<20) {
	  std::string sv3("secondaryVertex");
	  if (jet->hasTagInfo(sv3)) {
	    //	    std::cout << noJets << "   :   Tag info exist" << std::endl;
	    const reco::SecondaryVertexTagInfo* svti = jet->tagInfoSecondaryVertex(sv3);
	    if (svti != NULL) {
	      //	      std::cout << "n vertices = " << svti -> nVertices() << std::endl;
	      for ( unsigned int j = 0; j < svti -> nVertices() ; ++j )
		{
		  double fdSig =  svti->reco::SecondaryVertexTagInfo::flightDistance(j).significance();
		  if (fdSig>2.0) {
		    fdSv[nSv]    = float(svti->reco::SecondaryVertexTagInfo::flightDistance(j).value());
		    fdSigSv[nSv] = float(svti->reco::SecondaryVertexTagInfo::flightDistance(j).significance());
		    
		    const reco::Vertex & secondaryVertex = svti -> secondaryVertex(j);
		    massSv[nSv]  = secondaryVertex.p4().mass();	
		    nTrkSv[nSv]  = secondaryVertex.nTracks();
		    int nTrkCounter = 0;
		    for(reco::Vertex::trackRef_iterator iTrack  = secondaryVertex.tracks_begin(); iTrack != secondaryVertex.tracks_end();++iTrack) {
		      if (nTrkCounter<10) {
			pxTrkSv[nSv][nTrkCounter] =(*iTrack)->px();
			pyTrkSv[nSv][nTrkCounter] =(*iTrack)->py();
			pzTrkSv[nSv][nTrkCounter] =(*iTrack)->pz();
		      }
		      nTrkCounter++;
		    }
		    jetSv[nSv] = noJets;
		    nSv++;
		  }
		}
	    }
	  }
	}
	noJets++;
      }
    }
  }

  // *************************
  // accessing reco jets ====>
  // *************************
  edm::Handle <std::vector<reco::PFJet> > recoJetHandle;
  iEvent.getByLabel( _recoJetSrc, recoJetHandle );
  int numberOfRecoJets = recoJetHandle->size();
  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "Number of reco Jets = " << numberOfRecoJets << std::endl;
  }
	//  std::sort(correctedJetPairs.begin(), correctedJetPairs.end(), PFJetPairPtSorter);

  noRecoJets = 0;
  for (int iJ=0;iJ<numberOfRecoJets;++iJ) {
     // copy original jet
     const reco::PFJet& jet = (*recoJetHandle)[iJ];

     // remove jets which correspond to selected muons
     assert(selectedPair.size() == 2);
     if(reco::deltaR(*selectedPair[0], jet) < 0.4 || reco::deltaR(*selectedPair[1], jet) < 0.4)
       continue;

    if(jet.et() > _jetEnergyMinCut && noRecoJets < 50)
    {
      recoJetPt[noRecoJets] = jet.pt();
      recoJetEt[noRecoJets] = jet.et();
      recoJetEta[noRecoJets] = jet.eta();
      recoJetPhi[noRecoJets] = jet.phi();
      recoJetMass[noRecoJets] = jet.mass();
      recoJetArea[noRecoJets] = jet.jetArea();

      ret_loose.set(false);
      recoJetIDLoose[noRecoJets] = jetIDLooseFunctor(jet, ret_loose);

      ret_medium.set(false);
      recoJetIDMedium[noRecoJets] = jetIDMediumFunctor(jet, ret_medium);

      ret_tight.set(false);
      recoJetIDTight[noRecoJets] = jetIDTightFunctor(jet, ret_tight);

      ++noRecoJets;
    }
  }

  // ******************
  // accessing MET ===>
  // ******************

  double met = 0.0;
  edm::Handle<pat::METCollection> patMetHandle;
  iEvent.getByLabel(_patMetSrc,patMetHandle);
  int nPatMET = patMetHandle->size();
  if (nPatMET==0)
    return;
  pat::METRef patMet(patMetHandle,0);

  met = patMet->et();

  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "PAT MET  :   (Px,Py)=(" << patMet->px() << "," << patMet->py() << ")   MET=" << met << std::endl;
  }

  h_MET_->Fill(float(met));
  if (met<_metCut) {
    step = 7.0;
    h_cutFlow_->Fill(step);
  }
  else {
    return;
  }
  
  TLorentzVector METLorentz;
  METLorentz.SetPxPyPzE(patMet->px(),patMet->py(),0.0,patMet->pt());

  /*edm::Handle<pat::METCollection> patMetMVAHandle;
  iEvent.getByLabel(_patMetMVASrc, patMetMVAHandle);
  if(patMetMVAHandle->size() == 0) throw cms::Exception("No MVA MET available");
  pat::METRef patMetMVA(patMetMVAHandle,0);*/

  // ** accessing leptons -->
  
  
  reco::Candidate::LorentzVector negMuonLV = negMuon->p4();
  reco::Candidate::LorentzVector posMuonLV = posMuon->p4();

  double x1CollApprox = 0.;
  double x2CollApprox = 0.;

  double pxLeg1 = negMuonLV.Px();
  double pyLeg1 = negMuonLV.Py();

  double pxLeg2 = posMuonLV.Px();
  double pyLeg2 = posMuonLV.Py();

  double pxMEt  = patMet->px();
  double pyMEt  = patMet->py();

  compX1X2byCollinearApprox(x1CollApprox, x2CollApprox, 
			    pxLeg1,       pyLeg1, 
			    pxLeg2,       pyLeg2, 
			    pxMEt,        pyMEt);

  bool isX1WithinPhysRange = true;
  getPhysX(x1CollApprox, isX1WithinPhysRange);

  bool isX2WithinPhysRange = true;
  getPhysX(x2CollApprox, isX2WithinPhysRange);

  if (_printOut > 0 ) {
    std::cout << std::endl;
    if (isX1WithinPhysRange&&isX1WithinPhysRange) {
      reco::Candidate::LorentzVector diMuLV = negMuonLV + posMuonLV;
      double collDiTauMass = diMuLV.M()/TMath::Sqrt(x1CollApprox*x2CollApprox);
      std::cout << "TauAnalysisTools : X1 = " << x1CollApprox << "  X2 = " << x2CollApprox
		<< "  DiTauMass(CA) = " << collDiTauMass << std::endl;
    }
    else {
      std::cout << "TauAnalysisTools : Collinear approx. has no phys. solution... -> " 
		<< "  X1 = " << x1CollApprox 
		<< "  X2 = " << x2CollApprox << std::endl;
    }
    std::cout << std::endl;
    
  }

  //  ntpMassInfo.analyze(*posMuon, *negMuon, *(patMetMVA.get()) );

  

  // ** instantiating kinematics calculator --->
  HTauTauKinematicsCalculator kinematicsCalculator(PosLorentz,NegLorentz,METLorentz); // TODO: Switch to MVA MET here?


  // ** deltaPhi between muons --->
  double dileptonDPhi = kinematicsCalculator.diLeptonDeltaPhi();
  h_2LepDeltaPhi_->Fill(float(dileptonDPhi));
  
  if (_printOut > 0) {
    std::cout << std::endl;
    std::cout << "DeltaPhi 2Muon = " << dileptonDPhi << std::endl;
  }  
  if (dileptonDPhi>_diLepDPhiCut) {
    step = 8.0;
    h_cutFlow_->Fill(step); 
  }
  else {
    return;
  }

  // ** lepton isolation -->
  //  fillMuonMVAVariables(
  //		       *(posMuon.get()),
  //		       *(negMuon.get()),
  //		       primaryVtx,
  //		       *pflow,
  //		       iEvent, iSetup);

  _isoTrkNeg   = negMuon->trackIso();
  _isoECalNeg  = negMuon->ecalIso();
  _isoHCalNeg  = negMuon->hcalIso();
  _isoNeg      = float((_isoTrkNeg+_isoECalNeg+_isoHCalNeg)/(fmax(0.1,negMuon->pt())));

  _isoTrkPos   = posMuon->trackIso();
  _isoECalPos  = posMuon->ecalIso();
  _isoHCalPos  = posMuon->hcalIso();
  _isoPos      = float((_isoTrkPos+ _isoECalPos+_isoHCalPos)/(fmax(0.1,posMuon->pt())));

  // PFParticle based isolation ---->


  _isoChargedHadPFsNeg = negMuon->chargedHadronIso();
  _isoChargedHadPFsPos = posMuon->chargedHadronIso();
  
  _isoNeutralHadPFsNeg = negMuon->neutralHadronIso();
  _isoNeutralHadPFsPos = posMuon->neutralHadronIso();
  
  _isoPhotonsPFsNeg = negMuon->photonIso();
  _isoPhotonsPFsPos = posMuon->photonIso();

  _isoPFsNeg = _isoChargedHadPFsNeg + _isoNeutralHadPFsNeg + _isoPhotonsPFsNeg;
  _isoPFsPos = _isoChargedHadPFsPos + _isoNeutralHadPFsPos + _isoPhotonsPFsPos;
  

  if (_printOut > 0) {
    std::cout << std::endl;
    std::cout << "PF based isolation --> pos muon iso = " << _isoPFsPos << "  neg muon iso = " << _isoPFsNeg << std::endl;    
    std::cout << "Charged Had Iso    --> pos muon iso = " << _isoChargedHadPFsPos << "  neg muon iso = " << _isoChargedHadPFsNeg << std::endl;    
    std::cout << "Neutral Had Iso    --> pos muon iso = " << _isoNeutralHadPFsPos << "  neg muon iso = " << _isoNeutralHadPFsNeg << std::endl;    
    std::cout << "Photons Iso        --> pos muon iso = " << _isoPhotonsPFsPos << "  neg muon iso = " << _isoPhotonsPFsNeg << std::endl;    
    std::cout << std::endl;
  }

  if (_doMC) {
    _nFSRPos = 0;
    _nFSRNeg = 0;
    //  generated Z/gamma* and leptons ---->

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel(_genParticleSrc,genParticles);
    
    int nParticles = genParticles->size();
    
    bool posFound = false;
    bool negFound = false;
    //bool ZFound = false;
    //bool GammaFound = false;
    bool HiggsFound = false;
    
    VBFid = 0;
    isVBF = false;
    bbbarPhi = false;
    bool bfound =false;
    bool bbarfound = false;
    //bool parton1found = false;
    //bool parton2found = false;

    TLorentzVector posLep;
    TLorentzVector negLep;
    
    const Candidate * posMuonCand = NULL;
    const Candidate * negMuonCand = NULL;
    
    parton1pt = 0;
    parton2pt = 0;
    
    _bosonExist = false;
    for (int iP=0;iP<nParticles;++iP) {
      GenParticleRef part(genParticles,iP);
      int pdg = part->pdgId();
      
      // looking for b-jets ---->
      if (pdg==5){
	bfound=true;
	bPx   = part->px();
	bPy   = part->py();
	bPz   = part->pz();
	benergy=part->energy();
      }
      if (pdg==-5){
	bbarfound=true;
	bbarPx   = part->px();
	bbarPy   = part->py();
	bbarPz   = part->pz();
	bbarenergy=part->energy();
      }

      if(TMath::Abs(pdg)<=5){
	if (part->pt()>parton1pt){
	  
	  // assign variables for the second leading parton
	  parton2pt = parton1pt;

	  parton2pdg    = parton1pdg;
          parton2Px     = parton1Px;
          parton2Py     = parton1Py;
          parton2Pz     = parton1Pz;
          parton2energy = parton1energy;

	  parton1pt     = part->pt();
	  //parton1found  = true;
	  parton1pdg    = pdg;
	  parton1Px     = part->px();
	  parton1Py     = part->py();
	  parton1Pz     = part->pz();
	  parton1energy = part->energy();
	}
	else if (part->pt()>parton2pt) {
	  parton2pt     = part->pt();
	  //parton2found  = true;
	  parton2pdg    = pdg;
	  parton2Px     = part->px();
	  parton2Py     = part->py();
	  parton2Pz     = part->pz();
	  parton2energy = part->energy();

	}

      }

      if (pdg==23||pdg==22||pdg==25||pdg==35||pdg==36) { // Z or gamma* or Higgs
	int nD = part->numberOfDaughters();
	int decay = 0;
	/*if (pdg==23) 
	  ZFound = true;
	if (pdg==22)
	  GammaFound = true;*/
	if (pdg==25||pdg==35||pdg==36){
	  HiggsFound = true;

	  _HMass = part->mass();
	  _HPx   = part->px();
	  _HPy   = part->py();
	  _HPz   = part->pz();
	  
	  int nM = part->numberOfMothers();
	  for (int iM=0;iM<nM;++iM) {
	    const Candidate * mum = part->mother(iM);
	    if (mum->pdgId()==24 || mum->pdgId()==23)  
	      VBFid = mum->pdgId();
	  }
	}
	for (int iD=0;iD<nD;++iD) {
	  const Candidate * cand = part->daughter(iD);
	  if (cand->pdgId()==11||cand->pdgId()==13||cand->pdgId()==15) {
	    negFound = true;
	    negLep.SetPxPyPzE(cand->px(),cand->py(),cand->pz(),cand->energy());
	    decay = cand->pdgId();
	    negMuonCand = cand;
	  }
	  if (cand->pdgId()==-11||cand->pdgId()==-13||cand->pdgId()==-15) {
	    posFound = true;
	    posLep.SetPxPyPzE(cand->px(),cand->py(),cand->pz(),cand->energy());
	    posMuonCand = cand;
	  }
	}
	if (posFound&&negFound) {
	  _bosonPDG = pdg;
	  _decayPDG = TMath::Abs(decay);
	  
	  _posDecayMuonPx = float(posMuonCand->px());
	  _posDecayMuonPy = float(posMuonCand->py());
	  _posDecayMuonPz = float(posMuonCand->pz());
	  _posDecayMuonPdg = posMuonCand->pdgId();

	  _negDecayMuonPx = float(negMuonCand->px());
	  _negDecayMuonPy = float(negMuonCand->py());
	  _negDecayMuonPz = float(negMuonCand->pz());
	  _negDecayMuonPdg = negMuonCand->pdgId();

	  _ZMass = part->mass();
	  _ZPx   = part->px();
	  _ZPy   = part->py();
	  _ZPz   = part->pz();
	  
	  _nFSRPos = posMuonCand->numberOfDaughters();
	  for (int iD=0;iD<_nFSRPos;++iD) {
	    if (iD<100) {
	      const Candidate * daughter = posMuonCand->daughter(iD);
	      _FSRPosPdg[iD] = daughter->pdgId();
	      _FSRPosPx[iD] = daughter->px();
	      _FSRPosPy[iD] = daughter->py();
	      _FSRPosPz[iD] = daughter->pz();	      
	    }   
	  }	
	  
	  _nFSRNeg = negMuonCand->numberOfDaughters();
	  for (int iD=0;iD<_nFSRNeg;++iD) {
	    if (iD<100) {
	      const Candidate * daughter = negMuonCand->daughter(iD);
	      _FSRNegPdg[iD] = daughter->pdgId();
	      _FSRNegPx[iD] = daughter->px();
	      _FSRNegPy[iD] = daughter->py();
	      _FSRNegPz[iD] = daughter->pz();
	      
	    }   
	  }		  
	}
      }
      if (posFound&&negFound)
	break;
    }

    _bosonExist = posFound && negFound;

    bbbarPhi = bfound && bbarfound && HiggsFound;

    isVBF = (VBFid!=0) && HiggsFound;

     if (_printOut>0 && _bosonExist) { // debugging printout
       std::cout << std::endl;
      std::cout << "Generator info ---->   PDG of boson : " << _bosonPDG << "   boson mass = " << _ZMass << std::endl;
      std::cout << "Pos Lepton (Px,Py,Pz) = (" 
		<< posLep.Px() << "," 
		<< posLep.Py() << ","
		<< posLep.Pz() << ")" << std::endl;
      int nD = posMuonCand->numberOfDaughters();
      std::cout << "Status = " << posMuonCand->status() << "   Number of daughters = " << nD << "  :  " ;
      for (int iD=0; iD<nD; ++iD) {
	const Candidate * daughter = posMuonCand->daughter(iD);
	int pdgCode = daughter->pdgId();
	std::cout << "[" << pdgCode << "," << daughter->status() << "] ;  ";          
      }
      std::cout << std::endl;
      
      std::cout << "Neg Lepton (Px,Py,Pz) = (" 
		<< negLep.Px() << "," 
		<< negLep.Py() << ","
		<< negLep.Pz() << ")" << std::endl;
      nD = negMuonCand->numberOfDaughters();
      std::cout << "Status = " << negMuonCand->status() << "   Number of daughters = " << nD << "  :  " ;
      for (int iD=0; iD<nD; ++iD) {
	const Candidate * daughter = negMuonCand->daughter(iD);
	int pdgCode = daughter->pdgId();
	std::cout << "[" << pdgCode << "," << daughter->status() << "] ;   ";          
      }
      std::cout << std::endl;
    }

   /* // Embededd MC: Get decay information of original event
    if(!_origGenParticleSrc.label().empty())
    {
      iEvent.getByLabel(_origGenParticleSrc,genParticles);
      for (int iP=0;iP<nParticles;++iP) {
        GenParticleRef part(genParticles,iP);
        int pdg = part->pdgId();

        if (pdg==23||pdg==22||pdg==25||pdg==35||pdg==36) { // Z or gamma* or Higgs
	  int nD = part->numberOfDaughters();
 	  //int decay = 0;

          bool negFound = false, posFound = false;
          int origDecay = 0;
	  for (int iD=0;iD<nD;++iD) {
	    const Candidate * cand = part->daughter(iD);
	    if (cand->pdgId()==11||cand->pdgId()==13||cand->pdgId()==15) {
	      negFound = true;
	      origDecay = cand->pdgId();
	    }
	    if (cand->pdgId()==-11||cand->pdgId()==-13||cand->pdgId()==-15) {
	      posFound = true;
	    }
	  }

	  if (posFound&&negFound) {
	    _origBosonPDG = pdg;
	    _origDecayPDG = TMath::Abs(origDecay);
            break;
          }
        }
      }
    }*/

    // checking mothers and daughters in generator muons ---->

    const reco::GenParticle * negMuonGen = negMuon->genLepton();
    const reco::GenParticle * posMuonGen = posMuon->genLepton();

    // **************************************
    // Positive Generated Muon -------------> 
    // **************************************

    posMuonGenQ = 0.0;
    posMuonGenExist = false;
    PMmotherExist = false;
    PMGrmotherExist = false;
    posMuonGenpdg = 0;
    PMmotherpdg = 0;
    PMGrmotherpdg = 0;

    if (posMuonGen!=NULL) {

      posMuonGenExist = true;      
      posMuonGenpdg=posMuonGen->pdgId();
      posMuonGenQ = posMuonGen->charge();
      
      posMuonGenpx= posMuonGen->px();
      posMuonGenpy= posMuonGen->py();
      posMuonGenpz= posMuonGen->pz();
      posMuonGenpt= posMuonGen->pt();
      posMuonGeneta= posMuonGen->eta();
      posMuonGenphi= posMuonGen->phi();

      posMuonGenNoMothers= posMuonGen->numberOfMothers();

      if (posMuonGenNoMothers>0) {
	
	PMmotherExist = true;
	PMmomISmuon=true;
	const Candidate * PMmother=posMuonGen->mother(0);
	PMmotherpdg=PMmother->pdgId();
	
	if (PMmotherpdg!=-13 ) {
	  PMmomISmuon =false;
	  PMmotherpdg=PMmotherpdg;
	  PMmothereta=PMmother->eta();
	  PMmotherphi=PMmother->phi();
	  PMmotherpt=PMmother->pt();
	  PMmotherpx=PMmother->px();
	  PMmotherpy=PMmother->py();
	  PMmotherpz=PMmother->pz();
	}
	else  
	  PMmomISmuon=true;
	while( PMmomISmuon ) {
	  PMmotherNoMothers=PMmother->numberOfMothers();
	  if(PMmotherNoMothers>0){
	    const Candidate * grandmother=PMmother->mother(0);
	    PMmother=grandmother;
	    PMmotherpdg=PMmother->pdgId();	      
	    if (PMmotherpdg!=-13 ){
	      PMmomISmuon =false;
	      PMmothereta=PMmother->eta();
	      PMmotherphi=PMmother->phi();
	      PMmotherpt=PMmother->pt();
	      PMmotherpx=PMmother->px();
	      PMmotherpy=PMmother->py();
	      PMmotherpz=PMmother->pz();
	    }
	  }
	  else {
	    PMmotherExist = false;
	    PMmomISmuon = false;
	  }
	}
	
	bool isPosDmeson = isDmeson(PMmotherpdg);
	
	posMuonGenNoGrmothers=PMmother->numberOfMothers();

	if(posMuonGenNoGrmothers>0){
	  
	  PMGrmotherExist = true;
	  const Candidate * PMGrmother=PMmother->mother(0);
	  PMGrmotherpdg=PMGrmother->pdgId();
	  bool Looping = true;
	  
	  if (isPosDmeson) {
	    if (isBmeson(PMGrmotherpdg)) {
	      Looping = false;
	      PMGrmothereta=PMGrmother->eta();
	      PMGrmotherphi=PMGrmother->phi();
	      PMGrmotherpt=PMGrmother->pt();
	      PMGrmotherpx=PMGrmother->px();
	      PMGrmotherpy=PMGrmother->py();
	      PMGrmotherpz=PMGrmother->pz();
	    }
	  }
	  else {
	    if (PMGrmotherpdg != PMmotherpdg ){
	      Looping = false;
	      PMGrmothereta=PMGrmother->eta();
	      PMGrmotherphi=PMGrmother->phi();
	      PMGrmotherpt=PMGrmother->pt();
	      PMGrmotherpx=PMGrmother->px();
	      PMGrmotherpy=PMGrmother->py();
	      PMGrmotherpz=PMGrmother->pz();
	    }
	  }
	  
	  while( Looping ){
	    
	    PMGrmotherNoMothers=PMGrmother->numberOfMothers();
	    
	    if(PMGrmotherNoMothers>0){
	      
	      const Candidate * ancestor = PMGrmother->mother(0);
	      PMGrmother = ancestor;
	      PMGrmotherpdg = PMGrmother->pdgId();		

	      if ( isPosDmeson ) {
		if (isBmeson(PMGrmotherpdg)) {
		  Looping = false;
		  PMGrmothereta=PMGrmother->eta();
		  PMGrmotherphi=PMGrmother->phi();
		  PMGrmotherpt=PMGrmother->pt();
		  PMGrmotherpx=PMGrmother->px();
		  PMGrmotherpy=PMGrmother->py();
		  PMGrmotherpz=PMGrmother->pz();
		}
	      }
	      else {
		if ( PMGrmotherpdg != PMmotherpdg ) {
		  Looping = false;
		  PMGrmothereta=PMGrmother->eta();
		  PMGrmotherphi=PMGrmother->phi();
		  PMGrmotherpt=PMGrmother->pt();
		  PMGrmotherpx=PMGrmother->px();
		  PMGrmotherpy=PMGrmother->py();
		  PMGrmotherpz=PMGrmother->pz();
		}
	      }		  
	    }
	    else {		
	      Looping = false;
	      PMGrmotherExist = false;
	    }
	  }	    
	}
      }      
    }

    // **********************************
    // Negative Generated Muon --------->
    // **********************************

    NMmotherExist = false;
    NMGrmotherExist = false;
    negMuonGenExist = false;
    negMuonGenQ = 0.0;
    negMuonGenpdg = 0;
    NMmotherpdg = 0;
    NMGrmotherpdg = 0;
    if (negMuonGen!=NULL) {

      negMuonGenExist = true;
      negMuonGenQ = negMuonGen->charge();
      negMuonGenpdg=negMuonGen->pdgId();

      negMuonGenpx= negMuonGen->px();
      negMuonGenpy= negMuonGen->py();
      negMuonGenpz= negMuonGen->pz();
      negMuonGenpt= negMuonGen->pt();
      negMuonGeneta= negMuonGen->eta();
      negMuonGenphi= negMuonGen->phi();

      negMuonGenNoMothers = negMuonGen->numberOfMothers();

      if ( negMuonGenNoMothers>0) {

	NMmotherExist = true; 
	NMmomISmuon = true;
	const Candidate * NMmother=negMuonGen->mother(0);
	NMmotherpdg=NMmother->pdgId();	
	if (NMmotherpdg!=13 ) {
	  NMmomISmuon =false;
	  NMmotherpdg=NMmotherpdg;
	  NMmothereta=NMmother->eta();
	  NMmotherphi=NMmother->phi();
	  NMmotherpt=NMmother->pt();
	  NMmotherpx=NMmother->px();
	  NMmotherpy=NMmother->py();
	  NMmotherpz=NMmother->pz();
	}
	else  
	  NMmomISmuon=true;
	
	while ( NMmomISmuon ) {

	  NMmotherNoMothers=NMmother->numberOfMothers();
	  
	  if(NMmotherNoMothers>0){
	    const Candidate * grandmother=NMmother->mother(0);
	    NMmother=grandmother;
	    NMmotherpdg=NMmother->pdgId();	      
	    
	    if (NMmotherpdg!=13 ){
	      NMmomISmuon =false;
	      NMmothereta=NMmother->eta(); 
	      NMmotherphi=NMmother->phi();
	      NMmotherpt=NMmother->pt();
	      NMmotherpx=NMmother->px();
	      NMmotherpy=NMmother->py();
	      NMmotherpz=NMmother->pz();
	    }
	  }
	  else {
	    NMmotherExist = false;
	    NMmomISmuon = false;
	  }
	}

	negMuonGenNoGrmothers=NMmother->numberOfMothers();
	
	bool isNegDmeson = isDmeson(NMmotherpdg);

	if(negMuonGenNoGrmothers>0) {
	  
	  NMGrmotherExist = true;
	  const Candidate * NMGrmother = NMmother->mother(0);
	  NMGrmotherpdg = NMGrmother->pdgId();
	  bool Looping = true;	    	    
	  
	  if (isNegDmeson) {
	    if (isBmeson(NMGrmotherpdg)) {
	      Looping =false;
	      NMGrmothereta=NMGrmother->eta();
	      NMGrmotherphi=NMGrmother->phi();
	      NMGrmotherpt=NMGrmother->pt();
	      NMGrmotherpx=NMGrmother->px();
	      NMGrmotherpy=NMGrmother->py();
	      NMGrmotherpz=NMGrmother->pz();
	    }
	  }
	  else {
	    if (NMGrmotherpdg!=NMmotherpdg) {
	      Looping =false;
	      NMGrmothereta=NMGrmother->eta();
	      NMGrmotherphi=NMGrmother->phi();
	      NMGrmotherpt=NMGrmother->pt();
	      NMGrmotherpx=NMGrmother->px();
	      NMGrmotherpy=NMGrmother->py();
	      NMGrmotherpz=NMGrmother->pz();
	    }
	  }
	  
	  while( Looping ) {
	    
	    NMGrmotherNoMothers=NMGrmother->numberOfMothers();

	    if(NMGrmotherNoMothers>0) {
	      const Candidate * ancestor = NMGrmother->mother(0);
	      NMGrmother = ancestor;
	      NMGrmotherpdg = NMGrmother->pdgId();
	      
	      if (isNegDmeson) {
		if (isBmeson(NMGrmotherpdg)) {
		  Looping = false;
		  NMGrmothereta=NMGrmother->eta();
		  NMGrmotherphi=NMGrmother->phi();
		  NMGrmotherpt=NMGrmother->pt();
		  NMGrmotherpx=NMGrmother->px();
		  NMGrmotherpy=NMGrmother->py();
		  NMGrmotherpz=NMGrmother->pz();
		}
	      }
	      else {
		if (NMGrmotherpdg!=NMmotherpdg ){
		  Looping = false;
		  NMGrmothereta=NMGrmother->eta();
		  NMGrmotherphi=NMGrmother->phi();
		  NMGrmotherpt=NMGrmother->pt();
		  NMGrmotherpx=NMGrmother->px();
		  NMGrmotherpy=NMGrmother->py();
		  NMGrmotherpz=NMGrmother->pz();
		}
	      }

	    }
	    else {
	      NMGrmotherExist =false;
	      Looping = false;	      
	    }
	  }	    
	}
      }
    }

    if (NMmotherExist && PMmotherExist) {

      if (_printOut>0) {
	std::cout << std::endl;
	std::cout << "First muon  -->  q (reco) =  " << posMuon->charge() 
		  << "  q (gen)  = " << posMuonGenQ 
		  << "  mother pdg = " << PMmotherpdg 
		  << "  (pt,eta,phi)=(" << PMmotherpt 
		  << "," << PMmothereta 
		  << "," << PMmotherphi << ")";
	if (PMGrmotherExist)
	  std::cout << "  grandmother pdg = " << PMGrmotherpdg 
		    << "  (pt,eta,phi)=(" << PMGrmotherpt 
		    << "," << PMGrmothereta 
		    << "," << PMGrmotherphi << ")";
	else 
	  std::cout << "  no grandmother";      
	std::cout << std::endl;
	std::cout << "Second muon -->  q (reco) = " << negMuon->charge() 
		  << "  q (gen) = " << negMuonGenQ 
		  << "  mother pdg = " << NMmotherpdg 
		  << "  (pt,eta,phi)=(" << NMmotherpt 
		  << "," << NMmothereta 
		  << "," << NMmotherphi << ")";
	if (NMGrmotherExist)
	  std::cout << "  grandmother pdg = " << NMGrmotherpdg 
		    << "  (pt,eta,phi)=(" << NMGrmotherpt 
		    << "," << NMGrmothereta 
		    << "," << NMGrmotherphi << ")";
	else 
	  std::cout << "  no grandmother";      
	std::cout << std::endl;
	std::cout << std::endl;
	
      }

    }

  }

  _isPosGlobalMu = posMuon->isGlobalMuon();
  _isNegGlobalMu = negMuon->isGlobalMuon();

  _isPosTrackerMu = posMuon->isTrackerMuon();
  _isNegTrackerMu = negMuon->isTrackerMuon();

  _isPosPFMu = isPFMuon(posMuon, *pflow);
  _isNegPFMu = isPFMuon(negMuon, *pflow);

  // *** accessing global muon tracks

  reco::TrackRef globalTrackPos = posMuon->globalTrack();
  reco::TrackRef globalTrackNeg = negMuon->globalTrack();
  
  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "checking muon inner tracks ---> " << std::endl;
  }

  if (globalTrackPos.isNull()) {
    if (_printOut>0)
      std::cout << " inner positive track absent!" << std::endl;
    return;
  }
  else {
    if (_printOut>0)
      std::cout << " inner positive track found!" << std::endl;
  }

  if (globalTrackNeg.isNull()) {
    if (_printOut>0)
      std::cout << " inner negative track absent!" << std::endl;
    return;
  }
  else {
    if (_printOut>0)
      std::cout << " inner negative track found!" << std::endl;
  }
  
  _chi2PosMu = globalTrackPos->chi2();
  _chi2NegMu = globalTrackNeg->chi2();

  _ndofPosMu = globalTrackPos->ndof();
  _ndofNegMu = globalTrackNeg->ndof();

  const reco::HitPattern hitPatternGlobalPos = globalTrackPos->hitPattern();
  const reco::HitPattern hitPatternGlobalNeg = globalTrackNeg->hitPattern();

  _nMuonHitsPos = hitPatternGlobalPos.numberOfValidMuonHits();
  _nMuonHitsNeg = hitPatternGlobalNeg.numberOfValidMuonHits();

  // *** muon station Matches

  _nMuonStationsPos = posMuon->numberOfMatchedStations();
  _nMuonStationsNeg = negMuon->numberOfMatchedStations();


  // *** accessing inner muon tracks ---->
  
  reco::TrackRef trackPos = posMuon->innerTrack();
  reco::TrackRef trackNeg = negMuon->innerTrack();

  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "checking muon inner tracks ---> " << std::endl;
  }

  if (trackPos.isNull()) {
    if (_printOut>0)
      std::cout << " inner positive track absent!" << std::endl;
    return;
  }
  else {
    if (_printOut>0)
      std::cout << " inner positive track found!" << std::endl;
  }

  if (trackNeg.isNull()) {
    if (_printOut>0)
      std::cout << " inner negative track absent!" << std::endl;
    return;
  }
  else {
    if (_printOut>0)
      std::cout << " inner negative track found!" << std::endl;
  }

  const reco::HitPattern hitPatternPos = trackPos->hitPattern();
  const reco::HitPattern hitPatternNeg = trackNeg->hitPattern();

  _nPixelHitsPos = hitPatternPos.numberOfValidPixelHits();
  _nPixelHitsNeg = hitPatternNeg.numberOfValidPixelHits();

  _nTrackerHitsPos = hitPatternPos.numberOfValidHits();
  _nTrackerHitsNeg = hitPatternNeg.numberOfValidHits();

  if (_printOut > 0) {
    std::cout << std::endl;
    std::cout << "pos muon : n(pixel) = " << _nPixelHitsPos
	      << "  ; n(tracker) = " << _nTrackerHitsPos
	      << "  ; n(muon) = " << _nMuonHitsPos
	      << "  ; n(stations) = " << _nMuonStationsPos
	      << "  ; isGlobal = " << _isPosGlobalMu
	      << "  ; isTracker = " << _isPosTrackerMu
	      << "  ; chi2/ndof = " << _chi2PosMu/_ndofPosMu << std::endl;  
    std::cout << "neg muon : n(pixel) = " << _nPixelHitsNeg
	      << "  ; n(tracker) = " << _nTrackerHitsNeg
	      << "  ; n(muon) = " << _nMuonHitsNeg
	      << "  ; n(stations) = " << _nMuonStationsNeg
	      << "  ; isGlobal = " << _isNegGlobalMu
	      << "  ; isTracker = " << _isNegTrackerMu 
	      << "  : chi2/ndof = " << _chi2NegMu/_ndofNegMu << std::endl;

  }

  TransientTrack transientPosTrack = transientTrackBuilder->build(*trackPos);
  TransientTrack transientNegTrack = transientTrackBuilder->build(*trackNeg);

  TransientTrack * trPosTrkPtr = & transientPosTrack;
  TransientTrack * trNegTrkPtr = & transientNegTrack;

  FreeTrajectoryState posState = trPosTrkPtr->impactPointTSCP().theState();
  FreeTrajectoryState negState = trNegTrkPtr->impactPointTSCP().theState();

  _twoMuonDist3D  = -1.0;
  _twoMuonDist3DE = -1.0;
  
  _twoMuonDist2D  = -1.0;
  _twoMuonDist2DE = -1.0;
  
  _twoMuonDistRPhi3D  = -1.0;
  _twoMuonDistRPhi3DE = -1.0;
  
  _twoMuonDistRPhi2D  = -1.0;
  _twoMuonDistRPhi2DE = -1.0;
  

  if (trPosTrkPtr->impactPointTSCP().isValid() &&  trNegTrkPtr->impactPointTSCP().isValid()) {

    ClosestApproachInRPhi cApp;
    TwoTrackMinimumDistance minDist;

    typedef ROOT::Math::SVector<double, 3> SVector3;
    typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;	   

    minDist.calculate(posState,negState);
    if (minDist.status()) {

      float dist3D = minDist.distance();
      std::pair<GlobalPoint,GlobalPoint> pcaMuons = minDist.points();
      GlobalPoint posPCA = pcaMuons.first;
      GlobalPoint negPCA = pcaMuons.second;

      ParticleMass muon_mass = 0.105658;
      float muon_sigma = muon_mass*1.e-6;

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;

      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      RefCountedKinematicParticle posMuonPart = pFactory.particle(transientPosTrack,muon_mass,chi,ndf,muon_sigma); 
      RefCountedKinematicParticle negMuonPart = pFactory.particle(transientNegTrack,muon_mass,chi,ndf,muon_sigma); 

      SVector3 distanceVector(posPCA.x()-negPCA.x(),
			      posPCA.y()-negPCA.y(),
			      posPCA.z()-negPCA.z());

      _twoMuonDist3D = ROOT::Math::Mag(distanceVector);
      
      std::vector<float> vvv(6);

      vvv[0] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,0);
      vvv[1] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,1);	   
      vvv[2] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(1,1);
      vvv[3] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,2);
      vvv[4] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(1,2);	   
      vvv[5] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(2,2);
      
      SMatrixSym3D posPCACov(vvv.begin(),vvv.end());

      vvv[0] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,0);
      vvv[1] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,1);	   
      vvv[2] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(1,1);
      vvv[3] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,2);
      vvv[4] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(1,2);	   
      vvv[5] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(2,2);
      
      SMatrixSym3D negPCACov(vvv.begin(),vvv.end());

      SMatrixSym3D totCov = posPCACov + negPCACov;
      
      _twoMuonDist3DE = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/_twoMuonDist3D;

      if (_printOut>0)
	std::cout << "distance between muons : " << dist3D << "  vs.  3D distance sig = " 
		  << _twoMuonDist3D << "/" << _twoMuonDist3DE << " = " 
		  << _twoMuonDist3D/_twoMuonDist3DE << std::endl;

      distanceVector(2) = 0.0;
      _twoMuonDist2D = ROOT::Math::Mag(distanceVector);
      _twoMuonDist2DE = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/_twoMuonDist2D;


    }

    cApp.calculate(posState,negState);
    if (cApp.status()) {

      float dist3D = cApp.distance();
      std::pair<GlobalPoint,GlobalPoint> pcaMuons = cApp.points();
      GlobalPoint posPCA = pcaMuons.first;
      GlobalPoint negPCA = pcaMuons.second;

      ParticleMass muon_mass = 0.105658;
      float muon_sigma = muon_mass*1.e-6;

      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;

      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      RefCountedKinematicParticle posMuonPart = pFactory.particle(transientPosTrack,muon_mass,chi,ndf,muon_sigma); 
      RefCountedKinematicParticle negMuonPart = pFactory.particle(transientNegTrack,muon_mass,chi,ndf,muon_sigma); 

      SVector3 distanceVector(posPCA.x()-negPCA.x(),
			      posPCA.y()-negPCA.y(),
			      posPCA.z()-negPCA.z());

      _twoMuonDistRPhi3D = ROOT::Math::Mag(distanceVector);
      
      std::vector<float> vvv(6);

      vvv[0] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,0);
      vvv[1] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,1);	   
      vvv[2] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(1,1);
      vvv[3] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(0,2);
      vvv[4] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(1,2);	   
      vvv[5] = posMuonPart->stateAtPoint(posPCA).kinematicParametersError().matrix()(2,2);
      
      SMatrixSym3D posPCACov(vvv.begin(),vvv.end());

      vvv[0] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,0);
      vvv[1] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,1);	   
      vvv[2] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(1,1);
      vvv[3] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(0,2);
      vvv[4] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(1,2);	   
      vvv[5] = negMuonPart->stateAtPoint(negPCA).kinematicParametersError().matrix()(2,2);
      
      SMatrixSym3D negPCACov(vvv.begin(),vvv.end());

      SMatrixSym3D totCov = posPCACov + negPCACov;
      
      _twoMuonDistRPhi3DE = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/_twoMuonDistRPhi3D;

      if (_printOut>0)
	std::cout << "rphi distance between muons : " << dist3D << "  vs. rphi distance sig = " 
		  << _twoMuonDistRPhi3D << "/" << _twoMuonDistRPhi3DE << " = " 
		  << _twoMuonDistRPhi3D/_twoMuonDistRPhi3DE << std::endl;

      distanceVector(2) = 0.0;

      _twoMuonDistRPhi2D = ROOT::Math::Mag(distanceVector);

      _twoMuonDistRPhi2DE = sqrt(ROOT::Math::Similarity(totCov, distanceVector))/_twoMuonDistRPhi2D;


    }

  }


  for (size_t iPV=0;iPV<pVtxs.size();++iPV) {

    reco::Vertex currentVtx = pVtxs.at(iPV);

    if (_printOut>0)
      std::cout << "Vertex " << iPV << "   (x,y,z)=(" 
		<< currentVtx.x() << "," 
		<< currentVtx.y() << "," 
		<< currentVtx.z() << ")" << std::endl;

    IPMeasurement ipMeasPosTrk = ipTool.impactParameter(transientPosTrack,currentVtx);
    IPMeasurement ipMeasNegTrk = ipTool.impactParameter(transientNegTrack,currentVtx);

    if (_printOut>0)
      std::cout << "ip is measured for PV[" << iPV << "]" << std::endl;
      

    Measurement1D ipSigPosIP = ipMeasPosTrk.ip;
    Measurement1D ipSigNegIP = ipMeasNegTrk.ip;
  
    Measurement1D ipSigPosIPZ = ipMeasPosTrk.ipZ;
    Measurement1D ipSigNegIPZ = ipMeasNegTrk.ipZ;
  
    Measurement1D ipSigPosIP3D = ipMeasPosTrk.ip3D;
    Measurement1D ipSigNegIP3D = ipMeasNegTrk.ip3D;
  

    _ipSig2DPos[iPV] = ipSigPosIP.significance();
    _ipSig2DNeg[iPV] = ipSigNegIP.significance();
    _ip2DPos[iPV] = ipSigPosIP.value();
    _ip2DNeg[iPV] = ipSigNegIP.value();
  
    _ipSig3DPos[iPV] = ipSigPosIP3D.significance();
    _ipSig3DNeg[iPV] = ipSigNegIP3D.significance();
    _ip3DPos[iPV] = ipSigPosIP3D.value();
    _ip3DNeg[iPV] = ipSigNegIP3D.value();
  
    _ipSigZPos[iPV] = ipSigPosIPZ.significance();
    _ipSigZNeg[iPV] = ipSigNegIPZ.significance();
    _ipZPos[iPV] = ipSigPosIPZ.value();
    _ipZNeg[iPV] = ipSigNegIPZ.value();
  
    if (_printOut>0) {     
      std::cout << "       ipSig(r-phi,+) = " <<  _ipSig2DPos[iPV] << "  ipSig(r-phi,mu-) = " <<  _ipSig2DNeg[iPV] << std::endl;
    }

  }

  
  TLorentzVector diLepton = kinematicsCalculator.diLepton4P();
  
  float diLeptonPL = float(fabs(diLepton.Pz()));
  h_diLeptonPL_->Fill(diLeptonPL);
  
  float diLeptonEta = float(diLepton.Eta());
  h_diLeptonEta_->Fill(diLeptonEta);
  
  float diLeptonMass = float(kinematicsCalculator.diLeptonMass());

  if (_printOut > 0) {
    std::cout << std::endl;
    std::cout << "dilepton mass from kinematicsCalculator = " << diLeptonMass << std::endl;
  }

  
  h_diLepMass_->Fill(diLeptonMass);
  
  float cosNegDiLepton = float(kinematicsCalculator.cosNegDiLeptonRestFrame());
  h_cosNegDiLepton_->Fill(cosNegDiLepton);
  
  float mTDiLepton = float(kinematicsCalculator.mTdiLepton());
  h_mTDiLepton_->Fill(mTDiLepton);
  
  float mTDiLeptonMET = float(kinematicsCalculator.mTdiLeptonMET());
  h_mTDiLeptonMET_->Fill(mTDiLeptonMET);

  float posLepMETDPhi = float(kinematicsCalculator.posLepMETDeltaPhi());
  h_lepPosMETDPhi_->Fill(posLepMETDPhi);
  
  float negLepMETDPhi = float(kinematicsCalculator.negLepMETDeltaPhi());
  h_lepNegMETDPhi_->Fill(negLepMETDPhi);

  float diLepMETDPhi = float(kinematicsCalculator.diLeptonMETDeltaPhi());
  h_diLepMETDPhi_->Fill(diLepMETDPhi);

  float eDiff = float(posMuon->energy() - negMuon->energy());
  h_lepEDif_->Fill(eDiff);

  float ptDiff = float(posMuon->pt() - negMuon->pt());
  h_lepPtDif_->Fill(ptDiff);

  _diLepMassVar  = diLeptonMass;
  _diLepDPhi = float(kinematicsCalculator.diLeptonDeltaPhi());
  _diLepDEta = float(kinematicsCalculator.diLeptonDeltaEta());
  _diLepDR = float(kinematicsCalculator.diLeptonDeltaR());
  _posLepMETDPhi = posLepMETDPhi;
  _negLepMETDPhi = negLepMETDPhi;
  _met = float(met);

  _metPx = patMet->px();
  _metPy = patMet->py();
  _metPz = patMet->pz();
  _metEn = patMet->energy();
  /*_metMVAPx = patMetMVA->px();
  _metMVAPy = patMetMVA->py();
  _metMVAPz = patMetMVA->pz();
  _metMVAEn = patMetMVA->energy();*/

  /*TMatrixD covMET(2,2);
  covMET=patMet->getSignificanceMatrix();
  _metCovXX = float(covMET(0,0));
  _metCovXY = float(covMET(0,1));
  _metCovYX = float(covMET(1,0));
  _metCovYY = float(covMET(1,1));

  TMatrix covMETMVA(2,2);
  covMETMVA=patMetMVA->getSignificanceMatrix();
  _metMVACovXX = float(covMETMVA(0,0));
  _metMVACovXY = float(covMETMVA(0,1));
  _metMVACovYX = float(covMETMVA(1,0));
  _metMVACovYY = float(covMETMVA(1,1));*/


  // TODO:
_diLepEta = diLeptonEta;
_posLepPt = posMuon->pt();
_negLepPt = negMuon->pt();
_posLepPx = posMuon->px();
_negLepPx = negMuon->px();
_posLepPy = posMuon->py();
_negLepPy = negMuon->py();
_posLepPz = posMuon->pz();
_negLepPz = negMuon->pz();
_posLepEta = posMuon->eta();
_negLepEta = negMuon->eta();
_posLepPhi = posMuon->phi();
_negLepPhi = negMuon->phi();
_posLepEn = posMuon->energy();
_negLepEn = negMuon->energy();

  _posLepQ = posMuon->charge();
  _negLepQ = negMuon->charge();

  _diLepPhi = float(diLepton.Phi());
  _diLepPt  = float(diLepton.Pt());


  double ipPos = trackPos->dxy(primaryVtx.position());
  double ipNeg = trackNeg->dxy(primaryVtx.position());
  
  double chi2TrkPos = trackPos->chi2();
  int ndfTrkPos = int(trackPos->ndof());
  double probTrkPos = TMath::Prob(chi2TrkPos,ndfTrkPos);
  
  double chi2TrkNeg = trackNeg->chi2();
  int ndfTrkNeg = int(trackNeg->ndof());
  double probTrkNeg = TMath::Prob(chi2TrkNeg,ndfTrkNeg);
  
  int nHitsTrkPos = trackPos->numberOfValidHits();
  const reco::HitPattern hitPatternP = trackPos->hitPattern();
  int nPixelHitsTrkPos = hitPatternP.numberOfValidPixelHits();
  
  int nHitsTrkNeg = trackNeg->numberOfValidHits();
  const reco::HitPattern hitPatternN = trackNeg->hitPattern();
  int nPixelHitsTrkNeg = hitPatternN.numberOfValidPixelHits();
  
  double ipMax = ipPos;
  if (ipNeg>ipMax)
    ipMax = ipNeg;
  
  if (ipMax>1.0) {
    h_xVtxStrange_->Fill(primaryVtx.x());
    h_yVtxStrange_->Fill(primaryVtx.y());     
  }
  else {
    h_xVtxNormal_->Fill(primaryVtx.x());
    h_yVtxNormal_->Fill(primaryVtx.y());
  }
  
  
  if (ipPos>1.0) {
    h_nHitsStrange_->Fill(float(nHitsTrkPos));
    h_nVtxHitsStrange_->Fill(float(nPixelHitsTrkPos));
    h_probTrkStrange_->Fill(float(probTrkPos));
    ipPos = 0.99;
  }
  else {
    h_nHitsNormal_->Fill(float(nHitsTrkPos));
    h_nVtxHitsNormal_->Fill(float(nPixelHitsTrkPos));
    h_probTrkNormal_->Fill(float(probTrkPos));
  }
  
  if (ipNeg>1.0) {
    h_nHitsStrange_->Fill(float(nHitsTrkNeg));
    h_nVtxHitsStrange_->Fill(float(nPixelHitsTrkNeg));
    h_probTrkStrange_->Fill(float(probTrkNeg));
    ipNeg = 0.99;
  }
  else {
    h_nHitsNormal_->Fill(float(nHitsTrkNeg));
    h_nVtxHitsNormal_->Fill(float(nPixelHitsTrkNeg));
    h_probTrkNormal_->Fill(float(probTrkNeg));     
  }
  
  h_ipLeptonPosTrk_->Fill(fabs(float(ipPos)));
  h_ipLeptonNegTrk_->Fill(fabs(float(ipNeg)));

  _validDiTau = false;
  
  if ( _printOut > 0 ) {
    if (kinematicsCalculator.isValidDiTau())
      std::cout << "Collinear approximation solution is found..." << std::endl;
    else 
      std::cout << "Unphysical solution for collinear approximation..." << std::endl;
  }

  _diTauMass = -1.0;

  if (kinematicsCalculator.isValidDiTau()) {
    
    _validDiTau = true;

    step = 9.0;
    h_cutFlow_->Fill(step);
    
    _diTauMass = float(kinematicsCalculator.diTauMass());

    if (_printOut > 0 ) {
      std::cout << "DiTauMass(kinematicCalculator) = " << _diTauMass << std::endl;
      std::cout << std::endl;
    }

    h_diTauMass_->Fill(_diTauMass);
    
    float massDif = _diTauMass - diLeptonMass;
    h_diTaudiLepMass_->Fill(massDif);
    
    _posLepETauRest = float(kinematicsCalculator.posLepETauRestFrame());
    _negLepETauRest = float(kinematicsCalculator.negLepETauRestFrame());
    h_lepETauRestFrame_->Fill(_posLepETauRest);
    h_lepETauRestFrame_->Fill(_negLepETauRest);
    
    _posLepEdiTauRest = float(kinematicsCalculator.posLepEdiTauRestFrame());
    _negLepEdiTauRest = float(kinematicsCalculator.negLepEdiTauRestFrame());
    h_lepEDifDiTauRestFrame_->Fill(_negLepEdiTauRest-_posLepEdiTauRest);
    h_lepESumDiTauRestFrame_->Fill(_negLepEdiTauRest+_posLepEdiTauRest);
    
    _cosNegDiTau = float(kinematicsCalculator.cosNegTauRestFrame());
    h_cosNegDiTau_->Fill(_cosNegDiTau);
    
    _posLepTauERatio = float(kinematicsCalculator.posTauELepETauRatio());
    _negLepTauERatio = float(kinematicsCalculator.negTauELepETauRatio());
    h_lepTauERatio_->Fill(_posLepTauERatio);
    h_lepTauERatio_->Fill(_negLepTauERatio);
    
    TLorentzVector diTau = kinematicsCalculator.diTau4P();
    _diTauPx = float(fabs(diTau.Px()));
    _diTauPy = float(fabs(diTau.Py()));
    _diTauPz = float(fabs(diTau.Pz()));
    _diTauEta = float(diTau.Eta());
    _diTauPt = float(diTau.Pt());
    _diTauE = float(diTau.E());

    h_diTauPL_->Fill(_diTauPz);
    h_diTauEta_->Fill(_diTauEta);
    
  }

  if (_printOut>0) 
    std::cout << "about to fill TTree" << std::endl;

  _tree->Fill();

  if (_printOut>0)
    std::cout << "TTree is filled " << std::endl;


  _selEvent++;

  if (_printOut>0) {
    std::cout << std::endl;
    std::cout << "Event selected ----> " << std::endl;
    std::cout << "Number of selected events = " << _selEvent << std::endl;
  }
  
}


// ------------ method called once each job just before starting event loop  ------------
void HTauTauNMSSMAnalysis::beginJob() {

  _event = 0;
  _selEvent = 0;

  _firstEvent = true;



  edm::Service<TFileService> fs;

  // ** declaration of histograms --->

  h_cutFlow_=fs->make<TH1F>("h_cutFlow_","Cut flow",21,-0.5,20.5);

  h_numLep_=fs->make<TH1F>("h_numLep_","Number of muons",11,-0.5,10.5);
  h_ptLep_=fs->make<TH1F>("h_ptLep_","Muon pt",50,0.,100.);
  h_etaLep_=fs->make<TH1F>("h_etaLep_","Muon eta",50,-3.,3.);
  h_isoLep_=fs->make<TH1F>("h_isoLep_","Muon isolation",50,0.,2.);
  h_ptHard_=fs->make<TH1F>("h_ptHard_","Pt of hard Muon",50,0.,100.);
  h_jetEt_=fs->make<TH1F>("h_jetEt_","Jet Et",50,0.,150.);
  h_numJets_=fs->make<TH1F>("h_numJets_","number of jets",16,-0.5,15.5);
  h_2LepDeltaPhi_=fs->make<TH1F>("h_2LepDeltaPhi_","delta phi of muons",50,0.,TMath::Pi());
  h_MET_=fs->make<TH1F>("h_MET_","MET",50,0.,100.);

  h_diLepMass_=fs->make<TH1F>("h_diLepMass_","Dilepton mass",40,2.5,202.5);
  h_diLepMassEtaCut_=fs->make<TH1F>("h_diLepMassEtaCut_","Dilepton mass",40,2.5,202.5);

  h_cosNegDiLepton_=fs->make<TH1F>("h_cosNegDiLepton_","cos(negLep) in dilepton rest frame",20,-1.,1.);
  h_mTDiLepton_=fs->make<TH1F>("h_mTDiLepton_","Dilepton transverse mass",40,0.,200.);
  h_mTDiLeptonMET_=fs->make<TH1F>("h_mTDiLeptonMET_","DileptonMET transverse mass",40,0.,200.);

  h_lepPosMETDPhi_=fs->make<TH1F>("h_lepPosMETDPhi_","Pos Lep MET dPhi",50,0.,TMath::Pi());
  h_lepNegMETDPhi_=fs->make<TH1F>("h_lepNegMETDPhi_","Neg Lep MET dPhi",50,0.,TMath::Pi());

  h_diLepMETDPhi_=fs->make<TH1F>("h_diLepMETDPhi_","diLep MET dPhi",50,0.,TMath::Pi());

  h_lepEDif_=fs->make<TH1F>("h_lepEDif_","pos lep E - neg lep E",200,-100.,100.);
  h_lepPtDif_=fs->make<TH1F>("h_lepPtDif_","pos lep Pt - neg lep Pt",100,-50.,50.);
  
  h_diLeptonPL_=fs->make<TH1F>("h_diLeptonPL_","PL of dilepton",100,0.,500.);
  h_diLeptonEta_=fs->make<TH1F>("h_diLeptonEta_","Eta dilepton",100,-5.,5.);

  h_ipLeptonPosTrk_=fs->make<TH1F>("h_ipLeptonPosTrk_","pos lepton ip",250,0.,1.0);
  h_ipLeptonNegTrk_=fs->make<TH1F>("h_ipLeptonNegTrk_","neg lepton ip",250,0.,1.0);
  h_ipSigLeptonPosTrk_=fs->make<TH1F>("h_ipSigLeptonPosTrk_","IP sign pos lept",120,-2.,1.);
  h_ipSigLeptonNegTrk_=fs->make<TH1F>("h_ipSigLeptonNegTrk_","IP sign neg lept",120,-2.,1.);

  h_diTauMass_=fs->make<TH1F>("h_diTauMass_","di tau mass",40,0.,200.);
  h_diTaudiLepMass_=fs->make<TH1F>("h_diTaudiLepMass_","diTau diLep mass difference",50,-50.,200.);
  h_lepETauRestFrame_=fs->make<TH1F>("h_lepETauRestFrame_","lepton energy in tau rest frame",100,0.,2.5);
  h_lepEDifDiTauRestFrame_=fs->make<TH1F>("h_lepEDifDiTauRestFrame_","lep energy difference in tau rest frame",100,-100.,100.);
  h_lepESumDiTauRestFrame_=fs->make<TH1F>("h_lepESumDiTauRestFrame_","lep energy sum in tau rest frame",100,0.,200.);
  h_cosNegDiTau_=fs->make<TH1F>("h_cosNegDiTau_","cos(negtau) in ditau rest frame",50,-1.,1.);
  h_lepTauERatio_=fs->make<TH1F>("h_lepTauERatio_","lep to tau (Energy) ratio",50,0.,1.);
  h_diTauPL_=fs->make<TH1F>("h_diTauPL_","ditau PL",200,0.,1000.);
  h_diTauEta_=fs->make<TH1F>("h_diTauEta_","ditau Eta",100,-5.,5.);

  h_xVtxStrange_=fs->make<TH1F>("h_xVtxStrange_","x Vtx for big IP",200,-1.,1.);
  h_yVtxStrange_=fs->make<TH1F>("h_yVtxStrange_","y Vtx for big IP",200,-1.,1.);

  h_xVtxNormal_=fs->make<TH1F>("h_xVtxNormal_","x Vtx for norm IP",200,-1.,1.);
  h_yVtxNormal_=fs->make<TH1F>("h_yVtxNormal_","y Vtx for norm IP",200,-1.,1.);

  h_nHitsStrange_=fs->make<TH1F>("h_nHitsStrange_","n Hits for big IP",31,-0.5,30.5);
  h_nVtxHitsStrange_=fs->make<TH1F>("h_nVtxHitsStrange_","n Vtx Hits for big IP",21,-0.5,20.5);
  h_probTrkStrange_=fs->make<TH1F>("h_probTrkStrange_","chi2 trk prob for big IP",100,0.,1.);

  h_nHitsNormal_=fs->make<TH1F>("h_nHitsNormal_","n Hits for normal IP",31,-0.5,30.5);
  h_nVtxHitsNormal_=fs->make<TH1F>("h_nVtxHitsNormal_","n Vtx Hits for normal IP",21,-0.5,20.5);
  h_probTrkNormal_=fs->make<TH1F>("h_probTrkNormal_","chi2 trk prob for normal IP",100,0.,1.);

  h_pt_pfNeutral_ = fs->make<TH1F>("h_pt_pfNeutral_","",100,0.,10.);

  // ** declaration of trees ---->

  if (_doMC) {

    _genTree = fs->make<TTree>("Gen4Taus","Gen4Taus");

    _genTree->Branch("H2pt",&H2Pt,"H2pt/F");
    _genTree->Branch("H2px",&H2Px,"H2px/F");
    _genTree->Branch("H2py",&H2Py,"H2py/F");
    _genTree->Branch("H2pz",&H2Pz,"H2pz/F");
    _genTree->Branch("H2phi",&H2Phi,"H2phi/F");
    _genTree->Branch("H2eta",&H2Eta,"H2eta/F");
    _genTree->Branch("H2mass",&H2mass,"H2mass/F");

    _genTree->Branch("PtTau1",&tau1pt,"PtTau1/F");
    _genTree->Branch("PxTau1",&tau1px,"PxTau1/F");
    _genTree->Branch("PyTau1",&tau1py,"PyTau1/F");
    _genTree->Branch("PzTau1",&tau1pz,"PzTau1/F");
    _genTree->Branch("PhiTau1",&tau1phi,"PhiTau1/F");
    _genTree->Branch("EtaTau1",&tau1eta,"EtaTau1/F");
    _genTree->Branch("Tau1FinalStates",&FSO1,"Tau1FinalStates/I");
    _genTree->Branch("FinalStateIdTau1",finalstateIdTau1,"FinalStateIdTau1[Tau1FinalStates]/I");

    _genTree->Branch("PtTau2",&tau2pt,"PtTau2/F");
    _genTree->Branch("PxTau2",&tau2px,"PxTau2/F");
    _genTree->Branch("PyTau2",&tau2py,"PyTau2/F");
    _genTree->Branch("PzTau2",&tau2pz,"PzTau2/F");
    _genTree->Branch("PhiTau2",&tau2phi,"PhiTau2/F");
    _genTree->Branch("EtaTau2",&tau2eta,"EtaTau2/F");
    _genTree->Branch("Tau2FinalStates",&FSO2,"Tau2FinalStates/I");
    _genTree->Branch("FinalStateIdTau2",finalstateIdTau2,"FinalStateIdTau2[Tau2FinalStates]/I");

    _genTree->Branch("PtTau3",&tau3pt,"PtTau3/F");
    _genTree->Branch("PxTau3",&tau3px,"PxTau3/F");
    _genTree->Branch("PyTau3",&tau3py,"PyTau3/F");
    _genTree->Branch("PzTau3",&tau3pz,"PzTau3/F");
    _genTree->Branch("PhiTau3",&tau3phi,"PhiTau3/F");
    _genTree->Branch("EtaTau3",&tau3eta,"EtaTau3/F");
    _genTree->Branch("Tau3FinalStates",&FSO3,"Tau3FinalStates/I");
    _genTree->Branch("FinalStateIdTau3",finalstateIdTau3,"FinalStateIdTau3[Tau3FinalStates]/I");

    _genTree->Branch("PtTau4",&tau4pt,"PtTau4/F");
    _genTree->Branch("PxTau4",&tau4px,"PxTau4/F");
    _genTree->Branch("PyTau4",&tau4py,"PyTau4/F");
    _genTree->Branch("PzTau4",&tau4pz,"PzTau4/F");
    _genTree->Branch("PhiTau4",&tau4phi,"PhiTau4/F");
    _genTree->Branch("EtaTau4",&tau4eta,"EtaTau4/F");
    _genTree->Branch("Tau4FinalStates",&FSO4,"Tau4FinalStates/I");
    _genTree->Branch("FinalStateIdTau4",finalstateIdTau4,"FinalStateIdTau4[Tau4FinalStates]/I");
    
    _genTree->Branch("PtFirstH1",&ptFirstH1,"PtFirstH1/F");
    _genTree->Branch("PxFirstH1",&pxFirstH1,"PxFirstH1/F");
    _genTree->Branch("PyFirstH1",&pyFirstH1,"PyFirstH1/F");
    _genTree->Branch("PzFirstH1",&pzFirstH1,"PzFirstH1/F");
    _genTree->Branch("PhiFirstH1",&phiFirstH1,"PhiFirstH1/F");
    _genTree->Branch("EtaFirstH1",&etaFirstH1,"EtaFirstH1/F");
    _genTree->Branch("MassFirstH1",&massFirstH1,"MassFirstH1/F");

    _genTree->Branch("PtSecondH1",&ptSecondH1,"PtSecondH1/F");
    _genTree->Branch("PxSecondH1",&pxSecondH1,"PxSecondH1/F");
    _genTree->Branch("PySecondH1",&pySecondH1,"PySecondH1/F");
    _genTree->Branch("PzSecondH1",&pzSecondH1,"PzSecondH1/F");
    _genTree->Branch("PhiSecondH1",&phiSecondH1,"PhiSecondH1/F");
    _genTree->Branch("EtaSecondH1",&etaSecondH1,"EtaSecondH1/F");
    _genTree->Branch("MassSecondH1",&massSecondH1,"MassSecondH1/F");

    _genTree->Branch("NPartFirstH1",&nPartFirstH1,"NPartFirstH1/I");
    _genTree->Branch("IdPartFirstH1",pdgIdPartFirstH1,"IdPartFirstH1[NPartFirstH1]/I");
    _genTree->Branch("QPartFirstH1", qPartFirstH1, "QPartFirstH1[NPartFirstH1]/I");
    _genTree->Branch("PxPartFirstH1",pxPartFirstH1,"PxPartFirstH1[NPartFirstH1]/F");
    _genTree->Branch("PyPartFirstH1",pyPartFirstH1,"PyPartFirstH1[NPartFirstH1]/F");
    _genTree->Branch("PzPartFirstH1",pzPartFirstH1,"PzPartFirstH1[NPartFirstH1]/F");
    _genTree->Branch("PtPartFirstH1",ptPartFirstH1,"PtPartFirstH1[NPartFirstH1]/F");
    _genTree->Branch("EtaPartFirstH1",etaPartFirstH1,"EtaPartFirstH1[NPartFirstH1]/F");
    _genTree->Branch("PhiPartFirstH1",phiPartFirstH1,"PhiPartFirstH1[NPartFirstH1]/F");
    
    _genTree->Branch("NPartAroundMuFirstH1",&nPartAroundMuFirstH1,"NPartAroundMuFirstH1/I");
    _genTree->Branch("IdPartAroundMuFirstH1",pdgIdPartAroundMuFirstH1,"IdPartAroundMuFirstH1[NPartAroundMuFirstH1]/I");
    _genTree->Branch("QPartAroundMuFirstH1", qPartAroundMuFirstH1, "QPartAroundMuFirstH1[NPartAroundMuFirstH1]/I");
    _genTree->Branch("PxPartAroundMuFirstH1",pxPartAroundMuFirstH1,"PxPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");
    _genTree->Branch("PyPartAroundMuFirstH1",pyPartAroundMuFirstH1,"PyPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");
    _genTree->Branch("PzPartAroundMuFirstH1",pzPartAroundMuFirstH1,"PzPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");
    _genTree->Branch("PtPartAroundMuFirstH1",ptPartAroundMuFirstH1,"PtPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");
    _genTree->Branch("EtaPartAroundMuFirstH1",etaPartAroundMuFirstH1,"EtaPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");
    _genTree->Branch("PhiPartAroundMuFirstH1",phiPartAroundMuFirstH1,"PhiPartAroundMuFirstH1[NPartAroundMuFirstH1]/F");

    _genTree->Branch("NPartSecondH1",&nPartSecondH1,"NPartSecondH1/I");
    _genTree->Branch("IdPartSecondH1",pdgIdPartSecondH1,"IdPartSecondH1[NPartSecondH1]/I");
    _genTree->Branch("QPartSecondH1", qPartSecondH1, "QPartSecondH1[NPartSecondH1]/I");
    _genTree->Branch("PxPartSecondH1",pxPartSecondH1,"PxPartSecondH1[NPartSecondH1]/F");
    _genTree->Branch("PyPartSecondH1",pyPartSecondH1,"PyPartSecondH1[NPartSecondH1]/F");
    _genTree->Branch("PzPartSecondH1",pzPartSecondH1,"PzPartSecondH1[NPartSecondH1]/F");
    _genTree->Branch("PtPartSecondH1",ptPartSecondH1,"PtPartSecondH1[NPartSecondH1]/F");
    _genTree->Branch("EtaPartSecondH1",etaPartSecondH1,"EtaPartSecondH1[NPartSecondH1]/F");
    _genTree->Branch("PhiPartSecondH1",phiPartSecondH1,"PhiPartSecondH1[NPartSecondH1]/F");

    _genTree->Branch("NPartAroundMuSecondH1",&nPartAroundMuSecondH1,"NPartAroundMuSecondH1/I");
    _genTree->Branch("IdPartAroundMuSecondH1",pdgIdPartAroundMuSecondH1,"IdPartAroundMuSecondH1[NPartAroundMuSecondH1]/I");
    _genTree->Branch("QPartAroundMuSecondH1", qPartAroundMuSecondH1, "QPartAroundMuSecondH1[NPartAroundMuSecondH1]/I");
    _genTree->Branch("PxPartAroundMuSecondH1",pxPartAroundMuSecondH1,"PxPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");
    _genTree->Branch("PyPartAroundMuSecondH1",pyPartAroundMuSecondH1,"PyPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");
    _genTree->Branch("PzPartAroundMuSecondH1",pzPartAroundMuSecondH1,"PzPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");
    _genTree->Branch("PtPartAroundMuSecondH1",ptPartAroundMuSecondH1,"PtPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");
    _genTree->Branch("EtaPartAroundMuSecondH1",etaPartAroundMuSecondH1,"EtaPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");
    _genTree->Branch("PhiPartAroundMuSecondH1",phiPartAroundMuSecondH1,"PhiPartAroundMuSecondH1[NPartAroundMuSecondH1]/F");

  }

  _tree = fs->make<TTree>("HTauTau2Muons","HTauTau2Muons");

  if (_doMC) {

    _tree->Branch("H2pt",&H2Pt,"H2pt/F");
    _tree->Branch("H2px",&H2Px,"H2px/F");
    _tree->Branch("H2py",&H2Py,"H2py/F");
    _tree->Branch("H2pz",&H2Pz,"H2pz/F");
    _tree->Branch("H2phi",&H2Phi,"H2phi/F");
    _tree->Branch("H2eta",&H2Eta,"H2eta/F");
    _tree->Branch("H2mass",&H2mass,"H2mass/F");

    _tree->Branch("PtFirstH1",&ptFirstH1,"PtFirstH1/F");
    _tree->Branch("PxFirstH1",&pxFirstH1,"PxFirstH1/F");
    _tree->Branch("PyFirstH1",&pyFirstH1,"PyFirstH1/F");
    _tree->Branch("PzFirstH1",&pzFirstH1,"PzFirstH1/F");
    _tree->Branch("PhiFirstH1",&phiFirstH1,"PhiFirstH1/F");
    _tree->Branch("EtaFirstH1",&etaFirstH1,"EtaFirstH1/F");
    _tree->Branch("MassFirstH1",&massFirstH1,"MassFirstH1/F");

    _tree->Branch("PtSecondH1",&ptSecondH1,"PtSecondH1/F");
    _tree->Branch("PxSecondH1",&pxSecondH1,"PxSecondH1/F");
    _tree->Branch("PySecondH1",&pySecondH1,"PySecondH1/F");
    _tree->Branch("PzSecondH1",&pzSecondH1,"PzSecondH1/F");
    _tree->Branch("PhiSecondH1",&phiSecondH1,"PhiSecondH1/F");
    _tree->Branch("EtaSecondH1",&etaSecondH1,"EtaSecondH1/F");
    _tree->Branch("MassSecondH1",&massSecondH1,"MassSecondH1/F");

    _tree->Branch("NPartFirstH1",&nPartFirstH1,"NPartFirstH1/I");
    _tree->Branch("IdPartFirstH1",pdgIdPartFirstH1,"IdPartFirstH1[NPartFirstH1]/I");
    _tree->Branch("QPartFirstH1", qPartFirstH1, "QPartFirstH1[NPartFirstH1]/I");
    _tree->Branch("PxPartFirstH1",pxPartFirstH1,"PxPartFirstH1[NPartFirstH1]/F");
    _tree->Branch("PyPartFirstH1",pyPartFirstH1,"PyPartFirstH1[NPartFirstH1]/F");
    _tree->Branch("PzPartFirstH1",pzPartFirstH1,"PzPartFirstH1[NPartFirstH1]/F");
    _tree->Branch("PtPartFirstH1",ptPartFirstH1,"PtPartFirstH1[NPartFirstH1]/F");
    _tree->Branch("EtaPartFirstH1",etaPartFirstH1,"EtaPartFirstH1[NPartFirstH1]/F");
    _tree->Branch("PhiPartFirstH1",phiPartFirstH1,"PhiPartFirstH1[NPartFirstH1]/F");

    _tree->Branch("NPartSecondH1",&nPartSecondH1,"NPartSecondH1/I");
    _tree->Branch("IdPartSecondH1",pdgIdPartSecondH1,"IdPartSecondH1[NPartSecondH1]/I");
    _tree->Branch("QPartSecondH1", qPartSecondH1, "QPartSecondH1[NPartSecondH1]/I");
    _tree->Branch("PxPartSecondH1",pxPartSecondH1,"PxPartSecondH1[NPartSecondH1]/F");
    _tree->Branch("PyPartSecondH1",pyPartSecondH1,"PyPartSecondH1[NPartSecondH1]/F");
    _tree->Branch("PzPartSecondH1",pzPartSecondH1,"PzPartSecondH1[NPartSecondH1]/F");
    _tree->Branch("PtPartSecondH1",ptPartSecondH1,"PtPartSecondH1[NPartSecondH1]/F");
    _tree->Branch("EtaPartSecondH1",etaPartSecondH1,"EtaPartSecondH1[NPartSecondH1]/F");
    _tree->Branch("PhiPartSecondH1",phiPartSecondH1,"PhiPartSecondH1[NPartSecondH1]/F");

  }


  //ntpGeneralInfo.beginJob(_tree, fs);
  ntpTriggerInfo.beginJob(_tree, fs);
  ntpPileUpInfo.beginJob(_tree, fs);
  ntpVertexInfo.beginJob(_tree, fs);
  std::cout <<"Trees Ok?"<< std::endl;

  //Tracks from PV with highest pt2 sum
  _tree->Branch("NPVtracks",&nPVtracks,"NPVtracks/I");
  _tree->Branch("PtTrack",ptTrack,"PtTrack[NPVtracks]/F");
  _tree->Branch("PxTrack",pxTrack,"PxTrack[NPVtracks]/F");
  _tree->Branch("PyTrack",pyTrack,"PyTrack[NPVtracks]/F");
  _tree->Branch("PzTrack",pzTrack,"PzTrack[NPVtracks]/F");
  _tree->Branch("PhiTrack",phiTrack,"PhiPVTrack[NPVtracks]/F");
  _tree->Branch("EtaTrack",etaTrack,"EtaTrack[NPVtracks]/F");
  _tree->Branch("ChargeTrack",chargeTrack,"ChargeTrack[NPVtracks]/I");
  _tree->Branch("Chi2Track",chi2Track,"Chi2Track[NPVtracks]/D");

  _tree->Branch("NPVPostracks",&nPVPostracks,"NPVPostracks/I");
  _tree->Branch("PtPVTrackPos",ptPVTracksPos,"PtPVTrackPos[NPVPostracks]/F");
  _tree->Branch("PhiPVTrackPos",phiPVTracksPos,"PhiPVTrackPos[NPVPostracks]/F");
  _tree->Branch("EtaPVTrackPos",etaPVTracksPos,"EtaPVTrackPos[NPVPostracks]/F");
  _tree->Branch("PxPVTrackPos",pxPVTracksPos,"PxPVTrackPos[NPVPostracks]/F");
  _tree->Branch("PyPVTrackPos",pyPVTracksPos,"PyPVTrackPos[NPVPostracks]/F");
  _tree->Branch("PzPVTrackPos",pzPVTracksPos,"PzPVTrackPos[NPVPostracks]/F");
  _tree->Branch("ChargePVTracksPos",qPVTracksPos,"ChargePVTracksPos[NPVPostracks]/F");

  _tree->Branch("NPVNegtracks",&nPVNegtracks,"NPVNegtracks/I");
  _tree->Branch("PtPVTrackNeg",ptPVTracksNeg,"PtPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("PhiPVTrackNeg",phiPVTracksNeg,"PhiPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("EtaPVTrackNeg",etaPVTracksNeg,"EtaPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("PxPVTrackNeg",pxPVTracksNeg,"PxPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("PyPVTrackNeg",pyPVTracksNeg,"PyPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("PzPVTrackNeg",pzPVTracksNeg,"PzPVTrackNeg[NPVNegtracks]/F");
  _tree->Branch("ChargePVTracksNeg",qPVTracksNeg,"ChargePVTracksNeg[NPVNegtracks]/F");
  
   
  // Dilepton kinematics
  _tree->Branch("DiLeptonMass",&_diLepMassVar,"DiLeptonMass/F");
  _tree->Branch("MET",&_met,"MET/F");
  _tree->Branch("METPx",&_metPx,"METPx/F");
  _tree->Branch("METPy",&_metPy,"METPy/F");
  _tree->Branch("METPz",&_metPz,"METPz/F");
  _tree->Branch("METEn",&_metEn,"METEn/F");
  /*_tree->Branch("METCovXX",&_metCovXX,"METCovXX/F");
  _tree->Branch("METCovXY",&_metCovXY,"METCovXY/F");
  _tree->Branch("METCovYX",&_metCovYX,"METCovYX/F");
  _tree->Branch("METCovYY",&_metCovYY,"METCovYY/F");

  _tree->Branch("METMVAPx",&_metMVAPx,"METMVAPx/F");
  _tree->Branch("METMVAPy",&_metMVAPy,"METMVAPy/F");
  _tree->Branch("METMVAPz",&_metMVAPz,"METMVAPz/F");
  _tree->Branch("METMVAEn",&_metMVAEn,"METMVAEn/F");
  _tree->Branch("METMVACovXX",&_metMVACovXX,"METMVACovXX/F");
  _tree->Branch("METMVACovXY",&_metMVACovXY,"METMVACovXY/F");
  _tree->Branch("METMVACovYX",&_metMVACovYX,"METMVACovYX/F");
  _tree->Branch("METMVACovYY",&_metMVACovYY,"METMVACovYY/F");*/

  _tree->Branch("DiLepDPhi",&_diLepDPhi,"DiLepDPhi/F");
  _tree->Branch("DiLepDEta",&_diLepDEta,"DiLepDEta/F");
  _tree->Branch("DiLepDR",&_diLepDR,"DiLepDR/F");
  _tree->Branch("PosMETDPhi",&_posLepMETDPhi,"PosMETDPhi/F");
  _tree->Branch("NegMETDPhi",&_negLepMETDPhi,"NegMETDPhi/F");
  _tree->Branch("DiLepEta",&_diLepEta,"DiLepEta/F");
  _tree->Branch("NegLepPt",&_negLepPt,"NegLepPt/F");
  _tree->Branch("PosLepPt",&_posLepPt,"PosLepPt/F");
  _tree->Branch("NegLepPx",&_negLepPx,"NegLepPx/F");
  _tree->Branch("PosLepPx",&_posLepPx,"PosLepPx/F");
  _tree->Branch("NegLepPy",&_negLepPy,"NegLepPy/F");
  _tree->Branch("PosLepPy",&_posLepPy,"PosLepPy/F");
  _tree->Branch("NegLepPz",&_negLepPz,"NegLepPz/F");
  _tree->Branch("PosLepPz",&_posLepPz,"PosLepPz/F");
  _tree->Branch("NegLepEta",&_negLepEta,"NegLepEta/F");
  _tree->Branch("PosLepEta",&_posLepEta,"PosLepEta/F");
  _tree->Branch("NegLepPhi",&_negLepPhi,"NegLepPhi/F");
  _tree->Branch("PosLepPhi",&_posLepPhi,"PosLepPhi/F");
  _tree->Branch("NegLepEn",&_negLepEn,"NegLepEn/F");
  _tree->Branch("PosLepEn",&_posLepEn,"PosLepEn/F");
  _tree->Branch("NegLepQ",&_negLepQ,"NegLepQ/F");
  _tree->Branch("PosLepQ",&_posLepQ,"PosLepQ/F");

  _tree->Branch("NegLep_IsoMu24",&_negLep_IsoMu24,"NegLep_IsoMu24/O");
  _tree->Branch("PosLep_IsoMu24",&_posLep_IsoMu24,"PosLep_IsoMu24/O");
  _tree->Branch("NegLep_Mu20",&_negLep_Mu20,"NegLep_Mu20/O");
  _tree->Branch("PosLep_Mu20",&_posLep_Mu20,"PosLep_Mu20/O");
  _tree->Branch("NegLep_Mu30",&_negLep_Mu30,"NegLep_Mu30/O");
  _tree->Branch("PosLep_Mu30",&_posLep_Mu30,"PosLep_Mu30/O");
  _tree->Branch("NegLep_Mu40",&_negLep_Mu40,"NegLep_Mu40/O");
  _tree->Branch("PosLep_Mu40",&_posLep_Mu40,"PosLep_Mu40/O");
  _tree->Branch("NegLep_Mu17_Mu8",&_negLep_Mu17_Mu8,"NegLep_Mu17_Mu8/O");
  _tree->Branch("PosLep_Mu17_Mu8",&_posLep_Mu17_Mu8,"PosLep_Mu17_Mu8/O");
  _tree->Branch("NegLep_Mu17_TkMu8",&_negLep_Mu17_TkMu8,"NegLep_Mu17_TkMu8/O");
  _tree->Branch("PosLep_Mu17_TkMu8",&_posLep_Mu17_TkMu8,"PosLep_Mu17_TkMu8/O");

  _tree->Branch("DiLepPt",&_diLepPt,"DiLepPt/F");
  _tree->Branch("DiLepPhi",&_diLepPhi,"DiLepPhi/F");
 
  //PfNoPileUp objects 
  _tree->Branch("NpfNoPU",&npfNoPU,"NpfNoPU/I");
  _tree->Branch("PxpfNoPU",pxpfNoPU,"PxpfNoPU[NpfNoPU]/F");
  _tree->Branch("PypfNoPU",pypfNoPU,"PypfNoPU[NpfNoPU]/F");
  _tree->Branch("PzpfNoPU",pzpfNoPU,"PzpfNoPU[NpfNoPU]/F");
  _tree->Branch("PtpfNoPU",ptpfNoPU,"PtpfNoPU[NpfNoPU]/F");
  _tree->Branch("EnpfNoPU",enpfNoPU,"EnpfNoPU[NpfNoPU]/F");
  _tree->Branch("EtapfNoPU",etapfNoPU,"EtapfNoPU[NpfNoPU]/F");
  _tree->Branch("ChargepfNoPU",chargepfNoPU,"ChargepfNoPU[NpfNoPU]/F");
  _tree->Branch("IppfNoPU",ippfNoPU,"IppfNoPU[NpfNoPU]/F");
  _tree->Branch("IpSigpfNoPU",ipSigpfNoPU,"IpSigpfNoPU[NpfNoPU]/F");
  _tree->Branch("IpZpfNoPU",ipZpfNoPU,"IpZpfNoPU[NpfNoPU]/F");
  _tree->Branch("IpZSigpfNoPU",ipZSigpfNoPU,"IpZSigpfNoPU[NpfNoPU]/F");
  _tree->Branch("Ip3DpfNoPU",ip3DpfNoPU,"Ip3DpfNoPU[NpfNoPU]/F");
  _tree->Branch("Ip3DSigpfNoPU",ip3DSigpfNoPU,"Ip3DSigpfNoPU[NpfNoPU]/F");
  _tree->Branch("IdpfNoPU",idpfNoPU,"IdpfNoPU[NpfNoPU]/I");

  // PfNoPileUp objects around muon tracks
  _tree->Branch("NTracksPos",&ntracksPos,"NTracksPos/I");
  _tree->Branch("PxTrackPos",pxTrackPos ,"PxTrackPos[NTracksPos]/F");
  _tree->Branch("PyTrackPos",pyTrackPos ,"PyTrackPos[NTracksPos]/F");
  _tree->Branch("PzTrackPos",pzTrackPos ,"PzTrackPos[NTracksPos]/F");
  _tree->Branch("PtTrackPos",ptTrackPos ,"PtTrackPos[NTracksPos]/F");
  _tree->Branch("PhiTrackPos",phiTrackPos ,"PhiTrackPos[NTracksPos]/F");
  _tree->Branch("EnergyTrackPos",enTrackPos ,"EnergyTrackPos[NTracksPos]/F");
  _tree->Branch("EtaTrackPos",etaTrackPos ,"EtaTrackPos[NTracksPos]/F");
  _tree->Branch("QTrackPos", qTrackPos  , "QTrackPos[NTracksPos]/F");
  _tree->Branch("IdTrackPos",idTrackPos ,"IdTrackPos[NTracksPos]/I");
  _tree->Branch("IpTrackPos",ipTrackPos ,"IpTrackPos[NTracksPos]/F");
  _tree->Branch("IpSigTrackPos",ipSigTrackPos ,"IpSigTrackPos[NTracksPos]/F");
  _tree->Branch("IpZTrackPos",ipZTrackPos ,"IpZTrackPos[NTracksPos]/F");
  _tree->Branch("IpZSigTrackPos",ipZSigTrackPos ,"IpZSigTrackPos[NTracksPos]/F");
  _tree->Branch("Ip3DTrackPos",ip3DTrackPos ,"Ip3DTrackPos[NTracksPos]/F");
  _tree->Branch("Ip3DSigTrackPos",ip3DSigTrackPos ,"Ip3DSigTrackPos[NTracksPos]/F");

  _tree->Branch("NTracksNeg",&ntracksNeg,"NTracksNeg/I");
  _tree->Branch("PxTrackNeg",pxTrackNeg ,"PxTrackNeg[NTracksNeg]/F");
  _tree->Branch("PyTrackNeg",pyTrackNeg ,"PyTrackNeg[NTracksNeg]/F");
  _tree->Branch("PzTrackNeg",pzTrackNeg ,"PzTrackNeg[NTracksNeg]/F");
  _tree->Branch("PtTrackNeg",ptTrackNeg ,"PtTrackNeg[NTracksNeg]/F");
  _tree->Branch("PhiTrackNeg",phiTrackNeg ,"PhiTrackNeg[NTracksNeg]/F");
  _tree->Branch("EtaTrackNeg",etaTrackNeg ,"EtaTrackNeg[NTracksNeg]/F");
  _tree->Branch("EnergyTrackNeg",enTrackNeg ,"EnergyTrackNeg[NTracksNeg]/F");
  _tree->Branch("QTrackNeg",qTrackNeg   , "QTrackNeg[NTracksNeg]/F");
  _tree->Branch("IdTrackNeg",idTrackNeg ,"IdTrackNeg[NTracksNeg]/I");
  _tree->Branch("IpTrackNeg",ipTrackNeg ,"IpTrackNeg[NTracksNeg]/F");
  _tree->Branch("IpSigTrackNeg",ipSigTrackNeg ,"IpSigTrackNeg[NTracksNeg]/F");
  _tree->Branch("IpZTrackNeg",ipZTrackNeg ,"IpZTrackNeg[NTracksNeg]/F");
  _tree->Branch("IpZSigTrackNeg",ipZSigTrackNeg ,"IpZSigTrackNeg[NTracksNeg]/F");
  _tree->Branch("Ip3DTrackNeg",ip3DTrackNeg ,"Ip3DTrackNeg[NTracksNeg]/F");
  _tree->Branch("Ip3DSigTrackNeg",ip3DSigTrackNeg ,"Ip3DSigTrackNeg[NTracksNeg]/F");


  //PfPileUp objects around muon tracks
  _tree->Branch("NTracksPileUpPos",&ntracksPileUpPos,"NTracksPileUpPos/I");
  _tree->Branch("PxTrackPileUpPos",pxTrackPileUpPos ,"PxTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("PyTrackPileUpPos",pyTrackPileUpPos ,"PyTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("PzTrackPileUpPos",pzTrackPileUpPos ,"PzTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("PtTrackPileUpPos",ptTrackPileUpPos ,"PtTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("PhiTrackPileUpPos",phiTrackPileUpPos ,"PhiTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("EnergyTrackPileUpPos",enTrackPileUpPos ,"EnergyTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("EtaTrackPileUpPos",etaTrackPileUpPos ,"EtaTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("QTrackPileUpPos", qTrackPileUpPos  , "QTrackPileUpPos[NTracksPileUpPos]/F");
  _tree->Branch("IdTrackPileUpPos",idTrackPileUpPos ,"IdTrackPileUpPos[NTracksPileUpPos]/I");

  _tree->Branch("NTracksPileUpNeg",&ntracksPileUpNeg,"NTracksPileUpNeg/I");
  _tree->Branch("PxTrackPileUpNeg",pxTrackPileUpNeg ,"PxTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("PyTrackPileUpNeg",pyTrackPileUpNeg ,"PyTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("PzTrackPileUpNeg",pzTrackPileUpNeg ,"PzTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("PtTrackPileUpNeg",ptTrackPileUpNeg ,"PtTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("PhiTrackPileUpNeg",phiTrackPileUpNeg ,"PhiTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("EtaTrackPileUpNeg",etaTrackPileUpNeg ,"EtaTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("EnergyTrackPileUpNeg",enTrackPileUpNeg ,"EnergyTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("QTrackPileUpNeg",qTrackPileUpNeg   , "QTrackPileUpNeg[NTracksPileUpNeg]/F");
  _tree->Branch("IdTrackPileUpNeg",idTrackPileUpNeg ,"IdTrackPileUpNeg[NTracksPileUpNeg]/I");

  // Original muons in embedding
  _tree->Branch("EmbedWeight", &_embedWeight, "EmbedWeight/F");
  /*if(!_origMuonSrc.label().empty())
  {
    _tree->Branch("OrigDiLepMass", &_origDiLepMass, "OrigDiLepMass/F");
    _tree->Branch("OrigNegLepPt", &_origNegLepPt, "OrigNegLepPt/F");
    _tree->Branch("OrigPosLepPt", &_origPosLepPt, "OrigPosLepPt/F");
    _tree->Branch("OrigNegLepEta", &_origNegLepEta, "OrigNegLepEta/F");
    _tree->Branch("OrigPosLepEta", &_origPosLepEta, "OrigPosLepEta/F");
    _tree->Branch("OrigNegLepPhi", &_origNegLepPhi, "OrigNegLepPhi/F");
    _tree->Branch("OrigPosLepPhi", &_origPosLepPhi, "OrigPosLepPhi/F");
    _tree->Branch("OrigNegLepChargedHadronIso", &_origNegLepChargedHadronIso, "OrigNegLepChargedHadronIso/F");
    _tree->Branch("OrigPosLepChargedHadronIso", &_origPosLepChargedHadronIso, "OrigPosLepChargedHadronIso/F");
  }*/

  // Lepton isolation
  _tree->Branch("IsoNeg",&_isoNeg,"IsoNeg/F");
  _tree->Branch("IsoPos",&_isoPos,"IsoPos/F");
  _tree->Branch("IsoTrkNeg",&_isoTrkNeg,"IsoTrkNeg/F");
  _tree->Branch("IsoTrkPos",&_isoTrkPos,"IsoTrkPos/F");
  _tree->Branch("IsoECalNeg",&_isoECalNeg,"IsoECalNeg/F");
  _tree->Branch("IsoECalPos",&_isoECalPos,"IsoECalPos/F");
  _tree->Branch("IsoHCalNeg",&_isoHCalNeg,"IsoHCalNeg/F");
  _tree->Branch("IsoHCalPos",&_isoHCalPos,"IsoHCalPos/F");
  _tree->Branch("IsoPFsNeg",&_isoPFsNeg,"IsoPFsNeg/F");
  _tree->Branch("IsoPFsPos",&_isoPFsPos,"IsoPFsPos/F");

  _tree->Branch("IsoChargedHadPFsNeg",&_isoChargedHadPFsNeg,"IsoChargedHadPFsNeg/F");
  _tree->Branch("IsoChargedHadPFsPos",&_isoChargedHadPFsPos,"IsoChargedHadPFsPos/F");
  _tree->Branch("IsoNeutralHadPFsNeg",&_isoNeutralHadPFsNeg,"IsoNeutralHadPFsNeg/F");
  _tree->Branch("IsoNeutralHadPFsPos",&_isoNeutralHadPFsPos,"IsoNeutralHadPFsPos/F");
  _tree->Branch("IsoPhotonsPFsNeg",&_isoPhotonsPFsNeg,"IsoPhotonsPFsNeg/F");
  _tree->Branch("IsoPhotonsPFsPos",&_isoPhotonsPFsPos,"IsoPhotonsPFsPos/F");
  

  _tree->Branch("isGlobalMuPos",&_isPosGlobalMu,"isGlobalMuPos/O");
  _tree->Branch("isGlobalMuNeg",&_isNegGlobalMu,"isGlobalMuNeg/O");

  _tree->Branch("isTrackerMuPos",&_isPosTrackerMu,"isTrackerMuPos/O");
  _tree->Branch("isTrackerMuNeg",&_isNegTrackerMu,"isTrackerMuNeg/O");

  _tree->Branch("isPFMuPos",&_isPosPFMu,"isPFMuPos/O");
  _tree->Branch("isPFMuNeg",&_isNegPFMu,"isPFMuNeg/O");

  _tree->Branch("nPixelHitsPos",&_nPixelHitsPos,"nPixelHitsPos/I");
  _tree->Branch("nPixelHitsNeg",&_nPixelHitsNeg,"nPixelHitsNeg/I");

  _tree->Branch("nTrackerHitsPos",&_nTrackerHitsPos,"nTrackerHitsPos/I");
  _tree->Branch("nTrackerHitsNeg",&_nTrackerHitsNeg,"nTrackerHitsNeg/I");

  _tree->Branch("nMuonHitsPos",&_nMuonHitsPos,"nMuonHitsPos/I");
  _tree->Branch("nMuonHitsNeg",&_nMuonHitsNeg,"nMuonHitsNeg/I");

  _tree->Branch("nMuonStationsPos",&_nMuonStationsPos,"nMuonStationsPos/I");
  _tree->Branch("nMuonStationsNeg",&_nMuonStationsNeg,"nMuonStationsNeg/I");

  _tree->Branch("Chi2PosMu",&_chi2PosMu,"Chi2PosMu/F");
  _tree->Branch("Chi2NegMu",&_chi2NegMu,"Chi2NegMu/F");

  _tree->Branch("NdofPosMu",&_ndofPosMu,"NdofPosMu/F");
  _tree->Branch("NdofNegMu",&_ndofNegMu,"NdofNegMu/F");

  // impact parameter significance
  _tree->Branch("IpSig2DPos",_ipSig2DPos,"IpSig2DPos[nPV]/F");
  _tree->Branch("IpSig2DNeg",_ipSig2DNeg,"IpSig2DNeg[nPV]/F");
  _tree->Branch("IpSig3DPos",_ipSig3DPos,"IpSig3DPos[nPV]/F");
  _tree->Branch("IpSig3DNeg",_ipSig3DNeg,"IpSig3DNeg[nPV]/F");
  _tree->Branch("IpSigZPos",_ipSigZPos,"IpSigZPos[nPV]/F");
  _tree->Branch("IpSigZNeg",_ipSigZNeg,"IpSigZNeg[nPV]/F");  
  
  // impact parameters
  _tree->Branch("Ip2DPos",_ip2DPos,"Ip2DPos[nPV]/F");
  _tree->Branch("Ip2DNeg",_ip2DNeg,"Ip2DNeg[nPV]/F");
  _tree->Branch("Ip3DPos",_ip3DPos,"Ip3DPos[nPV]/F");
  _tree->Branch("Ip3DNeg",_ip3DNeg,"Ip3DNeg[nPV]/F");
  _tree->Branch("IpZPos",_ipZPos,"IpZPos[nPV]/F");
  _tree->Branch("IpZNeg",_ipZNeg,"IpZNeg[nPV]/F");  
  
  // distance between muon tracks
  _tree->Branch("TwoMuonDist3D",&_twoMuonDist3D,"TwoMuonDistance3D/F");
  _tree->Branch("TwoMuonDist3DE",&_twoMuonDist3DE,"TwoMuonDistance3DE/F");
  _tree->Branch("TwoMuonDist2D",&_twoMuonDist2D,"TwoMuonDistance2D/F");
  _tree->Branch("TwoMuonDist2DE",&_twoMuonDist2DE,"TwoMuonDistance2DE/F");
  _tree->Branch("TwoMuonDistRPhi3D",&_twoMuonDistRPhi3D,"TwoMuonDistanceRPhi3D/F");
  _tree->Branch("TwoMuonDistRPhi3DE",&_twoMuonDistRPhi3DE,"TwoMuonDistanceRPhi3DE/F");
  _tree->Branch("TwoMuonDistRPhi2D",&_twoMuonDistRPhi2D,"TwoMuonDistanceRPhi2D/F");
  _tree->Branch("TwoMuonDistRPhi2DE",&_twoMuonDistRPhi2DE,"TwoMuonDistanceRPhi2DE/F");

  // jets
  _tree->Branch("noJets",&noJets,"noJets/I"); 
  _tree->Branch("noRecoJets",&noRecoJets,"noRecoJets/I"); 
  _tree->Branch("jetIDLoose",jetIDLoose,"jetIDLoose[noJets]/O");
  _tree->Branch("jetIDMedium",jetIDMedium,"jetIDMedium[noJets]/O");
  _tree->Branch("jetIDTight",jetIDTight,"jetIDTight[noJets]/O");
  _tree->Branch("jetPt",jetPt,"jetPt[noJets]/F");
  _tree->Branch("jetEt",jetEt,"jetEt[noJets]/F");
  _tree->Branch("jetPhi",jetPhi,"jetPhi[noJets]/F");
  _tree->Branch("jetEta",jetEta,"jetEta[noJets]/F");
  _tree->Branch("jetPx",jetPx,"jetPx[noJets]/F");
  _tree->Branch("jetPy",jetPy,"jetPy[noJets]/F");
  _tree->Branch("jetPz",jetPz,"jetPz[noJets]/F");
  _tree->Branch("jetEn",jetEn,"jetEn[noJets]/F");
  _tree->Branch("jetMass",jetMass,"jetMass[noJets]/F");
  _tree->Branch("jetArea",jetArea,"jetArea[noJets]/F");
  _tree->Branch("jetCorrFactor",jetCorrFactor,"jetCorrFactor[noJets]/F");

  if (_doMC) {

    _tree->Branch("genJetExist",genJetExist,"genJetExist[noJets]/O");
    _tree->Branch("genJetPx",genJetPx,"genJetPx[noJets]/F");
    _tree->Branch("genJetPy",genJetPy,"genJetPy[noJets]/F");
    _tree->Branch("genJetPz",genJetPz,"genJetPz[noJets]/F");
    _tree->Branch("genJetEn",genJetEn,"genJetEn[noJets]/F");
    _tree->Branch("genJetCollisionId",genJetCollisionId,"genJetCollisionId[noJets]/I");

    _tree->Branch("genPartonExist",genPartonExist,"genPartonExist[noJets]/O");
    _tree->Branch("genPartonPdgId",genPartonPdg,"genPartonPdgId[noJets]/I");
    _tree->Branch("genPartonPx",genPartonPx,"genPartonPx[noJets]/F");
    _tree->Branch("genPartonPy",genPartonPy,"genPartonPy[noJets]/F");
    _tree->Branch("genPartonPz",genPartonPz,"genPartonPz[noJets]/F");
    _tree->Branch("genPartonEn",genPartonEn,"genPartonEn[noJets]/F");
  
  }

  _tree->Branch("jetBProbBJetTag",jetBProbBJetTag,"jetBProbBJetTag[noJets]/F");
  _tree->Branch("jetProbBJetTag",jetProbBJetTag,"jetProbBJetTag[noJets]/F");
  _tree->Branch("tcHPBJetTag",tcHPBJetTag,"tcHPBJetTag[noJets]/F");
  _tree->Branch("tcHEBJetTag",tcHEBJetTag,"tcHEBJetTag[noJets]/F");
  _tree->Branch("svHEBJetTag",svHEBJetTag,"svHEBJetTag[noJets]/F");
  _tree->Branch("svHPBJetTag",svHPBJetTag,"svHPBJetTag[noJets]/F");
  _tree->Branch("combSVBJetTag",combSVBJetTag,"combSVBJetTag[noJets]/F");
  _tree->Branch("combSVMVABJetTag",combSVMVABJetTag,"combSVMVABJetTag[noJets]/F");

  // jet constituents ===>
  _tree->Branch("noConstituents",noConstituents,"noConstituents[noJets]/I");
  _tree->Branch("neuEMEnergyFr",neuEMEnergyFr,"neuEMEnergyFr[noJets]/F");
  _tree->Branch("neuHadronEnergyFr",neuHadronEnergyFr,"neuHadronEnergyFr[noJets]/F");
  _tree->Branch("neuHadronMultiplicity",neuHadronMultiplicity,"neuHadronMultiplicity[noJets]/I");
  _tree->Branch("neuMultiplicity",neuMultiplicity,"neuMultiplicity[noJets]/I");
  _tree->Branch("photonEnergyFr",photonEnergyFr,"photonEnergyFr[noJets]/F");
  _tree->Branch("photonMultiplicity",photonMultiplicity,"photonMultiplicity[noJets]/I");
  _tree->Branch("hfEMEnergyFr",hfEMEnergyFr,"hfEMEnergyFr[noJets]/F");
  _tree->Branch("hfEMMultiplicity",hfEMMultiplicity,"hfEMMultiplicity[noJets]/I");
  _tree->Branch("hfHadronEnergyFr",hfHadronEnergyFr,"hfHadronEnergyFr[noJets]/F");
  _tree->Branch("hfHadronMultiplicity",hfHadronMultiplicity,"hfHadronMultiplicity[noJets]/I");
  _tree->Branch("noChargedConstituents",noChargedConstituents,"noChargedConstituents[noJets]/I");
  _tree->Branch("chargedEnergyFr",chargedEnergyFr,"chargedEnergyFr[noJets]/F");
  _tree->Branch("chargedMultiplicity",chargedMultiplicity,"chargedMultiplicity[noJets]/I");
  _tree->Branch("electronEnergyFr",electronEnergyFr,"electronEnergyFr[noJets]/F");
  _tree->Branch("electronMultiplicity",electronMultiplicity,"electronMultiplicity[noJets]/I");
  _tree->Branch("muonEnergyFr",muonEnergyFr,"muonEnergyFr[noJets]/F");
  _tree->Branch("muonMultiplicity",muonMultiplicity,"muonMultiplicity[noJets]/I");

  _tree->Branch("recoJetIDLoose",recoJetIDLoose,"recoJetIDLoose[noRecoJets]/O");
  _tree->Branch("recoJetIDMedium",recoJetIDMedium,"recoJetIDMedium[noRecoJets]/O");
  _tree->Branch("recoJetIDTight",recoJetIDTight,"recoJetIDTight[noRecoJets]/O");
  _tree->Branch("recoJetPt",recoJetPt,"recoJetPt[noRecoJets]/F");
  _tree->Branch("recoJetEt",recoJetEt,"recoJetEt[noRecoJets]/F");
  _tree->Branch("recoJetPhi",recoJetPhi,"recoJetPhi[noRecoJets]/F");
  _tree->Branch("recoJetEta",recoJetEta,"recoJetEta[noRecoJets]/F");
  _tree->Branch("recoJetMass",recoJetMass,"recoJetMass[noRecoJets]/F");
  _tree->Branch("recoJetArea",recoJetArea,"recoJetArea[noRecoJets]/F");


  // jet PU ID ===>
  /*_tree->Branch("puJetInputJetEta", puJetInputJetEta, "puJetInputJetEta[noJets]/F");
  _tree->Branch("puJetInputJetPt", puJetInputJetPt, "puJetInputJetPt[noJets]/F");
  _tree->Branch("puJetInputJetNCharged", puJetInputJetNCharged, "puJetInputJetNCharged[noJets]/F");
  _tree->Branch("puJetInputJetNNeutrals", puJetInputJetNNeutrals, "puJetInputJetNNeutrals[noJets]/F");
  _tree->Branch("puJetInputJetDZ", puJetInputJetDZ, "puJetInputJetDZ[noJets]/F");
  _tree->Branch("puJetInputJetNParticles", puJetInputJetNParticles, "puJetInputJetNParticles[noJets]/F");
  _tree->Branch("puJetInputJetDR2Mean", puJetInputJetDR2Mean, "puJetInputJetDR2Mean[noJets]/F");
  _tree->Branch("puJetInputJetDRMean", puJetInputJetDRMean, "puJetInputJetDRMean[noJets]/F");
  _tree->Branch("puJetInputJetFrac01", puJetInputJetFrac01, "puJetInputJetFrac01[noJets]/F");
  _tree->Branch("puJetInputJetFrac02", puJetInputJetFrac02, "puJetInputJetFrac02[noJets]/F");
  _tree->Branch("puJetInputJetFrac03", puJetInputJetFrac03, "puJetInputJetFrac03[noJets]/F");
  _tree->Branch("puJetInputJetFrac04", puJetInputJetFrac04, "puJetInputJetFrac04[noJets]/F");
  _tree->Branch("puJetInputJetFrac05", puJetInputJetFrac05, "puJetInputJetFrac05[noJets]/F");
  _tree->Branch("puJetInputJetFrac06", puJetInputJetFrac06, "puJetInputJetFrac06[noJets]/F");
  _tree->Branch("puJetInputJetFrac07", puJetInputJetFrac07, "puJetInputJetFrac07[noJets]/F");
  _tree->Branch("puJetInputJetBeta", puJetInputJetBeta, "puJetInputJetBeta[noJets]/F");
  _tree->Branch("puJetInputJetBetaStar", puJetInputJetBetaStar, "puJetInputJetBetaStar[noJets]/F");
  _tree->Branch("puJetInputJetBetaClassic", puJetInputJetBetaClassic, "puJetInputJetBetaClassic[noJets]/F");
  _tree->Branch("puJetInputJetBetaStarClassic", puJetInputJetBetaStarClassic, "puJetInputJetBetaStarClassic[noJets]/F");
  _tree->Branch("puJetInputJetPtD", puJetInputJetPtD, "puJetInputJetPtD[noJets]/F");
  _tree->Branch("puJetInputJetNvtx", puJetInputJetNvtx, "puJetInputJetNvtx[noJets]/F");
  _tree->Branch("puJetFullMVA", puJetFullMVA, "puJetFullMVA[noJets]/F");
  //_tree->Branch("puJetFullLoose", puJetFullLoose, "puJetFullLoose[noJets]/O");
  _tree->Branch("puJetFullMedium", puJetFullMedium, "puJetFullMedium[noJets]/O");
  _tree->Branch("puJetFullTight", puJetFullTight, "puJetFullTight[noJets]/O");*/

  // secondary Verices
  _tree->Branch("NumberOfSV",&nSv,"NumberOfSV/I");
  _tree->Branch("FlightDistanceSV",fdSv,"FlightDistanceSV[NumberOfSV]/F");
  _tree->Branch("FlightDistanceSigSV",fdSigSv,"FlightDistanceSigSV[NumberOfSV]/F");
  _tree->Branch("MassSV",massSv,"MassSV[NumberOfSV]/F");
  _tree->Branch("nTrkSv",nTrkSv,"nTrkSv[NumberOfSV]/I");
  _tree->Branch("PxTrkSv",pxTrkSv,"pxTrkSv[NumberOfSV][10]/F");
  _tree->Branch("PyTrkSv",pyTrkSv,"pyTrkSv[NumberOfSV][10]/F");
  _tree->Branch("PzTrkSv",pzTrkSv,"pzTrkSv[NumberOfSV][10]/F");
  _tree->Branch("JetSv",jetSv,"JetSv[NumberOfSV]/I");


  // collinear approximation =====>
  _tree->Branch("ValidDiTau",&_validDiTau,"ValidDiTau/O");
  _tree->Branch("DiTauMass",&_diTauMass,"DiTauMass/F");
  _tree->Branch("PosLepETauRF",&_posLepETauRest,"PosLepETauRF/F");
  _tree->Branch("NegLepETauRF",&_negLepETauRest,"NegLepETauRF/F");
  _tree->Branch("PosLepEDiTauRF",&_posLepEdiTauRest,"PosLepEDiTauRF/F");
  _tree->Branch("NegLepEDiTauRF",&_posLepEdiTauRest,"NegLepEDiTauRF/F");
  _tree->Branch("PosLepTauERatio",&_posLepTauERatio,"PosLepTauERatio/F");
  _tree->Branch("NegLepTauERatio",&_negLepTauERatio,"NegLepTauERatio/F");
  _tree->Branch("CosNegDiTau",&_cosNegDiTau,"CosNegDiTau/F");

  // di-tau kinematics from collinear approximation
  _tree->Branch("DiTauPx",&_diTauPx,"DiTauPx/F");
  _tree->Branch("DiTauPy",&_diTauPy,"DiTauPy/F");
  _tree->Branch("DiTauPz",&_diTauPz,"DiTauPz/F");
  _tree->Branch("DiTauEta",&_diTauEta,"DiTauEta/F");
  _tree->Branch("DiTauPt",&_diTauPt,"DiTauPt/F");
  _tree->Branch("DiTauE",&_diTauE,"DiTauE/F");

  // Generator information ---->
  if (_doMC) {

    _tree->Branch("Boson",&_bosonPDG,"Boson/I");
    _tree->Branch("Decay",&_decayPDG,"DecayPDG/I");

    /*if(!_origGenParticleSrc.label().empty())
    {
      _tree->Branch("OrigBoson",&_origBosonPDG,"OrigBoson/I");
      _tree->Branch("OrigDecay",&_origDecayPDG,"OrigDecayPDG/I");
    }*/

    _tree->Branch("ZExist",&_bosonExist,"ZExist/O"); 
    _tree->Branch("ZMass",&_ZMass,"ZMass/F");
    _tree->Branch("ZPx",&_ZPx,"ZPx/F");
    _tree->Branch("ZPy",&_ZPy,"ZPy/F");
    _tree->Branch("ZPz",&_ZPz,"ZPz/F");
    
    _tree->Branch("HMass",&_HMass,"HMass/F");
    _tree->Branch("HPx",&_HPx,"HPx/F");
    _tree->Branch("HPy",&_HPy,"HPy/F");
    _tree->Branch("HPz",&_HPz,"HPz/F");
    
    _tree->Branch("bbbarPhi",&bbbarPhi,"bbbarPhi/O");
    _tree->Branch("isVBF",&isVBF,"isVBF/O");
    _tree->Branch("VBFid",&VBFid,"VBFid/I");

    _tree->Branch("bPx",&bPx,"bPx/F");
    _tree->Branch("bPy",&bPy,"bPy/F");
    _tree->Branch("bPz",&bPz,"bPz/F");
    _tree->Branch("benergy",&benergy,"benergy/F");

    _tree->Branch("bbarPx",&bbarPx,"bbarPx/F");
    _tree->Branch("bbarPy",&bbarPy,"bbarPy/F");
    _tree->Branch("bbarPz",&bbarPz,"bbarPz/F");
    _tree->Branch("bbarenergy",&bbarenergy,"bbarenergy/F");

    _tree->Branch("parton1pdg",&parton1pdg,"parton1pdg/I");
    _tree->Branch("parton1Px",&parton1Px,"parton1Px/F");
    _tree->Branch("parton1Py",&parton1Py,"parton1Py/F");
    _tree->Branch("parton1Pz",&parton1Pz,"parton1Pz/F");
    _tree->Branch("parton1energy",&parton1energy,"parton1energy/F");

    _tree->Branch("parton2pdg",&parton2pdg,"parton2pdg/I");
    _tree->Branch("parton2Px",&parton2Px,"parton2Px/F");
    _tree->Branch("parton2Py",&parton2Py,"parton2Py/F");
    _tree->Branch("parton2Pz",&parton2Pz,"parton2Pz/F");
    _tree->Branch("parton2energy",&parton2energy,"parton2energy/F");

    _tree->Branch("PosDecayLeptonPx",&_posDecayMuonPx,"PosDecayLeptonPx/F");
    _tree->Branch("PosDecayLeptonPy",&_posDecayMuonPy,"PosDecayLeptonPy/F");
    _tree->Branch("PosDecayLeptonPz",&_posDecayMuonPz,"PosDecayLeptonPz/F");
    _tree->Branch("PosDecayLeptonPdg",&_posDecayMuonPdg,"PosDecayLeptonPdg/I");

    _tree->Branch("NegDecayLeptonPx",&_negDecayMuonPx,"NegDecayLeptonPx/F");
    _tree->Branch("NegDecayLeptonPy",&_negDecayMuonPy,"NegDecayLeptonPy/F");
    _tree->Branch("NegDecayLeptonPz",&_negDecayMuonPz,"NegDecayLeptonPz/F");
    _tree->Branch("NegDecayLeptonPdg",&_negDecayMuonPdg,"NegDecayLeptonPdg/I");

    _tree->Branch("NFsrPos",&_nFSRPos,"NFsrPos/I");
    _tree->Branch("FsrPosPdg",_FSRPosPdg,"FsrPosPdg[NFsrPos]/I");
    _tree->Branch("FsrPosPx",_FSRPosPx,"FsrPosPx[NFsrPos]/F");
    _tree->Branch("FsrPosPy",_FSRPosPy,"FsrPosPy[NFsrPos]/F");
    _tree->Branch("FsrPosPz",_FSRPosPz,"FsrPosPz[NFsrPos]/F");

    _tree->Branch("NFsrNeg",&_nFSRNeg,"NFsrNeg/I");
    _tree->Branch("FsrNegPdg",_FSRNegPdg,"FsrNegPdg[NFsrNeg]/I");
    _tree->Branch("FsrNegPx",_FSRNegPx,"FsrNegPx[NFsrNeg]/F");
    _tree->Branch("FsrNegPy",_FSRNegPy,"FsrNegPy[NFsrNeg]/F");
    _tree->Branch("FsrNegPz",_FSRNegPz,"FsrNegPz[NFsrNeg]/F");


    // positive generated muon
    _tree->Branch("posMuonGenExist",&posMuonGenExist,"posMuonGenExist/O");
    _tree->Branch("posMuonGenpdg",&posMuonGenpdg,"posMuonGenpdg/I");
    _tree->Branch("posMuonGenQ",&posMuonGenQ,"posMuonGenQ/F");
    _tree->Branch("posMuonGenpx",&posMuonGenpx,"posMuonGenpx/F");
    _tree->Branch("posMuonGenpy",&posMuonGenpy,"posMuonGenpy/F");
    _tree->Branch("posMuonGenpz",&posMuonGenpz,"posMuonGenpz/F");
    _tree->Branch("posMuonGenpt",&posMuonGenpt,"posMuonGenpt/F");
    _tree->Branch("posMuonGenphi",&posMuonGenphi,"posMuonGenphi/F");
    _tree->Branch("posMuonGeneta",&posMuonGeneta,"posMuonGeneta/F");

    _tree->Branch("PMmotherExist",&PMmotherExist,"PMmotherExist/O");    
    _tree->Branch("PMmotherpdg",&PMmotherpdg,"PMmotherpdg/I");
    _tree->Branch("PMmotherphi",&PMmotherphi,"PMmotherphi/F");
    _tree->Branch("PMmothereta",&PMmothereta,"PMmothereta/F");
    _tree->Branch("PMmotherpx",&PMmotherpx,"PMmotherpx/F");
    _tree->Branch("PMmotherpy",&PMmotherpy,"PMmotherpy/F");
    _tree->Branch("PMmotherpz",&PMmotherpz,"PMmotherpz/F");
    _tree->Branch("PMmotherpt",&PMmotherpt,"PMmotherpt/F");
    
    _tree->Branch("PMGrmotherExist",&PMGrmotherExist,"PMGrmotherExist/O");
    _tree->Branch("PMGrmotherpdg",&PMGrmotherpdg,"PMGrmotherpdg/I");
    _tree->Branch("PMGrmotherphi",&PMGrmotherphi,"PMGrmotherphi/F");
    _tree->Branch("PMGrmothereta",&PMGrmothereta,"PMGrmothereta/F");
    _tree->Branch("PMGrmotherpx",&PMGrmotherpx,"PMGrmotherpx/F");
    _tree->Branch("PMGrmotherpy",&PMGrmotherpy,"PMGrmotherpy/F");
    _tree->Branch("PMGrmotherpz",&PMGrmotherpz,"PMGrmotherpz/F");
    _tree->Branch("PMGrmotherpt",&PMGrmotherpt,"PMGrmotherpt/F");
    

    // negative generated muon
    _tree->Branch("negMuonGenExist",&negMuonGenExist,"negMuonGenExist/O");
    _tree->Branch("negMuonGenpdg",&negMuonGenpdg,"negMuonGenpdg/I");
    _tree->Branch("negMuonGenQ",&negMuonGenQ,"negMuonGenQ/F");
    _tree->Branch("negMuonGenpx",&negMuonGenpx,"negMuonGenpx/F");
    _tree->Branch("negMuonGenpy",&negMuonGenpy,"negMuonGenpy/F");
    _tree->Branch("negMuonGenpz",&negMuonGenpz,"negMuonGenpz/F");
    _tree->Branch("negMuonGenpt",&negMuonGenpt,"negMuonGenpt/F");
    _tree->Branch("negMuonGenphi",&negMuonGenphi,"negMuonGenphi/F");
    _tree->Branch("negMuonGeneta",&negMuonGeneta,"negMuonGeneta/F");

    _tree->Branch("NMmotherExist",&NMmotherExist,"NMmotherExist/O");  
    _tree->Branch("NMmotherpdg",&NMmotherpdg,"NMmotherpdg/I");
    _tree->Branch("NMmotherphi",&NMmotherphi,"NMmotherphi/F");
    _tree->Branch("NMmothereta",&NMmothereta,"NMmothereta/F");
    _tree->Branch("NMmotherpx",&NMmotherpx,"NMmotherpx/F");
    _tree->Branch("NMmotherpy",&NMmotherpy,"NMmotherpy/F");
    _tree->Branch("NMmotherpz",&NMmotherpz,"NMmotherpz/F");
    _tree->Branch("NMmotherpt",&NMmotherpt,"NMmotherpt/F");
    
    _tree->Branch("NMGrmotherExist",&NMGrmotherExist,"NMGrmotherExist/O");
    _tree->Branch("NMGrmotherpdg",&NMGrmotherpdg,"NMGrmotherpdg/I");
    _tree->Branch("NMGrmotherphi",&NMGrmotherphi,"NMGrmotherphi/F");
    _tree->Branch("NMGrmothereta",&NMGrmothereta,"NMGrmothereta/F");
    _tree->Branch("NMGrmotherpx",&NMGrmotherpx,"NMGrmotherpx/F");
    _tree->Branch("NMGrmotherpy",&NMGrmotherpy,"NMGrmotherpy/F");
    _tree->Branch("NMGrmotherpz",&NMGrmotherpz,"NMGrmotherpz/F");
    _tree->Branch("NMGrmotherpt",&NMGrmotherpt,"NMGrmotherpt/F");
  }

  std::cout << " HTauTauNMSSMAnalysis::beginJob() --> TTree is declared " << std::endl; 

}

void HTauTauNMSSMAnalysis::beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &setup)
{
	//ntpGeneralInfo.beginLuminosityBlock(lumiBlock, setup);
}

// ------------ method called once each job just after ending the event loop  ------------
void HTauTauNMSSMAnalysis::endJob() {
  std::cout << std::endl;
  std::cout << " HTauTauNMSSMAnalysis::endJob() --> " << std::endl;
  std::cout << " Input events = " << _event << std::endl;
  std::cout << " Selected events = " << _selEvent << std::endl;
  std::cout << std::endl;
  
  _tree->GetDirectory()->cd();
  _tree->Write();

	//ntpGeneralInfo.endJob();
}

// void HTauTauNMSSMAnalysis::fillMuonMVAVariables(
// 		const reco::Muon& posMuon,
// 		const reco::Muon& negMuon,
// 		const reco::Vertex& primaryVertex,
// 		const std::vector<reco::PFCandidate> pfCandidates,
// 		const edm::Event& iEvent, const edm::EventSetup& iSetup)
// {
// 	edm::Handle<double> hRho;
// 	edm::InputTag tag("kt6PFJets","rho");
// 	iEvent.getByLabel(tag,hRho);

// 	_mvaIDMuPos = -1;
// 	_mvaIDMuNeg = -1;
// 	_mvaIsoMuPos = -1;
// 	_mvaIsoMuNeg = -1;

// 	std::vector<reco::Muon> IdentifiedMuons;
// 	std::vector<reco::GsfElectron> IdentifiedElectrons;

// 	edm::Handle<std::vector<reco::GsfElectron> > inElectrons;
// 	if (!iEvent.getByLabel("gsfElectrons", inElectrons))
// 	{
// 		std::cout << "gsfElectrons not found" << std::endl;
// 		return;
// 	}

// 	edm::Handle<edm::View<reco::Muon> > inMuons;
// 	if (!iEvent.getByLabel(_muonSrc, inMuons))
// 	{
// 		std::cout << _muonSrc << " not found" << std::endl;
// 		assert(false);
// 		return;
// 	}


// 	for (std::vector<reco::GsfElectron>::const_iterator iE = inElectrons->begin(); iE != inElectrons->end(); ++iE)
// 	{
// 		double electronTrackZ = 0;
// 		if (iE->gsfTrack().isNonnull())
// 		{
// 			electronTrackZ = iE->gsfTrack()->dz(primaryVertex.position());
// 		}
// 		else if (iE->closestCtfTrackRef().isNonnull())
// 		{
// 			electronTrackZ = iE->closestCtfTrackRef()->dz(primaryVertex.position());
// 		}
// 		if(fabs(electronTrackZ) > 0.2)
// 			continue;

// 		if(fabs(iE->superCluster()->eta())<1.479)
// 		{
// 			if(iE->pt() > 20)
// 			{
// 				if(iE->sigmaIetaIeta() > 0.01) continue;
// 				if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
// 				if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
// 				if(iE->hadronicOverEm() > 0.15) continue;
// 			}
// 			else
// 			{
// 				if(iE->sigmaIetaIeta() > 0.012) continue;
// 				if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
// 				if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
// 				if(iE->hadronicOverEm() > 0.15) continue;
// 			}
// 		}
// 		else
// 		{
// 			if(iE->pt() > 20)
// 			{
// 				if(iE->sigmaIetaIeta() > 0.03) continue;
// 				if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
// 				if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
// 			}
// 			else
// 			{
// 				if(iE->sigmaIetaIeta() > 0.032) continue;
// 				if(fabs(iE->deltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
// 				if(fabs(iE->deltaPhiSuperClusterTrackAtVtx()) > 0.8) continue;
// 			}
// 		}
// 		IdentifiedElectrons.push_back(*iE);
// 	}

// 	for (edm::View<reco::Muon>::const_iterator iM = inMuons->begin(); iM != inMuons->end(); ++iM)
// 	{

// 		if(!(iM->innerTrack().isNonnull()))
// 			continue;

// 		if(!(iM->isGlobalMuon() || iM->isTrackerMuon()))
// 			continue;

// 		if(iM->innerTrack()->numberOfValidHits() < 11 )
// 			continue;

// 		IdentifiedMuons.push_back(*iM);
// 	}

// 	MuonEffectiveArea::MuonEffectiveAreaTarget currentAreaTarget = MuonEffectiveArea::kMuEANoCorr;

// 	if (_year == 2011)
// 	{
// 		if (iEvent.isRealData())
// 			currentAreaTarget = MuonEffectiveArea::kMuEAData2011;
// 		else
// 		{
// 			if (_period == "Summer11")
// 				currentAreaTarget = MuonEffectiveArea::kMuEASummer11MC;
// 			if (_period == "Fall11")
// 				currentAreaTarget = MuonEffectiveArea::kMuEAFall11MC;
// 		}
// 	}
// 	if (_year == 2012)
// 	{
// 		if (iEvent.isRealData())
// 			currentAreaTarget = MuonEffectiveArea::kMuEAData2012;
// 	}

// 	_mvaIsoMuPos = fMuonIsoMVA.mvaValue(posMuon, primaryVertex, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);
// 	_mvaIsoMuNeg = fMuonIsoMVA.mvaValue(negMuon, primaryVertex, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);

// 	_mvaIDMuPos = fMuonIDMVA.mvaValue(posMuon, primaryVertex, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);
// 	_mvaIDMuNeg = fMuonIDMVA.mvaValue(negMuon, primaryVertex, pfCandidates, *hRho, currentAreaTarget, IdentifiedElectrons, IdentifiedMuons);

// }

bool HTauTauNMSSMAnalysis::isDmeson(int pdg) {

  int absPdg = TMath::Abs(pdg);
  bool yes = (absPdg >= 411) && (absPdg <= 439);
  return yes;

}

bool HTauTauNMSSMAnalysis::isBmeson(int pdg) {

  int absPdg = TMath::Abs(pdg);
  bool yes = (absPdg >= 511) && (absPdg <= 549);
  return yes;

}

bool HTauTauNMSSMAnalysis::isJPsi(int pdg) {

  int absPdg = TMath::Abs(pdg);
  bool yes = (absPdg >= 441) && (absPdg <= 499);
  return yes;

}

bool HTauTauNMSSMAnalysis::isPFMuon(const pat::MuonRef& muon, const std::vector<reco::PFCandidate> pfCandidates)
{
	bool isPFMuon = false;
	for(std::vector<reco::PFCandidate>::const_iterator iter = pfCandidates.begin(); iter != pfCandidates.end() && !isPFMuon; ++iter)
		if(abs(iter->pdgId()) == 13 && iter->muonRef().get() == muon->originalObject())
			isPFMuon = true;
	return isPFMuon;
}
bool HTauTauNMSSMAnalysis::isYpsilon(int pdg) {

  int absPdg = TMath::Abs(pdg);
  bool yes = (absPdg >= 551) && (absPdg <= 599);
  return yes;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HTauTauNMSSMAnalysis);
