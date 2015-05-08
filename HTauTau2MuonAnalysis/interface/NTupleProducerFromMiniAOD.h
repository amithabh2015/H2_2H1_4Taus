#ifndef NTupleProducerFromMiniAOD_h
#define NTupleProducerFromMiniAOD_h
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include <Math/Vector3D.h>
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "CondFormats/L1TObjects/interface/L1GtPrescaleFactors.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsAlgoTrigRcd.h"
#include "CondFormats/DataRecord/interface/L1GtPrescaleFactorsTechTrigRcd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h" 
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h" 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "RecoEgamma/EgammaTools/interface/ConversionInfo.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"

using namespace std;
using namespace reco;

#define M_trackmaxcount 1000
#define M_primvertexmaxcount 1000
#define M_muonmaxcount 1000
#define M_jetmaxcount 1000
#define M_genparticlesmaxcount 1000
#define M_trigobjectmaxcount 1000
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

bool doDebug = false;
class NTupleProducerFromMiniAOD : public edm::EDAnalyzer{ 
 public:
  explicit NTupleProducerFromMiniAOD( const edm::ParameterSet& iConfig );
  ~NTupleProducerFromMiniAOD();

 private:
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
  virtual void beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup);
  virtual void analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup );

  bool AddGenParticles(const edm::Event& iEvent);
  unsigned int AddMuons(const edm::Event& iEvent);
  unsigned int AddPackedPFCand(const edm::Event& iEvent);
  unsigned int AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  unsigned int AddTriggerObjects(const edm::Event& iEvent);
  bool foundCompatibleInnerHits(const reco::HitPattern& hitPatA, const reco::HitPattern& hitPatB);
  
  UInt_t GenParticleInfo(const GenParticle* particle);
  bool GetL1ExtraTriggerMatch(const l1extra::L1JetParticleCollection* l1jets,  const l1extra::L1JetParticleCollection* l1taus, const LeafCandidate& leg2);
  Int_t HasAnyMother(const GenParticle* particle, int id);
  math::XYZPoint PositionOnECalSurface(reco::TransientTrack&);
  double ComputeDiTauMass(LorentzVector leg1, LorentzVector leg2, LorentzVector met, TMatrixD cov);
  //  TLorentzVector RecoilCorrectedMET(LorentzVector pfMet_, LorentzVector Leg1p4_, LorentzVector Leg2p4_, const reco::GenParticle *boson_, string sampleName_, int nJets_);

  struct DCA {
    float dca2d;
    float dca2dErr;
    float dca3d;
    float dca3dErr;
  };
  DCA calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2);

  TTree* tree;
  TTree* lumitree;
  TTree* runtree;
  TH1D*  nEvents;

  //Configuration (steering cards)

  bool cdata;

  vector<string> cHLTriggerPaths;
  string cTriggerProcess;

  double cMuPtMin;
  double cMuEtaMax;
  vector<string> cMuHLTriggerMatching;
  int cMuNum;

  double cTrackPtMin;
  double cTrackEtaMax;
  int cTrackNum;

  double cJetPtMin;
  double cJetEtaMax;
  vector<string> cBtagDiscriminators;
  vector<string> cJetHLTriggerMatching;
  int cJetNum;

  edm::InputTag MuonCollectionTag_;
  edm::InputTag JetCollectionTag_;
  edm::InputTag MetCollectionTag_;
  edm::InputTag TrackCollectionTag_;
  edm::InputTag GenParticleCollectionTag_;
  edm::InputTag TriggerObjectCollectionTag_;
  edm::InputTag BeamSpotTag_;
  edm::InputTag PVTag_;
  std::string sampleName;

  //Variables
  edm::ESHandle<TransientTrackBuilder>  TTrackBuilder         ;
  edm::ESHandle<MagneticField>          magneticField         ; 
  Cylinder::ConstCylinderPointer        ecalBarrel            ;
  Plane::ConstPlanePointer              ecalNegativeEtaEndcap ;
  Plane::ConstPlanePointer              ecalPositiveEtaEndcap ;
  //  MEtRecoilCorrection *                 metRecCorr; 
  //  RecoilCorrector *                     corrector_ ;  
  
  HLTConfigProvider HLTConfiguration;
  edm::Handle<edm::TriggerResults> HLTrigger;
  edm::Handle<trigger::TriggerEvent> HLTriggerEvent;
  //edm::Handle<l1extra::L1MuonParticleCollection> L1Muons;
  //edm::Handle<l1extra::L1EmParticleCollection> L1Electrons;
  //edm::Handle<l1extra::L1EmParticleCollection> L1ElectronsIso;
  vector<std::pair<string,string> > muontriggers;
  vector<std::pair<string,string> > jettriggers;

  vector<int> HLTriggerIndexSelection;
  vector<int> tauIndexSelection;
  vector<int> DiTauIndex;

  math::XYZPoint pv_position;
  Vertex primvertex;

  //Data		
  Int_t njets4RC;

  UInt_t errors;
  UInt_t event_nr;
  UInt_t event_luminosityblock;
  UInt_t event_run;
  UInt_t event_timeunix;
  UInt_t event_timemicrosec;
  UChar_t trigger_level1bits[8];
  UChar_t trigger_level1[128];
  UChar_t trigger_HLT[128];

  // beam spot   
  Float_t beamspot_x;
  Float_t beamspot_y;
  Float_t beamspot_z;
  Float_t beamspot_xwidth;
  Float_t beamspot_ywidth;
  Float_t beamspot_zsigma;
  Float_t beamspot_cov[6];

  // primary vertex
  UInt_t  primvertex_count;
  Float_t primvertex_x;
  Float_t primvertex_y;
  Float_t primvertex_z;
  Float_t primvertex_chi2;
  Float_t primvertex_ndof;
  Float_t primvertex_ptq;
  Int_t   primvertex_ntracks;
  Float_t primvertex_cov[6];

  //tracks
  UInt_t track_count; 
  Float_t track_px[M_trackmaxcount];
  Float_t track_py[M_trackmaxcount];
  Float_t track_pz[M_trackmaxcount];
  Float_t track_pt[M_trackmaxcount];
  Float_t track_eta[M_trackmaxcount];
  Float_t track_phi[M_trackmaxcount];
  Float_t track_mass[M_trackmaxcount];
  Int_t track_charge[M_trackmaxcount];
  /*Float_t track_outerx[M_trackmaxcount];
  Float_t track_outery[M_trackmaxcount];
  Float_t track_outerz[M_trackmaxcount];
  Float_t track_closestpointx[M_trackmaxcount];
  Float_t track_closestpointy[M_trackmaxcount];
  Float_t track_closestpointz[M_trackmaxcount];
  Float_t track_chi2[M_trackmaxcount];
  Float_t track_ndof[M_trackmaxcount];*/
  Float_t track_dxy[M_trackmaxcount];
  Float_t track_dxyerr[M_trackmaxcount];
  Float_t track_dz[M_trackmaxcount];
  Float_t track_dzerr[M_trackmaxcount];
  Int_t track_ID[M_trackmaxcount];
  /*Float_t track_dedxharmonic2[M_trackmaxcount];
  Int_t track_charge[M_trackmaxcount];
  UChar_t track_nhits[M_trackmaxcount];
  UChar_t track_nmissinghits[M_trackmaxcount];
  UChar_t track_npixelhits[M_trackmaxcount];
  UChar_t track_npixellayers[M_trackmaxcount];
  UChar_t track_nstriplayers[M_trackmaxcount];*/

  //electrons from tracks
  UInt_t track_e_pos_count; 
  Float_t track_e_pos_px[M_trackmaxcount];
  Float_t track_e_pos_py[M_trackmaxcount];
  Float_t track_e_pos_pz[M_trackmaxcount];
  Float_t track_e_pos_pt[M_trackmaxcount];
  Float_t track_e_pos_eta[M_trackmaxcount];
  Float_t track_e_pos_phi[M_trackmaxcount];
  Float_t track_e_pos_mass[M_trackmaxcount];
  Int_t track_e_pos_charge[M_trackmaxcount];
  Float_t track_e_pos_dxy[M_trackmaxcount];
  Float_t track_e_pos_dxyerr[M_trackmaxcount];
  Float_t track_e_pos_dz[M_trackmaxcount];
  Float_t track_e_pos_dzerr[M_trackmaxcount];
  UInt_t track_e_neg_count;
  Float_t track_e_neg_px[M_trackmaxcount];
  Float_t track_e_neg_py[M_trackmaxcount];
  Float_t track_e_neg_pz[M_trackmaxcount];
  Float_t track_e_neg_pt[M_trackmaxcount];
  Float_t track_e_neg_eta[M_trackmaxcount];
  Float_t track_e_neg_phi[M_trackmaxcount];
  Float_t track_e_neg_mass[M_trackmaxcount];
  Int_t track_e_neg_charge[M_trackmaxcount];
  Float_t track_e_neg_dxy[M_trackmaxcount];
  Float_t track_e_neg_dxyerr[M_trackmaxcount];
  Float_t track_e_neg_dz[M_trackmaxcount];
  Float_t track_e_neg_dzerr[M_trackmaxcount];
 
  //muons from tracks
  UInt_t track_m_pos_count; 
  Float_t track_m_pos_px[M_trackmaxcount];
  Float_t track_m_pos_py[M_trackmaxcount];
  Float_t track_m_pos_pz[M_trackmaxcount];
  Float_t track_m_pos_pt[M_trackmaxcount];
  Float_t track_m_pos_eta[M_trackmaxcount];
  Float_t track_m_pos_phi[M_trackmaxcount];
  Float_t track_m_pos_mass[M_trackmaxcount];
  Int_t track_m_pos_charge[M_trackmaxcount];
  Float_t track_m_pos_dxy[M_trackmaxcount];
  Float_t track_m_pos_dxyerr[M_trackmaxcount];
  Float_t track_m_pos_dz[M_trackmaxcount];
  Float_t track_m_pos_dzerr[M_trackmaxcount];
  UInt_t track_m_neg_count;
  Float_t track_m_neg_px[M_trackmaxcount];
  Float_t track_m_neg_py[M_trackmaxcount];
  Float_t track_m_neg_pz[M_trackmaxcount];
  Float_t track_m_neg_pt[M_trackmaxcount];
  Float_t track_m_neg_eta[M_trackmaxcount];
  Float_t track_m_neg_phi[M_trackmaxcount];
  Float_t track_m_neg_mass[M_trackmaxcount];
  Int_t track_m_neg_charge[M_trackmaxcount];
  Float_t track_m_neg_dxy[M_trackmaxcount];
  Float_t track_m_neg_dxyerr[M_trackmaxcount];
  Float_t track_m_neg_dz[M_trackmaxcount];
  Float_t track_m_neg_dzerr[M_trackmaxcount];
  
  //pions from tracks
  UInt_t track_pi_pos_count; 
  Float_t track_pi_pos_px[M_trackmaxcount];
  Float_t track_pi_pos_py[M_trackmaxcount];
  Float_t track_pi_pos_pz[M_trackmaxcount];
  Float_t track_pi_pos_pt[M_trackmaxcount];
  Float_t track_pi_pos_eta[M_trackmaxcount];
  Float_t track_pi_pos_phi[M_trackmaxcount];
  Float_t track_pi_pos_mass[M_trackmaxcount];
  Int_t track_pi_pos_charge[M_trackmaxcount];
  Float_t track_pi_pos_dxy[M_trackmaxcount];
  Float_t track_pi_pos_dxyerr[M_trackmaxcount];
  Float_t track_pi_pos_dz[M_trackmaxcount];
  Float_t track_pi_pos_dzerr[M_trackmaxcount];
  UInt_t track_pi_neg_count;
  Float_t track_pi_neg_px[M_trackmaxcount];
  Float_t track_pi_neg_py[M_trackmaxcount];
  Float_t track_pi_neg_pz[M_trackmaxcount];
  Float_t track_pi_neg_pt[M_trackmaxcount];
  Float_t track_pi_neg_eta[M_trackmaxcount];
  Float_t track_pi_neg_phi[M_trackmaxcount];
  Float_t track_pi_neg_mass[M_trackmaxcount];
  Int_t track_pi_neg_charge[M_trackmaxcount];
  Float_t track_pi_neg_dxy[M_trackmaxcount];
  Float_t track_pi_neg_dxyerr[M_trackmaxcount];
  Float_t track_pi_neg_dz[M_trackmaxcount];
  Float_t track_pi_neg_dzerr[M_trackmaxcount];


 
  // muons
  UInt_t muon_count;
  Float_t muon_px[M_muonmaxcount];
  Float_t muon_py[M_muonmaxcount];
  Float_t muon_pz[M_muonmaxcount];
  Float_t muon_pt[M_muonmaxcount];
  Float_t muon_eta[M_muonmaxcount];
  Float_t muon_phi[M_muonmaxcount];
  Float_t muon_pterror[M_muonmaxcount];
  Float_t muon_chi2[M_muonmaxcount];
  Float_t muon_ndof[M_muonmaxcount];
  Float_t muon_charge[M_muonmaxcount];
  UInt_t muon_nMuonStations[M_muonmaxcount];
  UInt_t muon_nMuonHits[M_muonmaxcount];
  UInt_t muon_nPixelHits[M_muonmaxcount];
  UInt_t muon_nTrackerHits[M_muonmaxcount];

  Float_t muon_dz[M_muonmaxcount];
  Float_t muon_dzerr[M_muonmaxcount];
  Float_t muon_dxy[M_muonmaxcount];
  Float_t muon_dxyerr[M_muonmaxcount];

  Float_t muon_chargedHadIso[M_muonmaxcount];
  Float_t muon_neutralHadIso[M_muonmaxcount];
  Float_t muon_photonIso[M_muonmaxcount];
  Float_t muon_puIso[M_muonmaxcount];

  Bool_t muon_isPF[M_muonmaxcount];
  Bool_t muon_isGlobal[M_muonmaxcount];
  Bool_t muon_isTracker[M_muonmaxcount];
  Bool_t muon_isTight[M_muonmaxcount];
  Bool_t muon_isLoose[M_muonmaxcount];

  Bool_t muon_globalTrack[M_muonmaxcount];
  Bool_t muon_innerTrack[M_muonmaxcount];

  // pf jets 
  UInt_t pfjet_count;
  Float_t pfjet_e[M_jetmaxcount];
  Float_t pfjet_px[M_jetmaxcount];
  Float_t pfjet_py[M_jetmaxcount];
  Float_t pfjet_pz[M_jetmaxcount];
  Float_t pfjet_pt[M_jetmaxcount];
  Float_t pfjet_eta[M_jetmaxcount];
  Float_t pfjet_phi[M_jetmaxcount];
  Float_t pfjet_hadronicenergy[M_jetmaxcount];
  Float_t pfjet_chargedhadronicenergy[M_jetmaxcount];
  Float_t pfjet_emenergy[M_jetmaxcount];
  Float_t pfjet_chargedemenergy[M_jetmaxcount];
  UInt_t pfjet_chargedmulti[M_jetmaxcount];	
  UInt_t pfjet_neutralmulti[M_jetmaxcount];	
  Float_t pfjet_energycorr[M_jetmaxcount];
  Float_t pfjet_energycorr_l1fastjet[M_jetmaxcount];
  Float_t pfjet_energycorr_l2relative[M_jetmaxcount];
  Float_t pfjet_energycorr_l3absolute[M_jetmaxcount];
  Float_t pfjet_energycorr_l2l3residual[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_cut_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_cut_mva[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_simple_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_simple_mva[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_loose[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_medium[M_jetmaxcount];
  Bool_t pfjet_pu_jet_full_tight[M_jetmaxcount];
  Float_t pfjet_pu_jet_full_mva[M_jetmaxcount];
  Int_t pfjet_flavour[M_jetmaxcount];
  Float_t pfjet_btag[M_jetmaxcount][10];

  
  // rho neutral
  Float_t rhoNeutral;

  // met
  Float_t pfmet_ex;
  Float_t pfmet_ey;

  Float_t pfmet_sigxx;
  Float_t pfmet_sigxy;
  Float_t pfmet_sigyx;
  Float_t pfmet_sigyy;

  Float_t genmet_ex;
  Float_t genmet_ey;

  // Electrons
     //Generator Information
  Float_t genweight;
  Float_t genid1;
  Float_t genx1;
  Float_t genid2;
  Float_t genx2;
  Float_t genScale;

  Int_t numpileupinteractionsminus;
  Int_t numpileupinteractions;
  Int_t numpileupinteractionsplus;
  Float_t numtruepileupinteractions;
  Int_t hepNUP_;
  
  UInt_t genparticles_count;
  Float_t genparticles_e[M_genparticlesmaxcount];
  Float_t genparticles_px[M_genparticlesmaxcount];
  Float_t genparticles_py[M_genparticlesmaxcount];
  Float_t genparticles_pz[M_genparticlesmaxcount];
  Float_t genparticles_vx[M_genparticlesmaxcount];
  Float_t genparticles_vy[M_genparticlesmaxcount];
  Float_t genparticles_vz[M_genparticlesmaxcount];
  Int_t genparticles_pdgid[M_genparticlesmaxcount];
  Int_t genparticles_status[M_genparticlesmaxcount];
  UInt_t genparticles_info[M_genparticlesmaxcount];
  string genparticles_mother[M_genparticlesmaxcount];

  // trigger objects
  UInt_t trigobject_count;
  Float_t trigobject_px[M_trigobjectmaxcount];
  Float_t trigobject_py[M_trigobjectmaxcount];
  Float_t trigobject_pz[M_trigobjectmaxcount];
  Float_t trigobject_pt[M_trigobjectmaxcount];
  Float_t trigobject_eta[M_trigobjectmaxcount];
  Float_t trigobject_phi[M_trigobjectmaxcount];
  Bool_t  trigobject_filters[M_trigobjectmaxcount][50];

  //lumitree
  UInt_t lumi_run;
  UInt_t lumi_block;
  Float_t lumi_value;
  Float_t lumi_valueerr;
  Float_t lumi_livefrac;
  Float_t lumi_deadfrac;
  UInt_t lumi_quality;
  UInt_t lumi_eventsprocessed;
  UInt_t lumi_eventsfiltered;
  UInt_t lumi_hltprescaletable;
  UInt_t lumi_l1algoprescaletable;
  UInt_t lumi_l1techprescaletable;

  //runtree
  UInt_t run_number;
  UInt_t run_hltcount;
  vector<string> run_hltnames;
  vector<string> run_hltfilters;
  vector<string> run_btagdiscriminators;
  UInt_t run_hltprescaletablescount;
  UInt_t run_hltprescaletables[10000];
  UInt_t run_l1algocount;
  UInt_t run_l1algoprescaletablescount;
  UInt_t run_l1algoprescaletables[10000];
  UInt_t run_l1techcount;
  UInt_t run_l1techprescaletablescount;
  UInt_t run_l1techprescaletables[10000];		

  std::map<std::string, int>* hltriggerresults_;
  std::map<std::string, int>* hltriggerprescales_;
  std::vector<std::string>hltriggerresultsV_;

};

DEFINE_FWK_MODULE(NTupleProducerFromMiniAOD);

#endif

