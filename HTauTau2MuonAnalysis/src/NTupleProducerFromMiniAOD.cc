#include "H2to2H1to4Taus/HTauTau2MuonAnalysis/interface/NTupleProducerFromMiniAOD.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include <DataFormats/RecoCandidate/interface/IsoDepositVetos.h>
#include <DataFormats/METReco/interface/GenMET.h>
#include "RecoBTag/BTagTools/interface/SignedImpactParameter3D.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


#include <TString.h>

using namespace reco;

//typedef std::vector<NSVfitEventHypothesisByIntegration> NSVfitEventHypothesisByIntegrationCollection;
typedef ROOT::Math::XYZVector Vector;

// http://root.cern.ch/root/html/ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___.html#ROOT__Math__SMatrix_double_3_3_-p1MatRepSym_double_3___:_SMatrix_double_3_3_ROOT__Math__MatRepSym_double_3___
// http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SVectorDoc.html
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;  // Standard Matrix representation for a general 3x3 matrix of type double.
typedef ROOT::Math::SVector<double, 3> SVector3; //// SVector: vector of size 3 

// static const unsigned int SKIM_MUTAUTAU   = (1 << 0);     //1  : mu+tau+tau
// static const unsigned int SKIM_ETAUTAU    = (1 << 1);     //2  : e+tau+tau  
// static const unsigned int SKIM_MUMU       = (1 << 2);     //4  : mu+mu
// static const unsigned int SKIM_EE         = (1 << 3);     //8  : e+e
// static const unsigned int SKIM_ETAU       = (1 << 4);     //16 : e+tau
// static const unsigned int SKIM_ALL        = (1 << 5);     //32 : all
// static const unsigned int SKIM_EMU        = (1 << 6);     //64 : e+mu
// static const unsigned int SKIM_TAUTAU     = (1 << 7);     //128: tau+tau


//to set the values from parameter
NTupleProducerFromMiniAOD::NTupleProducerFromMiniAOD(const edm::ParameterSet& iConfig) :  
  // 
  cdata(iConfig.getUntrackedParameter<bool>("IsData",false)),
  // triggers
  cHLTriggerPaths(iConfig.getUntrackedParameter<vector<string> >("HLTriggerPaths")),
  cTriggerProcess(iConfig.getUntrackedParameter<string>("TriggerProcess", "HLT")),
  // muons
  cMuPtMin(iConfig.getUntrackedParameter<double>("RecMuonPtMin", 10.)),
  cMuEtaMax(iConfig.getUntrackedParameter<double>("RecMuonEtaMax", 2.4)),
  cMuHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecMuonHLTriggerMatching")),
  cMuNum(iConfig.getUntrackedParameter<int>("RecMuonNum", 0)),
  // tracks
  cTrackPtMin(iConfig.getUntrackedParameter<double>("RecTrackPtMin", 0.5)),
  cTrackEtaMax(iConfig.getUntrackedParameter<double>("RecTrackEtaMax", 2.4)),
  cTrackNum(iConfig.getUntrackedParameter<int>("RecTrackNum", 0)),
  // jets
  cJetPtMin(iConfig.getUntrackedParameter<double>("RecJetPtMin", 30.)),
  cJetEtaMax(iConfig.getUntrackedParameter<double>("RecJetEtaMax", 4.5)),
  cBtagDiscriminators(iConfig.getUntrackedParameter<vector<string> >("RecJetBtagDiscriminators")),
  cJetHLTriggerMatching(iConfig.getUntrackedParameter<vector<string> >("RecJetHLTriggerMatching")),
  cJetNum(iConfig.getUntrackedParameter<int>("RecJetNum", 0)),
  // collections
  MuonCollectionTag_(iConfig.getParameter<edm::InputTag>("MuonCollectionTag")),
  JetCollectionTag_(iConfig.getParameter<edm::InputTag>("JetCollectionTag")),
  MetCollectionTag_(iConfig.getParameter<edm::InputTag>("MetCollectionTag")),
  TrackCollectionTag_(iConfig.getParameter<edm::InputTag>("TrackCollectionTag")),
  GenParticleCollectionTag_(iConfig.getParameter<edm::InputTag>("GenParticleCollectionTag")),
  TriggerObjectCollectionTag_(iConfig.getParameter<edm::InputTag>("TriggerObjectCollectionTag")),
  BeamSpotTag_(iConfig.getParameter<edm::InputTag>("BeamSpotCollectionTag")),
  PVTag_(iConfig.getParameter<edm::InputTag>("PVCollectionTag")),
  sampleName(iConfig.getUntrackedParameter<std::string>("SampleName", "MC"))
{
  
  double barrelRadius = 129.;  //p81, p50, ECAL TDR
  double endcapZ      = 320.5; // fig 3.26, p81, ECAL TDR
  Surface::RotationType rot;
  ecalBarrel         = Cylinder::build(Surface::PositionType(0, 0, 0), rot, barrelRadius);
  ecalNegativeEtaEndcap = Plane::build(Surface::PositionType(0, 0, -endcapZ), rot);
  ecalPositiveEtaEndcap = Plane::build(Surface::PositionType(0, 0, endcapZ), rot);
  
  const char* cmssw_base = getenv("CMSSW_BASE");
  if(!cmssw_base) throw cms::Exception("No CMSSW runtime environment found");
  std::string prefix = std::string(cmssw_base) + "/src/";

}//NTupleProducerFromMiniAOD::NTupleProducerFromMiniAOD(const edm::ParameterSet& iConfig)


//destructor
NTupleProducerFromMiniAOD::~NTupleProducerFromMiniAOD(){
  
}


void NTupleProducerFromMiniAOD::beginJob(){
  edm::Service<TFileService> FS;
  tree = FS->make<TTree>("AC1B", "AC1B", 1);
  nEvents = FS->make<TH1D>("nEvents", "nEvents", 2, -0.5, +1.5);
  
  tree->Branch("errors", &errors, "errors/i");
  tree->Branch("event_nr", &event_nr, "event_nr/i");
  tree->Branch("event_run", &event_run, "event_run/i");
  tree->Branch("event_timeunix", &event_timeunix, "event_timeunix/i");
  tree->Branch("event_timemicrosec", &event_timemicrosec, "event_timemicrosec/i");
  tree->Branch("event_luminosityblock", &event_luminosityblock, "event_luminosityblock/i");
  tree->Branch("trigger_level1bits", &trigger_level1bits, "trigger_level1bits[8]/b");
  tree->Branch("trigger_level1", &trigger_level1, "trigger_level1[128]/b");
  tree->Branch("trigger_HLT", &trigger_HLT, "trigger_HLT[128]/b");
  
  tree->Branch("beamspot_x", &beamspot_x, "beamspot_x/F");
  tree->Branch("beamspot_y", &beamspot_y, "beamspot_y/F");
  tree->Branch("beamspot_z", &beamspot_z, "beamspot_z/F");
  tree->Branch("beamspot_xwidth", &beamspot_xwidth, "beamspot_xwidth/F");
  tree->Branch("beamspot_ywidth", &beamspot_ywidth, "beamspot_ywidth/F");
  tree->Branch("beamspot_zsigma", &beamspot_zsigma, "beamspot_zsigma/F");
  tree->Branch("beamspot_cov", &beamspot_cov, "beamspot_cov[6]/F");

  // tracks (added)
  tree->Branch("track_count", &track_count, "track_count/i"); 
  tree->Branch("track_px", track_px, "track_px[track_count]/F");
  tree->Branch("track_py", track_py, "track_py[track_count]/F");
  tree->Branch("track_pz", track_pz, "track_pz[track_count]/F");
  tree->Branch("track_pt", track_pt, "track_pt[track_count]/F");
  tree->Branch("track_eta", track_eta, "track_eta[track_count]/F");
  tree->Branch("track_phi", track_phi, "track_phi[track_count]/F");
  tree->Branch("track_charge", track_charge, "track_charge[track_count]/F");
  tree->Branch("track_mass", track_mass, "track_mass[track_count]/F");
  tree->Branch("track_dxy", track_dxy, "track_dxy[track_count]/F");
  tree->Branch("track_dxyerr", track_dxyerr, "track_dxyerr[track_count]/F");
  tree->Branch("track_dz", track_dz, "track_dz[track_count]/F");
  tree->Branch("track_dzerr", track_dzerr, "track_dzerr[track_count]/F");
  tree->Branch("track_ID", track_ID, "track_ID[track_count]/I");
    //tree->Branch("track_muon_count", &track_muon_count, "track_muon_count/i");
  //tree->Branch("track_pion_count", &track_pion_count, "track_pion_count/i");
  // tracks(To be added)
  /*tree->Branch("track_outerx", track_outerx, "track_outerx[track_count]/F");
  tree->Branch("track_outery", track_outery, "track_outery[track_count]/F");
  tree->Branch("track_outerz", track_outerz, "track_outerz[track_count]/F");
  tree->Branch("track_closestpointx", track_closestpointx, "track_closestpointx[track_count]/F");
  tree->Branch("track_closestpointy", track_closestpointy, "track_closestpointy[track_count]/F");
  tree->Branch("track_closestpointz", track_closestpointz, "track_closestpointz[track_count]/F");
  tree->Branch("track_chi2", track_chi2, "track_chi2[track_count]/F");
  tree->Branch("track_ndof", track_ndof, "track_ndof[track_count]/F");
    tree->Branch("track_dedxharmonic2", track_dedxharmonic2, "track_dedxharmonic2[track_count]/F");
  tree->Branch("track_nhits", track_nhits, "track_nhits[track_count]/b");
  tree->Branch("track_npixelhits", track_npixelhits, "track_npixelhits[track_count]/b");
  tree->Branch("track_nmissinghits", track_nmissinghits, "track_nmissinghits[track_count]/b");
  tree->Branch("track_npixellayers", track_npixellayers, "track_npixellayers[track_count]/b");
  tree->Branch("track_nstriplayers", track_nstriplayers, "track_nstriplayers[track_count]/b");*/
  
  // primary vertex
  tree->Branch("primvertex_count", &primvertex_count, "primvertex_count/i"); 
  tree->Branch("primvertex_x", &primvertex_x, "primvertex_x/F");
  tree->Branch("primvertex_y", &primvertex_y, "primvertex_y/F");
  tree->Branch("primvertex_z", &primvertex_z, "primvertex_z/F");
  tree->Branch("primvertex_chi2", &primvertex_chi2, "primvertex_chi2/F");
  tree->Branch("primvertex_ndof", &primvertex_ndof, "primvertex_ndof/F");
  tree->Branch("primvertex_ptq", &primvertex_ptq, "primvertex_pdf/F");
  tree->Branch("primvertex_ntracks", &primvertex_ntracks, "primvertex_ntracks/I");
  tree->Branch("primvertex_cov", primvertex_cov, "primvertex_cov[6]/F");

  // muons
  tree->Branch("muon_count", &muon_count, "muon_count/i");
  tree->Branch("muon_px", muon_px, "muon_px[muon_count]/F");
  tree->Branch("muon_py", muon_py, "muon_py[muon_count]/F");
  tree->Branch("muon_pz", muon_pz, "muon_pz[muon_count]/F");
  tree->Branch("muon_pt", muon_pt, "muon_pt[muon_count]/F");
  tree->Branch("muon_eta", muon_eta, "muon_eta[muon_count]/F");
  tree->Branch("muon_phi", muon_phi, "muon_phi[muon_count]/F");
  tree->Branch("muon_pterror", muon_pterror, "muon_pterror[muon_count]/F");
  tree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_count]/F");
  tree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_count]/F");
  tree->Branch("muon_charge", muon_charge, "muon_charge[muon_count]/F");
  tree->Branch("muon_nMuonStations", muon_nMuonStations,"muon_nMuonStations[muon_count]/i");
  tree->Branch("muon_nMuonHits", muon_nMuonHits,"muon_nMuonHits[muon_count]/i");
  tree->Branch("muon_nPixelHits", muon_nPixelHits,"muon_nPixelHits[muon_count]/i");
  tree->Branch("muon_nTrackerHits", muon_nTrackerHits,"muon_nTrackerHits[muon_count]/i");
  tree->Branch("muon_dxy",muon_dxy,"muon_dxy[muon_count]/F");
  tree->Branch("muon_dxyerr",muon_dxyerr,"muon_dxyerr[muon_count]/F");
  tree->Branch("muon_dz",muon_dz,"muon_dz[muon_count]/F");
  tree->Branch("muon_dzerr",muon_dzerr,"muon_dzerr[muon_count]/F");
  tree->Branch("muon_chargedHadIso",muon_chargedHadIso,"muon_chargedHadIso[muon_count]/F");
  tree->Branch("muon_neutralHadIso",muon_neutralHadIso,"muon_neutralHadIso[muon_count]/F");
  tree->Branch("muon_photonIso",muon_photonIso,"muon_photonIso[muon_count]/F");
  tree->Branch("muon_puIso",muon_puIso,"muon_puIso[muon_count]/F");
  tree->Branch("muon_isPF",muon_isPF,"muon_isPF[muon_count]/O");
  tree->Branch("muon_isGlobal",muon_isGlobal,"muon_isGlobal[muon_count]/O");
  tree->Branch("muon_isTracker",muon_isTracker,"muon_isTracker[muon_count]/O");
  tree->Branch("muon_isTight",muon_isTight,"muon_isTight[muon_count]/O");
  tree->Branch("muon_isLoose",muon_isLoose,"muon_isLoose[muon_count]/O");

  // pf jets
  tree->Branch("pfjet_count", &pfjet_count, "pfjet_count/i");
  tree->Branch("pfjet_e", pfjet_e, "pfjet_e[pfjet_count]/F");
  tree->Branch("pfjet_px", pfjet_px, "pfjet_px[pfjet_count]/F");
  tree->Branch("pfjet_py", pfjet_py, "pfjet_py[pfjet_count]/F");
  tree->Branch("pfjet_pz", pfjet_pz, "pfjet_pz[pfjet_count]/F");
  tree->Branch("pfjet_pt", pfjet_pt, "pfjet_pt[pfjet_count]/F");
  tree->Branch("pfjet_eta", pfjet_eta, "pfjet_eta[pfjet_count]/F");
  tree->Branch("pfjet_phi", pfjet_phi, "pfjet_phi[pfjet_count]/F");
  tree->Branch("pfjet_hadronicenergy", pfjet_hadronicenergy, "pfjet_hadronicenergy[pfjet_count]/F");
  tree->Branch("pfjet_chargedhadronicenergy", pfjet_chargedhadronicenergy, "pfjet_chargedhadronicenergy[pfjet_count]/F");
  tree->Branch("pfjet_emenergy", pfjet_emenergy, "pfjet_emenergy[pfjet_count]/F");
  tree->Branch("pfjet_chargedemenergy", pfjet_chargedemenergy, "pfjet_chargedemenergy[pfjet_count]/F");
  tree->Branch("pfjet_chargedmulti", pfjet_chargedmulti, "pfjet_chargedmulti[pfjet_count]/i");	
  tree->Branch("pfjet_neutralmulti", pfjet_neutralmulti, "pfjet_neutralmulti[pfjet_count]/i");	
  tree->Branch("pfjet_energycorr", pfjet_energycorr, "pfjet_energycorr[pfjet_count]/F");
  tree->Branch("pfjet_energycorr_l1fastjet", pfjet_energycorr_l1fastjet, "pfjet_energycorr_l1fastjet[pfjet_count]/F");
  tree->Branch("pfjet_energycorr_l2relative", pfjet_energycorr_l2relative, "pfjet_energycorr_l2relative[pfjet_count]/F");
  tree->Branch("pfjet_energycorr_l3absolute", pfjet_energycorr_l3absolute, "pfjet_energycorr_l3absolute[pfjet_count]/F");
  // tree->Branch("pfjet_energycorr_l2l3residual", pfjet_energycorr_l2l3residual, "pfjet_energycorr_l2l3residual[pfjet_count]/F");
  // tree->Branch("pfjet_pu_jet_cut_loose", pfjet_pu_jet_cut_loose, "pfjet_pu_jet_cut_loose[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_cut_medium", pfjet_pu_jet_cut_medium, "pfjet_pu_jet_cut_medium[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_cut_tight", pfjet_pu_jet_cut_tight, "pfjet_pu_jet_cut_tight[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_cut_mva", pfjet_pu_jet_cut_mva, "pfjet_pu_jet_cut_mva[pfjet_count]/F");
  // tree->Branch("pfjet_pu_jet_simple_loose", pfjet_pu_jet_simple_loose, "pfjet_pu_jet_simple_loose[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_simple_medium", pfjet_pu_jet_simple_medium, "pfjet_pu_jet_simple_medium[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_simple_tight", pfjet_pu_jet_simple_tight, "pfjet_pu_jet_simple_tight[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_simple_mva", pfjet_pu_jet_simple_mva, "pfjet_pu_jet_simple_mva[pfjet_count]/F");
  // tree->Branch("pfjet_pu_jet_full_loose", pfjet_pu_jet_full_loose, "pfjet_pu_jet_full_loose[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_full_medium", pfjet_pu_jet_full_medium, "pfjet_pu_jet_full_medium[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_full_tight", pfjet_pu_jet_full_tight, "pfjet_pu_jet_full_tight[pfjet_count]/O");
  // tree->Branch("pfjet_pu_jet_full_mva", pfjet_pu_jet_full_mva, "pfjet_pu_jet_full_mva[pfjet_count]/F");
  tree->Branch("pfjet_flavour", pfjet_flavour, "pfjet_flavour[pfjet_count]/I");
  tree->Branch("pfjet_btag", pfjet_btag,"pfjet_btag[pfjet_count][10]/F");

  // Met
  tree->Branch("pfmet_ex", &pfmet_ex, "pfmet_ex/F");
  tree->Branch("pfmet_ey", &pfmet_ey, "pfmet_ey/F");
  tree->Branch("pfmet_sigxx", &pfmet_sigxx, "pfmet_sigxx/F");
  tree->Branch("pfmet_sigxy", &pfmet_sigxy, "pfmet_sigxy/F");
  tree->Branch("pfmet_sigyx", &pfmet_sigyx, "pfmet_sigyx/F");
  tree->Branch("pfmet_sigyy", &pfmet_sigyy, "pfmet_sigyy/F");
  
  // Electrons
    if (!cdata) {
    // gen met
    tree->Branch("genmet_ex", &genmet_ex, "genmet_ex/F");
    tree->Branch("genmet_ey", &genmet_ey, "genmet_ey/F");

    // generator info
    tree->Branch("genweight", &genweight, "genweight/F");
    tree->Branch("genid1", &genid1, "genid1/F");
    tree->Branch("genx1", &genx1, "genx1/F");
    tree->Branch("genid2", &genid2, "genid2/F");
    tree->Branch("genx2", &genx2, "genx2/F");
    tree->Branch("genScale", &genScale, "genScale/F");
    
    tree->Branch("numpileupinteractionsminus", &numpileupinteractionsminus, "numpileupinteractionsminus/I");
    tree->Branch("numpileupinteractions", &numpileupinteractions, "numpileupinteractions/I");
    tree->Branch("numpileupinteractionsplus", &numpileupinteractionsplus, "numpileupinteractionsplus/I");
    tree->Branch("numtruepileupinteractions", &numtruepileupinteractions, "numtruepileupinteractions/F");

    // generated particles
    tree->Branch("genparticles_count", &genparticles_count, "genparticles_count/i");
    tree->Branch("genparticles_e", genparticles_e, "genparticles_e[genparticles_count]/F");
    tree->Branch("genparticles_px", genparticles_px, "genparticles_px[genparticles_count]/F");
    tree->Branch("genparticles_py", genparticles_py, "genparticles_py[genparticles_count]/F");
    tree->Branch("genparticles_pz", genparticles_pz, "genparticles_pz[genparticles_count]/F");
    tree->Branch("genparticles_vx", genparticles_vx, "genparticles_vx[genparticles_count]/F");
    tree->Branch("genparticles_vy", genparticles_vy, "genparticles_vy[genparticles_count]/F");
    tree->Branch("genparticles_vz", genparticles_vz, "genparticles_vz[genparticles_count]/F");
    tree->Branch("genparticles_pdgid", genparticles_pdgid, "genparticles_pdgid[genparticles_count]/I");
    tree->Branch("genparticles_status", genparticles_status, "genparticles_status[genparticles_count]/I");
    tree->Branch("genparticles_info", genparticles_info, "genparticles_info[genparticles_count]/i");
    tree->Branch("genparticles_mother", genparticles_mother, "genparticles_mother[genparticles_count]/C");
  }    

  // trigger objects
  tree->Branch("trigobject_count",&trigobject_count,"trigobject_count/i");
  tree->Branch("trigobject_px",trigobject_px,"trigobject_px[trigobject_count]/F");
  tree->Branch("trigobject_py",trigobject_py,"trigobject_py[trigobject_count]/F");
  tree->Branch("trigobject_pz",trigobject_pz,"trigobject_pz[trigobject_count]/F");
  tree->Branch("trigobject_pt",trigobject_pt,"trigobject_pt[trigobject_count]/F");
  tree->Branch("trigobject_eta",trigobject_eta,"trigobject_eta[trigobject_count]/F");
  tree->Branch("trigobject_phi",trigobject_phi,"trigobject_phi[trigobject_count]/F");
  tree->Branch("trigobject_filters",trigobject_filters,"trigobject_filters[trigobject_count][50]/O");

  //add these branches to main tree as well
  tree->Branch("run_hltnames", "std::vector<std::string>", &run_hltnames);
  tree->Branch("run_hltfilters", "std::vector<std::string>",&run_hltfilters);
  tree->Branch("run_btagdiscriminators", "std::vector<std::string>", &run_btagdiscriminators);
  hltriggerresults_ = new std::map<std::string, int>();
  hltriggerprescales_ = new std::map<std::string, int>(); 
  tree->Branch("hltriggerresults", "std::map<std::string, int>", &hltriggerresults_);
  tree->Branch("hltriggerprescales", "std::map<std::string, int>", &hltriggerprescales_);
  tree->Branch("hltriggerresultsV", "std::vector<std::string>", &hltriggerresultsV_);
  
  lumitree = FS->make<TTree>("AC1Blumi", "AC1Blumi", 1);
  lumitree->Branch("lumi_run", &lumi_run, "lumi_run/i");
  lumitree->Branch("lumi_block", &lumi_block, "lumi_block/i");
  
  if(cdata)
    {
      lumitree->Branch("lumi_value", &lumi_value, "lumi_value/F");
      lumitree->Branch("lumi_valueerr", &lumi_valueerr, "lumi_valueerr/F");
      lumitree->Branch("lumi_livefrac", &lumi_livefrac, "lumi_livefrac/F");
      lumitree->Branch("lumi_deadfrac", &lumi_deadfrac, "lumi_deadfrac/F");
      lumitree->Branch("lumi_quality", &lumi_quality, "lumi_quality/i");
      lumitree->Branch("lumi_eventsprocessed", &lumi_eventsprocessed, "lumi_eventsprocessed/i");
      lumitree->Branch("lumi_eventsfiltered", &lumi_eventsfiltered, "lumi_eventsfiltered/i");
      lumitree->Branch("lumi_hltprescaletable", &lumi_hltprescaletable, "lumi_hltprescaletable/i");
      lumitree->Branch("lumi_l1algoprescaletable", &lumi_l1algoprescaletable, "lumi_l1algoprescaletable/i");
      lumitree->Branch("lumi_l1techprescaletable", &lumi_l1techprescaletable, "lumi_l1techprescaletable/i");
    }
  
  runtree = FS->make<TTree>("AC1Brun", "AC1Brun", 1);
  runtree->Branch("run_number", &run_number, "run_number/i");
  runtree->Branch("run_hltcount", &run_hltcount, "run_hltcount/i");
  runtree->Branch("run_hltnames", "std::vector<std::string>", &run_hltnames);
  runtree->Branch("run_hltfilters", "std::vector<std::string>",&run_hltfilters);
  runtree->Branch("run_btagdiscriminators", "std::vector<std::string>", &run_btagdiscriminators);
  runtree->Branch("run_hltprescaletablescount", &run_hltprescaletablescount, "run_hltprescaletablescount/i");
  runtree->Branch("run_hltprescaletables", run_hltprescaletables, "run_hltprescaletables[run_hltprescaletablescount]/i");
  runtree->Branch("run_l1algocount", &run_l1algocount, "run_l1algocount/i");
  runtree->Branch("run_l1algoprescaletablescount", &run_l1algoprescaletablescount, "run_l1algoprescaletablescount/i");
  runtree->Branch("run_l1algoprescaletables", run_l1algoprescaletables, "run_l1algoprescaletables[run_l1algoprescaletablescount]/i");
  runtree->Branch("run_l1techcount", &run_l1techcount, "run_l1techcount/i");
  runtree->Branch("run_l1techprescaletablescount", &run_l1techprescaletablescount, "run_l1techprescaletablescount/i");
  runtree->Branch("run_l1techprescaletables", run_l1techprescaletables, "run_l1techprescaletables[run_l1techprescaletablescount]/i");
} //void NTupleProducerFromMiniAOD::beginJob()



//http://www.boost.org/doc/libs/1_55_0/libs/regex/doc/html/boost_regex/introduction_and_overview.html
void AddTriggerList(const edm::Run& run, 
		    const HLTConfigProvider& hltConfig,
                    const std::vector<std::pair<boost::regex, std::string> >& regexes,
                    std::vector<std::pair<std::string, std::string> >& foundTriggers,
                    std::string& triggerNames   )
{
  for(std::size_t i = 0; i < hltConfig.size(); ++i) {
    
    // In some early 2011 runs saveTagsModules() does not give the
    // modules which actually save tags, but an empty list. So check all
    // the modules instead.
    
    //const std::vector<std::string> saveTagsModules = hltConfig.saveTagsModules(i);
    const std::vector<std::string> saveTagsModules = hltConfig.moduleLabels(i);
    
    for(std::size_t j = 0; j < regexes.size(); ++j) {
      boost::cmatch what;
      
      if(boost::regex_match(hltConfig.triggerName(i), regexes[j].first))
	//if(boost::regex_match(HLTConfiguration.triggerName(i).c_str(), what, muonregexes[j].first) && muontriggers.size() < 32)
	{
	  // Check for filter
	  if(regexes[j].second.size() != 0)
	      {
		// First, check for leg
		std::string::size_type legpos = regexes[j].second.find('|');
		std::string leg, filters;
		if(legpos != std::string::npos)
		  {
		    leg     = regexes[j].second.substr(0, legpos); 
		    filters = regexes[j].second.substr(legpos + 1); 
		  }
		else
		  {
		    filters = regexes[j].second;
		  }
		
		std::vector<std::string> strs;
		boost::split(strs, filters, boost::is_any_of(","));
		bool foundFilter = false;
		std::string filter;
		
		for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		  {
		    for(std::size_t l = 0; l < strs.size() && !foundFilter; ++l)
		      {
			if(saveTagsModules[k] == strs[l])
			  {
			    filter = saveTagsModules[k];
			    foundFilter = true;
			  }
		      }
		  }//for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		
		if(!foundFilter)
		  {
		    for(std::size_t k = 0; k < saveTagsModules.size() && !foundFilter; ++k)
		      std::cout << saveTagsModules[k] << std::endl;
		    throw cms::Exception("NTupleProducerFromMiniAOD") << "Did not find filter for trigger " << hltConfig.triggerName(i) << " in run " << run.run() << std::endl;
		  }
		
		if(leg.empty())
		  triggerNames += hltConfig.triggerName(i) + string(":") + filter + string(" ");
		else
		  triggerNames += hltConfig.triggerName(i) + string("|") + leg + string(":") + filter + string(" ");
		
		foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), filter));
	      }//if(regexes[j].second.size() != 0)
	  
	  else
	    {
	      triggerNames += hltConfig.triggerName(i) + string(" ");
	      foundTriggers.push_back(std::make_pair(hltConfig.triggerName(i), saveTagsModules.back()));
	    }
	}
    }//for(std::size_t j = 0; j < regexes.size(); ++j)
  }//for(std::size_t i = 0; i < hltConfig.size(); ++i)
}// std::string& triggerNames


void NTupleProducerFromMiniAOD::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", TTrackBuilder);
  
  run_number = iRun.run();
  //L1 prescales
  edm::ESHandle<L1GtPrescaleFactors> l1GtPfAlgo;
  iSetup.get<L1GtPrescaleFactorsAlgoTrigRcd>().get(l1GtPfAlgo);
  
  unsigned numl1algo = (l1GtPfAlgo.product()->gtPrescaleFactors())[0].size();
  unsigned numl1algotables = (l1GtPfAlgo.product()->gtPrescaleFactors()).size();
  
  run_l1algoprescaletablescount = numl1algo*numl1algotables;
  run_l1algocount = numl1algo;
  if(l1GtPfAlgo.isValid())
    {
      for(unsigned i = 0 ; i < numl1algotables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1algo ; j++)
	    {
	      run_l1algoprescaletables[j+numl1algo*i] = (l1GtPfAlgo.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }
  
  edm::ESHandle<L1GtPrescaleFactors> l1GtPfTech;
  iSetup.get<L1GtPrescaleFactorsTechTrigRcd>().get(l1GtPfTech);
  
  unsigned numl1tech = (l1GtPfTech.product()->gtPrescaleFactors())[0].size();
  unsigned numl1techtables = (l1GtPfTech.product()->gtPrescaleFactors()).size();
  
  run_l1techprescaletablescount = numl1tech*numl1techtables;
  run_l1techcount = numl1tech;
  if(l1GtPfTech.isValid())
    {
      for(unsigned i = 0 ; i < numl1techtables ; i++)
	{
	  for(unsigned j = 0 ; j < numl1tech ; j++)
	    {
	      run_l1techprescaletables[j+numl1tech*i] = (l1GtPfTech.product()->gtPrescaleFactors())[i][j];
	    }
	}
    }
  
  //HLT names and prescales
  muontriggers.clear();
  jettriggers.clear();
  
  bool changed = true;
  HLTConfiguration.init(iRun, iSetup, cTriggerProcess, changed);
  
  vector<pair<boost::regex, string> > muonregexes;
  for(unsigned i = 0 ; i < cMuHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cMuHLTriggerMatching[i], boost::is_any_of(":"));
      if(strs.size() == 1) strs.push_back(string(""));
      muonregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }

  vector<pair<boost::regex, string> > jetregexes;
  for(unsigned i = 0 ; i < cJetHLTriggerMatching.size() ; i++)
    {
      vector<string> strs;
      boost::split(strs, cJetHLTriggerMatching[i], boost::is_any_of(":"));	
      if(strs.size() == 1) strs.push_back(string(""));
      jetregexes.push_back(pair<boost::regex, string>(boost::regex(strs[0].c_str()), strs[1]));
    }
  
  run_hltcount = HLTConfiguration.size();
  string allnames;
  string allmuonnames;
  string alljetnames;
  
  AddTriggerList(iRun, HLTConfiguration, muonregexes,     muontriggers,     allmuonnames);
  AddTriggerList(iRun, HLTConfiguration, jetregexes,      jettriggers,      alljetnames);
  
  run_hltnames.clear();
  for (unsigned int i=0; i<cHLTriggerPaths.size(); ++i)
    run_hltnames.push_back(cHLTriggerPaths.at(i));
	
  if(muontriggers.size() > 32) throw cms::Exception("NTupleProducerFromMiniAOD") << "Too many muon triggers!" << std::endl;
  if(jettriggers.size() > 32) throw cms::Exception("NTupleProducerFromMiniAOD") << "Too many jet triggers!" << std::endl;

  // adding all filters
  run_hltfilters.clear();
  for (unsigned int i = 0; i<muontriggers.size(); ++i) { 
    run_hltfilters.push_back(muontriggers.at(i).second);
  }
  for (unsigned int i = 0; i<jettriggers.size(); ++i) {
    run_hltfilters.push_back(jettriggers.at(i).second);
  }
  
  if (run_hltfilters.size()>50) throw cms::Exception("NTupleProducerFromMiniAOD") << "Too many HLT filters!" << std::endl;

  run_hltprescaletablescount = HLTConfiguration.prescaleSize()*HLTConfiguration.size();

  // adding btag discriminators
  run_btagdiscriminators.clear();
  for (unsigned int i=0; i<cBtagDiscriminators.size(); ++i)
    run_btagdiscriminators.push_back(cBtagDiscriminators.at(i));

  if (run_btagdiscriminators.size()>10) throw cms::Exception("NTupleProducerFromMiniAOD") << "Too many btag discriminators!" << std::endl;

  for(unsigned j = 0 ; j < HLTConfiguration.prescaleSize() ; j++)
    {
      for(unsigned i = 0 ; i < HLTConfiguration.size() ; i++)
	{
	  run_hltprescaletables[i+HLTConfiguration.size()*j] = HLTConfiguration.prescaleValue(j, HLTConfiguration.triggerName(i));
	}
    }	
  runtree->Fill();
}//void NTupleProducerFromMiniAOD::beginRun

void NTupleProducerFromMiniAOD::beginLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumi_run = iLumiBlock.run();
  lumi_block = iLumiBlock.luminosityBlock();
  
  if(cdata)
    {
      //edm::Handle<LumiSummary> lumiSummary;
      //iLumiBlock.getByLabel(edm::InputTag("lumiProducer"), lumiSummary);
      lumi_value = 0;//lumiSummary->avgInsDelLumi();
      lumi_valueerr = 0;//lumiSummary->avgInsDelLumiErr();
      lumi_livefrac = 0;//lumiSummary->lumiSecQual();
      lumi_deadfrac = 0;//lumiSummary->deadFrac();
      lumi_quality = 0;//lumiSummary->liveFrac();
      lumi_eventsprocessed = 0;
      lumi_eventsfiltered = 0;
    }
}//void NTupleProducerFromMiniAOD::beginLuminosityBlock

void NTupleProducerFromMiniAOD::endLuminosityBlock(const edm::LuminosityBlock& iLumiBlock, const edm::EventSetup& iSetup)
{
  lumitree->Fill();
}


void NTupleProducerFromMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if(doDebug)  cout<<"inside the analyze function"<< endl; 

  track_count = 0;
  primvertex_count = 0;
  muon_count = 0;
  pfjet_count = 0;
  genparticles_count = 0;
  errors = 0;
  trigobject_count = 0;
  //bool takeevent = true;

  nEvents->Fill(0);
  pv_position = math::XYZPoint(0.,0.,0.);
  
  lumi_eventsprocessed++;
  
  event_nr      = iEvent.id().event();
  event_run      = iEvent.id().run();
  event_timeunix = iEvent.time().unixTime();
  event_timemicrosec = iEvent.time().microsecondOffset();
  event_luminosityblock = iEvent.getLuminosityBlock().luminosityBlock();

  // L1TriggerBits
  // https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h?v=CMSSW_6_2_0_SLHC2#04
  edm::Handle<L1GlobalTriggerReadoutRecord> L1trigger;
  iEvent.getByLabel(edm::InputTag("gtDigis"), L1trigger);
  assert(L1trigger.isValid());
  
  const TechnicalTriggerWord& L1triggerbits = L1trigger->technicalTriggerWord();
  for(int i  = 0  ; i < 8 ; i++) trigger_level1bits[i] = 0;
  
  for(unsigned i = 0 ; i < min(unsigned(L1triggerbits.size()), unsigned(64)) ; i++)
    trigger_level1bits[i/8] |= (Byte_t)L1triggerbits[i] << (i % 8);  // bitwise OR -> | 
  
  //trigger_level1bits[i/8] = trigger_level1bits[i/8]  | (Byte_t) L1triggerbits[i] << (i % 8);  // bitwise OR -> | 
  
  
  //L1TriggerAlgos
  const DecisionWord& L1triggeralgos = L1trigger->decisionWord();
  for(int i = 0  ; i < 128 ; i++){trigger_level1[i] = 0;}
  for(unsigned i = 0 ; i < min(unsigned(L1triggeralgos.size()), unsigned(1024)) ; i++)
    {
      trigger_level1[i/8] |= (Byte_t)L1triggeralgos[i] << (i%8);
    }
  lumi_l1techprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexTech();
  lumi_l1algoprescaletable = (L1trigger->gtFdlWord()).gtPrescaleFactorIndexAlgo();	
  lumi_hltprescaletable = -1;
  
  //HLTriggerResults
  iEvent.getByLabel(edm::InputTag("TriggerResults", "", cTriggerProcess), HLTrigger);
  assert(HLTrigger.isValid());
  for(int i = 0  ; i < 128 ; i++){trigger_HLT[i] = 0;}
  

  //store trigger bits for selected trigger paths
  hltriggerresults_->clear();
  hltriggerprescales_->clear();
  hltriggerresultsV_.clear();
  const edm::TriggerNames& TrigNames_ = iEvent.triggerNames(*HLTrigger);
  for(unsigned i = 0 ; i < HLTrigger->size(); i++)
    {
      if(!HLTrigger->wasrun(i) )continue;
      std::string trigName=TrigNames_.triggerName(i);
      if(cHLTriggerPaths.size() > 0){
  	for(size_t ip = 0; ip < cHLTriggerPaths.size(); ip++){
  	  if(trigName.find(cHLTriggerPaths[ip]) != string::npos){
	    hltriggerprescales_->insert(std::pair<string, int>(trigName, HLTConfiguration.prescaleValue(iEvent,iSetup,trigName)));
  	    hltriggerresults_->insert(std::pair<string, int>(trigName, HLTrigger->accept(i)));
	    TString TriggerName(trigName);
	    //	 std::cout << trigName << " : " 
	    //   << HLTrigger->accept(i) << " ; prescale : " 
	    //   << HLTConfiguration.prescaleValue(iEvent,iSetup,trigName) << std::endl;
  	    if(HLTrigger->accept(i)) hltriggerresultsV_.push_back(trigName);
  	  }
  	}
      }
    }

  // beam spot  
  edm::Handle<BeamSpot> TheBeamSpot;
  iEvent.getByLabel(BeamSpotTag_, TheBeamSpot);
  if(TheBeamSpot.isValid())
    {
      beamspot_x = TheBeamSpot->x0();
      beamspot_y = TheBeamSpot->y0();
      beamspot_z = TheBeamSpot->z0();
      beamspot_xwidth = TheBeamSpot->BeamWidthX();
      beamspot_ywidth = TheBeamSpot->BeamWidthY();
      beamspot_zsigma = TheBeamSpot->sigmaZ();
      beamspot_cov[0] = TheBeamSpot->covariance(0,0); 
      beamspot_cov[1] = TheBeamSpot->covariance(0,1); 
      beamspot_cov[2] = TheBeamSpot->covariance(0,2); 
      beamspot_cov[3] = TheBeamSpot->covariance(1,1); 
      beamspot_cov[4] = TheBeamSpot->covariance(1,2); 
      beamspot_cov[5] = TheBeamSpot->covariance(2,2); 
      pv_position = math::XYZPoint(TheBeamSpot->x0(), TheBeamSpot->y0(), TheBeamSpot->z0());
    }
  else
    {
      beamspot_x = 0.;
      beamspot_y = 0.;
      beamspot_z = 0.;
      beamspot_xwidth = 0.;
      beamspot_ywidth = 0.;
      beamspot_zsigma = 0.;
      beamspot_cov[0] = 0.;
      beamspot_cov[1] = 0.;
      beamspot_cov[2] = 0.;
      beamspot_cov[3] = 0.;
      beamspot_cov[4] = 0.;
      beamspot_cov[5] = 0.;
    }

  // primary vertex
  edm::Handle<VertexCollection> Vertex;
  iEvent.getByLabel(PVTag_, Vertex);
  if(Vertex.isValid())
    {
      for(unsigned i = 0 ; i < Vertex->size(); i++)
	{
	  if((*Vertex)[i].isValid() && !(*Vertex)[i].isFake())
	    {
	      if((*Vertex)[i].ndof() >= 4 && (*Vertex)[i].z() > -24 && (*Vertex)[i].z() < 24 && (*Vertex)[i].position().Rho() < 2.){
		if(primvertex_count == 0)
		  {
		    primvertex_x = (*Vertex)[i].x();
		    primvertex_y = (*Vertex)[i].y();
		    primvertex_z = (*Vertex)[i].z();
		    primvertex_chi2 = (*Vertex)[i].chi2();
		    primvertex_ndof = (*Vertex)[i].ndof();
		    primvertex_ntracks = (*Vertex)[i].tracksSize();
		    primvertex_cov[0] = (*Vertex)[i].covariance(0,0); // xError()
		    primvertex_cov[1] = (*Vertex)[i].covariance(0,1); 
		    primvertex_cov[2] = (*Vertex)[i].covariance(0,2);
		    primvertex_cov[3] = (*Vertex)[i].covariance(1,1); // yError()
		    primvertex_cov[4] = (*Vertex)[i].covariance(1,2);
		    primvertex_cov[5] = (*Vertex)[i].covariance(2,2); // zError()
		    Float_t ptq = 0.;
		    for(Vertex::trackRef_iterator it = (*Vertex)[i].tracks_begin() ; it != (*Vertex)[i].tracks_end() ; ++it)
		      {
			ptq += (*it)->pt() * (*it)->pt();
		      }
		    primvertex_ptq = ptq;
		    
		    pv_position = (*Vertex)[i].position();
		    primvertex = (*Vertex)[i];
		  }
		
		primvertex_count++;
	      }
	    }
	}
    }

  // muons
  int numberOfMuons = int(AddMuons(iEvent));
  if (numberOfMuons<cMuNum) return;

  // Packed PF Candidates
  int numberOfTracks = int(AddPackedPFCand(iEvent));
  if (numberOfTracks<cTrackNum) return;

  // pf jets
  int numberOfJets = int(AddPFJets(iEvent,iSetup));
  if (numberOfJets<cJetNum) return;

  // pf met 
  edm::Handle<pat::METCollection> patMet;
  iEvent.getByLabel(MetCollectionTag_, patMet);
  
  assert(patMet->size() > 0);
  pfmet_ex = (*patMet)[0].px();
  pfmet_ey = (*patMet)[0].py();
  
  pfmet_sigxx = (*patMet)[0].getSignificanceMatrix()(0,0);
  pfmet_sigxy = (*patMet)[0].getSignificanceMatrix()(0,1);
  pfmet_sigyx = (*patMet)[0].getSignificanceMatrix()(1,0);
  pfmet_sigyy = (*patMet)[0].getSignificanceMatrix()(1,1);
  
  if (!cdata) {
    const reco::GenMET * genMET = (*patMet)[0].genMET();
    genmet_ex = genMET->px();
    genmet_ey = genMET->py();
  }

  // rho neutral
  edm::Handle<double> rho;
  iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral"), rho);
  assert(rho.isValid());
  rhoNeutral = *rho;

  //  std::cout << "rhoNeutral = " << rhoNeutral << std::endl;

  genweight = 1.;
  numpileupinteractionsminus = -1;
  numpileupinteractions      = -1;
  numpileupinteractionsplus  = -1;
  numtruepileupinteractions  = -1.0f;
  hepNUP_ = -1;

  // generator info and generated particles 
  if(!cdata)
    {
      //bool haveGenParticles = AddGenParticles(iEvent);

      edm::Handle<GenEventInfoProduct> HEPMC;
      iEvent.getByLabel(edm::InputTag("generator"), HEPMC);
      if(HEPMC.isValid())
	{
	  genweight = HEPMC->weight();
	  genid1 = HEPMC->pdf()->id.first;
	  genx1 = HEPMC->pdf()->x.second;
	  genid2 = HEPMC->pdf()->id.first;
	  genx2 = HEPMC->pdf()->x.second;
	  genScale = HEPMC->qScale();
	}

      edm::Handle<vector<PileupSummaryInfo> > PUInfo;
      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PUInfo);
      if(PUInfo.isValid())
	{
	  for(vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI)
	    {
	      int BX = PVI->getBunchCrossing();
	      if(BX == -1)
		{ 
		  numpileupinteractionsminus = PVI->getPU_NumInteractions();
		}
	      else if(BX == 0)
		{ 
		  numpileupinteractions = PVI->getPU_NumInteractions();
		}
	      else if(BX == 1)
		{ 
		  numpileupinteractionsplus = PVI->getPU_NumInteractions();
		}
	      
	      numtruepileupinteractions = PVI->getTrueNumInteractions();
	    }
	}
    } // cgen
    
  
  // trigger objects
  //int numberOfTriggerObjects = int(AddTriggerObjects(iEvent));

  // tracks to be added

  tree->Fill();

} //void NTupleProducerFromMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


Int_t NTupleProducerFromMiniAOD::HasAnyMother(const GenParticle* particle, int id)
{
  vector<unsigned> bknummother;
  vector<const GenParticle*> bkparticle;
  bknummother.reserve(10);
  bkparticle.reserve(10);
  int level = 0;
  bkparticle.push_back(particle);
  bknummother.push_back(0);
  
  unsigned j = 0;
  while(true)
    {
      if(j == bkparticle[level]->numberOfMothers())
	{
	  level--;
	  if(level == -1){return(0);}
	  j = bknummother[level];
	  bkparticle.resize(level+1);
	  bknummother.resize(level+1);
	  continue;
	}
      
      if(bkparticle[level]->mother(j)->pdgId() == id) return(2);
      if(abs(bkparticle[level]->mother(j)->pdgId()) == abs(id)) return(1);
      
      if(bkparticle[level]->mother(j)->numberOfMothers() > 0)
	{
	  bknummother[level] = j+1;
	  bkparticle.push_back(dynamic_cast<const GenParticle*>(bkparticle[level]->mother(j)));
	  bknummother.push_back(0);
	  j = 0;
	  level++;
	  continue;
	}
      j++;
    }
  return(0);
} // Int_t NTupleProducerFromMiniAOD::HasAnyMother(const GenParticle* particle, int id)

void NTupleProducerFromMiniAOD::endJob()
{
}

bool NTupleProducerFromMiniAOD::AddGenParticles(const edm::Event& iEvent) {

  edm::Handle<reco::GenParticleCollection> GenParticles;
  iEvent.getByLabel(GenParticleCollectionTag_, GenParticles);

  bool passed = false;
  
  if(GenParticles.isValid())
    {
      passed = true;
      for(unsigned i = 0 ; i < GenParticles->size() ; i++)
	{
	  bool fill = false;
	  UInt_t info = 0;
	  string mother = "";
	  if(abs((*GenParticles)[i].pdgId()) == 13 && (*GenParticles)[i].pt() > 8. && (*GenParticles)[i].status()==1) // gen muons
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother="Z";}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother="W";}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother="tau";}
	      //	      std::cout << "GenMuon : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta() 
	      //			<< "   phi = " << (*GenParticles)[i].phi() 
	      //			<< "   status = " << (*GenParticles)[i].status() 
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 11 && (*GenParticles)[i].pt() > 8. && (*GenParticles)[i].status()==1) // gen electrons
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother="Z";}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother="W";}
	      if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2; mother="tau";}
	      //	      std::cout << "GenElectron : " << (*GenParticles)[i].pdgId()
	      //			<< "   pt = " << (*GenParticles)[i].pt()
	      //			<< "   eta = " << (*GenParticles)[i].eta() 
	      //			<< "   phi = " << (*GenParticles)[i].phi() 
	      //			<< "   status = " << (*GenParticles)[i].status() 
	      //			<< "   mother = " << mother << std::endl;
	    }
	  else if(abs((*GenParticles)[i].pdgId()) == 15 && (*GenParticles)[i].pt() > 10.) // gen taus
	    {
	      fill = true;
	      if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0; mother="Z";}
	      if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1; mother="W";}
	    }
	  else if( (abs((*GenParticles)[i].pdgId()) == 16 || abs((*GenParticles)[i].pdgId()) == 14 || abs((*GenParticles)[i].pdgId()) == 12)) // gen v's
	    {
	      if ((*GenParticles)[i].status()==1) {
		fill = true;
		if(HasAnyMother(&(*GenParticles)[i], 23) > 0 || HasAnyMother(&(*GenParticles)[i], 22) > 0) {info |= 1<<0;mother="Z";}
		if(HasAnyMother(&(*GenParticles)[i], 24) > 0) {info |= 1<<1;mother="W";}
		if(HasAnyMother(&(*GenParticles)[i], 15) > 0) {info |= 1<<2;mother="tau";}
		//		std::cout << "GenNeutrino : " << (*GenParticles)[i].pdgId() 
		//		 	  << "   pt = " << (*GenParticles)[i].pt() 
		//		 	  << "   eta = " << (*GenParticles)[i].eta() 
		//		 	  << "   phi = " << (*GenParticles)[i].phi()
		//			  << "   status = " << (*GenParticles)[i].status() 
		//		 	  << "   mother =  " << mother << std::endl;

	      }
	    }
	  // Save all partons (quarks)
	  else if(abs((*GenParticles)[i].pdgId()) < 6 && (*GenParticles)[i].pt() > 20.) // quarks
	    {
	      fill = true;
	    }
	  // Save all tops from Madgraph
	  else if ( abs((*GenParticles)[i].pdgId()) == 6 && (*GenParticles)[i].status()==62 )  // top-quarks
	    {
	      fill = true;
	      //	      	std::cout << "GenTop : " << (*GenParticles)[i].pdgId() 
	      //			  << "   pt = " << (*GenParticles)[i].pt() 
	      //		 	  << "   eta = " << (*GenParticles)[i].eta()
	      //		 	  << "   phi = " << (*GenParticles)[i].phi() 
	      //		 	  << "   status = " << (*GenParticles)[i].status() << std::endl;
	    }
	  // Save all partons (gluons)
	  else if(abs((*GenParticles)[i].pdgId()) == 21 && (*GenParticles)[i].pt() > 20.) // gluons
	    {
	      fill = true;
	    }
	  // Save all W/Z bosons from Madgraph
	  else if(abs((*GenParticles)[i].pdgId()) == 23 || abs((*GenParticles)[i].pdgId()) == 24 )
	    {
	      if ( (*GenParticles)[i].status()==62 || (*GenParticles)[i].status()==52 ) {
		fill = true;
		//		std::cout << "GenBoson : " << (*GenParticles)[i].pdgId() 
		//			  << "   pt = " << (*GenParticles)[i].pt() 
		//		 	  << "   eta = " << (*GenParticles)[i].eta()
		//		 	  << "   phi = " << (*GenParticles)[i].phi() 
		//		 	  << "   status = " << (*GenParticles)[i].status() << std::endl;
	      }
	    }
	  if(fill)
	    {
	      genparticles_e[genparticles_count] = (*GenParticles)[i].energy();
	      genparticles_px[genparticles_count] = (*GenParticles)[i].px();
	      genparticles_py[genparticles_count] = (*GenParticles)[i].py();
	      genparticles_pz[genparticles_count] = (*GenParticles)[i].pz();
	      genparticles_vx[genparticles_count] = (*GenParticles)[i].vx();
	      genparticles_vy[genparticles_count] = (*GenParticles)[i].vy();
	      genparticles_vz[genparticles_count] = (*GenParticles)[i].vz();
	      genparticles_pdgid[genparticles_count] = (*GenParticles)[i].pdgId();
	      genparticles_status[genparticles_count] = (*GenParticles)[i].status();
	      genparticles_info[genparticles_count] = info;
	      genparticles_mother[genparticles_count] = mother;
              std::cout << "Gen_mother : " << genparticles_mother[genparticles_count] << std::endl;
	      genparticles_count++;
	    }
	} // for(unsigned i = 0 ; i < GenParticles->size() ; i++)
    } // if(GenParticles.isValid())
  
  return passed;

} // bool NTupleProducerFromMiniAOD::AddGenParticles(const edm::Event& iEvent) 

unsigned int NTupleProducerFromMiniAOD::AddMuons(const edm::Event& iEvent)
{

  edm::Handle<pat::MuonCollection> Muons;
  //	iEvent.getByLabel(edm::InputTag("muons"), Muons);
  iEvent.getByLabel(MuonCollectionTag_, Muons);
  
  if(Muons.isValid())
    {
      for(unsigned i = 0 ; i < Muons->size() ; i++){
	
	if ((*Muons)[i].pt() < cMuPtMin) continue;
	if (fabs(((*Muons)[i].eta()))>cMuEtaMax) continue;

	//	std::cout << "Selected pat::Muon " << i << std::endl;
	
	muon_px[muon_count] = (*Muons)[i].px();
	muon_py[muon_count] = (*Muons)[i].py();
	muon_pz[muon_count] = (*Muons)[i].pz();
	muon_pt[muon_count] = (*Muons)[i].pt();
	muon_eta[muon_count] = (*Muons)[i].eta();
	muon_phi[muon_count] = (*Muons)[i].phi();
	muon_charge[muon_count] = (*Muons)[i].charge();

	if((*Muons)[i].globalTrack().isNonnull())
	  {
	    muon_globalTrack[muon_count] = true;
	    muon_pterror[muon_count] = (*Muons)[i].globalTrack()->ptError();
	    muon_chi2[muon_count] = (*Muons)[i].globalTrack()->chi2();
	    muon_ndof[muon_count] = (*Muons)[i].globalTrack()->ndof();
	    muon_nMuonHits[muon_count] = (*Muons)[i].globalTrack()->hitPattern().numberOfValidMuonHits();
	  }
	else
	  {
	    muon_globalTrack[muon_count] = false;
	    muon_pterror[muon_count] = -1.;
	    muon_chi2[muon_count] = -1.;
	    muon_ndof[muon_count] = 0;
	    muon_nMuonHits[muon_count] = 0;
	  }

	//	std::cout << "  chi2 = " << muon_chi2[muon_count] << "  ndof = " << muon_ndof[muon_count] << std::endl;

	muon_nMuonStations[muon_count] = (*Muons)[i].numberOfMatchedStations();
	
	muon_isTracker[muon_count] = (*Muons)[i].isTrackerMuon();
	muon_isPF[muon_count] = (*Muons)[i].isPFMuon();
	muon_isTight[muon_count] = (*Muons)[i].isTightMuon(primvertex); 
	muon_isLoose[muon_count] = (*Muons)[i].isLooseMuon();
	  
	muon_chargedHadIso[muon_count] = (*Muons)[i].chargedHadronIso();
	muon_neutralHadIso[muon_count] = (*Muons)[i].neutralHadronIso();
	muon_photonIso[muon_count] = (*Muons)[i].photonIso();
	muon_puIso[muon_count] = (*Muons)[i].puChargedHadronIso();

	TrackRef innertrack = (*Muons)[i].innerTrack();

	if(innertrack.isNonnull())
	  {
	    TransientTrack TTrack = TTrackBuilder->build(innertrack);
	    TrajectoryStateClosestToPoint TTrackState = TTrack.trajectoryStateClosestToPoint(GlobalPoint(pv_position.x(), pv_position.y(), pv_position.z()));
	    muon_innerTrack[muon_count] = true;
	    muon_dxy[muon_count]    = TTrackState.perigeeParameters().transverseImpactParameter();
	    muon_dxyerr[muon_count] = TTrackState.perigeeError().transverseImpactParameterError();
	    muon_dz[muon_count]     = TTrackState.perigeeParameters().longitudinalImpactParameter();
	    muon_dzerr[muon_count]  = TTrackState.perigeeError().longitudinalImpactParameterError();
	    muon_nPixelHits[muon_count] = innertrack->hitPattern().numberOfValidPixelHits();
	    muon_nTrackerHits[muon_count] = innertrack->hitPattern().trackerLayersWithMeasurement();
	  }
	else 
	  {
	    muon_innerTrack[muon_count] = false;
	    muon_dxy[muon_count] = -9999;
	    muon_dxyerr[muon_count] = -9999;
	    muon_dxy[muon_count] = -9999;
	    muon_dxyerr[muon_count] = -9999;
	    muon_nPixelHits[muon_count] = 0; 
	    muon_nTrackerHits[muon_count] = 0;
	  }
	
	muon_count++;
	
	if (muon_count==M_muonmaxcount) {
	  cerr << "number of muons > M_muonmaxcount. They are missing." << endl; errors |= 1<<1; 
	  break;
	}

      }

    }
  return muon_count;
}//unsigned int NTupleProducerFromMiniAOD::AddMuons

unsigned int NTupleProducerFromMiniAOD::AddPackedPFCand(const edm::Event& iEvent)
{

  edm::Handle<pat::PackedCandidateCollection> Tracks;
  iEvent.getByLabel(TrackCollectionTag_, Tracks);
  if(Tracks.isValid())
    {
      for(unsigned i = 0 ; i < Tracks->size() ; i++){
	if ((*Tracks)[i].pt() < cTrackPtMin) continue;
	if (fabs(((*Tracks)[i].eta())) > cTrackEtaMax) continue;
        if (fabs(((*Tracks)[i].charge())) < 0.5) continue;
       	track_px[track_count] = (*Tracks)[i].px();
	track_py[track_count] = (*Tracks)[i].py();
	track_pz[track_count] = (*Tracks)[i].pz();
	track_pt[track_count] = (*Tracks)[i].pt();
	track_eta[track_count] = (*Tracks)[i].eta();
	track_phi[track_count] = (*Tracks)[i].phi();
	track_charge[track_count] = (*Tracks)[i].charge();
        track_mass[track_count] = (*Tracks)[i].mass();
        track_dxy[track_count] = (*Tracks)[i].dxy();
        track_dz[track_count] = (*Tracks)[i].dz();
        track_dxyerr[track_count] = (*Tracks)[i].dxyError();
        track_dzerr[track_count] = (*Tracks)[i].dzError();
        track_ID[track_count] = (*Tracks)[i].pdgId();
  	track_count++;
	
	if (track_count==M_trackmaxcount) {
	  cerr << "number of tracks > M_trackmaxcount. They are missing." << endl; errors |= 1<<1; 
	  break;
	}
					   }
		
    }
  return track_count;
}

bool NTupleProducerFromMiniAOD::GetL1ExtraTriggerMatch(const l1extra::L1JetParticleCollection* l1jets,  
						    const l1extra::L1JetParticleCollection* l1taus, 
						    const LeafCandidate& leg2) 
{
  bool matched = false;
  //check matching to l1tau 44 or l1jet 64
  if(l1taus)
    {
      //check matching with l1tau Pt>44 |eta|<2.172
      matched = false;
      for(unsigned int i=0; i<l1taus->size(); ++i)
	{
	  if( (*l1taus)[i].pt() < 44 || fabs((*l1taus)[i].eta() ) > 2.172 ) continue;
	  if( ROOT::Math::VectorUtil::DeltaR( (*l1taus)[i].p4(), leg2.p4() )  < 0.5 )
	    {
	      matched = true;
	      break;
	    }
	}// for(unsigned int i=0; i<l1taus->size(); ++i) 
      
      if(!matched)
	{ 
	  if(l1jets){//check matching with l1jet Pt>64 |eta|<2.172
	    for(unsigned int i=0; i < l1jets->size(); ++i)
	      {
		if( (*l1jets)[i].pt() < 64 || fabs((*l1jets)[i].eta() ) > 2.172 ) continue;
		if( ROOT::Math::VectorUtil::DeltaR((*l1jets)[i].p4(), leg2.p4() ) < 0.5 ) {
		  matched = true;
		  break;
		}
	      }//for(unsigned int i=0; i<l1jets->size(); ++i)
	  }
	}
    } //if(l1taus)
  return matched;
}//bool NTupleProducerFromMiniAOD::GetL1ExtraTriggerMatch

unsigned int NTupleProducerFromMiniAOD::AddTriggerObjects(const edm::Event& iEvent) {

  // trigger objects
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByLabel(TriggerObjectCollectionTag_, triggerObjects);
  assert(triggerObjects.isValid());
  
  for (unsigned int iTO=0; iTO<triggerObjects->size(); ++iTO) {
    vector<string> filterLabels = (*triggerObjects)[iTO].filterLabels();
    bool matchFound = false;
    std::vector<bool> passedFilters; passedFilters.clear();
    std::vector<std::string> matchedFilters; matchedFilters.clear();
    for (unsigned int n=0; n < run_hltfilters.size(); ++n) {
      TString HltFilter(run_hltfilters.at(n));
      bool thisMatched = false;
      for (unsigned int i=0; i < filterLabels.size(); ++i) {
	TString FilterName(filterLabels.at(i));
	if (HltFilter==FilterName) {
	  matchFound = true;
	  thisMatched = true;
	  matchedFilters.push_back(filterLabels.at(i));
	  break;
	}
      }
      passedFilters.push_back(thisMatched);
    } 
    if (matchFound) {
      // std::cout << "   trigger object " << iTO 
      // 		<< "   pt = " << (*triggerObjects)[iTO].pt() 
      // 		<< "   eta = " << (*triggerObjects)[iTO].eta() 
      // 		<< "   phi = " << (*triggerObjects)[iTO].phi() << std::endl;
      // for (unsigned int ifilter=0; ifilter<matchedFilters.size(); ++ifilter)
      // 	std::cout << "    " << matchedFilters[ifilter] << std::endl;
      for (unsigned int n=0; n < 50; ++n) {
	if (n<passedFilters.size())
	  trigobject_filters[trigobject_count][n] = passedFilters.at(n);
	else
	  trigobject_filters[trigobject_count][n] = false;
      }
      trigobject_px[trigobject_count] = (*triggerObjects)[iTO].px();
      trigobject_py[trigobject_count] = (*triggerObjects)[iTO].py();
      trigobject_pz[trigobject_count] = (*triggerObjects)[iTO].pz();
      trigobject_pt[trigobject_count] = (*triggerObjects)[iTO].pt();
      trigobject_eta[trigobject_count] = (*triggerObjects)[iTO].eta();
      trigobject_phi[trigobject_count] = (*triggerObjects)[iTO].phi();
      trigobject_count++;
      if (trigobject_count==M_trigobjectmaxcount) {
	 cerr << "number of trigger objects > M_trigobjectmaxcount. They are missing." << endl; 
	 errors |= 1<<5; 
	 break;
      }
    }
  }
  return trigobject_count;
}//unsigned int NTupleProducerFromMiniAOD::AddTriggerObjects

NTupleProducerFromMiniAOD::DCA NTupleProducerFromMiniAOD::calculateDCA(const pat::Tau& tau1, const pat::Tau& tau2)
{
        // TODO: Use the reconstructed decay mode, and then only the mass of the decay products?
	float tauMass = 1.777;
	// TODO: What is this sigma?
	float tauSigma = tauMass*1e-6;

	DCA dca = { -1.0, -1.0f, -1.0f, -1.0f };
	reco::TrackRef track1 = tau1.leadTrack();
	reco::TrackRef track2 = tau2.leadTrack();
	if(track1.isNonnull() && track2.isNonnull())
	{
		reco::TransientTrack transientTrack1 = TTrackBuilder->build(*track1);
		reco::TransientTrack transientTrack2 = TTrackBuilder->build(*track2);
		if(transientTrack1.impactPointTSCP().isValid() && transientTrack2.impactPointTSCP().isValid())
		{
			FreeTrajectoryState state1 = transientTrack1.impactPointTSCP().theState();
			FreeTrajectoryState state2 = transientTrack2.impactPointTSCP().theState();
			TwoTrackMinimumDistance minDist;
			minDist.calculate(state1, state2);
			if(minDist.status())
			{
				const float dist3D = minDist.distance();
				std::pair<GlobalPoint,GlobalPoint> pcas = minDist.points();
				GlobalPoint pca1 = pcas.first;
				GlobalPoint pca2 = pcas.second;

				ROOT::Math::SVector<double, 3> distanceVector(pca1.x()-pca2.x(), pca1.y()-pca2.y(), pca1.z()-pca2.z());
				const float twoTauDist3D = ROOT::Math::Mag(distanceVector);
				assert(fabs(dist3D - twoTauDist3D) < 1e-3);

				float chi2 = 0.0f, ndf = 0.0f;
				KinematicParticleFactoryFromTransientTrack pFactory;
				RefCountedKinematicParticle tau1Particle = pFactory.particle(transientTrack1, tauMass, chi2, ndf, tauSigma);
				RefCountedKinematicParticle tau2Particle = pFactory.particle(transientTrack2, tauMass, chi2, ndf, tauSigma);

				float sig[6];
				sig[0] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,0);
				sig[1] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,1);
				sig[2] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,1);
				sig[3] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(0,2);
				sig[4] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(1,2);
				sig[5] = tau1Particle->stateAtPoint(pca1).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca1Cov(sig, sig+6);

				sig[0] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,0);
				sig[1] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,1);
				sig[2] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,1);
				sig[3] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(0,2);
				sig[4] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(1,2);
				sig[5] = tau2Particle->stateAtPoint(pca2).kinematicParametersError().matrix()(2,2);
				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > pca2Cov(sig, sig+6);

				ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > totCov = pca1Cov + pca2Cov;
				const float twoTauDist3DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist3D;

				distanceVector(2) = 0.0;
				const float twoTauDist2D = ROOT::Math::Mag(distanceVector);
				const float twoTauDist2DErr = TMath::Sqrt(ROOT::Math::Similarity(totCov, distanceVector)) / twoTauDist2D;

				dca.dca2d = twoTauDist2D;
				dca.dca2dErr = twoTauDist2DErr;
				dca.dca3d = twoTauDist3D;
				dca.dca3dErr = twoTauDist3DErr;
			}
		}
	}
	
	return dca;
}//NTupleProducerFromMiniAOD::DCA NTupleProducerFromMiniAOD::calculateDCA

unsigned int NTupleProducerFromMiniAOD::AddPFJets(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<pat::JetCollection> pfjets;
  iEvent.getByLabel(JetCollectionTag_, pfjets);
  
  //	edm::Handle<std::vector<reco::SecondaryVertexTagInfo> > svInfos;
  //	iEvent.getByLabel(edm::InputTag("secondaryVertexTagInfosEI"), svInfos);
  //	assert(svInfos.isValid());
  
  //	edm::Handle<edm::ValueMap<float> > puJetIdMVAFull;
  //	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullDiscriminant"), puJetIdMVAFull);
  
  //	edm::Handle<edm::ValueMap<int> > puJetIdFlagFull;
  //	iEvent.getByLabel(edm::InputTag("pileupJetIdProducer","fullId"), puJetIdFlagFull);

  if(pfjets.isValid())
    {
      // if (pfjets->size()>0) {
      // 	const std::vector< std::pair< std::string, float > > pairDiscriVector = (*pfjets)[0].getPairDiscri();
      // 	int nDiscri = pairDiscriVector.size();
      // 	std::cout << "Number of discriminators = " << nDiscri << std::endl;
      // 	for (int iD=0;iD<nDiscri;++iD) {
      // 	  std::pair<std::string, float> pairDiscri = pairDiscriVector[iD];
      // 	  std::cout << "Dicsr = " << pairDiscriVector[iD].first << std::endl;
      // 	}
      // }
      // std::cout << std::endl;
      
      
      for(unsigned i = 0 ; i < pfjets->size() ; i++)
	{
	  if((*pfjets)[i].pt() < cJetPtMin) continue;
	  if(fabs((*pfjets)[i].eta()) > cJetEtaMax) continue;
	  //	  std::cout << "Jet  " << i <<  ", pT=" <<  (*pfjets)[i].pt() << std::endl;
	  
	  pfjet_e[pfjet_count] = (*pfjets)[i].energy();
	  pfjet_px[pfjet_count] = (*pfjets)[i].px();
	  pfjet_py[pfjet_count] = (*pfjets)[i].py();
	  pfjet_pz[pfjet_count] = (*pfjets)[i].pz();
          pfjet_pt[pfjet_count] = (*pfjets)[i].pt();
          pfjet_eta[pfjet_count] = (*pfjets)[i].eta();
          pfjet_phi[pfjet_count] = (*pfjets)[i].phi();
	  pfjet_hadronicenergy[pfjet_count] = (*pfjets)[i].chargedHadronEnergy() + (*pfjets)[i].neutralHadronEnergy();
	  pfjet_chargedhadronicenergy[pfjet_count] = (*pfjets)[i].chargedHadronEnergy();
	  pfjet_emenergy[pfjet_count] = (*pfjets)[i].chargedEmEnergy() + (*pfjets)[i].neutralEmEnergy();
	  pfjet_chargedemenergy[pfjet_count] = (*pfjets)[i].chargedEmEnergy();
	  pfjet_chargedmulti[pfjet_count] = (*pfjets)[i].chargedMultiplicity();
	  pfjet_neutralmulti[pfjet_count] = (*pfjets)[i].neutralMultiplicity();
	  
	  pfjet_energycorr[pfjet_count] = -1.;
	  pfjet_energycorr_l1fastjet[pfjet_count] = -1.;
	  pfjet_energycorr_l2relative[pfjet_count] = -1.;
	  pfjet_energycorr_l3absolute[pfjet_count] = -1.;
	  pfjet_energycorr_l2l3residual[pfjet_count] = -1.;
		
	  if((*pfjets)[i].jecSetsAvailable())
	    {
	      pfjet_energycorr[pfjet_count] = (*pfjets)[i].jecFactor("Uncorrected");
	      pfjet_energycorr_l1fastjet[pfjet_count] = (*pfjets)[i].jecFactor("L1FastJet");
	      pfjet_energycorr_l2relative[pfjet_count] = (*pfjets)[i].jecFactor("L2Relative");
	      pfjet_energycorr_l3absolute[pfjet_count] = (*pfjets)[i].jecFactor("L3Absolute");
	      if (cdata) pfjet_energycorr_l2l3residual[pfjet_count] = (*pfjets)[i].jecFactor("L2L3Residual");
	    }
		    
	  // std::cout << "Jet Energy corrections : " << std::endl;
	  // std::cout << "    L1FastJet    = " << pfjet_energycorr_l1fastjet[pfjet_count] << std::endl;
	  // std::cout << "    L2Relative   = " << pfjet_energycorr_l2relative[pfjet_count] << std::endl;
	  // std::cout << "    L3Absolute   = " << pfjet_energycorr_l3absolute[pfjet_count] << std::endl;
	  // std::cout << "    L2L3Residual = " << pfjet_energycorr_l2l3residual[pfjet_count] << std::endl;
	  // std::cout << "    Total (Uncor)= " << pfjet_energycorr[pfjet_count] << std::endl;
	  
	  
	  pfjet_pu_jet_simple_loose[pfjet_count] = false;
	  pfjet_pu_jet_simple_medium[pfjet_count] = false;
	  pfjet_pu_jet_simple_tight[pfjet_count] = false;
	  pfjet_pu_jet_simple_mva[pfjet_count] = -1.0f;
	  
	  
	  //		    pfjet_pu_jet_full_loose[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kLoose);
	  //		    pfjet_pu_jet_full_medium[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kMedium);
	  //		    pfjet_pu_jet_full_tight[pfjet_count] = PileupJetIdentifier::passJetId( (*puJetIdFlagFull)[(*pfjets)[i].originalObjectRef()], PileupJetIdentifier::kTight);
	  //		    pfjet_pu_jet_full_mva[pfjet_count] = (*puJetIdMVAFull)[(*pfjets)[i].originalObjectRef()];
	  pfjet_flavour[pfjet_count] = (*pfjets)[i].partonFlavour();
		
	  for(unsigned n = 0 ; n < cBtagDiscriminators.size() ; n++)
	    {
	      pfjet_btag[pfjet_count][n] = -1000;
	      if(cBtagDiscriminators[n] != "F"){
		//		std::cout << " " << cBtagDiscriminators.at(n) << "  : " <<  (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) << std::endl;
		pfjet_btag[pfjet_count][n] = (*pfjets)[i].bDiscriminator(cBtagDiscriminators[n]) ;
	      }
	    }
	  pfjet_count++;
	  if(pfjet_count == M_jetmaxcount){
	    cerr << "number of pfjet > M_jetmaxcount. They are missing." << endl; 
	    errors |= 1<<4; 
	    break;
	  }
	}
    }
  
  return  pfjet_count;
}//unsigned int NTupleProducerFromMiniAOD::AddPFJets





