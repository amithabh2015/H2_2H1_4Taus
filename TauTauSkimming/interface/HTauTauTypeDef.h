#ifndef DesyHTauTau_TauTauSkimming_HTauTauTypeDef
#define DesyHTauTau_TauTauSkimming_HTauTauTypeDef

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h>

#include <vector>

typedef std::vector<edm::InputTag> InputTagVector;
typedef std::vector<reco::GsfElectronRef> GsfElectronVector; 
typedef std::vector<reco::PFTauRef> PFTauVector;
typedef std::vector<reco::GenParticleRef> GenParticleVector;
typedef std::pair<reco::PFTauRef,reco::GenParticleRef> pfTauGenParticleRef;
typedef std::pair<int,reco::GenParticleRef> pfTauIndexGenParticleRef; 
typedef std::vector<pfTauGenParticleRef> pfTauGenParticleRefVector;
typedef std::vector<pfTauIndexGenParticleRef> pfTauIndexGenParticleRefVector; 

typedef std::vector<pat::TauRef> patTauVector;
typedef std::vector<pat::ElectronRef> patElectronVector;
typedef std::vector<pat::MuonRef> patMuonVector;
typedef std::vector<reco::TrackRef> recoTrackVector;
typedef std::vector<pat::PFParticleRef> patPFParticleVector;
typedef std::vector<reco::PileUpPFCandidateRef> puPFCandidateVector;
typedef std::vector<reco::PFCandidateRef> PFCandidateVector;

//typedef std::map<std::string,double> MapOfElecId;
//typedef std::vector<MapOfElecId> ElecIdVector;
//typedef std::map<pat::ElectronRef,ElecIdVector> patElectronIdMap;
//typedef std::vector<patElectronIdMap> patElectronIdMapVector;

#endif
