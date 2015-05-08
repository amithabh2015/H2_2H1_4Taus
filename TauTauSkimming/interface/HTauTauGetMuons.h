#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauGetMuons
#define H2to2H1to4Taus_TauTauSkimming_HTauTauGetMuons

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauGetMuons {

 public:
  HTauTauGetMuons(edm::InputTag tag);
  ~HTauTauGetMuons();
  void getMuons(const edm::Event& iEvent, patMuonVector & output);
  double getIsolation(pat::MuonRef & muon, patPFParticleVector allPfs, double deltaRMax, double ptThreshold);
  double getIsolation(pat::MuonRef & muon, puPFCandidateVector allPFs, double deltaRMax, double ptThreshold);
  double getIsolation(pat::MuonRef & muon, PFCandidateVector allPFs, double deltaRMax, double ptThreshold);

 private:
  edm::InputTag _muonColName; 

};

#endif
