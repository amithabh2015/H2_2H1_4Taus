#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauGetTracks
#define H2to2H1to4Taus_TauTauSkimming_HTauTauGetTracks

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauGetTracks {

 public:
  HTauTauGetTracks(edm::InputTag tag);
  ~HTauTauGetTracks();
  void getTracks(const edm::Event& iEvent, recoTrackVector & output);
  double getIsolation(reco::TrackRef & track, recoTrackVector allTracks,  double deltaRMax, double ptThreshold);
  double getIsolation(reco::TrackRef & track, patPFParticleVector allPfs, double deltaRMax, double ptThreshold);
 

 private:
  edm::InputTag _trkColName; 

};

#endif
