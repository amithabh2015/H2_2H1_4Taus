#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauMuonTrackSelector
#define H2to2H1to4Taus_TauTauSkimming_HTauTauMuonTrackSelector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauMuonTrackSelector {

 public:
  HTauTauMuonTrackSelector(double ipCut, double probCut, 
			   int nValidHits, int nValidPixelHits, 
			   double PtMin, double ptHardMin, 
			   double etaMax, reco::Vertex & vtx);
  ~HTauTauMuonTrackSelector();
  void selectMuons(patMuonVector & input, patMuonVector & output);
  void selectTracks(recoTrackVector & input, recoTrackVector & output);
  void setResonanceMassAndWindow(double mass, double massWindow);
  double getTrackIsolation(reco::TrackRef & track, recoTrackVector allTracks, double deltaRMax, double ptThreshold); 
  int muonTrkAssociatorByMass(pat::MuonRef & muon, recoTrackVector & trkVector, recoTrackVector & track);

 private:

  double _probCut;
  int _nValidHits;
  int _nValidPixelHits;
  double _ipCut;

  double _centralMass;
  double _massWindow;

  double _lowerMass;
  double _upperMass;

  double _muonMassPDG;

  double _ptMin;
  double _ptHardMin;
  double _etaMax;

  reco::Vertex _vtx;

};

#endif
