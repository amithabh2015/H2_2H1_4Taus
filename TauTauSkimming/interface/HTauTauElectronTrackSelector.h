#ifndef DesyHTauTau_TauTauSkimming_HTauTauElectronTrackSelector
#define DesyHTauTau_TauTauSkimming_HTauTauElectronTrackSelector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauElectronTrackSelector {

 public:
  HTauTauElectronTrackSelector(double ipCut, double probCut, int nValidHits, int nValidPixelHits);
  ~HTauTauElectronTrackSelector();
  void selectElectrons(reco::Vertex & vtx, patElectronVector & input, patElectronVector & output);  
  
 private:

  double _probCut;
  int _nValidHits;
  int _nValidPixelHits;
  double _ipCut;

};

#endif
