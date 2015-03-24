#ifndef DesyHTauTau_TauTauSkimming_HtauTauVertexSelector
#define DesyHTauTau_TauTauSkimming_HTauTauVertexSelector


#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class HTauTauVertexSelector {

 public:
  HTauTauVertexSelector(edm::InputTag tag, double zOffset, double dOffset, double probCut);
  ~HTauTauVertexSelector();
  bool primaryVertexFound(const edm::Event& ev, reco::Vertex & vertex, std::vector<reco::Vertex> & pVtxs);

 private:   

  edm::InputTag _vtxColName;
  double _zOffset;
  double _probCut;
  double _dOffset;

};


#endif
