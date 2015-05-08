#ifndef DesyHTauTau_TauTauSkimming_HTauTauGetPFParticles
#define DesyHTauTau_TauTauSkimming_HTauTauGetPFParticles

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauGetPFParticles {

 public:
  HTauTauGetPFParticles(edm::InputTag tag);
  ~HTauTauGetPFParticles();
  void getPFParticles(const edm::Event& iEvent, patPFParticleVector & output);

 private:
  edm::InputTag _pfPartColName; 

};

#endif
