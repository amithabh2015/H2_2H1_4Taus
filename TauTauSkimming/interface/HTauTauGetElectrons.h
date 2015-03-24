#ifndef DesyHTauTau_TauTauSkimming_HTauTauGetElectrons
#define DesyHTauTau_TauTauSkimming_HTauTauGetElectrons

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauGetElectrons {

 public:
  HTauTauGetElectrons(edm::InputTag tag);
  ~HTauTauGetElectrons();
  void getElectrons(const edm::Event& iEvent, patElectronVector & output);
  double getIsolation(pat::ElectronRef & elec, puPFCandidateVector allPfs, double deltaRMax, double ptThreshold); 
  double getIsolation(pat::ElectronRef & elec, PFCandidateVector allPfs, double deltaRMax, double ptThreshold); 

 private:
  edm::InputTag _elecColName; 

};

#endif
