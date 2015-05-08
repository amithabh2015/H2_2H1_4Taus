#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauAssociateDecayMode
#define H2to2H1to4Taus_TauTauSkimming_HTauTauAssociateDecayMode 

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauAssociateDecayMode {

 public:
  HTauTauAssociateDecayMode(double);
  ~HTauTauAssociateDecayMode();
  void selectTausByDecay(int mode, 
			 GenParticleVector & genParticles, 
			 PFTauVector & pfTaus,
			 pfTauIndexGenParticleRefVector & output);

 private:
  
  void selectGeneratedTaus(GenParticleVector & input,
			   GenParticleVector & output);

  bool isTau2PiNu(reco::GenParticleRef & tau);
  bool isTau2RhoNu(reco::GenParticleRef & tau);
  bool isTau2Ele(reco::GenParticleRef & tau);
  bool isTau2Mu(reco::GenParticleRef & tau);
  bool isTauDecay(int, reco::GenParticleRef & tau);

  double _deltaRCut;

};

#endif
