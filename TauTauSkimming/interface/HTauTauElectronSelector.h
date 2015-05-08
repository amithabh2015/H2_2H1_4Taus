#ifndef DesyHTauTau_TauTauSkimming_HTauTauElectronSelector
#define DesyHTauTau_TauTauSkimming_HTauTauElectronSelector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"
#include <string.h>
//#include <pair.h>
//#include <vector.h>

typedef std::pair<pat::ElectronRef,pat::ElectronRef> PairElectrons;
typedef std::vector<PairElectrons> PairElectronsVector;

class HTauTauElectronSelector {

 public:
  HTauTauElectronSelector(double ptMin, double etaMax, double ptMinDiElec, bool sameSign);
  ~HTauTauElectronSelector();
  void selectElectrons(const std::string label, double IdThreshold, patElectronVector & input, patElectronVector & output);
  void selectElectronsByKinematics(patElectronVector & input, patElectronVector & output);
  void selectElectronsById(bool mode, const std::string label, double IdThreshold, patElectronVector & input, patElectronVector & output);
  double selectPairElectrons(patElectronVector & input, patElectronVector & output);


 private:
  double _ptMin;
  double _etaMax;  
  double _ptMinDiElec;
  bool _sameSign;


};


#endif
