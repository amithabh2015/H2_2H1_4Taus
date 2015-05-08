#ifndef DesyHTauTau_TauTauSkimming_HTauTauElectronIsolator
#define DesyHTauTau_TauTauSkimming_HTauTauElectronIsolator

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"


class HTauTauElectronIsolator {

 public:
  HTauTauElectronIsolator(double,double,double,
			  double,double,double,
			  double,double);
  ~HTauTauElectronIsolator();
  void isolateElectronsAbs(bool mode, patElectronVector &, patElectronVector &);
  void isolateElectronsRel(bool mode, patElectronVector &, patElectronVector &);

 private:

  double _isoTrkCutEndCap;  
  double _isoECalCutEndCap;
  double _isoHCalCutEndCap;
  double _isoTrkCutBarrel;  
  double _isoECalCutBarrel;
  double _isoHCalCutBarrel;
  double _relIsoCut;

  double _cone;

};

#endif
