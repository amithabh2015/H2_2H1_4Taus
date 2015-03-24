#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauMuonSelector
#define H2to2H1to4Taus_TauTauSkimming_HTauTauMuonSelector

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"
#include <string.h>
//#include <pair.h>
//#include <vector.h>

typedef std::pair<pat::MuonRef,pat::MuonRef> PairMuons;
typedef std::vector<PairMuons> PairMuonsVector;

class HTauTauMuonSelector {

 public:
  HTauTauMuonSelector(double ptMin, double etaMax, double ptMinDiMuon, bool sameSign);
  ~HTauTauMuonSelector();
  void selectMuons(patMuonVector & input, patMuonVector & output);
  void selectMuonsByKinematics(patMuonVector & input, patMuonVector & output);
  double selectPairMuons(patMuonVector & input, patMuonVector & output);
  void selectMuonsByResonanceMass(patMuonVector & input, patMuonVector & output);
  void setResonanceMassAndWindow(double mass, double massWindow);

 private:
  double _ptMin;
  double _etaMax;  
  double _ptMinDiMuon;
  bool _sameSign;

  double _centralMass;
  double _massWindow;

  double _lowerMass;
  double _upperMass;

  double _muonMassPDG;

};


#endif
