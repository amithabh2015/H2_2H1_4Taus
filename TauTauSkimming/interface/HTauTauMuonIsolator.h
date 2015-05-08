#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauMuonIsolator
#define H2to2H1to4Taus_TauTauSkimming_HTauTauMuonIsolator

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"


class HTauTauMuonIsolator {

 public:
  HTauTauMuonIsolator(double);
  ~HTauTauMuonIsolator();
  void isolateMuonsRel(bool mode, patMuonVector &, patMuonVector &);

 private:

  double _relIsoCut;
};

#endif
