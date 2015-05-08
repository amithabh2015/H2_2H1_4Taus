#ifndef DesyHTauTau_TauTauSkimming_HTauTauElectronMuonSelector
#define DesyHTauTau_TauTauSkimming_HTauTauElectronMuonSelector

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"

class HTauTauElectronMuonSelector
{
 public:
  HTauTauElectronMuonSelector(double const ptMin, double const etaMax, bool const sameSign, patElectronVector& inputElectrons, patMuonVector& inputMuons);
  ~HTauTauElectronMuonSelector();
	
	bool OnePairPassesSelection();
	
	pat::ElectronRef GetSelectedElectron();
	pat::MuonRef     GetSelectedMuon();
	
 private:
  double _ptMin;
  double _etaMax;
  bool   _sameSign;
	
	bool _pairFound;
	
	pat::ElectronRef _selectedElectron;
	pat::MuonRef     _selectedMuon;
};


#endif
