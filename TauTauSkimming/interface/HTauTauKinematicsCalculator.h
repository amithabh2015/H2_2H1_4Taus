#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauKinematicsCalculator
#define H2to2H1to4Taus_TauTauSkimming_HTauTauKinematicsCalculator

#include <TLorentzVector.h>
#include <TVector3.h>
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

class HTauTauKinematicsCalculator {

 public:
  HTauTauKinematicsCalculator(TLorentzVector posLepton, 
			      TLorentzVector negLepton,
			      TLorentzVector met);
  ~HTauTauKinematicsCalculator();

  double diLeptonMass() { return _diLeptonMass; } ;
  double diLeptonDeltaPhi() {return _diLeptonDeltaPhi; } ;
  double diLeptonDeltaR() { return _diLeptonDeltaR; } ;
  double diLeptonDeltaEta() { return _diLeptonDeltaEta; } ;
  TLorentzVector & diLepton4P() { return _diLepton; } ;

  double posLepMETDeltaPhi() { return _posLepMETDeltaPhi; } ;
  double negLepMETDeltaPhi() { return _negLepMETDeltaPhi; } ;  
  double diLeptonMETDeltaPhi() { return _diLepMETDeltaPhi; };

  double mTdiLepton() { return _mTdiLepton; } ;
  double mTdiLeptonMET() { return _mTdiLeptonMET; } ;

  bool isValidDiTau() {return _isValidDiTau; } ;

  double diTauMass() { return _diTauMass; } ;
  TLorentzVector & diTau4P() { return _diTau; } ;
  
  double posTauELepETauRatio() { return _posTauELepETauRatio; };
  double negTauELepETauRatio() { return _negTauELepETauRatio; };

  double posLepETauRestFrame() { return _posLepETauRestFrame; };
  double negLepETauRestFrame() { return _negLepETauRestFrame; };

  double posLepEdiTauRestFrame() { return _posLepEdiTauRestFrame; };
  double negLepEdiTauRestFrame() { return _negLepEdiTauRestFrame; };

  TLorentzVector & posTau4P() { return _posTau; } ;
  TLorentzVector & negTau4P() { return _negTau; } ;

  double cosNegDiLeptonRestFrame() { return _cosNegDiLeptonRestFrame; };
  double cosNegTauRestFrame() { return _cosNegTauRestFrame; }

  double cosRestFrame(TLorentzVector boost, TLorentzVector vect);
  double energyRestFrame(TLorentzVector boost, TLorentzVector vect);

  double deltaR(double eta1, double phi1,
		double eta2, double phi2); 

 private:
  

  TLorentzVector _posLepton;
  TLorentzVector _negLepton;
  TLorentzVector _met;

  TLorentzVector _diLepton;
  TLorentzVector _diTau;

  TLorentzVector _posTau;
  TLorentzVector _negTau;

  bool _isValidDiTau;

  double _tauMass;

  double _diLeptonMass;
  double _diLeptonDeltaPhi;
  double _diLeptonDeltaEta;
  double _diLeptonDeltaR;

  double _posLepMETDeltaPhi;
  double _negLepMETDeltaPhi;
  double _diLepMETDeltaPhi;

  double _mTposLepMET;
  double _mTnegLepMET;
  double _mTdiLepton;
  double _mTdiLeptonMET;

  double _diTauMass;

  double _cosNegDiLeptonRestFrame;
  double _cosNegTauRestFrame;

  double _posTauELepETauRatio;
  double _negTauELepETauRatio;

  double _posLepETauRestFrame;
  double _negLepETauRestFrame;

  double _posLepEdiTauRestFrame;
  double _negLepEdiTauRestFrame;

};

#endif
