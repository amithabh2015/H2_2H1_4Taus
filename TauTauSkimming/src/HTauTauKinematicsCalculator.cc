#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauKinematicsCalculator.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TMath.h>

HTauTauKinematicsCalculator::HTauTauKinematicsCalculator(TLorentzVector posLepton, 
							 TLorentzVector negLepton,
							 TLorentzVector met) {
  //  Calculation of all the variables --->

  _tauMass = 1.77684;

  _posLepton = posLepton;
  _negLepton = negLepton;
  _met = met;
  _diLepton = _posLepton + _negLepton;

  _diLeptonMass = _diLepton.M();

  double posLepMod = TMath::Sqrt(_posLepton.Px()*_posLepton.Px()+
				 _posLepton.Py()*_posLepton.Py());

  double negLepMod = TMath::Sqrt(_negLepton.Px()*_negLepton.Px()+
				 _negLepton.Py()*_negLepton.Py());

  double diLepMod = TMath::Sqrt(_diLepton.Px()*_diLepton.Px()+
				_diLepton.Py()*_diLepton.Py());

  double product = 
    _negLepton.Px()*_posLepton.Px()+
    _negLepton.Py()*_posLepton.Py();

  double cosinus = product/(posLepMod*negLepMod);

  _diLeptonDeltaPhi = TMath::ACos(cosinus);

	// protect against the case of one particle having zero pT (where ROOT throws an exception)
	if (_posLepton.Pt() != 0.0 && _negLepton.Pt() != 0.0)
		_diLeptonDeltaEta = _posLepton.Eta() - _negLepton.Eta();
	else
	{
		_diLeptonDeltaEta = 0.0;
	}

  _diLeptonDeltaR = TMath::Sqrt(_diLeptonDeltaPhi*_diLeptonDeltaPhi+
				_diLeptonDeltaEta*_diLeptonDeltaEta);


  double prodPosLepMET = _posLepton.Px()*_met.Px()+_posLepton.Py()*_met.Py();
  double prodNegLepMET = _negLepton.Px()*_met.Px()+_negLepton.Py()*_met.Py();
  double prodDiLepMET = _diLepton.Px()*_met.Px()+_diLepton.Py()*_met.Py();
  double METMod = TMath::Sqrt(_met.Px()*_met.Px()+_met.Py()*_met.Py());

  _posLepMETDeltaPhi = TMath::ACos(prodPosLepMET/(posLepMod*METMod));
  _negLepMETDeltaPhi = TMath::ACos(prodNegLepMET/(negLepMod*METMod));
  _diLepMETDeltaPhi = TMath::ACos(prodDiLepMET/(diLepMod*METMod));

  double pxDiLepMET = _diLepton.Px() + _met.Px();
  double pyDiLepMET = _diLepton.Py() + _met.Py();

  _mTdiLepton = TMath::Sqrt(_diLepton.Px()*_diLepton.Px()+_diLepton.Py()*_diLepton.Py());
  _mTdiLeptonMET = TMath::Sqrt(pxDiLepMET*pxDiLepMET+pyDiLepMET*pyDiLepMET);

  double cosPos = TMath::Cos(_posLepton.Phi());
  double cosNeg = TMath::Cos(_negLepton.Phi());

  double sinPos = TMath::Sin(_posLepton.Phi());
  double sinNeg = TMath::Sin(_negLepton.Phi());

  double det = cosPos*sinNeg-cosNeg*sinPos;

  double detP = _met.Px()*sinNeg - cosNeg*_met.Py();
  double detN = cosPos*_met.Py() - _met.Px()*sinPos;

  double EPosNeu = detP/det;
  double ENegNeu = detN/det;

  _cosNegDiLeptonRestFrame = cosRestFrame(_diLepton,_negLepton);
    

  _isValidDiTau = true;

  if (EPosNeu<0.0)
    _isValidDiTau = false;

  if (ENegNeu<0.0)
    _isValidDiTau = false;

  if (_isValidDiTau) {
    double PxPosTau = _posLepton.Px() + EPosNeu*cosPos;
    double PyPosTau = _posLepton.Py() + EPosNeu*sinPos;
    _posTauELepETauRatio = _posLepton.Px()/PxPosTau;
    double PzPosTau = _posLepton.Pz()/_posTauELepETauRatio;
    double PPosTau = TMath::Sqrt(PxPosTau*PxPosTau+PyPosTau*PyPosTau+PzPosTau*PzPosTau);
    double EPosTau = TMath::Sqrt(PPosTau*PPosTau+_tauMass*_tauMass);

    _posTau.SetPx(PxPosTau);
    _posTau.SetPy(PyPosTau);
    _posTau.SetPz(PzPosTau);
    _posTau.SetE(EPosTau);

    double PxNegTau = _negLepton.Px() + ENegNeu*cosNeg;
    double PyNegTau = _negLepton.Py() + ENegNeu*sinNeg;
    _negTauELepETauRatio = _negLepton.Px()/PxNegTau;
    double PzNegTau = _negLepton.Pz()/_negTauELepETauRatio;
    double PNegTau = TMath::Sqrt(PxNegTau*PxNegTau+PyNegTau*PyNegTau+PzNegTau*PzNegTau);
    double ENegTau = TMath::Sqrt(PNegTau*PNegTau+_tauMass*_tauMass);

    _negTau.SetPx(PxNegTau);
    _negTau.SetPy(PyNegTau);
    _negTau.SetPz(PzNegTau);
    _negTau.SetE(ENegTau);

    _diTau = _posTau + _negTau;

    _diTauMass = _diTau.M();

    _cosNegTauRestFrame = cosRestFrame(_diTau,_negTau);

    _posLepETauRestFrame = energyRestFrame(_posTau,_posLepton);
    _negLepETauRestFrame = energyRestFrame(_negTau,_negLepton);

    _posLepEdiTauRestFrame = energyRestFrame(_posTau,_diLepton);
    _negLepEdiTauRestFrame = energyRestFrame(_negTau,_diLepton);

  }


}

HTauTauKinematicsCalculator::~HTauTauKinematicsCalculator() {}


double HTauTauKinematicsCalculator::energyRestFrame(TLorentzVector boost, TLorentzVector dummy) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  TLorentzVector vect = dummy;

  vect.Boost(bx,by,bz);
  
  return vect.E();

}

double HTauTauKinematicsCalculator::cosRestFrame(TLorentzVector boost, TLorentzVector dummy) {

  double bx = -boost.Px()/boost.E();
  double by = -boost.Py()/boost.E();
  double bz = -boost.Pz()/boost.E();

  TLorentzVector vect = dummy;

  vect.Boost(bx,by,bz);
  double prod = -vect.Px()*bx-vect.Py()*by-vect.Pz()*bz;
  double modBeta = TMath::Sqrt(bx*bx+by*by+bz*bz); 
  double modVect = TMath::Sqrt(vect.Px()*vect.Px()+vect.Py()*vect.Py()+vect.Pz()*vect.Pz());
  
  double cosinus = prod/(modBeta*modVect);

  return cosinus;

}

double HTauTauKinematicsCalculator::deltaR(double eta1, 
					   double phi1,
					   double eta2,
					   double phi2) {

  double p1x = TMath::Cos(phi1);
  double p2x = TMath::Cos(phi2);
  
  double p1y = TMath::Sin(phi1);
  double p2y = TMath::Sin(phi2);

  double cosDeltaPhi = p1x*p2x + p1y*p2y;
  double deltaPhi = TMath::ACos(cosDeltaPhi);

  double deltaEta = eta1 - eta2;

  double dR = TMath::Sqrt(deltaPhi*deltaPhi+deltaEta*deltaEta);

  return dR;

}

