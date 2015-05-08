#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonSelector.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" 
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <DataFormats/CaloRecHit/interface/CaloCluster.h>
#include <math.h>
#include <TLorentzVector.h>

HTauTauMuonSelector::HTauTauMuonSelector(double ptMin, double etaMax, double ptMinDiMuon, bool sameSign) {

  _ptMin = ptMin;
  _etaMax = etaMax;
  _ptMinDiMuon = ptMinDiMuon;
  _sameSign = sameSign;
  _centralMass = 0;
  _massWindow = 0;
  _lowerMass = 0;
  _upperMass = 0;
  _muonMassPDG = 0.105658;

}

HTauTauMuonSelector::~HTauTauMuonSelector() {}


void HTauTauMuonSelector::setResonanceMassAndWindow(double mass, double massWindow) {
  
  _centralMass = mass;
  _massWindow = massWindow;

  _lowerMass = _centralMass - _massWindow;
  _upperMass = _centralMass - _massWindow;
  
}


void HTauTauMuonSelector::selectMuons(patMuonVector & input, patMuonVector & output) {

  using namespace edm;
  using namespace pat;

  output.clear();

  int numberOfMuons = input.size();
  
  for (int iE=0;iE<numberOfMuons;++iE) {
    MuonRef muonRef = input[iE];

    if ( (muonRef->pt()>_ptMin) && (fabs(muonRef->eta()) < _etaMax) ) 
      output.push_back( muonRef );
  }

}


void HTauTauMuonSelector::selectMuonsByKinematics(patMuonVector & input, patMuonVector & output) {

  using namespace edm;
  using namespace pat;

  output.clear();
  int numberOfMuons = input.size();
  for (int iE=0;iE<numberOfMuons;++iE) {
    MuonRef muonRef = input[iE];
    if (muonRef->pt() > _ptMin && fabs(muonRef->eta()) < _etaMax)
       output.push_back( muonRef );
  }  

}

double HTauTauMuonSelector::selectPairMuons(patMuonVector & input, 
					    patMuonVector & output) {

  // returns max pt of leptons
  // first muon "-"
  // second muons "+"


  using namespace pat;
  using namespace edm;

  output.clear();

  int nMuont = input.size();

  double ptSumMax = 0.;
  double ptMax = 0.;

  PairMuons pairMuons;
  
  bool pairMuonsFound = false;

  for (int iMuont=0;iMuont<nMuont-1;++iMuont) {
    MuonRef firstMuont = input[iMuont];
    for (int jMuont=iMuont+1;iMuont<nMuont;++iMuont) {
      MuonRef secondMuont = input[jMuont];

      bool pairCharge = true;
      if (_sameSign) {
	pairCharge = firstMuont->charge()*secondMuont->charge()>0.;
      }
      else {
	pairCharge = firstMuont->charge()*secondMuont->charge()<0.;
      }

      if (pairCharge) {
	ptMax = fmax(firstMuont->pt(),secondMuont->pt()); 
	double ptSum = firstMuont->pt()+secondMuont->pt();
	pairMuonsFound = true;
	if (ptSum>ptSumMax) {
	  ptSumMax = ptSum;
	  ptMax = fmax(firstMuont->pt(),secondMuont->pt()); 
	  if (firstMuont->charge()<0) {
	    pairMuons.first = firstMuont;
	    pairMuons.second = secondMuont;
	  }
	  else {
	    pairMuons.first = secondMuont;
	    pairMuons.second = firstMuont;
	  }
	}
      }
    }
  }

  if (pairMuonsFound&&ptMax>_ptMinDiMuon) {
    output.push_back(pairMuons.first);
    output.push_back(pairMuons.second);
  }

  


  return ptMax;

}

void HTauTauMuonSelector::selectMuonsByResonanceMass(patMuonVector & input, patMuonVector & output) {

  using namespace pat;
  using namespace edm;

  output.clear();

  int nMuont = input.size();
  
  for (int iMuont=0;iMuont<nMuont;++iMuont) {
    MuonRef firstMuont = input[iMuont];
    float chargeFirst = firstMuont->charge();
    bool selectMuon = false;
    TLorentzVector firstMuonTL(0.,0.,0.,0.);
    firstMuonTL.SetXYZM(firstMuont->px(),firstMuont->py(),firstMuont->pz(),_muonMassPDG);
    for (int jMuont=0;iMuont<nMuont;++iMuont) {
      if (iMuont==jMuont) continue;
      MuonRef secondMuont = input[jMuont];
      float chargeSecond = secondMuont->charge();
      if (chargeFirst*chargeSecond>0.) continue;
      TLorentzVector secondMuonTL(0.,0.,0.,0.);
      secondMuonTL.SetXYZM(secondMuont->px(),secondMuont->py(),secondMuont->pz(),_muonMassPDG);
      TLorentzVector resonance = firstMuonTL + secondMuonTL;
      double resonanceMass = resonance.M();
      double massDiff = double(TMath::Abs(resonanceMass-_centralMass));
      if (massDiff<_massWindow) {
	selectMuon = true;
	break;
      }

    }
    if (selectMuon)
      output.push_back(firstMuont);
  }

}
