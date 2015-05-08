#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonIsolator.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h" 
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"

#include <math.h>



HTauTauMuonIsolator::HTauTauMuonIsolator(double relIsoCut) {

  _relIsoCut = relIsoCut;

}

HTauTauMuonIsolator::~HTauTauMuonIsolator() {}

void HTauTauMuonIsolator::isolateMuonsRel(bool mode, patMuonVector & input, patMuonVector & output) {
  // mode = true - direct isolator
  // mode = false - inverse isolator

  output.clear(); 

  int nMuons = input.size();

  for (int iE=0;iE<nMuons;++iE) {

    pat::MuonRef muonRef = input[iE];

    double muonTrackIso = muonRef->ecalIso ();
    double muonEcalIso = muonRef->hcalIso();
    double muonHcalIso = muonRef->trackIso();

    bool muonIso = ((muonTrackIso+muonEcalIso+muonHcalIso)/(fmax(0.1,muonRef->pt()))) < _relIsoCut;


    if (mode) {
      if (muonIso)
	output.push_back(muonRef);
    }
    else {
      if (!muonIso)
	output.push_back(muonRef);
    }

  }

}
