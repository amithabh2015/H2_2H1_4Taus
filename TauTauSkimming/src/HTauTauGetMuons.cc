#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetMuons.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"

HTauTauGetMuons::HTauTauGetMuons(edm::InputTag tag) 
{
  _muonColName = tag;
}
HTauTauGetMuons::~HTauTauGetMuons() {}
void HTauTauGetMuons::getMuons(const edm::Event& iEvent, patMuonVector & output) {
  
  using namespace edm;
  using namespace pat;

  output.clear();

  Handle<MuonCollection> muonHandle;
  iEvent.getByLabel(_muonColName,muonHandle);
  
  for (unsigned int j=0;j<muonHandle->size();j++) {
    MuonRef thisMuonRef(muonHandle,j);
    output.push_back(thisMuonRef);
  }

}

double HTauTauGetMuons::getIsolation(pat::MuonRef & muon, patPFParticleVector allPfs, double deltaRMax, double ptThreshold) {

  double iso = 0;
  int nOfPfs = int(allPfs.size());

  for (int iT=0; iT<nOfPfs; ++iT) {
    pat::PFParticleRef currentPf = allPfs[iT];
    if ( currentPf->pt()<ptThreshold ) continue;
    double dR = deltaR(*muon,*currentPf);
    if (dR>deltaRMax) continue;
    iso += currentPf->pt();
  }

  return iso;

}

double HTauTauGetMuons::getIsolation(pat::MuonRef & muon, puPFCandidateVector allPfs, double deltaRMax, double ptThreshold) {

  double iso = 0;
  int nOfPfs = int(allPfs.size());

  TLorentzVector muonTL(muon->px(),muon->py(),muon->pz(),muon->energy());

  for (int iT=0; iT<nOfPfs; ++iT) {
    reco::PileUpPFCandidateRef currentPf = allPfs[iT];
    TLorentzVector pfTL(currentPf->px(),currentPf->py(),currentPf->pz(),currentPf->energy());
    TLorentzVector diffTL = muonTL - pfTL;
    double deltaP = diffTL.P();
    if (deltaP<0.01) continue;
    if ( currentPf->pt()<ptThreshold ) continue;
    double dR = deltaR(*muon,*currentPf);
    if (dR>deltaRMax) continue;
    iso += currentPf->pt();
  }

  return iso;

}

double HTauTauGetMuons::getIsolation(pat::MuonRef & muon, PFCandidateVector allPfs, double deltaRMax, double ptThreshold) {

  double iso = 0;
  int nOfPfs = int(allPfs.size());

  TLorentzVector muonTL(muon->px(),muon->py(),muon->pz(),muon->energy());

  for (int iT=0; iT<nOfPfs; ++iT) {
    reco::PFCandidateRef currentPf = allPfs[iT];
    TLorentzVector pfTL(currentPf->px(),currentPf->py(),currentPf->pz(),currentPf->energy());
    TLorentzVector diffTL = muonTL - pfTL;
    double deltaP = diffTL.P();
    if (deltaP<0.01) continue;
    if ( currentPf->pt()<ptThreshold ) continue;
    double dR = deltaR(*muon,*currentPf);
    if (dR>deltaRMax) continue;
    iso += currentPf->pt();
  }

  return iso;

}

