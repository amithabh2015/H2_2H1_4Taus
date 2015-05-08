#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetTracks.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"

HTauTauGetTracks::HTauTauGetTracks(edm::InputTag tag) 
{
  _trkColName = tag;
}
HTauTauGetTracks::~HTauTauGetTracks() {}
void HTauTauGetTracks::getTracks(const edm::Event& iEvent, recoTrackVector & output) {
  
  using namespace edm;
  using namespace reco;

  output.clear();

  Handle<TrackCollection> trkHandle;
  iEvent.getByLabel(_trkColName,trkHandle);
  
  for (unsigned int j=0;j<trkHandle->size();j++) {
    TrackRef thisTrackRef(trkHandle,j);
    output.push_back(thisTrackRef);
  }

}

double HTauTauGetTracks::getIsolation(reco::TrackRef & track, recoTrackVector allTracks, double deltaRMax, double ptThreshold) {

  double iso = 0;
  int nOfTracks = int(allTracks.size());

  for (int iT=0; iT<nOfTracks; ++iT) {
    reco::TrackRef currentTrk = allTracks[iT];
    if ( currentTrk->pt()<ptThreshold ) continue;
    double dR = deltaR(*track,*currentTrk);
    if (dR>deltaRMax) continue;
    iso += currentTrk->pt();
  }

  iso -= track->pt();

  return iso;

}

double HTauTauGetTracks::getIsolation(reco::TrackRef & track, patPFParticleVector allPfs, double deltaRMax, double ptThreshold) {

  double iso = 0;
  int nOfPfs = int(allPfs.size());

  TLorentzVector trackTL;
  trackTL.SetXYZM(track->px(),track->py(),track->pz(),0.0);

  for (int iT=0; iT<nOfPfs; ++iT) {
    pat::PFParticleRef currentPf = allPfs[iT];
    if ( currentPf->pt()<ptThreshold ) continue;
    TLorentzVector pfTL(currentPf->px(),currentPf->py(),currentPf->pz(),currentPf->energy());
    TLorentzVector diff = trackTL - pfTL;
    double dP = diff.P();
    if (dP<1e-3) continue;
    double dR = deltaR(*track,*currentPf);
    if (dR>deltaRMax) continue;
    iso += currentPf->pt();
  }

  return iso;

}

