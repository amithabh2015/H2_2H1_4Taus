#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauMuonTrackSelector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <DataFormats/TrackReco/interface/TrackBase.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>
#include <DataFormats/TrackReco/interface/HitPattern.h>
#include <TMath.h>
#include <TLorentzVector.h>

HTauTauMuonTrackSelector::HTauTauMuonTrackSelector(double ipCut, double probCut, 
						   int nValidHits, int nValidPixelHits, 
						   double ptMin, double ptHardMin, double etaMax,
						   reco::Vertex & vtx) {

  _probCut = probCut;
  _ipCut = ipCut;
  _nValidHits = nValidHits;
  _nValidPixelHits = nValidPixelHits;

  _centralMass = 0;
  _massWindow = 0;

  _ptMin = ptMin;
  _ptHardMin = ptHardMin;

  _muonMassPDG = 0.105658;

  _etaMax = etaMax;

  _vtx = vtx;

}

HTauTauMuonTrackSelector::~HTauTauMuonTrackSelector() {}

void HTauTauMuonTrackSelector::setResonanceMassAndWindow(double mass, double window) {
  
  _centralMass = mass;
  _massWindow = window;

  _lowerMass = _centralMass - _massWindow;
  _upperMass = _centralMass - _massWindow;

}

void HTauTauMuonTrackSelector::selectMuons(patMuonVector & input, patMuonVector & output) {

  using namespace edm;
  using namespace pat;

  output.clear();

  int numberOfMuons = input.size();

  for (int iE=0;iE<numberOfMuons;++iE) {
    MuonRef muonRef = input[iE];
    reco::TrackRef tkRef = muonRef->innerTrack();
    double chi2 = tkRef->chi2();
    int ndf = int(tkRef->ndof());
    double prob = TMath::Prob(chi2,ndf);
    if (prob>_probCut) {
      reco::TrackBase::Point point = _vtx.position();
      double ip = tkRef->dxy( point );
      if (fabs(ip)<_ipCut) {
	int nValidHits = tkRef->numberOfValidHits();
	const reco::HitPattern hitPattern = tkRef->hitPattern();
	int nValidPixelHits = hitPattern.numberOfValidPixelHits();
	if (nValidHits >= _nValidHits && nValidPixelHits >= _nValidPixelHits) {
	  output.push_back(muonRef);
	}
      }
    }
  }


}  

void HTauTauMuonTrackSelector::selectTracks(recoTrackVector & input, recoTrackVector & output) {

  using namespace edm;
  using namespace reco;

  output.clear();

  int numberOfTrks = input.size();

  for (int iE=0;iE<numberOfTrks;++iE) {
    reco::TrackRef tkRef = input[iE];
    double chi2 = tkRef->chi2();
    int ndf = int(tkRef->ndof());
    double prob = TMath::Prob(chi2,ndf);
    if (prob>_probCut) {
      reco::TrackBase::Point point = _vtx.position();
      double ip = tkRef->dxy( point );
      if (fabs(ip)<_ipCut) {
	int nValidHits = tkRef->numberOfValidHits();
	const reco::HitPattern hitPattern = tkRef->hitPattern();
	int nValidPixelHits = hitPattern.numberOfValidPixelHits();
	float eta = tkRef->eta();
	if (abs(eta)<_etaMax) {
	  if (nValidHits >= _nValidHits && nValidPixelHits >= _nValidPixelHits) 
	    output.push_back(tkRef);
	}
      }
    }
  }


}  


int HTauTauMuonTrackSelector::muonTrkAssociatorByMass(pat::MuonRef & muon, 
						      recoTrackVector & tracks,
						      recoTrackVector & selectedTracks) {
  // this member function returns number of tracks ;

  using namespace reco;
  using namespace pat;

  selectedTracks.clear();
  int numberOfTracks = int(tracks.size());
  float muonCharge   = muon->charge();

  TLorentzVector muonTL(0.,0.,0.,0.);
  muonTL.SetXYZM(muon->px(),muon->py(),muon->pz(),_muonMassPDG);

//  std::cout << " selector # tracks = " << numberOfTracks << std::endl;	

  for (int iT=0;iT<numberOfTracks;++iT) {
    reco::TrackRef trk = tracks[iT];
    float trkCharge = trk->charge();
//    std::cout << "   selector : " << iT << "   track (px,py,pz)=(" << trk->px() << "," << trk->py() << "," << trk->pz() << ")     q = " << trkCharge << std::endl;  
    if (trkCharge*muonCharge>0.) continue;
    TLorentzVector trkTL(0.,0.,0.,0.);
    trkTL.SetXYZM(trk->px(),trk->py(),trk->pz(),_muonMassPDG);
    TLorentzVector resonance = muonTL + trkTL;
    double resonanceMass = resonance.M();
//    std::cout << "                selector mass = " << resonanceMass << std::endl;
    double massDiff = double(TMath::Abs(float(resonanceMass)-float(_centralMass)));
    if (massDiff<_massWindow)
      selectedTracks.push_back(trk);
  }

  int numberOfSelectedTracks = int(selectedTracks.size());
  return numberOfSelectedTracks;

}

