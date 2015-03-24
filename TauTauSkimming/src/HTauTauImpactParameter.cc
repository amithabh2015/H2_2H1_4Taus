#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauImpactParameter.h"
#include <RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h>
#include <RecoBTag/BTagTools/interface/SignedImpactParameter3D.h>
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


HTauTauImpactParameter::HTauTauImpactParameter() {}
HTauTauImpactParameter::~HTauTauImpactParameter() {}

IPMeasurement HTauTauImpactParameter::impactParameter(const reco::TransientTrack& transientTrack, const reco::Vertex & primaryVertex) {

  reco::Track track = transientTrack.track();
  GlobalVector direction(track.px(),track.py(),track.pz());
  return impactParameter(transientTrack,direction,primaryVertex);
}

IPMeasurement HTauTauImpactParameter::impactParameter(const reco::TransientTrack& transientTrack, const GlobalVector& direction, const reco::Vertex & primaryVertex) {

  SignedTransverseImpactParameter stip;
  Measurement1D ip  = stip.apply(transientTrack,direction,primaryVertex).second;
  Measurement1D ipZ = stip.zImpactParameter(transientTrack,direction,primaryVertex).second;

  SignedImpactParameter3D signed_ip3D;
  Measurement1D ip3D = signed_ip3D.apply(transientTrack,direction,primaryVertex).second;
  
  IPMeasurement ipMeas;
  
  ipMeas.ip = ip;
  ipMeas.ipZ = ipZ;
  ipMeas.ip3D = ip3D;

  return ipMeas;

}

