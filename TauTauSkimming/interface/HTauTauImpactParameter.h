#ifndef DesyHTauTau_TauTauSkimming_HTauTauImpactParameter
#define DesyHTauTau_TauTauSkimming_HTauTauImpactParameter 

#include <DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h>
#include <TrackingTools/TransientTrack/interface/TransientTrack.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

struct IPMeasurement {
  Measurement1D ip;
  Measurement1D ipZ;
  Measurement1D ip3D;
}; 

class HTauTauImpactParameter {

 public:
  HTauTauImpactParameter();
  ~HTauTauImpactParameter();

  IPMeasurement impactParameter(const reco::TransientTrack& transientTrack, const reco::Vertex & primaryVertex);
  IPMeasurement impactParameter(const reco::TransientTrack& transientTrack, const GlobalVector& direction, const reco::Vertex & primaryVertex);
    
};


#endif
