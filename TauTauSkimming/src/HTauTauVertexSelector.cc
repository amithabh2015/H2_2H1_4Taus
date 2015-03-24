#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauVertexSelector.h"
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h" 
#include <math.h>
#include <TMath.h>

HTauTauVertexSelector::HTauTauVertexSelector(edm::InputTag vtxColName, double zOffset, double dOffset, double probCut) {
  _vtxColName = vtxColName;
  _zOffset = zOffset;
  _dOffset = dOffset;
  _probCut = probCut;
}

HTauTauVertexSelector::~HTauTauVertexSelector() {}

bool HTauTauVertexSelector::primaryVertexFound(const edm::Event& iEvent, reco::Vertex & primaryVertex, std::vector<reco::Vertex> & pVtxs) {

  bool primaryVertexFound = false;
  pVtxs.clear();

  using namespace edm;
  using namespace reco;

  Handle<VertexCollection> vertexCol;
  try{
    iEvent.getByLabel(_vtxColName,vertexCol);
  }catch(const edm::Exception &e) {;}

  if(vertexCol.isValid()){
    int numberOfVertices = vertexCol->size();
    if(numberOfVertices > 0){
      double ptmax = 0;
      for(int iVertex = 0; iVertex<numberOfVertices; iVertex++){
	Vertex vtx = vertexCol->at(iVertex);
	pVtxs.push_back(vtx);
	//cout << "vertex x,y,z " << iVertex.x() << " "
	//                        << iVertex.y() << " "
	//                        << iVertex.z() << endl;
	double chi2 =  vtx.chi2();
	int ndf = int(vtx.ndof());
	double prob = TMath::Prob(chi2,ndf);
	double xVtx = vtx.x();
	double yVtx = vtx.y();
	double d = sqrt(xVtx*xVtx+yVtx*yVtx);
	if (fabs(vtx.z())<_zOffset && d < _dOffset && prob > _probCut) {
	  double ptsum = 0;
	  Vertex::trackRef_iterator iTrack;
	  for(iTrack  = vtx.tracks_begin();
	      iTrack != vtx.tracks_end();++iTrack){
	    ptsum += (*iTrack)->pt();
	  }	
	  if(ptsum > ptmax){
	    ptmax = ptsum;
	    primaryVertex = vtx;
	    primaryVertexFound = true;
	  }
	}
      }
    }
  }
  return primaryVertexFound;
}
