#ifndef DesyHTauTau_TauTauSkimming_HTauTauVertexInfo
#define DesyHTauTau_TauTauSkimming_HTauTauVertexInfo

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TMath.h"
#include "TString.h"
#include "TTree.h"

#include <vector>

class HTauTauVertexInfo
{
public:
	HTauTauVertexInfo(const edm::ParameterSet&);
	~HTauTauVertexInfo();

	void beginJob(TTree* tree, edm::Service<TFileService>& fs);
	void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
	void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void endJob();
	
	reco::Vertex getPrimaryVtx() const { return primaryVtx; };
	bool primaryVertexFound() const { return primaryVtxFound; };
	std::vector<reco::Vertex> getVertices() const { return vertices; };
	int getNumberOfVertices() const { return nPV; };

private:
	TTree* _tree;
	TH1F * h_vtxZ_;
	TH1F * h_vtxD_;
	TH1F * h_vtxProb_;

	int _verbosity;

	edm::InputTag _vtxSrc;
	double _zVtxCut;
	double _dVtxCut;
	double _probVtxCut;

	int nPV;
	float probPV[100];
	float ndofPV[100];
	int   nTrkPV[100];
	float chi2PV[100];
	float xPV[100];
	float yPV[100];
	float zPV[100];
	float sumPtPV[100];

	reco::Vertex primaryVtx;
	std::vector<reco::Vertex> vertices;
	bool primaryVtxFound;
};

#endif 
