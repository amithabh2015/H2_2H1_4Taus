#include "../interface/HTauTauVertexInfo.h"

HTauTauVertexInfo::HTauTauVertexInfo(const edm::ParameterSet& iConfig)
{
	_verbosity = iConfig.getParameter<int>("PrintOut");

	_vtxSrc         = iConfig.getParameter<edm::InputTag>  (  "PrimaryVertexSource" );
	_zVtxCut         = iConfig.getParameter<double>         (  "ZVertexCut" );
	_dVtxCut         = iConfig.getParameter<double>         (  "DVertexCut" );
	_probVtxCut      = iConfig.getParameter<double>         (  "ProbVertexCut" );

}

HTauTauVertexInfo::~HTauTauVertexInfo()
{

}

void HTauTauVertexInfo::beginJob(TTree* tree, edm::Service<TFileService>& fs)
{
	_tree = tree;

  // Primary vertices
  _tree->Branch("nPV",&nPV,"nPV/I");
  _tree->Branch("chi2PV",chi2PV,"chi2PV[nPV]/f");
  _tree->Branch("ndofPV",ndofPV,"ndofPV[nPV]/f");
  _tree->Branch("probPV",probPV,"probPV[nPV]/f");
  _tree->Branch("xPV",xPV,"xPV[nPV]/f");
  _tree->Branch("yPV",yPV,"yPV[nPV]/f");
  _tree->Branch("zPV",zPV,"zPV[nPV]/f");
  _tree->Branch("nTrkPV",nTrkPV,"nTrkPV[nPV]/i");
  _tree->Branch("sumPtPV",sumPtPV,"sumPtPV[nPV]/f");


  h_vtxZ_=fs->make<TH1F>("h_vtxZ_","Z Vertex",100,-100.,100.);
  h_vtxD_=fs->make<TH1F>("h_vtxD_","d Vertex",100,0.,10.);
  h_vtxProb_=fs->make<TH1F>("h_vtxProb_","Vertex Probability",50,0.,1.);

}

void HTauTauVertexInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

}

void HTauTauVertexInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	vertices.clear();
	nPV = 0;
	primaryVtxFound = false;
	edm::Handle< std::vector<reco::Vertex> > vertexCol;
	if (iEvent.getByLabel(_vtxSrc,vertexCol) && vertexCol.isValid())
	{
		//primaryVtx = vertexCol->at(0);
		//primaryVtxFound = true;
		//nPV = std::max(100, vertexCol.size());

		int numberOfVertices = vertexCol->size();
		double ptmax = 0;
		for (int iVertex = 0; iVertex< std::min(100,numberOfVertices); iVertex++)
		{
			reco::Vertex tmpVertex = vertexCol->at(iVertex);
			vertices.push_back(tmpVertex);
			//cout << "vertex x,y,z " << iVertex.x() << " "
			//                        << iVertex.y() << " "
			//                        << iVertex.z() << endl;

			double ptsum = 0;
			for(reco::Vertex::trackRef_iterator iTrack  = tmpVertex.tracks_begin(); iTrack != tmpVertex.tracks_end();++iTrack)
				ptsum += (*iTrack)->pt();


			xPV[nPV] = tmpVertex.x();
			yPV[nPV] = tmpVertex.y();
			zPV[nPV] = tmpVertex.z();

			chi2PV[nPV] = tmpVertex.chi2();
			ndofPV[nPV] = tmpVertex.ndof();
			probPV[nPV] = float(TMath::Prob(double(chi2PV[nPV]),int(ndofPV[nPV])));

			nTrkPV[nPV] = tmpVertex.tracksSize();
			sumPtPV[nPV] = ptsum;

			if(ptsum > ptmax)
			{
				// primary vertex?
				double chi2 = tmpVertex.chi2();
				int ndf = int(tmpVertex.ndof());
				double prob = TMath::Prob(chi2,ndf);
				double xVtx = tmpVertex.x();
				double yVtx = tmpVertex.y();
				double d = sqrt(xVtx*xVtx+yVtx*yVtx);
				if (fabs(tmpVertex.z())<_zVtxCut && d < _dVtxCut && prob > _probVtxCut)
				{
					ptmax = ptsum;
					primaryVtx = tmpVertex;
					primaryVtxFound = true;
				}
			}

			nPV++;
		}
	}


	if (primaryVtxFound)
	{
		h_vtxZ_->Fill(primaryVtx.z());
		h_vtxD_->Fill(sqrt(primaryVtx.x()*primaryVtx.x()+primaryVtx.y()*primaryVtx.y()));
		h_vtxProb_->Fill(TMath::Prob(primaryVtx.chi2(),int(primaryVtx.ndof())));
	}

	assert((unsigned int)nPV == vertices.size());
}

void HTauTauVertexInfo::endJob()
{

}

