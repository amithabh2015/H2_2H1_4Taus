#include "../interface/HTauTauPileUpInfo.h"

HTauTauPileUpInfo::HTauTauPileUpInfo(const edm::ParameterSet& iConfig)
{
	_verbosity = iConfig.getParameter<int>("PrintOut");
//	_rhoSrc         = iConfig.getParameter<edm::InputTag>("RhoSource");
//	_rhoNeutralSrc  = iConfig.getParameter<edm::InputTag>("RhoNeutralSource");

	_pfParticleSrc  = iConfig.getParameter<edm::InputTag>("PFParticleSource");
	_pfPileUpSrc    = iConfig.getParameter<edm::InputTag>("pfPileUpSource");
        _pfNoPileUpSrc  = iConfig.getParameter<edm::InputTag>("pfNoPileUpSource");

	// load default values
	_nPtBinsPFSum = 15;
	float const ptThresholds[15] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0};

	for (size_t i = 0; i != _nPtBinsPFSum; ++i)
	    _ptThresholds[i] = ptThresholds[i];
}

HTauTauPileUpInfo::~HTauTauPileUpInfo()
{

}

HTauTauPileUpInfo& HTauTauPileUpInfo::SetPtThresholds(std::vector<float> ptThresholds)
{
	if (ptThresholds.size() > 30) {
	  ReportNewPtThresholdVectorIsTooLong();
	  return *this;
	}

	_nPtBinsPFSum = ptThresholds.size();

	for (size_t i = 0; i != _nPtBinsPFSum; ++i ) {
	  _ptThresholds[i] = ptThresholds[i];
	}

	return *this;
}

void HTauTauPileUpInfo::ReportNewPtThresholdVectorIsTooLong()
{
  std::cout << "HTauTauPileUpInfo error: SetPtThresholds accepts a vector of lower edges of up to 30 floats, ignoring" << std::endl;
}

std::vector<float> HTauTauPileUpInfo::GetPtThresholds(){

	std::vector<float> returnedVector;

	for (size_t iPtBin=0; iPtBin < _nPtBinsPFSum; ++iPtBin){
	    returnedVector.push_back(_ptThresholds[iPtBin]);
	}

	return returnedVector;
}

void HTauTauPileUpInfo::beginJob(TTree* tree, edm::Service<TFileService>& fs)
{
	_tree = tree;
//	_tree->Branch("Rho",&_rho,"Rho/f");
//	_tree->Branch("RhoNeutral",&_rhoNeutral,"RhoNeutral/f");

	_tree->Branch("nPUTruth",&nPUTruth,"nPUTruth/F");
	_tree->Branch("nPUI",&nPUI,"nPUI/I");
	_tree->Branch("nPUIM1",&nPUIM1,"nPUIM1/I");
	_tree->Branch("nPUIP1",&nPUIP1,"nPUIP1/I");


	_tree->Branch("NPtBinsPFSum",&_nPtBinsPFSum,"NPtBinsPFSum/I");

	_tree->Branch("PFChargedSum",_pfChargedSum,"PFChargedSum[NPtBinsPFSum]/f");
	_tree->Branch("PFNeutralSum",_pfNeutralSum,"PFNeutralSum[NPtBinsPFSum]/f");
	_tree->Branch("PFChargedNoPileUpSum",_pfChargedNoPileUpSum,"PFChargedNoPileUpSum[NPtBinsPFSum]/f");
	_tree->Branch("PFNeutralNoPileUpSum",_pfNeutralNoPileUpSum,"PFNeutralNoPileUpSum[NPtBinsPFSum]/f");
	_tree->Branch("PFChargedPileUpSum",_pfChargedPileUpSum,"PFChargedPileUpSum[NPtBinsPFSum]/f");

	_tree->Branch("NPFCharged",_nPfCharged,"NPFCharged[NPtBinsPFSum]/I");
	_tree->Branch("NPFNeutral",_nPfNeutral,"NPFNeutral[NPtBinsPFSum]/I");
	_tree->Branch("NPFChargedNoPileUp",_nPfChargedNoPileUp,"NPFChargedNoPileUp[NPtBinsPFSum]/I");
	_tree->Branch("NPFNeutralNoPileUp",_nPfNeutralNoPileUp,"NPFNeutralNoPileUp[NPtBinsPFSum]/I");
	_tree->Branch("NPFChargedPileUp",_nPfChargedPileUp,"NPFChargedPileUp[NPtBinsPFSum]/I");

	h_nPUI = fs->make<TH1F>("h_nPUI_", "n PileUp interactions", 51, -0.5, 50.5);
}

void HTauTauPileUpInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

}

void HTauTauPileUpInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	// store information on rho
	/*_rho = 0;
	_rhoNeutral = 0;

	edm::Handle<double> rhoPtr;

	if (iEvent.getByLabel(_rhoSrc,rhoPtr))
		_rho = *rhoPtr;
	else
		std::cout << "Collection " <<  _rhoSrc << " is absent in an Event Container" << std::endl;

	if (iEvent.getByLabel(_rhoNeutralSrc,rhoPtr))
		_rhoNeutral = *rhoPtr;
	else
		std::cout << "Collection " <<  _rhoNeutralSrc << " is absent in an Event Container" << std::endl;

	if (_verbosity > 0)
	{
		std::cout << std::endl;
		std::cout << "Rho value = " << _rho << std::endl;
		std::cout << "Rho (neutral) value = " << _rhoNeutral << std::endl;
	}*/


	// store information on the pile-up truth
	nPUTruth = -1.;
	nPUI = -1;
	nPUIM1 = -1;
	nPUIP1 = -1;

	edm::Handle<std::vector< PileupSummaryInfo > > hPUInfo;

	if (iEvent.getByLabel(edm::InputTag("addPileupInfo","","HLT"), hPUInfo))
	{
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		for (PVI = hPUInfo->begin(); PVI != hPUInfo->end(); ++PVI)
		{
			int BX = PVI->getBunchCrossing();

			if(BX == 0)
			{
				nPUI = PVI->getPU_NumInteractions();

				// the following statement is available beginning with 4.2.8
				// this information is not available in Summer11 (the number is then
				// filled randomly, in most cases with values around 1e-43)
				nPUTruth = PVI->getTrueNumInteractions();
			}

			if(BX == -1)
				nPUIM1 = PVI->getPU_NumInteractions();
			if(BX == 1)
				nPUIP1 = PVI->getPU_NumInteractions();
		}
	}
	else
	{
		if (_verbosity>0)
		{
			std::cout << std::endl;
			std::cout << "Collection addPileupInfo is not found" << std::endl;
		}
	}

	// Fill n PU interactions nPUI
	h_nPUI->Fill(nPUI);


	//*****************
	// PU pT spectrum *
	//*****************

	// reset counters
	for (size_t iPtBins=0; iPtBins < _nPtBinsPFSum; ++iPtBins){
	    _pfChargedSum[iPtBins] = 0.0;
	    _pfNeutralSum[iPtBins] = 0.0;
	    _pfChargedPileUpSum[iPtBins] = 0.0;
	    _pfChargedNoPileUpSum[iPtBins] = 0.0;
	    _pfNeutralNoPileUpSum[iPtBins] = 0.0;

	    _nPfCharged[iPtBins] = 0;
	    _nPfNeutral[iPtBins] = 0;
	    _nPfChargedPileUp[iPtBins] = 0;
	    _nPfChargedNoPileUp[iPtBins] = 0;
	    _nPfNeutralNoPileUp[iPtBins] = 0;
	}

	edm::Handle<reco::PFCandidateCollection> pflow;


	if (iEvent.getByLabel(_pfParticleSrc,pflow))
	{
		int Npflow 		= pflow->size();
		if (_verbosity > 0 )
			std::cout << "   N(pflow) = " << Npflow << "\n";
		for (int iP = 0; iP < Npflow; ++iP)
		{
			reco::PFCandidateRef pfCand(pflow, iP);
			float ptPfCand = pfCand->pt();
			float chargePfCand = fabs(pfCand->charge());
			for (size_t iPtBin = 0; iPtBin < _nPtBinsPFSum; ++iPtBin)
			{
				if (ptPfCand>_ptThresholds[iPtBin])
				{
					if (chargePfCand>0.1) {
						_pfChargedSum[iPtBin] += ptPfCand;
						_nPfCharged[iPtBin]++;
					}
					else {
						_pfNeutralSum[iPtBin] += ptPfCand;
						_nPfNeutral[iPtBin]++;
					}
				}
			}
		}
	}
	else
	{
		std::cout << "pflow not found in HTauTauPileUpInfo: " << _pfParticleSrc << "\n";
	}
	
	
	edm::Handle<reco::PFCandidateCollection> pfNoPileUp;
	
	if (iEvent.getByLabel(_pfNoPileUpSrc,pfNoPileUp))
	{
		int NpfNoPileUp = pfNoPileUp->size();
		if (_verbosity > 0 )
			std::cout << "   N(pfNoPileUp) = " << NpfNoPileUp << "\n";
		for (int iP = 0; iP < NpfNoPileUp; ++iP)
		{
			reco::PFCandidateRef pfCand(pfNoPileUp, iP);
			float ptPfCand = pfCand->pt();
			float chargePfCand = fabs(pfCand->charge());
			for (size_t iPtBin=0; iPtBin<_nPtBinsPFSum; ++iPtBin)
			{
				if (ptPfCand>_ptThresholds[iPtBin])
				{
					if (chargePfCand>0.1)
					{
						_pfChargedNoPileUpSum[iPtBin] += ptPfCand;
						_nPfChargedNoPileUp[iPtBin]++;
					}
					else
					{
						_pfNeutralNoPileUpSum[iPtBin] += ptPfCand;
						_nPfNeutralNoPileUp[iPtBin]++;
					}
				}
			}
		}
	}
	else
	{
		std::cout << "pfNoPileUp not found in HTauTauPileUpInfo: " << _pfNoPileUpSrc << "\n";
	}
	
	
	edm::Handle<reco::PFCandidateCollection> pfPileUp;
	
	if (iEvent.getByLabel(_pfPileUpSrc,pfPileUp))
	{
		int NpfPileUp 	= pfPileUp->size();
		if (_verbosity > 0 )
			std::cout << "   N(pfPileUp) = " << NpfPileUp << "\n";
		for (int iP = 0; iP < NpfPileUp; ++iP)
		{
			reco::PFCandidateRef pfCand(pfPileUp, iP);
			float ptPfCand = pfCand->pt();
			//	float pfCandCharge = pfCand->charge();

			//	std::cout << "PileUp Candidate q : " << pfCandCharge << " pt = " << ptPfCand;
			//	if (fabs(pfCandCharge)<0.1) std::cout << "  *** ";
			//	std::cout << std::endl;

			for (size_t iPtBin=0; iPtBin<_nPtBinsPFSum; ++iPtBin) {
				if (ptPfCand>_ptThresholds[iPtBin]) {
					_pfChargedPileUpSum[iPtBin] += ptPfCand;
					_nPfChargedPileUp[iPtBin]++;
				}
			}
		}
	}
	else
	{
		std::cout << "pfPileUp not found in HTauTauPileUpInfo: " << _pfPileUpSrc << "\n";
	}

}

void HTauTauPileUpInfo::endJob()
{

}

