#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauPileUpInfo
#define H2to2H1to4Taus_TauTauSkimming_HTauTauPileUpInfo

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TString.h"
#include "TTree.h"

#include <vector>
#include <iostream>

#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauTypeDef.h"
//#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetElectrons.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauGetMuons.h"

class HTauTauPileUpInfo
{
public:
	HTauTauPileUpInfo(const edm::ParameterSet&);
	~HTauTauPileUpInfo();

	void beginJob(TTree* tree, edm::Service<TFileService>& fs);
	void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
	void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void endJob();

	HTauTauPileUpInfo& SetPtThresholds(std::vector<float> ptThresholds);
	void ReportNewPtThresholdVectorIsTooLong();
	
	std::vector<float> GetPtThresholds();

private:
	TTree* _tree;
	TH1F * h_nPUI;

	int _verbosity;

	//edm::InputTag _rhoSrc;
	//edm::InputTag _rhoNeutralSrc;

	//float _rho;
	//float _rhoNeutral;

	float nPUTruth;
	int nPUI;
	int nPUIM1;	// bx = -1
	int nPUIP1;	// bx = +1


	edm::InputTag _pfParticleSrc;
	edm::InputTag _pfPileUpSrc;
	edm::InputTag _pfNoPileUpSrc;

	size_t	_nPtBinsPFSum;
	float	_ptThresholds[30];

	float	_pfChargedSum[30];
	float	_pfNeutralSum[30];
	float	_pfChargedPileUpSum[30];
	float	_pfChargedNoPileUpSum[30];
	float	_pfNeutralNoPileUpSum[30];

	int	_nPfCharged[30];
	int	_nPfNeutral[30];
	int	_nPfChargedPileUp[30];
	int	_nPfChargedNoPileUp[30];
	int	_nPfNeutralNoPileUp[30];

};

#endif 
