#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauGeneralInfo
#define H2to2H1to4Taus_TauTauSkimming_HTauTauGeneralInfo

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandFlat.h"

#include "TString.h"
#include "TTree.h"

#include <vector>

class HTauTauGeneralInfo
{
public:
	HTauTauGeneralInfo(const edm::ParameterSet&);
	~HTauTauGeneralInfo();

	void beginJob(TTree* tree, edm::Service<TFileService>& fs);
	void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
	void beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &setup);
	void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void endJob();

private:
	TTree* _tree;
	TTree* _treeGenInfo;
	TTree* _treeLumi;
	TH1F* hDataset;

	int _run; // run number
	int _event; // event number
	int _lumi; // lumi section


	CLHEP::RandFlat* rndFlat;

	int _verbosity;
	bool _isFirstEvent;

	int _datasetID;
	std::string _datasetNickname;

	// ** generator info (x-sections)
	float _internalXsec;
	float _internalXsecE;
	float _externalXsecLO;
	float _externalXsecLOE;
	float _externalXsecNLO;
	float _externalXsecNLOE;
	float _filterEfficiency;

	float _randomNumber;
};
#endif
