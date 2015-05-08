#ifndef H2to2H1to4Taus_TauTauSkimming_HTauTauTriggerInfo
#define H2to2H1to4Taus_TauTauSkimming_HTauTauTriggerInfo

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>

#include "TString.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <vector>

class HTauTauTriggerInfo
{
public:
	HTauTauTriggerInfo(const edm::ParameterSet&);
	~HTauTauTriggerInfo();

	void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
	void analyzeObject(const reco::Candidate::LorentzVector& p4, const char* name, const edm::Event& iEvent, const edm::EventSetup& iSetup, float deltaR);
	void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);
	void beginJob(TTree* tree, edm::Service<TFileService>& fs);
	void endJob();

	// hasMatchingTriggerObject() returns a vector with booleans, one
	// for each candidate, indicating if the given trigger module
	// provides a matching trigger object
	std::vector<bool> hasMatchingTriggerObject(/*std::vector<reco::Candidate> candidates*/reco::Candidate::LorentzVector CandP4, std::string triggerModuleName, const edm::Event& iEvent, const edm::EventSetup& iSetup, float deltaR = 0.5);
	// return for a given HLT name a list of trigger modules which
	// produce trigger objects
	std::vector<std::string> getTriggerModules(std::string hltName, const edm::Event& iEvent, const edm::EventSetup& iSetup);
	// same as previous one, but only the last one
	std::string getLastTriggerModule(std::string hltName, const edm::Event& iEvent, const edm::EventSetup& iSetup);
	// input:  HLT_DoubleEle33_CaloIdL
	// output: HLT_DoubleEle33_CaloIdL_v9
	std::string getTriggerName(std::string hltName, const edm::Event& iEvent);
private:
	TTree* _tree;

	edm::InputTag _triggerSrc;
	//edm::InputTag _origTriggerSrc;

	int _verbosity;
	bool _firstTriggerEvent;

	std::vector<std::string> _selTrigger;

	std::map<std::string, bool> _hlt;
	std::map<std::string, int> _hlt_prescale;
	//std::map<std::string, bool> _origHLT;
	//std::map<std::string, int> _origHLT_prescale;

	// For trigger matching:
	struct TriggerMatchFilter {
		std::vector<std::string> modules;
		bool allowNonExisting;
		std::map<std::string, bool> match; // object -> match
	};

	// triggername -> filternick -> TriggerMatchFilter
	std::map<std::string, std::map<std::string, TriggerMatchFilter> > _triggerMatchFilters;

	// orig HLT config to read orig HLT prescale
	//HLTConfigProvider origHLTConfig;

	HLTConfigProvider hltConfigProvider;
	std::map<std::string, std::string> hltMapping;
};

#endif 
