#include "../interface/HTauTauTriggerInfo.h"

HTauTauTriggerInfo::HTauTauTriggerInfo(const edm::ParameterSet& iConfig)
{
	_verbosity = iConfig.getParameter<int>("PrintOut");

	_triggerSrc     = iConfig.getParameter<edm::InputTag>("TriggerSource");
	//_origTriggerSrc = iConfig.getParameter<edm::InputTag>("OrigTriggerSource");

	_selTrigger = iConfig.getParameter<std::vector<std::string> >("Triggers");

	const edm::ParameterSet& trigMatchPset = iConfig.getParameter<edm::ParameterSet>("TriggerMatching");
	const std::vector<std::string> trigMatchObjects = trigMatchPset.getParameter<std::vector<std::string> >("Objects");
	const std::vector<edm::ParameterSet> trigMatchTriggers = trigMatchPset.getParameter<std::vector<edm::ParameterSet> >("Triggers");
	for(std::vector<edm::ParameterSet>::const_iterator iter = trigMatchTriggers.begin(); iter != trigMatchTriggers.end(); ++iter)
	{
		const std::string triggerName = iter->getParameter<std::string>("TriggerName");
		std::map<std::string, TriggerMatchFilter>& nickFilterMap = _triggerMatchFilters[triggerName];

		const std::vector<edm::ParameterSet> trigMatchFilterModules = iter->getParameter<std::vector<edm::ParameterSet> >("FilterModules");
		for(std::vector<edm::ParameterSet>::const_iterator iter = trigMatchFilterModules.begin(); iter != trigMatchFilterModules.end(); ++iter) // note shadowing
		{
			const std::string filterNick = iter->getParameter<std::string>("Nick");
			TriggerMatchFilter& filter = nickFilterMap[filterNick];

			filter.allowNonExisting = iter->getUntrackedParameter<bool>("AllowNonExisting", false);
			filter.modules = iter->getParameter<std::vector<std::string> >("Modules");
			for(std::vector<std::string>::const_iterator iter = trigMatchObjects.begin(); iter != trigMatchObjects.end(); ++iter)
				filter.match[*iter] = false; // default value
		}
	}
}

HTauTauTriggerInfo::~HTauTauTriggerInfo()
{

}

void HTauTauTriggerInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
	/*if(!_origTriggerSrc.label().empty())
	{
		bool hltSetupChanged = false;
		if(!origHLTConfig.init(iRun, iSetup, _origTriggerSrc.process(), hltSetupChanged))
			throw cms::Exception("Failed to init orig HLT config");
	}*/

	if(_triggerSrc.process().empty())
	{
		std::cout << "process name is not specified for triggerSrc -> determine the process name automatically" << std::endl;

		std::cout << "the following processes are available:" << std::endl;

		edm::Handle<trigger::TriggerEvent> tmpTriggerEventHLT;

		const edm::ProcessHistory& processHistory(iRun.processHistory());
		for (edm::ProcessHistory::const_iterator it = processHistory.begin(); it != processHistory.end(); ++it)
		{
			std::cout << "\t" << it->processName();
			edm::ProcessConfiguration processConfiguration;
			if (processHistory.getConfigurationForProcess(it->processName(), processConfiguration))
			{
				edm::ParameterSet processPSet;
				if (edm::pset::Registry::instance()->getMapped(processConfiguration.parameterSetID(), processPSet))
				{
					if (processPSet.exists("hltTriggerSummaryAOD"))
					{
						_triggerSrc = edm::InputTag(_triggerSrc.label(), "", it->processName());
						std::cout << "*";
					}
				}
			}
			std::cout << std::endl;
		}
		std::cout << "* process with hltTriggerSummaryAOD" << std::endl;
		std::cout << "selected process: " << _triggerSrc.process() << std::endl;
	}

	bool hltSetupChanged = false;
	if(!hltConfigProvider.init(iRun, iSetup, _triggerSrc.process(), hltSetupChanged))
		throw cms::Exception("Failed to init the main HLT config provider");

	hltMapping.clear();
}

void HTauTauTriggerInfo::beginJob(TTree* tree, edm::Service<TFileService>& fs)
{
	_tree = tree;

	_firstTriggerEvent = true;

	for (std::vector<std::string>::const_iterator it = _selTrigger.begin(); it != _selTrigger.end(); ++it)
	{
		std::string tmpString = *it;
		if(_verbosity > 0)
			std::cout << "register branches for HLT_" << tmpString << "\n";
		_tree->Branch(("HLT_"+tmpString).c_str(), &_hlt[*it], ("HLT_"+tmpString+"/O").c_str());
		_tree->Branch(("HLT_PRESCALE_"+tmpString).c_str(), &_hlt_prescale[*it], ("HLT_PRESCALE_"+tmpString+"/I").c_str());

		/*if(!_origTriggerSrc.label().empty())
		{
			_tree->Branch(("ORIG_HLT_"+tmpString).c_str(), &_origHLT[*it], ("ORIG_HLT_"+tmpString+"/O").c_str());
			_tree->Branch(("ORIG_HLT_PRESCALE_"+tmpString).c_str(), &_origHLT_prescale[*it], ("ORIG_HLT_PRESCALE_"+tmpString+"/I").c_str());
		}*/
	}

	for(std::map<std::string, std::map<std::string, TriggerMatchFilter> >::iterator trig_iter = _triggerMatchFilters.begin(); trig_iter != _triggerMatchFilters.end(); ++trig_iter)
	{
		std::map<std::string, TriggerMatchFilter>& nickFilterMap = trig_iter->second;
		for(std::map<std::string, TriggerMatchFilter>::iterator nick_iter = nickFilterMap.begin(); nick_iter != nickFilterMap.end(); ++nick_iter)
		{
			TriggerMatchFilter& filter = nick_iter->second;
			for(std::map<std::string, bool>::iterator obj_iter = filter.match.begin(); obj_iter != filter.match.end(); ++obj_iter)
			{
				// Branch name: HLT_<Name>_TrigMatch_<FilterNick>_<Object>
				const std::string branchName = "HLT_" + trig_iter->first + "_TrigMatch_" + nick_iter->first + "_" + obj_iter->first;
				_tree->Branch(branchName.c_str(), &obj_iter->second, (branchName + "/O").c_str());
			}
		}
	}
}

void HTauTauTriggerInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	// initialise variables for triggers
	for (std::vector<std::string>::const_iterator it = _selTrigger.begin(); it != _selTrigger.end(); ++it)
	{
		_hlt[*it] = false;
		_hlt_prescale[*it] = 0;
//		_origHLT[*it] = false;
//		_origHLT_prescale[*it] = 0;
	}

	// Trigger Information
	edm::Handle<edm::TriggerResults> triggerResults;

	if (iEvent.getByLabel(_triggerSrc, triggerResults) && triggerResults.isValid())
	{
		edm::TriggerNames const & triggerNames = iEvent.triggerNames(*triggerResults);
		int ntrigs = triggerResults->size();

		if(!hltConfigProvider.inited())
			throw cms::Exception("HLT config was not initialized yet. This should have happened in beginRun!");
		assert(ntrigs == (int)hltConfigProvider.size());

		if (_firstTriggerEvent)
		{
			for (int it=0;it<ntrigs;++it)
			{
				TString trigName = triggerNames.triggerName(it);
				bool accept = triggerResults->accept(it);
				if (_verbosity > 0)
					std::cout << it << "	 " << trigName << "	" << accept << std::endl;
			}
			_firstTriggerEvent = false;
		}

		for (int it=0;it<ntrigs;++it)
		{

			std::string trigName = triggerNames.triggerName(it);
			const bool accept = triggerResults->accept(it);

			for (std::vector<std::string>::const_iterator it = _selTrigger.begin(); it != _selTrigger.end(); ++it)
			{
				std::string tmpString = "HLT_"+(*it)+"_v";
				if (TString(trigName).Contains(tmpString)||trigName==TString("HLT_"+(*it)))
				{
				  //					std::pair<int, int> prescale = hltConfigProvider.prescaleValues(iEvent, iSetup, trigName);

					_hlt[*it] = accept;
					//					if(prescale.first >= 0 && prescale.second >= 0)
					//						_hlt_prescale[*it] = prescale.first * prescale.second;
					//					else
					//
					_hlt_prescale[*it] = 1;

					if (_verbosity>1)
						std::cout << trigName << " : " << accept << "	prescale : " << _hlt_prescale[*it] << std::endl;
				}
			}
		}
		if (_verbosity>1)
			std::cout << std::endl;
	}

	/*if(!_origTriggerSrc.label().empty())
	{
		edm::Handle<edm::TriggerResults> origTriggerResults;
		iEvent.getByLabel(_origTriggerSrc, origTriggerResults);

		const edm::TriggerNames& origTriggerNames = iEvent.triggerNames(*origTriggerResults);
		int nOrigTrigs = origTriggerResults->size();

		if(!origHLTConfig.inited())
			throw cms::Exception("orig HLT config was not initialized yet. This should have happened in beginRun!");
		assert(nOrigTrigs == (int)origHLTConfig.size());

		for (int i=0;i<nOrigTrigs;++i)
		{
			std::string trigName = origTriggerNames.triggerName(i);
			for (std::vector<std::string>::const_iterator it = _selTrigger.begin(); it != _selTrigger.end(); ++it)
			{
				std::string tmpString = "HLT_"+(*it)+"_v";
				if (TString(trigName).Contains(tmpString)||trigName==TString("HLT_"+(*it)))
				{
					std::pair<int, int> prescale = origHLTConfig.prescaleValues(iEvent, iSetup, trigName);

					_origHLT[*it] = origTriggerResults->accept(i);
					if(prescale.first >= 0 && prescale.second >= 0)
						_origHLT_prescale[*it] = prescale.first * prescale.second;
					else
						_origHLT_prescale[*it] = 0;
					break;
				}
			}
		}
	}*/
}

void HTauTauTriggerInfo::analyzeObject(const reco::Candidate::LorentzVector& p4, const char* name, const edm::Event& iEvent, const edm::EventSetup& iSetup, float deltaR)
{
	edm::Handle<edm::TriggerResults> triggerResults;
	iEvent.getByLabel(_triggerSrc, triggerResults);

	edm::Handle<trigger::TriggerEvent> triggerEventHandle;
	iEvent.getByLabel("hltTriggerSummaryAOD",triggerEventHandle);

	for(std::map<std::string, std::map<std::string, TriggerMatchFilter> >::iterator trig_iter = _triggerMatchFilters.begin(); trig_iter != _triggerMatchFilters.end(); ++trig_iter)
	{
		const std::string& triggerName = trig_iter->first;
		std::map<std::string, TriggerMatchFilter>& nickFilterMap = trig_iter->second;
		for(std::map<std::string, TriggerMatchFilter>::iterator nick_iter = nickFilterMap.begin(); nick_iter != nickFilterMap.end(); ++nick_iter)
		{
			TriggerMatchFilter& triggerMatchFilter = nick_iter->second;

			std::map<std::string, bool>::iterator obj_iter = triggerMatchFilter.match.find(name);
			if(obj_iter == triggerMatchFilter.match.end())
				throw cms::Exception("HTauTauTriggerInfo") << "No such triggermatch object registered: " << name;

			// Note we could use hasMatchingTriggerObject() here, but this is more efficient since we need to
			// traverse the filter list only once for all the different modules whose presence we check.
			std::size_t filterIndex = 0;
			bool foundFilter = false;
			for(std::size_t iF = 0; iF < triggerEventHandle->sizeFilters() && !foundFilter; ++iF)
			{
				const std::string& filterName = triggerEventHandle->filterTag(iF).label();
				for(std::vector<std::string>::const_iterator module_iter = triggerMatchFilter.modules.begin(); module_iter != triggerMatchFilter.modules.end() && !foundFilter; ++module_iter)
				{
					if(filterName == *module_iter)
					{
						filterIndex = iF;
						foundFilter = true;
					}
				}
			}

			if(!foundFilter)
			{
				// We did not find one of the filters we are looking for. In principle this is an error and should not happen.
				// However, it can happen if the trigger is not available in the trigger menu, so let's check this.
				const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults);
				bool foundTrigger = false;
				//bool triggerWasRun = false;
				bool triggerAccept = false;
				std::string realTriggerName;
				for (size_t j = 0; j < triggerNames.size() && !foundTrigger; ++j)
				{
					const std::string& trigName = triggerNames.triggerName(j);
					if(trigName.find("HLT_" + triggerName + "_v") == 0 || trigName == TString("HLT_" + triggerName))
					{
						foundTrigger = true;
						realTriggerName = trigName;
						//triggerWasRun = triggerResults->wasrun(j);
						triggerAccept = triggerResults->accept(j);
					}
				}

				// Even if the trigger is available, sometimes the filter does not show up. I am not exactly sure
				// why this is, maybe the filter was not run at all in such an event? Or if it has not fired?
				// What we want to make sure here is that there is at least one module in the trigger path that we are looking for.
				bool foundTriggerModule = false;
				if(foundTrigger)
				{
					// If the filter was not found, require that the trigger has not fired.
					// If the trigger has fired, the filter must be there, so we can query it
					// for trigger objects.
					assert(!triggerAccept || triggerMatchFilter.allowNonExisting);

					std::vector<std::string> triggerModules = hltConfigProvider.saveTagsModules(realTriggerName);
					for(std::vector<std::string>::const_iterator iter = triggerModules.begin(); iter != triggerModules.end() && !foundTriggerModule; ++iter)
						if(std::find(triggerMatchFilter.modules.begin(), triggerMatchFilter.modules.end(), *iter) != triggerMatchFilter.modules.end())
							foundTriggerModule = true;
				}

				// If there is no such module, we consider it as an error -- there is no chance to match something to a
				// nonexisting trigger object. The most likely cause of this is that the configuration is wrong, and
				// that a filter module is incorrectly specified or missing.
				if(!foundTrigger || foundTriggerModule || triggerMatchFilter.allowNonExisting)
				{
					// OK, fair enough, trigger does not exist, so we cannot match anything against it.
					obj_iter->second = false;
				}
				else
				{
					throw cms::Exception("HTauTauTriggerInfo") << "No filter module for trigger match " << nick_iter->first << " and trigger " << triggerName << " found in run " << iEvent.run() << "!\n"
					                                           << "Please check the trigger menu for that run and add the filter module to the configuration.";
				}
			}
			else
			{
				const trigger::Keys& keys = triggerEventHandle->filterKeys(filterIndex);
				for(std::size_t iK = 0; iK < keys.size(); ++iK)
				{
					const trigger::TriggerObject triggerObject(triggerEventHandle->getObjects().at(keys[iK]));
					const reco::Particle::PolarLorentzVector trigP4(triggerObject.pt(), triggerObject.eta(), triggerObject.phi(), triggerObject.mass());
					if(ROOT::Math::VectorUtil::DeltaR(p4, trigP4) < deltaR)
						obj_iter->second = true;
					else
						obj_iter->second = false;
				}
			}
		}
	}
}

std::vector<bool> HTauTauTriggerInfo::hasMatchingTriggerObject(/*std::vector<reco::Candidate> candidates*/reco::Candidate::LorentzVector CandP4, std::string triggerModuleName, const edm::Event& iEvent, const edm::EventSetup& iSetup, float deltaR)
{
   edm::Handle<trigger::TriggerEvent> triggerEventHandle;
   iEvent.getByLabel("hltTriggerSummaryAOD",triggerEventHandle);

	std::vector<bool> res(1/*candidates.size()*/, false);


	for (size_t iF = 0; iF < triggerEventHandle->sizeFilters(); ++iF)
	{

	  const std::string filterName(triggerEventHandle->filterTag(iF).label());
	  if (triggerModuleName != /*"TriggerObject_" +*/ filterName)
	    continue;
	  
	 
	  
	  const trigger::Keys & keys = triggerEventHandle->filterKeys(iF);
	  for (size_t iK = 0; iK < keys.size(); ++iK)
	    {
	      trigger::TriggerObject triggerObject(triggerEventHandle->getObjects().at(keys[iK]));
	      reco::Particle::PolarLorentzVector tmpTriggerObject(triggerObject.pt(), triggerObject.eta(), triggerObject.phi(), triggerObject.mass());
	      for (size_t idx = 0; idx < 1/*candidates.size()*/; ++idx)
		{
		 
		  if (ROOT::Math::VectorUtil::DeltaR(tmpTriggerObject, /*candidates.at(idx).p4()*/CandP4) < deltaR)
		    res[idx] = true;
		}
	    }
	}

	return res;
}

// input:  HLT_DoubleEle33_CaloIdL
// output: HLT_DoubleEle33_CaloIdL_v9
std::string HTauTauTriggerInfo::getTriggerName(std::string hltName, const edm::Event& iEvent)
{
	edm::Handle<edm::TriggerResults> triggerResults;
	if (iEvent.getByLabel(_triggerSrc, triggerResults) && triggerResults.isValid())
	{
		edm::TriggerNames const & triggerNames = iEvent.triggerNames(*triggerResults);

		if (hltMapping.find(hltName) == hltMapping.end())
		{
			for (size_t j = 0; j < hltConfigProvider.size(); ++j)
			{
				std::string tmpString =hltName+"_v";
				TString trigName = triggerNames.triggerName(j);
				if (trigName.Contains(tmpString)||trigName==TString(hltName))
				{
					hltMapping[hltName] = hltConfigProvider.triggerName(j);
					return hltMapping[hltName];
				}
			}
		}
		else
			return hltMapping[hltName];
	}
	return "";
}

std::vector<std::string> HTauTauTriggerInfo::getTriggerModules(std::string hltName, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if(!hltConfigProvider.inited())
		throw cms::Exception("HLT config was not initialized yet. This should have happened in beginRun!");

	std::vector<std::string> triggerObjects;
	hltName = getTriggerName(hltName, iEvent);

	std::vector<std::string> modules = hltConfigProvider.saveTagsModules(hltName);
	for (std::vector<std::string>::const_iterator it2 = modules.begin(); it2 != modules.end(); ++it2)
	{
		if (std::find(triggerObjects.begin(), triggerObjects.end(), *it2) == triggerObjects.end())
		{
			triggerObjects.push_back(*it2);
		}
	}
	return triggerObjects;
}

std::string HTauTauTriggerInfo::getLastTriggerModule(std::string hltName, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	return getTriggerModules(hltName, iEvent, iSetup).back();
}
