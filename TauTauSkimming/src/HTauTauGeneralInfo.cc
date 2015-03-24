#include "../interface/HTauTauGeneralInfo.h"

HTauTauGeneralInfo::HTauTauGeneralInfo(const edm::ParameterSet& iConfig)
{
	_verbosity = iConfig.getParameter<int>("PrintOut");
	_datasetID = iConfig.getParameter<int>("datasetID");
	_datasetNickname = iConfig.getParameter<std::string>("datasetNickname");
	_isFirstEvent = true;
	_treeGenInfo = 0;

	// random number
	/*edm::Service<edm::RandomNumberGenerator> rng;
	if ( ! rng.isAvailable())
	{
		throw cms::Exception("Configuration")
			<< "EcalTBMCInfoProducer requires the RandomNumberGeneratorService\n"
			"which is not present in the configuration file.  You must add the service\n"
			"in the configuration file or remove the modules that require it.";
	}
	CLHEP::HepRandomEngine& engine = rng->getEngine();
	rndFlat = new CLHEP::RandFlat(engine);*/

}

HTauTauGeneralInfo::~HTauTauGeneralInfo()
{

}

void HTauTauGeneralInfo::beginJob(TTree* tree, edm::Service<TFileService>& fs)
{
	_tree = tree;

	hDataset = new TH1F("dataset", "", 1, -0.5, 0.5);
	hDataset->GetXaxis()->SetBinLabel(1, _datasetNickname.c_str());
	hDataset->Fill(0.);

	_tree->Branch("Run",&_run,"Run/I"); //name, pointer to variable (internal to class), name of leaf
	_tree->Branch("Event",&_event,"Event/I");
	_tree->Branch("Lumi",&_lumi,"Lumi/I"); //luminosity section

	_tree->Branch("datasetID",&_datasetID,"datasetID/I");

	_tree->Branch("internalXsec",&_internalXsec,"internalXsec/F");
	_tree->Branch("internalXsecE",&_internalXsecE,"internalXsecE/F");

	_tree->Branch("externalXsecLO",&_externalXsecLO,"externalXsecLO/F");
	_tree->Branch("externalXsecLOE",&_externalXsecLOE,"externalXsecLOE/F");

	_tree->Branch("externalXsecNLO",&_externalXsecNLO,"externalXsecNLO/F");
	_tree->Branch("externalXsecNLOE",&_externalXsecNLOE,"externalXsecNLOE/F");

	_tree->Branch("filterEfficiency",&_filterEfficiency,"filterEfficiency/F");

	_tree->Branch("randomNumber",&_randomNumber,"randomNumber/F");


	_treeGenInfo=fs->make<TTree>("GenInfo","GenInfo");

	_treeGenInfo->Branch("internalXsec",&_internalXsec,"internalXsec/F");
	_treeGenInfo->Branch("internalXsecE",&_internalXsecE,"internalXsecE/F");

	_treeGenInfo->Branch("externalXsecLO",&_externalXsecLO,"externalXsecLO/F");
	_treeGenInfo->Branch("externalXsecLOE",&_externalXsecLOE,"externalXsecLOE/F");

	_treeGenInfo->Branch("externalXsecNLO",&_externalXsecNLO,"externalXsecNLO/F");
	_treeGenInfo->Branch("externalXsecNLOE",&_externalXsecNLOE,"externalXsecNLOE/F");

	_treeGenInfo->Branch("filterEfficiency",&_filterEfficiency,"filterEfficiency/F");

	_treeLumi =fs->make<TTree>("LumiInfo","LumiInfo");
	_treeLumi->Branch("Run",&_run,"Run/I");
	_treeLumi->Branch("Lumi",&_lumi,"Lumi/I");

}

void HTauTauGeneralInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

}

void HTauTauGeneralInfo::beginLuminosityBlock(const edm::LuminosityBlock &lumiBlock, const edm::EventSetup &setup)
{
	_run = lumiBlock.run();
	_lumi = lumiBlock.luminosityBlock();
	_treeLumi->Fill();
}

void HTauTauGeneralInfo::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	_internalXsec = 0.;
	_internalXsecE = 0.;

	_externalXsecLO = 0.;
	_externalXsecLOE = 0.;

	_externalXsecNLO = 0.;
	_externalXsecNLOE = 0.;

	_filterEfficiency = 0.;

	_run = iEvent.id().run();
	_event = iEvent.id().event();
	_lumi = iEvent.getLuminosityBlock().luminosityBlock();

	if (_verbosity > 0)
	{
		std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
		std::cout << std::endl;
		std::cout << "Run = " << _run << "  event = " << _event << "  lumi = " << _lumi << std::endl;
	}

        edm::Handle<GenRunInfoProduct> genInfoHandle;
        std::vector<edm::Handle<GenRunInfoProduct>> vec;
        vec.push_back(genInfoHandle);
        iEvent.getRun().getManyByType(vec);
	//	if (iEvent.getRun().getByType(genInfoHandle))
	//	{
		_internalXsec = float(genInfoHandle->internalXSec().value());
		_internalXsecE = float(genInfoHandle->internalXSec().error());

		_externalXsecLO = float(genInfoHandle->externalXSecLO().value());
		_externalXsecLOE = float(genInfoHandle->externalXSecLO().error());

		_externalXsecNLO = float(genInfoHandle->externalXSecNLO().value());
		_externalXsecNLOE = float(genInfoHandle->externalXSecNLO().error());

		_filterEfficiency = float(genInfoHandle->filterEfficiency());

		if (_isFirstEvent)
		{
			_treeGenInfo->Fill();
			_isFirstEvent = false;
		}
		//	}

	_randomNumber = rndFlat->fire(0., 1.);
}

void HTauTauGeneralInfo::endJob()
{
	_tree->GetDirectory()->cd();
	hDataset->Write();
	if (_treeGenInfo)
	{
		_treeGenInfo->GetDirectory()->cd();
		_treeGenInfo->Write();
	}
	if (_treeLumi)
	{
		_treeGenInfo->GetDirectory()->cd();
		_treeGenInfo->Write();
	}
}

