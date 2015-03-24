#ifndef DesyHTauTau_TauTauSkimming_HTauTauMassInfo
#define DesyHTauTau_TauTauSkimming_HTauTauMassInfo

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/METReco/interface/MET.h"

#include "DataFormats/Common/interface/Ref.h"

#include "TString.h"
#include "TTree.h"

#include <vector>

template <typename T1, typename T2>
class HTauTauMassInfo
{
public:
	HTauTauMassInfo(const edm::ParameterSet&);
	~HTauTauMassInfo();

	void beginJob(TTree* tree, edm::Service<TFileService>& fs);
	void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);

	void analyze(const reco::Candidate& lep1, const reco::Candidate& lep2, const reco::MET& met);
	void endJob();

private:
	TTree* _tree;

	int _verbosity;

	// by integration --->
	bool _validSVFitMassInt;
	float _svFitMassInt;
	float _svFitMassErrUpInt;
	float _svFitMassErrDownInt;
	float _svFitPtInt;
	float _svFitEtaInt;
	float _svFitPhiInt;
	float _svFitPxInt;
	float _svFitPyInt;
	float _svFitPzInt;
	float _svFitEnergyInt;

	float _svFitMassCorr;
	float _svFitMassErrUpCorr;
	float _svFitMassErrDownCorr;
	float _svFitPtCorr;
	float _svFitEtaCorr;
	float _svFitPhiCorr;
	float _svFitPxCorr;
	float _svFitPyCorr;
	float _svFitPzCorr;
	float _svFitEnergyCorr;

	// standalone
	bool _validSVFitMassStandalone;
	float _svFitMassStandalone;
	float _svFitMassErrUpStandalone;
	float _svFitMassErrDownStandalone;

	float _svFitPtStandalone;
	float _svFitEtaStandalone;
	float _svFitPhiStandalone;

	float _svFitPxStandalone;
	float _svFitPyStandalone;
	float _svFitPzStandalone;
	float _svFitEnergyStandalone;

	float _svFitTauPosPxStandalone;
	float _svFitTauPosPyStandalone;
	float _svFitTauPosPzStandalone;
	float _svFitTauPosPtStandalone;
	float _svFitTauPosEtaStandalone;
	float _svFitTauPosPhiStandalone;
	float _svFitTauPosMassStandalone;
	float _svFitTauPosEnergyStandalone;

	float _svFitTauNegPxStandalone;
	float _svFitTauNegPyStandalone;
	float _svFitTauNegPzStandalone;
	float _svFitTauNegPtStandalone;
	float _svFitTauNegEtaStandalone;
	float _svFitTauNegPhiStandalone;
	float _svFitTauNegMassStandalone;
	float _svFitTauNegEnergyStandalone;
};

#endif 
