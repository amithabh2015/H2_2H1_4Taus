#ifndef DesyHTauTau_TauTauSkimming_HTauTauLeptonPUIsolationInfo
#define DesyHTauTau_TauTauSkimming_HTauTauLeptonPUIsolationInfo

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
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

#include <vector>
#include <iostream>

#include "DesyHTauTau/TauTauSkimming/interface/HTauTauTypeDef.h"
#include "DesyHTauTau/TauTauSkimming/interface/HTauTauGetElectrons.h"
#include "DesyHTauTau/TauTauSkimming/interface/HTauTauGetMuons.h"

// #include "DesyHTauTau/TauTauSkimming/interface/HTauTauElectronMuonSelector.h"
// #include "DesyHTauTau/TauTauSkimming/interface/HTauTauElectronSelector.h"
// #include "DesyHTauTau/TauTauSkimming/interface/HTauTauMuonSelector.h"



class HTauTauLeptonPUIsolationInfoBase;
template <class T1, class T2> class HTauTauLeptonPUIsolationInfo;
template <> class HTauTauLeptonPUIsolationInfo <pat::MuonRef, pat::MuonRef>;
template <> class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::MuonRef>;
template <> class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::ElectronRef>;


///**************************************
///* Base class for Lepton-PU isolation *
///**************************************
class HTauTauLeptonPUIsolationInfoBase
{
  friend class HTauTauLeptonPUIsolationInfo <pat::MuonRef, pat::MuonRef>;
  friend class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::MuonRef>;
  friend class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::ElectronRef>;
  
public:
	virtual ~HTauTauLeptonPUIsolationInfoBase() { };
	
private:
	
	HTauTauLeptonPUIsolationInfoBase(const edm::ParameterSet&);
	
	void CalcPUPfCandidateVector(const edm::Event&, const edm::EventSetup&);
	void ReportPtThresholdVectorIsTooLong();
	
	// member data
	TTree* _tree;
	
	edm::InputTag		_pfPileUpSrc;
	puPFCandidateVector	_pileUpChargedPFs;
	
	bool	_sameSignLeptons;
	
	double		_etaMax;
	double		_ptMin;
	double		_ptMinHard;
	double		_deltaR;
	
	int		_nPtBinsPFIsolation;
	double		_ptThresholds[30];
};


///***************************************
///* Derived class for mumu-PU isolation *
///***************************************
template <> class HTauTauLeptonPUIsolationInfo <pat::MuonRef, pat::MuonRef> : public HTauTauLeptonPUIsolationInfoBase
{
public:
	HTauTauLeptonPUIsolationInfo(const edm::ParameterSet&);
	~HTauTauLeptonPUIsolationInfo() { };
	
	void beginJob(TTree*, edm::Service<TFileService>&);
	void analyze(pat::MuonRef&, pat::MuonRef&, const edm::Event&, const edm::EventSetup&);
	void endJob();
	
private:
	
	edm::InputTag 	_muonSrc;
	
	float _isoPUChargedHadPFsNeg[30];
	float _isoPUChargedHadPFsPos[30];
};


///*************************************
///* Derived class for eMu-PU isolation *
///*************************************
template <> class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::MuonRef> : public HTauTauLeptonPUIsolationInfoBase
{
public:
	HTauTauLeptonPUIsolationInfo(const edm::ParameterSet&);
	~HTauTauLeptonPUIsolationInfo() { };
	
	void beginJob(TTree*, edm::Service<TFileService>&);
	void analyze(pat::ElectronRef&, pat::MuonRef&, const edm::Event&, const edm::EventSetup&);
	void endJob();
	
private:
	edm::InputTag	_muonSrc;
	edm::InputTag	_electronSrc;
	
	float _isoPUChargedHadPFsElec[30];
	float _isoPUChargedHadPFsMuon[30];
	
};


///*************************************
///* Derived class for ee-PU isolation *
///*************************************
template <> class HTauTauLeptonPUIsolationInfo <pat::ElectronRef, pat::ElectronRef> : public HTauTauLeptonPUIsolationInfoBase
{
public:
	HTauTauLeptonPUIsolationInfo(const edm::ParameterSet&);
	~HTauTauLeptonPUIsolationInfo() { };
	
	void beginJob(TTree*, edm::Service<TFileService>&);
	void analyze(pat::ElectronRef&, pat::ElectronRef&, const edm::Event&, const edm::EventSetup&);
	void endJob();
	
private:
	edm::InputTag	_electronSrc;
	
	float _isoPUChargedHadPFsElecPos[30];
	float _isoPUChargedHadPFsElecNeg[30];
};


#endif
