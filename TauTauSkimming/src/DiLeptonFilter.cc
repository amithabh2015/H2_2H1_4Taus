// -*- C++ -*-
//
// Package:    DiLeptonFilter
// Class:      DiLeptonFilter
// 
/**\class DiLeptonFilter DiLeptonFilter.cc Test/DiLeptonFilter/src/DiLeptonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Manuel Zeise, 8/22, 7243
//         Created:  Fri Apr 27 11:18:52 CEST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Common/interface/View.h"
//
// class declaration
//

class DiLeptonFilter : public edm::EDFilter
{
public:
	explicit DiLeptonFilter(const edm::ParameterSet&);
	~DiLeptonFilter();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob() ;
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;

	virtual bool beginRun(edm::Run&, edm::EventSetup const&);
	virtual bool endRun(edm::Run&, edm::EventSetup const&);
	virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

	// ----------member data ---------------------------
	bool requireDiMuon;
	bool requireDiElectron;
	bool requireMuonElectron;

	edm::InputTag srcMuons;
	edm::InputTag srcElectrons;

	float minPtMuon1;
	float minPtMuon2;
	float maxEtaMuon;

	float minPtElectron1;
	float minPtElectron2;
	float maxEtaElectron;
};

//
// constructors and destructor
//
DiLeptonFilter::DiLeptonFilter(const edm::ParameterSet& iConfig):
	requireDiMuon(false),
	requireDiElectron(false),
	requireMuonElectron(false),
	minPtMuon1(0.0),
	minPtMuon2(0.0),
	maxEtaMuon(100.0),
	minPtElectron1(0.0),
	minPtElectron2(0.0),
	maxEtaElectron(100.0)
{
	//now do what ever initialization is needed
	srcMuons = iConfig.getParameter<edm::InputTag>("srcMuons");
	srcElectrons = iConfig.getParameter<edm::InputTag>("srcElectrons");

	if (iConfig.exists("minPtMuon"))
	{
		minPtMuon1 = iConfig.getParameter<double>("minPtMuon");
		minPtMuon2 = iConfig.getParameter<double>("minPtMuon");
	}
	else
	{
		minPtMuon1 = iConfig.getParameter<double>("minPtMuon1");
		minPtMuon2 = iConfig.getParameter<double>("minPtMuon2");
	}
	maxEtaMuon = iConfig.getParameter<double>("maxEtaMuon");

	if (iConfig.exists("minPtElectron"))
	{
		minPtElectron1 = iConfig.getParameter<double>("minPtElectron");
		minPtElectron2 = iConfig.getParameter<double>("minPtElectron");
	}
	else
	{
		minPtElectron1 = iConfig.getParameter<double>("minPtElectron1");
		minPtElectron2 = iConfig.getParameter<double>("minPtElectron2");
	}
	maxEtaElectron = iConfig.getParameter<double>("maxEtaElectron");

	requireDiMuon = iConfig.getParameter<bool>("requireDiMuon");
	requireDiElectron = iConfig.getParameter<bool>("requireDiElectron");
	requireMuonElectron = iConfig.getParameter<bool>("requireMuonElectron");

	assert(minPtMuon1 >= minPtMuon2);
	assert(minPtElectron1 >= minPtElectron2);
}


DiLeptonFilter::~DiLeptonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	bool diMuonResult = false;
	bool diElectronResult = false;
	bool muonElectronResult = false;

	unsigned int cntMuon1 = 0;	// tighter pt cut
	unsigned int cntMuon2 = 0;
	unsigned int cntElectron1 = 0;	// tighter pt cut
	unsigned int cntElectron2 = 0;

	if (requireDiMuon||requireMuonElectron)
	{
		edm::Handle< edm::View<reco::Candidate> > muons;

		if (iEvent.getByLabel(srcMuons, muons))
		{
			for (edm::View<reco::Candidate>::const_iterator it = muons->begin(); it != muons->end(); ++it)
			{
				if (it->pt() > minPtMuon1 && std::abs(it->eta())<maxEtaMuon && it->isGlobalMuon() && it->isTrackerMuon())
					++cntMuon1;
				if (it->pt() > minPtMuon2 && std::abs(it->eta())<maxEtaMuon && it->isGlobalMuon() && it->isTrackerMuon())
					++cntMuon2;

			}
		}
		else
			std::cout << srcMuons << " not found" << std::endl;
	}

	if (requireDiElectron||requireMuonElectron)
	{
		edm::Handle< edm::View<reco::Candidate> > electrons;

		if (iEvent.getByLabel(srcElectrons, electrons))
		{
			for (edm::View<reco::Candidate>::const_iterator it = electrons->begin(); it != electrons->end(); ++it)
			{
				if (it->pt() > minPtElectron1 && std::abs(it->eta())<maxEtaElectron)
					++cntElectron1;
				if (it->pt() > minPtElectron2 && std::abs(it->eta())<maxEtaElectron)
					++cntElectron2;
			}
		}
		else
			std::cout << srcElectrons << " not found" << std::endl;
	}

	if (requireDiMuon && cntMuon1 >= 1 && cntMuon2 >= 2)
		diMuonResult = true;

	if (requireDiElectron && cntElectron1 >=1 && cntElectron2 >= 2)
		diElectronResult = true;

	if (requireMuonElectron && cntMuon1 >= 1 && cntElectron1 >= 1)
		muonElectronResult = true;

	return diMuonResult || diElectronResult || muonElectronResult;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DiLeptonFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiLeptonFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
DiLeptonFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
DiLeptonFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
DiLeptonFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
DiLeptonFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiLeptonFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiLeptonFilter);
