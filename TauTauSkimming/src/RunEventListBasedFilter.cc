// -*- C++ -*-
//
// Package:    RunEventListBasedFilter
// Class:      RunEventListBasedFilter
// 
/**\class RunEventListBasedFilter RunEventListBasedFilter.cc Test/RunEventListBasedFilter/src/RunEventListBasedFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luigi Calligaris,,,DESY
//         Created:  Fri Jun 14 16:12:13 CEST 2013
// $Id$
//
//


#include <memory>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <functional>
#include <boost/functional/hash.hpp>
#include <boost/tokenizer.hpp>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace desyhtautau
{
struct RunEventId
{
	RunEventId() : run(0), event(0) {}
	RunEventId(unsigned int const setRun, unsigned int const setEvent) : run(setRun), event(setEvent) {}
	RunEventId(RunEventId const& rhs) : run(rhs.run), event(rhs.event) {}
	~RunEventId() {}
	RunEventId& operator=(RunEventId const& rhs) {run = rhs.run; event = rhs.event; return *this;}
	RunEventId& operator=(RunEventId&& rhs) {run = rhs.run; event = rhs.event; return *this;}
	
	bool operator==(RunEventId const& rhs) const
	{
		return (run == rhs.run) && (event == rhs.event);
	}
	
	bool operator<(RunEventId const& rhs) const
	{
		return (run < rhs.run) || (run == rhs.run && event < rhs.event);
	}
	
	bool operator>(RunEventId const& rhs) const
	{
		return (run > rhs.run) || (run == rhs.run && event > rhs.event);
	}
	
	unsigned int run;
	unsigned int event;
};
}

namespace std
{
	template <> struct hash<desyhtautau::RunEventId>
	{
		size_t operator()(const desyhtautau::RunEventId& rhs) const
		{
			std::size_t seed = 0;
			boost::hash_combine(seed, rhs.run);
			boost::hash_combine(seed, rhs.event);
			return seed;
		}
	};
}



inline std::vector< desyhtautau::RunEventId > GetRunEventIdVectorFromEventList(std::string const eventListFileName)
{
	using namespace std;
	
	std::vector< desyhtautau::RunEventId > outVec;
	outVec.reserve(70000);
	
	edm::FileInPath fp(eventListFileName);
	string const eventListFullPath = fp.fullPath();
	
	ifstream eventList(eventListFullPath.c_str());
	
	while (!eventList.eof())
	{
		char cLine[256];
		eventList.getline(cLine, 256);
		string line(cLine);
		
		// check if line is empty or is a comment
		if (line.size() == 0 || line[0]=='#')
			continue;
		
		// convert everything else but numbers to whitespace
		for (char & c : line)
		{
			switch (c)
			{
				case '0': break;
				case '1': break;
				case '2': break;
				case '3': break;
				case '4': break;
				case '5': break;
				case '6': break;
				case '7': break;
				case '8': break;
				case '9': break;
				default: c = ' '; break;
			}
		}
		
		boost::char_separator<char> sep(" ");
		boost::tokenizer< boost::char_separator<char> > tokens(line, sep);
		
		unsigned int run   = 0;
		unsigned int event = 0;
		
		int fields = 0;
		for (const string tok : tokens)
		{
			if (tok == "" || tok == " ")
				continue;
			
			++fields;
			
			if (fields == 1) {run   = atoi(tok.c_str());} // first numeric field is the run
			if (fields == 2) {event = atoi(tok.c_str()); break;} // second numeric field is the event
		}
		
		if (fields == 2)
			outVec.push_back( desyhtautau::RunEventId(run, event) );
	}
	
	return outVec;
}








class RunEventListBasedFilter : public edm::EDFilter
{
public:
	explicit RunEventListBasedFilter(const edm::ParameterSet&);
	~RunEventListBasedFilter();
	
private:
	virtual void beginJob() ;
	virtual bool beginRun(edm::Run&, edm::EventSetup const&);
	virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual bool filter(edm::Event&, const edm::EventSetup&);
	virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
	virtual bool endRun(edm::Run&, edm::EventSetup const&);
	virtual void endJob() ;
	
	std::string                 _eventsToAcceptListFile;
	std::unordered_set<desyhtautau::RunEventId> _eventsToAccept;
	int _verbosity;
};


RunEventListBasedFilter::RunEventListBasedFilter(const edm::ParameterSet& iConfig) :
	_eventsToAcceptListFile(iConfig.getParameter<std::string> ("EventsToAcceptListFile")),
	_verbosity(iConfig.getParameter<int> ("Verbosity"))
{
	std::vector< desyhtautau::RunEventId > eventsToBeAcceptedVector = GetRunEventIdVectorFromEventList(_eventsToAcceptListFile);
	
	std::cout << "- HTauTauElecMuonTTreeProducer will accept a maximum total number of events: " << eventsToBeAcceptedVector.size() << std::endl;
	std::cout << "-- listed in file: " << _eventsToAcceptListFile << std::endl;
	
	if (_verbosity)
	{
		std::cout << "- Will accept the following events, shown as (run,event) pairs:" << std::endl;
		for (desyhtautau::RunEventId eid : eventsToBeAcceptedVector)
			std::cout << "(" << eid.run << "," << eid.event << "),";
		std::cout << std::endl;
	}
	
	for (desyhtautau::RunEventId eid : eventsToBeAcceptedVector)
	{
		_eventsToAccept.insert(std::move(eid));
	}
}


RunEventListBasedFilter::~RunEventListBasedFilter()
{
	
}


bool RunEventListBasedFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	if (_eventsToAccept.find( desyhtautau::RunEventId(iEvent.id().run(), iEvent.id().event()) ) != _eventsToAccept.end())
		return true;
	
	return false;
}

void RunEventListBasedFilter::beginJob()
{
}

void RunEventListBasedFilter::endJob()
{
}

bool RunEventListBasedFilter::beginRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

bool RunEventListBasedFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

bool RunEventListBasedFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

bool RunEventListBasedFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

DEFINE_FWK_MODULE(RunEventListBasedFilter);
