#ifndef DesyHTauTau_TauTauSkimming_ElectronCutBasedIDGetter
#define DesyHTauTau_TauTauSkimming_ElectronCutBasedIDGetter
#endif

#include <vector>
#include <string>

#include "DataFormats/PatCandidates/interface/Electron.h"


class ElectronCutBasedIDGetter
{
public:
  ElectronCutBasedIDGetter();
  ElectronCutBasedIDGetter(std::vector<std::string> const& customElectronIDs);
  ~ElectronCutBasedIDGetter();
	
	std::vector<std::string> GetIdNames();
	std::vector<float>       GetIdValues(pat::ElectronRef const electron);
	float                    GetId(std::string idName, pat::ElectronRef const electron);
	
private:
	void CheckElectronIsValid(pat::ElectronRef const electron);
	
	void ProcessElectron(pat::ElectronRef const electron);
	
  std::vector<std::string> _electronIdNames;
  std::vector<float>       _electronIdValues;
};
