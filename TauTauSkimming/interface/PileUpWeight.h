#ifndef DesyHTauTau_TauTauSkimming_PileUpWeight
#define DesyHTauTau_TauTauSkimming_PileUpWeight

#include <string>
#include <TH1.h>

class PileUpWeight
{
public:
	explicit PileUpWeight(std::string const weightFileName);
	explicit PileUpWeight(std::string const weightFileNameData, std::string const weightFileNameMC);
	explicit PileUpWeight(std::string const weightFileName, std::string const weightHistogramNameData, std::string const weightHistogramNameMC);
	explicit PileUpWeight(std::string const weightFileNameData, std::string const weightFileNameMC, std::string const weightHistogramNameData, std::string const weightHistogramNameMC);
	~PileUpWeight();
	
	float GetWeight(float const nPU) const;
private:
	TH1* _h_pu_data;
	TH1* _h_pu_mc;
};

#endif