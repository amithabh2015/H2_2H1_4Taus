#include <H2to2H1to4Taus/TauTauSkimming/interface/PileUpWeight.h>

#include <stdexcept>
#include <string>
#include <iostream>

#include <TObject.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>



PileUpWeight::PileUpWeight(std::string const weightFileName)
{
	using namespace std;
	
	static const char* src_h_data_name = "scaleFactors/dataMc/pileUp/data";
	static const char* src_h_mc_name   = "scaleFactors/dataMc/pileUp/mc";
	
	TFile weightFile(weightFileName.c_str(), "READ");
	
	TH1* src_h_data = static_cast<TH1*>(weightFile.Get(src_h_data_name));
	TH1* src_h_mc   = static_cast<TH1*>(weightFile.Get(src_h_mc_name));
	
	if (!src_h_data)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileName + ":" + src_h_data_name).c_str());
		throw ex;
	}
	if (!src_h_mc)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileName + ":" + src_h_mc_name).c_str());
		throw ex;
	}
	
	_h_pu_data = static_cast<TH1*>(src_h_data->Clone());
	_h_pu_mc   = static_cast<TH1*>(src_h_mc->Clone());
	
	_h_pu_data->SetDirectory(0);
	_h_pu_mc->SetDirectory(0);
	
	weightFile.Close();
	
	if (_h_pu_data->GetNbinsX() != _h_pu_mc->GetNbinsX())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different number of bins " << endl;
	
	if (_h_pu_data->GetXaxis()->GetXmax() != _h_pu_mc->GetXaxis()->GetXmax() || _h_pu_data->GetXaxis()->GetXmin() != _h_pu_mc->GetXaxis()->GetXmin())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different ranges" << endl;
	
	_h_pu_data->Sumw2();
	_h_pu_mc  ->Sumw2();
	
	_h_pu_data->Scale(1.0/_h_pu_data->Integral());
	_h_pu_mc  ->Scale(1.0/_h_pu_mc->Integral());
}

PileUpWeight::PileUpWeight(std::string const weightFileNameData, std::string const weightFileNameMC)
{
	using namespace std;
	
	static const char* src_h_data_name = "scaleFactors/dataMc/pileUp/data";
	static const char* src_h_mc_name   = "scaleFactors/dataMc/pileUp/mc";
	
	
	TFile weightFileData(weightFileNameData.c_str(), "READ");
	TH1* src_h_data = static_cast<TH1*>(weightFileData.Get(src_h_data_name));
	
	if (!src_h_data)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileNameData + ":" + src_h_data_name).c_str());
		throw ex;
	}
	
	_h_pu_data = static_cast<TH1*>(src_h_data->Clone());
	_h_pu_data->SetDirectory(0);
	
	weightFileData.Close();
	
	
	
	TFile weightFileMC(weightFileNameMC.c_str(), "READ");
	TH1* src_h_mc   = static_cast<TH1*>(weightFileMC.Get(src_h_mc_name));
	
	if (!src_h_mc)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileNameMC + ":" + src_h_mc_name).c_str());
		throw ex;
	}
	
	_h_pu_mc   = static_cast<TH1*>(src_h_mc->Clone());
	_h_pu_mc  ->SetDirectory(0);
	
	weightFileMC.Close();
	
	
	
	if (_h_pu_data->GetNbinsX() != _h_pu_mc->GetNbinsX())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different number of bins " << endl;
	
	if (_h_pu_data->GetXaxis()->GetXmax() != _h_pu_mc->GetXaxis()->GetXmax() || _h_pu_data->GetXaxis()->GetXmin() != _h_pu_mc->GetXaxis()->GetXmin())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different ranges" << endl;
	
	_h_pu_data->Sumw2();
	_h_pu_mc  ->Sumw2();
	
	_h_pu_data->Scale(1.0/_h_pu_data->Integral());
	_h_pu_mc  ->Scale(1.0/_h_pu_mc->Integral());
}



PileUpWeight::PileUpWeight(std::string const weightFileName, std::string const weightHistogramNameData, std::string const weightHistogramNameMC)
{
	using namespace std;
	
	TFile weightFile(weightFileName.c_str(), "READ");
	
	TH1* src_h_data = static_cast<TH1*>(weightFile.Get(weightHistogramNameData.c_str()));
	TH1* src_h_mc   = static_cast<TH1*>(weightFile.Get(weightHistogramNameMC.c_str()));
	
	if (!src_h_data)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileName + ":" + weightHistogramNameData).c_str());
		throw ex;
	}
	if (!src_h_mc)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileName + ":" + weightHistogramNameMC).c_str());
		throw ex;
	}
	
	_h_pu_data = static_cast<TH1*>(src_h_data->Clone());
	_h_pu_mc   = static_cast<TH1*>(src_h_mc->Clone());
	
	_h_pu_data->SetDirectory(0);
	_h_pu_mc->SetDirectory(0);
	
	weightFile.Close();
	
	if (_h_pu_data->GetNbinsX() != _h_pu_mc->GetNbinsX())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different number of bins " << endl;
	
	if (_h_pu_data->GetXaxis()->GetXmax() != _h_pu_mc->GetXaxis()->GetXmax() || _h_pu_data->GetXaxis()->GetXmin() != _h_pu_mc->GetXaxis()->GetXmin())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different ranges" << endl;
	
	_h_pu_data->Sumw2();
	_h_pu_mc  ->Sumw2();
	
	_h_pu_data->Scale(1.0/_h_pu_data->Integral());
	_h_pu_mc  ->Scale(1.0/_h_pu_mc->Integral());
}

PileUpWeight::PileUpWeight(std::string const weightFileNameData, std::string const weightFileNameMC, std::string const weightHistogramNameData, std::string const weightHistogramNameMC)
{
	using namespace std;
	
	TFile weightFileData(weightFileNameData.c_str(), "READ");
	TH1* src_h_data = static_cast<TH1*>(weightFileData.Get(weightHistogramNameData.c_str()));
	
	if (!src_h_data)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileNameData + ":" + weightHistogramNameData).c_str());
		throw ex;
	}
	
	_h_pu_data = static_cast<TH1*>(src_h_data->Clone());
	_h_pu_data->SetDirectory(0);
	
	weightFileData.Close();
	
	
	
	TFile weightFileMC(weightFileNameMC.c_str(), "READ");
	TH1* src_h_mc   = static_cast<TH1*>(weightFileMC.Get(weightHistogramNameMC.c_str()));
	
	if (!src_h_mc)
	{
		std::runtime_error ex((string("Cannot get ") + weightFileNameMC + ":" + weightHistogramNameMC).c_str());
		throw ex;
	}
	
	_h_pu_mc   = static_cast<TH1*>(src_h_mc->Clone());
	_h_pu_mc  ->SetDirectory(0);
	
	weightFileMC.Close();
	
	
	
	if (_h_pu_data->GetNbinsX() != _h_pu_mc->GetNbinsX())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different number of bins " << endl;
	
	if (_h_pu_data->GetXaxis()->GetXmax() != _h_pu_mc->GetXaxis()->GetXmax() || _h_pu_data->GetXaxis()->GetXmin() != _h_pu_mc->GetXaxis()->GetXmin())
		cerr << "WARNING in PileUpWeight constructor: data and MC histograms have different ranges" << endl;
	
	_h_pu_data->Sumw2();
	_h_pu_mc  ->Sumw2();
	
	_h_pu_data->Scale(1.0/_h_pu_data->Integral());
	_h_pu_mc  ->Scale(1.0/_h_pu_mc->Integral());
}



PileUpWeight::~PileUpWeight()
{
	delete _h_pu_data;
	delete _h_pu_mc;
}


float PileUpWeight::GetWeight(float const nPU) const
{
	using namespace std;
	
	int const nBinData = _h_pu_data->FindBin(nPU);
	int const nBinMC   = _h_pu_mc->FindBin(nPU);
	
	if ( _h_pu_data->IsBinOverflow(nBinData))
	{
		cerr << "WARNING in PileUpWeight::GetWeight(nPU) :  nPU=" << nPU <<" overflows the data histogram" << endl;
		return 0.0;
	}
	else if (_h_pu_data->IsBinUnderflow(nBinData))
	{
		cerr << "WARNING in PileUpWeight::GetWeight(nPU) :  nPU=" << nPU <<" underflows the data histogram" << endl;
		return 0.0;
	}
	
	if (_h_pu_mc->IsBinOverflow(nBinMC))
	{
		cerr << "WARNING in PileUpWeight::GetWeight(nPU) :  nPU=" << nPU <<" overflows the MC histogram" << endl;
		return 0.0;
	}
	else if (_h_pu_mc->IsBinUnderflow(nBinMC))
	{
		cerr << "WARNING in PileUpWeight::GetWeight(nPU) :  nPU=" << nPU <<" underflows the MC histogram" << endl;
		return 0.0;
	}
	
	float const evtDensityData = _h_pu_data->GetBinContent(nBinData) / _h_pu_data->GetBinWidth(nBinData);
	float const evtDensityMC   = _h_pu_mc->GetBinContent(nBinMC)     / _h_pu_mc->GetBinWidth(nBinMC);
	
	if (evtDensityMC == 0.0)
	{
		cerr << "WARNING in PileUpWeight::GetWeight(nPU) :  MC histogram iz zero-valued for nPU=" << nPU << endl;
		return 0.0;
	}
	
	return evtDensityData/evtDensityMC;
}
