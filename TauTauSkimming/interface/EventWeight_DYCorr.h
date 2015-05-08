#ifndef EventWeight_DYCorr_h
#define EventWeight_DYCorr_h

#include <iostream>
#include <map>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "TH2F.h"

#include "DesyHTauTau/GenericTools/interface/StringTools.h"

class EventWeight_DYCorr
{
private:
	TH1F * bdtBins;
	std::vector<TH2F*> calibInclHistos;
	std::vector<std::string> histoNames;
	std::map<std::string, TH2F*> calibHistos;

public:
	EventWeight_DYCorr(std::string inputFilename = "", std::string inputDirectoryName = "Run2011_HCP12")
	{
		TFile * inputFile = new TFile(inputFilename.c_str());
		assert(inputFile->IsOpen());
		assert(!inputFile->IsZombie());

		if (inputDirectoryName == ".") {
			TList* histogramKeysList = inputFile->GetListOfKeys();
			
			histoNames.clear();
			int nCategoryNumbers = 10;
			for(int categoryNumber = 0; categoryNumber < nCategoryNumbers; ++categoryNumber) histoNames.push_back("");
			
			for(int i = 0; i < histogramKeysList->GetSize(); ++i) {
				std::string histogramName = histogramKeysList->At(i)->GetName();
				
				for(int categoryNumber = 0; categoryNumber < nCategoryNumbers; ++categoryNumber) {
					if(StringTools::contains(histogramName, boost::lexical_cast<string>(categoryNumber))) histoNames[categoryNumber] = histogramName;
				}
			}
		}
		else if (inputDirectoryName == "Zmm_SUSY13_2011") {
			histoNames.push_back("InclLowPt");
			histoNames.push_back("InclHighPt");
			histoNames.push_back("NoBTag");
			histoNames.push_back("BTag");
		}
		else if (inputDirectoryName == "Zmm_SUSY13_2012") {
			histoNames.push_back("InclLowPt");
			histoNames.push_back("InclHighPt");
			histoNames.push_back("NoBTag");
			histoNames.push_back("BTag");
		}
		else if (inputDirectoryName == "TTj_SUSY13_2011") {
		  	histoNames.push_back("NoBTag");
			histoNames.push_back("BTag");
		}
		else if (inputDirectoryName == "TTj_SUSY13_2012") {
		  	histoNames.push_back("NoBTag");
			histoNames.push_back("BTag");
		}
		else if (inputDirectoryName == "Run2012_HCP12") {
			histoNames.push_back("InclLowPt");
			histoNames.push_back("InclHighPt");
			histoNames.push_back("LowPt_0j");
			histoNames.push_back("HighPt_0j");
			histoNames.push_back("LowPt_1j");
			histoNames.push_back("HighPt_1j");
			histoNames.push_back("Vbf");
			histoNames.push_back("LowPt_btag");
			histoNames.push_back("HighPt_btag");
		}
		else if (inputDirectoryName == "Run2011_HCP12") {
			histoNames.push_back("InclLowPt");
			histoNames.push_back("InclHighPt");
			histoNames.push_back("LowPt_0j");
			histoNames.push_back("HighPt_0j");
			histoNames.push_back("LowPt_1j");
			histoNames.push_back("HighPt_1j");
			histoNames.push_back("Vbf");
			histoNames.push_back("LowPt_btag");
			histoNames.push_back("HighPt_btag");
		}
		else if (inputDirectoryName == "Run2012_Moriond13") {
			histoNames.push_back("InclLowPt");
			histoNames.push_back("InclHighPt");
			histoNames.push_back("LowPt_0j");
			histoNames.push_back("HighPt_0j");
			histoNames.push_back("LowPt_1j");
			histoNames.push_back("HighPt_1j");
			histoNames.push_back("Vbf");
			histoNames.push_back("LowPt_btag");
			histoNames.push_back("HighPt_btag");
		}
		else if (inputDirectoryName == "HCP12_EE") {
			histoNames.push_back("DYcorr0Jet");
			histoNames.push_back("DYcorr1Jet");
			histoNames.push_back("DYcorrVbf");
		}
		else if (inputDirectoryName == "Moriond13_EE") {
			histoNames.push_back("DYcorr0Jet");
			histoNames.push_back("DYcorr1Jet");
			histoNames.push_back("DYcorrVbf");
		}
		else if (inputDirectoryName == "DY_EE_Summer13") {
			histoNames.push_back("ZeroJetDYCor");
			histoNames.push_back("OneJetDYCor");
			histoNames.push_back("VbfDYCor");
		}
		else if (inputDirectoryName == "DY_MuMu_Summer13") {
			histoNames.push_back("ZeroJetDYCor");
			histoNames.push_back("OneJetDYCor");
			histoNames.push_back("VbfDYCor");
		}
		else if (inputDirectoryName == "Summer13MSSM_EE") {
		  	histoNames.push_back("BTagCatDYCor");
			histoNames.push_back("NoBTagCatDYCor");
		
		}
		else 
		{
			calibInclHistos.push_back((TH2F*)inputFile->Get((inputDirectoryName+"/calib0").c_str()));
			calibInclHistos.push_back((TH2F*)inputFile->Get((inputDirectoryName+"/calib1").c_str()));
			calibInclHistos.push_back((TH2F*)inputFile->Get((inputDirectoryName+"/calib2").c_str()));
			bdtBins = (TH1F*)inputFile->Get((inputDirectoryName+"/bdtBins").c_str());
		}

		std::cout << inputFilename << std::endl;
		for (std::vector<std::string>::const_iterator itName = histoNames.begin(); itName != histoNames.end(); ++itName)
		{
			if((calibHistos.count(*itName) < 1) && ((*itName) != "")) {
				
				std::cout << *itName << ": ";
				if((inputDirectoryName == ".") || (inputDirectoryName == "Summer13MSSM_EE")) {
					calibHistos[*itName] = (TH2F*)inputFile->Get((*itName).c_str());
				}
				else {
					calibHistos[*itName] = (TH2F*)inputFile->Get((inputDirectoryName+"/"+*itName).c_str());
				}
				std::cout << calibHistos[*itName] << std::endl;

				int nBinsX = calibHistos[*itName]->GetNbinsX();
				int nBinsY = calibHistos[*itName]->GetNbinsY();
				std::cout << "NBinsX = " << nBinsX << " NBinsY = " << nBinsY << std::endl;
				
				for (int iBX = 0; iBX<nBinsX; ++iBX) {
					for (int iBY = 0; iBY<nBinsY; ++iBY) {
						std::cout << "[" << iBX << "," << iBY << "]" << "  :  " << calibHistos[*itName]->GetBinContent(iBX+1,iBY+1) << std::endl;
					}
				}
			}
		}
		//		exit(-1);
	}
	
	float getWeight3D(float bdt, float mvis, float mtt) {
		int bin = bdtBins->FindBin(bdt) - 1;
		if (bin<0) bin=0;
		if (bin>2) bin=2;
		TH2F * hist = calibInclHistos.at(bin);
		float weight = hist->GetBinContent(hist->FindBin(mvis,mtt));
		return weight;
	}

	float getWeight(int category, float firstVar, float secondVar, float systShift = 0.0, bool useBinErrorForSystShift = true)
	{
		//assert(category >= 0);
		//assert(category < (int)histoNames.size());
		if((category < 0) || (category > (int)(histoNames.size())) || calibHistos.count(histoNames[category]) < 1) return 1.0;

		int bin = calibHistos[histoNames[category]]->FindBin(firstVar, secondVar);
		
		float weight = calibHistos[histoNames[category]]->GetBinContent(bin);
		
		if(systShift != 0.0) {
			if(useBinErrorForSystShift) { // take shift from bin error
				weight += (systShift * calibHistos[histoNames[category]]->GetBinError(bin));
			}
			else { // shift down: uncorrected; shift up: double shift from correction
				if(systShift < 0.0) weight = 1.0;
				else weight += (weight - 1.0);
			}
		}
		
		if(weight < 0.0) return 0.0;
		else return weight;
	}
	float getWeightErr(int category, float firstVar, float secondVar)
	{
		assert(category >= 0);
		assert(category < (int)histoNames.size());

		int bin = calibHistos[histoNames[category]]->FindBin(firstVar, secondVar);
		
		float weightErr = calibHistos[histoNames[category]]->GetBinError(bin);
		return weightErr;
	}
};

#endif
