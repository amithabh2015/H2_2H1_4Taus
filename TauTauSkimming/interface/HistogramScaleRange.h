#include <math.h>
#include <string.h>
#include <sstream>
#include "TH1.h"
#include "TDirectory.h"

#ifndef HISTOGRAMSCALERANGE_H
#define HISTOGRAMSCALERANGE_H

/*
This class manages a number of histograms all storing data from the same source but multiplied with a different factor.
It does so by having one array for storing the histogramms and factors each, passing methods to the histograms and scaling input data with the factors

Note for further additions:
all methods working with standard histogram methods should have the same name in order to provide a uniform interface
*/

class TH1FScaleRange {
	public:
	// constructors
		TH1FScaleRange(TDirectory * outputDir, std::string listname, float min, float max);
		TH1FScaleRange(TDirectory * outputDir, std::string listname, float min, float max, float stepsize);
		TH1FScaleRange(TDirectory * outputDir, std::string listname, std::string titlename, float min, float max, float stepsize);
		TH1FScaleRange(TDirectory * outputDir, std::string listname, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax);
		TH1FScaleRange(TDirectory * outputDir, std::string listname, std::string listsuffix, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax);
		~TH1FScaleRange();

	// histogram manipulation methods - use the same interface as primary histogram methods
		void Fill(float value, float weight);
		void Fill(float value);
		void SetDirectory(TDirectory* dir);

	private:
		void initialize(TDirectory * outputDir, std::string listname, std::string listsuffix, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax);
	// data content objects
		TH1F** histostore;
		float* histofactor;
		int histocount;
		float maxbin;
/*
TH1F * metBeforeCutPlus01H;
metBeforeCutPlus01H = new TH1F("metBeforeCutPlus01H","",100,0.,200.);
metBeforeCutPlus01H->Fill(TMath::Min(float(1.01*_met),float(199.5)),weight);
*/
};

// method called by all constructors
void TH1FScaleRange::initialize(TDirectory * outputDir, std::string listname, std::string listsuffix, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax){
	TDirectory * tmpDir = outputDir->mkdir(listname.c_str());
	tmpDir->cd();
	histocount = std::floor((max - min)/stepsize); // this is actually number of histograms-1
	histofactor = new float[histocount+1];
	histostore = new TH1F*[histocount+1];
	maxbin = Xbinmax - (Xbinmax-Xbinmin)/(nbins*2);
	for(int i = 0; i <= histocount; i++){
		float factor = min + stepsize * i;
		std::stringstream histoname;
		if(factor < 1 )
			histoname << listname << "Minus" << setprecision(2) <<  setfill('0') << setw(2) << std::abs((factor-1)*100) << listsuffix;
		else if(factor > 1 )
			histoname << listname << "Plus" << setprecision(2) << setfill('0') << setw(2) << ((factor-1)*100) << listsuffix;
		else
			histoname << listname << listsuffix;
		//std::string displayname = titlename + " ( x " + std::to_string(factor) + ")";
		TString displayname = titlename + " ( x ";
		displayname += factor;
		displayname += ")";

		histofactor[i] = factor;
		//histostore[i] = new TH1F(histoname.str().c_str(),displayname.c_str(),100,0.,200.);
		histostore[i] = new TH1F(histoname.str().c_str(),displayname,100,0.,200.);
	}
}

// constructors
TH1FScaleRange::TH1FScaleRange (TDirectory * outputDir, std::string listname, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax) {
	initialize(outputDir, listname,"",titlename,min,max,stepsize,nbins,Xbinmin,Xbinmax);
}
TH1FScaleRange::TH1FScaleRange (TDirectory * outputDir, std::string listname, std::string titlename, float min, float max, float stepsize) {
	initialize(outputDir, listname,"",titlename,min,max,stepsize,100,0.,200.);
}
TH1FScaleRange::TH1FScaleRange (TDirectory * outputDir, std::string listname, float min, float max, float stepsize) {
	initialize(outputDir, listname,"",listname,min,max,stepsize,100,0.,200.);
}
TH1FScaleRange::TH1FScaleRange (TDirectory * outputDir, std::string listname, float min, float max) {
	initialize(outputDir, listname,"",listname,min,max,0.1,100,0.,200.);
}
TH1FScaleRange::TH1FScaleRange (TDirectory * outputDir, std::string listname, std::string listsuffix, std::string titlename, float min, float max, float stepsize, int nbins, float Xbinmin, float Xbinmax) {
	initialize(outputDir, listname,listsuffix,titlename,min,max,stepsize,nbins,Xbinmin,Xbinmax);
}

// histogram methods
void TH1FScaleRange::Fill(float value, float weight){
	for(int i = 0; i <= histocount; i++){
		histostore[i]->Fill(std::min(float(histofactor[i] * value),maxbin), weight);
	}
}
void TH1FScaleRange::Fill(float value){
	for(int i = 0; i <= histocount; i++){
		histostore[i]->Fill(std::min(float(histofactor[i] * value),maxbin));
	}
}
void TH1FScaleRange::SetDirectory(TDirectory* dir){
	for(int i = 0; i <= histocount; i++){
		histostore[i]->SetDirectory(dir);
	}
}


#endif
