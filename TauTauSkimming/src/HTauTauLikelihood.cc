#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauLikelihood.h"

HTauTauLikelihood::HTauTauLikelihood(StringVector classFileNames, StringVector histoNames) {

  _numberOfClasses = classFileNames.size();
  _numberOfHistos  = histoNames.size();

  _xMinHisto.clear();
  _xMaxHisto.clear();
  _nBinsHisto.clear();
  _binWidthHisto.clear();

  _pdfVector.clear();
  for (unsigned int iFile=0; iFile<_numberOfClasses; iFile++) {
    TFile * file = new TFile(TString(classFileNames[iFile]));
    HistVector histVector;
    histVector.clear();
    for (unsigned int iHisto=0; iHisto<_numberOfHistos;++iHisto) {
      TH1F * histo = (TH1F*)file->Get(TString("analysis/"+TString(histoNames[iHisto])));
      int nBins = histo->GetNbinsX();
      float sum = 0.0;
      FloatVector histoContent;
      histoContent.clear();
      for (int iBin=0;iBin<nBins;++iBin) {
	float content = histo->GetBinContent(iBin+1);
	if (content<1.0) 
	  content = 1.0;
	histoContent.push_back(content);
	sum += content;
      }
      FloatVector histoContentNorm;
      histoContentNorm.clear();
      for (int iBin=0;iBin<nBins;++iBin) {
	float content = histoContent[iBin]/sum;
	histoContentNorm.push_back(content);
      }
      if (iFile==0) { // fill histo info;
	_nBinsHisto.push_back(nBins);
	_xMinHisto.push_back(histo->GetBinLowEdge(1));
	_xMaxHisto.push_back(histo->GetBinLowEdge(nBins));
	_binWidthHisto.push_back(histo->GetBinWidth(1));
      }
      histVector.push_back(histoContentNorm);
    }    
    _pdfVector.push_back(histVector);
    file->Close();
    delete file;
  }


}
HTauTauLikelihood::~HTauTauLikelihood() {}

float HTauTauLikelihood::getLikelihood(FloatVector variables) {
  
  float likelihood = 1.0;

  HistVector Prob;
  Prob.clear();
  for (unsigned int iClass=0;iClass<_numberOfClasses;++iClass) {
    FloatVector pdfs;
    pdfs.clear();
    for (unsigned int iHisto=0;iHisto<_numberOfHistos;++iHisto) {
      float sum = 0.0;
      unsigned int bin = getBinNumber(iHisto,variables[iHisto]);
      for (unsigned int iCl=0;iCl<_numberOfClasses;++iCl) {
	HistVector histVect = _pdfVector[iCl];
	FloatVector floatVect = histVect[iHisto];
	float value = floatVect[bin];
	sum += value;
      }
      HistVector histVector = _pdfVector[iClass];
      FloatVector floatVector = histVector[iHisto];
      float val = floatVector[bin];
      pdfs.push_back(val/sum);
    }    
    Prob.push_back(pdfs);
  }

  FloatVector classProb;
  classProb.clear();
  float sum = 0;
  for (unsigned int iClass=0;iClass<_numberOfClasses;++iClass) {
    float prod = 1.0;
    FloatVector floatVector = Prob[iClass];    
    for (unsigned int iHisto=0;iHisto<_numberOfHistos;++iHisto) {
      prod *= floatVector[iHisto];
    }
    classProb.push_back(prod);
    sum += prod;
  }

  likelihood = classProb[0]/sum;

  return likelihood;

}
unsigned int HTauTauLikelihood::getBinNumber(unsigned int iHisto, float variable) {

  unsigned int bin = 0;
  if (variable<_xMinHisto[iHisto]) {
    bin = 0;
    return bin;
  }
  if (variable>_xMaxHisto[iHisto]) {
    bin = _nBinsHisto[iHisto]-1;
    return bin;
  }

  int nbin = int((variable-_xMinHisto[iHisto])/_binWidthHisto[iHisto]);
  bin = unsigned(nbin);
  return bin;

}
