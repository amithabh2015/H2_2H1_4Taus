#ifndef DesyHTauTau_TauTauSkimming_HTauTauLikelihood
#define DesyHTauTau_TauTauSkimming_HTauTauLikelihood


#include <string.h>
//#include <vector.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>

typedef std::vector<float> FloatVector;
typedef std::vector<int> IntVector;
typedef std::vector<FloatVector> HistVector;
typedef std::vector<HistVector> ClassHistVector;
typedef std::vector<std::string> StringVector;

class HTauTauLikelihood {

 public:
  HTauTauLikelihood(StringVector classFileNames, StringVector histoNames);
  ~HTauTauLikelihood();
  float getLikelihood(FloatVector);

 private:

  unsigned int getBinNumber(unsigned int iHisto, float variable);

  unsigned int _numberOfClasses;
  unsigned int _numberOfHistos;

  StringVector _histoNames;
  FloatVector _xMinHisto;
  FloatVector _xMaxHisto;
  IntVector _nBinsHisto;
  FloatVector _binWidthHisto;
  StringVector _classFileNames;
  
  ClassHistVector _pdfVector;



};


#endif
