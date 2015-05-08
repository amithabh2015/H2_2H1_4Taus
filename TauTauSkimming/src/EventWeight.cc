#include "H2to2H1to4Taus/TauTauSkimming/interface/EventWeight.h"
#include "H2to2H1to4Taus/TauTauSkimming/interface/HTauTauUtils.h"

Double_t FitFuncExp(Double_t * x, Double_t * par) {

  Double_t a = (x[0]-par[1])/par[2];
  Double_t b = (x[0]-par[1])/par[3];
  Double_t c = (x[0]-par[1])/par[4];

  Double_t aLeft = (x[0]-par[1])/par[5];
  Double_t bLeft = (x[0]-par[1])/par[6];
  Double_t cLeft = (x[0]-par[1])/par[7];
  

  Double_t aGauss = par[8]*TMath::Exp(-0.5*a*a);
  Double_t bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*b*b);
  Double_t cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-TMath::Abs(c));

  if (x[0]<par[1]) {
    aGauss = par[8]*TMath::Exp(-0.5*aLeft*aLeft);
    bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*bLeft*bLeft);
    cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-TMath::Abs(cLeft));
  }

  Double_t res = par[0]*(aGauss+bGauss+cGauss);

  return res;


}


Double_t FitFuncGauss(Double_t * x, Double_t * par) {

  Double_t a = (x[0]-par[1])/par[2];
  Double_t b = (x[0]-par[1])/par[3];
  Double_t c = (x[0]-par[1])/par[4];

  Double_t aLeft = (x[0]-par[1])/par[5];
  Double_t bLeft = (x[0]-par[1])/par[6];
  Double_t cLeft = (x[0]-par[1])/par[7];
  

  Double_t aGauss = par[8]*TMath::Exp(-0.5*a*a);
  Double_t bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*b*b);
  Double_t cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-0.5*c*c);

  if (x[0]<par[1]) {
    aGauss = par[8]*TMath::Exp(-0.5*aLeft*aLeft);
    bGauss = (1-par[8])*par[9]*TMath::Exp(-0.5*bLeft*bLeft);
    cGauss = (1-par[8])*(1-par[9])*TMath::Exp(-0.5*cLeft*cLeft);
  }

  Double_t res = par[0]*(aGauss+bGauss+cGauss);

  return res;

}

TF1 * EventWeight::getFuncRecoil(TF1 * initFunc, bool left, bool isGauss) {

  double xminD;
  double xmaxD;

  initFunc->GetRange(xminD,xmaxD);

  float xmin = float(xminD);
  float xmax = float(xmaxD);

  TF1 * func = NULL;
  if (isGauss) 
    func = new TF1("func",FitFuncGauss,xmin,xmax,10);
  else
    func = new TF1("func",FitFuncExp,xmin,xmax,10);

  func->SetParameter(0,initFunc->GetParameter(0));
  func->SetParameter(1,initFunc->GetParameter(1));
  if (left) {
    func->SetParameter(2,initFunc->GetParameter(5));
    func->SetParameter(3,initFunc->GetParameter(6));
    func->SetParameter(4,initFunc->GetParameter(7));
    func->SetParameter(5,initFunc->GetParameter(5));
    func->SetParameter(6,initFunc->GetParameter(6));
    func->SetParameter(7,initFunc->GetParameter(7));
  }
  else {
     func->SetParameter(2,initFunc->GetParameter(2));
     func->SetParameter(3,initFunc->GetParameter(3));
     func->SetParameter(4,initFunc->GetParameter(4));
     func->SetParameter(5,initFunc->GetParameter(2));
     func->SetParameter(6,initFunc->GetParameter(3));
     func->SetParameter(7,initFunc->GetParameter(4));
  }
  func->SetParameter(8,initFunc->GetParameter(8));
  func->SetParameter(9,initFunc->GetParameter(9));
  

  return func;

}


EventWeight::EventWeight(TString baseDir, int dataset) : ggHWeightHqT_histo(0), ggHWeightFeHiPro_histo(0) {

  _baseDir = baseDir;

  _weightMin = 0.05;
  _weightMax = 20;

  _dataset = dataset;

  _puFall11MC = 0;
  _puFall11Data = 0;
  _puSummer12MC = 0;
  _puSummer12Data = 0;
}

float EventWeight::MuEff2012(float pt, float eta) {

  float absEta = TMath::Abs(eta);

  float eff = 1;

  if (absEta<1.48) {
    if (pt>10&&pt<=15)
      eff = 0.90;
    else if (pt>15&&pt<=20)
      eff = 0.92;
    else 
      eff = 0.95;
  }
  else {
    if (pt>10&&pt<=15)
      eff = 1.03;
    else if (pt>15&&pt<=20)
      eff = 0.96;
    else 
      eff = 0.99;
  }

  return eff;

}


EventWeight::~EventWeight() {

}


float EventWeight::PUWeightS3(int nPV) {

  float ratio = 1;

  if (_dataset<2) { // Run2011A
    if (nPV==0){ ratio = 1.0;}
    else if (nPV==1){ ratio=0.127483;}
    else if (nPV==2){ ratio=0.494128;}
    else if (nPV==3){ ratio=1.06568;}
    else if (nPV==4){ ratio=1.43145;}
    else if (nPV==5){ ratio=1.56026;}
    else if (nPV==6){ ratio=1.52698;}
    else if (nPV==7){ ratio=1.37306;}
    else if (nPV==8){ ratio=1.14338;}
    else if (nPV==9){ ratio=0.911818;}
    else if (nPV==10){ ratio=0.747833;}
    else if (nPV==11){ ratio=0.62478;}
    else if (nPV==12){ ratio=0.559537;}
    else if (nPV==13){ ratio=0.483121;}
    else if (nPV==14){ ratio=0.442508;}
    else if (nPV==15){ ratio=0.46228;}
    else if (nPV==16){ ratio=0.472298;}
    else if (nPV==17){ ratio=0.516508;}
    else if (nPV==18){ ratio=0.616794;}
    else if (nPV==19){ ratio=0.48273;}
    else {ratio=0.47;}

  }

  else { // Run2011B

    if (nPV==0){ ratio = 1.0;}
    else if (nPV==1){ ratio=0.00659921;}
    else if (nPV==2){ ratio=0.0415862;}
    else if (nPV==3){ ratio=0.143151;}
    else if (nPV==4){ ratio=0.305751;}
    else if (nPV==5){ ratio=0.513532;}
    else if (nPV==6){ ratio=0.772757;}
    else if (nPV==7){ ratio=1.05185;}
    else if (nPV==8){ ratio=1.30912;}
    else if (nPV==9){ ratio=1.55455;}
    else if (nPV==10){ ratio=1.87055;}
    else if (nPV==11){ ratio=2.3045;}
    else if (nPV==12){ ratio=2.98664;}
    else if (nPV==13){ ratio=3.78916;}
    else if (nPV==14){ ratio=4.9661;}
    else if (nPV==15){ ratio=7.3026;}
    else if (nPV==16){ ratio=9.80483;}
    else if (nPV==17){ ratio=14.5751;}
    else if (nPV==18){ ratio=23.3065;}
    else if (nPV==19){ ratio=24.0932;}
    else {ratio=1;}


  }
  return ratio;

}

float EventWeight::PUWeightS4(int nPV) {
  
  float ratio = 1;

  if (_dataset==0) { // Run2011A Single Mu

    if (nPV==0){ ratio = 1.0;}
    else if (nPV==1){ ratio=0.128697;}
    else if (nPV==2){ ratio=0.725498;}
    else if (nPV==3){ ratio=1.36142;}
    else if (nPV==4){ ratio=1.85477;}
    else if (nPV==5){ ratio=1.94621;}
    else if (nPV==6){ ratio=1.81656;}
    else if (nPV==7){ ratio=1.44193;}
    else if (nPV==8){ ratio=1.04542;}
    else if (nPV==9){ ratio=0.70794;}
    else if (nPV==10){ ratio=0.470454;}
    else if (nPV==11){ ratio=0.322615;}
    else if (nPV==12){ ratio=0.198063;}
    else if (nPV==13){ ratio=0.1426;}
    else if (nPV==14){ ratio=0.0656408;}
    else if (nPV==15){ ratio=0.0375877;}
    else if (nPV==16){ ratio=0.0384816;}
    else if (nPV==17){ ratio=0.0119606;}
    else if (nPV==18){ ratio=0.0145203;}
    else if (nPV==19){ ratio=0.0145572;}
    else { ratio=0.0;}

  }
  else if (_dataset==1)  { // Run2011A DoubleMu 
    if (nPV==0){ ratio = 1.0;}
    else if (nPV==1){ ratio=0.0721942;}
    else if (nPV==2){ ratio=0.439623;}
    else if (nPV==3){ ratio=0.877597;}
    else if (nPV==4){ ratio=1.30548;}
    else if (nPV==5){ ratio=1.61185;}
    else if (nPV==6){ ratio=1.70677;}
    else if (nPV==7){ ratio=1.63155;}
    else if (nPV==8){ ratio=1.47442;}
    else if (nPV==9){ ratio=1.26223;}
    else if (nPV==10){ ratio=1.02304;}
    else if (nPV==11){ ratio=0.814332;}
    else if (nPV==12){ ratio=0.644816;}
    else if (nPV==13){ ratio=0.500549;}
    else if (nPV==14){ ratio=0.374304;}
    else if (nPV==15){ ratio=0.287225;}
    else if (nPV==16){ ratio=0.228961;}
    else if (nPV==17){ ratio=0.179353;}
    else if (nPV==18){ ratio=0.135111;}
    else if (nPV==19){ ratio=0.132036;}
    else if (nPV==20){ ratio=0.102249;}
    else if (nPV==21){ ratio=0.0920071;}
    else if (nPV==22){ ratio=0.0622184;}
    else if (nPV==23){ ratio=0.0471337;}
    else if (nPV==24){ ratio=0.0795056;}
    else if (nPV==25){ ratio=0.0488271;}
    else { ratio = 0.0;}
  }
  else { // Run2011B sample

    if (nPV==0){ ratio = 1.0;}
    else if (nPV==1){ ratio=0.00494198;}
    else if (nPV==2){ ratio=0.0472267;}
    else if (nPV==3){ ratio=0.144856;}
    else if (nPV==4){ ratio=0.326938;}
    else if (nPV==5){ ratio=0.579248;}
    else if (nPV==6){ ratio=0.879075;}
    else if (nPV==7){ ratio=1.18919;}
    else if (nPV==8){ ratio=1.48624;}
    else if (nPV==9){ ratio=1.76241;}
    else if (nPV==10){ ratio=2.00601;}
    else if (nPV==11){ ratio=2.25499;}
    else if (nPV==12){ ratio=2.46866;}
    else if (nPV==13){ ratio=2.68827;}
    else if (nPV==14){ ratio=2.86616;}
    else if (nPV==15){ ratio=3.01842;}
    else if (nPV==16){ ratio=3.1753;}
    else if (nPV==17){ ratio=3.35453;}
    else if (nPV==18){ ratio=3.50421;}
    else if (nPV==19){ ratio=3.79908;}
    else if (nPV==20){ ratio=4.18774;}
    else if (nPV==21){ ratio=4.57301;}
    else if (nPV==22){ ratio=5.07736;}
    else if (nPV==23){ ratio=5.6735;}
    else if (nPV==24){ ratio=6.494;}
    else if (nPV==25){ ratio=7.91473;}
    else if (nPV==26){ ratio=8.0829;}
    else if (nPV==27){ ratio=10.2266;}
    else if (nPV==28){ ratio=10.4023;}
    else if (nPV==29){ ratio=10.2846;}
    else if (nPV==30){ ratio=11.0965;}
    else if (nPV==31){ ratio=9.50741;}
    else if (nPV==32){ ratio=9.19359;}
    else if (nPV==33){ ratio=8.0088;}
    else if (nPV==34){ ratio=7.0783;}
    else if (nPV==35){ ratio=5.49925;}
    else if (nPV==36){ ratio=4.92651;}
    else if (nPV==37){ ratio=4.83972;}
    else { ratio = 0.0;}

  }

  return ratio;

}

float EventWeight::PUWeightS7(int nPV) {

  float ratio = 1.0;

  if (nPV==0){ ratio = 1.0;}
  else if (nPV==1){ ratio=0.447307;}
  else if (nPV==2){ ratio=0.982477;}
  else if (nPV==3){ ratio=2.37712;}
  else if (nPV==4){ ratio=3.87165;}
  else if (nPV==5){ ratio=4.96636;}
  else if (nPV==6){ ratio=5.53855;}
  else if (nPV==7){ ratio=5.66279;}
  else if (nPV==8){ ratio=6.27496;}
  else if (nPV==9){ ratio=5.90384;}
  else if (nPV==10){ ratio=5.21282;}
  else if (nPV==11){ ratio=4.38568;}
  else if (nPV==12){ ratio=3.40988;}
  else if (nPV==13){ ratio=2.81099;}
  else if (nPV==14){ ratio=2.22016;}
  else if (nPV==15){ ratio=1.7313;}
  else if (nPV==16){ ratio=1.46289;}
  else if (nPV==17){ ratio=1.06598;}
  else if (nPV==18){ ratio=0.800481;}
  else if (nPV==19){ ratio=0.629233;}
  else if (nPV==20){ ratio=0.453126;}
  else if (nPV==21){ ratio=0.329058;}
  else if (nPV==22){ ratio=0.239505;}
  else if (nPV==23){ ratio=0.183793;}
  else if (nPV==24){ ratio=0.119271;}
  else if (nPV==25){ ratio=0.0839921;}
  else if (nPV==26){ ratio=0.0625807;}
  else if (nPV==27){ ratio=0.0508364;}
  else if (nPV==28){ ratio=0.0356859;}
  else if (nPV==29){ ratio=0.0258227;}
  else if (nPV==30){ ratio=0.0171003;}
  else if (nPV==31){ ratio=0.0125707;}
  else if (nPV==32){ ratio=0.00605623;}
  else if (nPV==33){ ratio=0.00512061;}
  else if (nPV==34){ ratio=0.00575931;}
  else if (nPV==35){ ratio=0.00402526;}
  else if (nPV==36){ ratio=0.0013942;}
  else if (nPV==37){ ratio=0.00207246;}
  else {ratio=0;}

  return ratio;

}


float EventWeight::PUWeightS6(int nPV ) {

  float ratio = 1;

  if (nPV==0){ ratio = 1.0;}
  else if (nPV==1){ ratio=0.499055;}
  else if (nPV==2){ ratio=0.848472;}
  else if (nPV==3){ ratio=1.11498;}
  else if (nPV==4){ ratio=1.28453;}
  else if (nPV==5){ ratio=1.37683;}
  else if (nPV==6){ ratio=1.40231;}
  else if (nPV==7){ ratio=1.38086;}
  else if (nPV==8){ ratio=1.32043;}
  else if (nPV==9){ ratio=1.24751;}
  else if (nPV==10){ ratio=1.14454;}
  else if (nPV==11){ ratio=1.03332;}
  else if (nPV==12){ ratio=0.895612;}
  else if (nPV==13){ ratio=0.772444;}
  else if (nPV==14){ ratio=0.64545;}
  else if (nPV==15){ ratio=0.51325;}
  else if (nPV==16){ ratio=0.402694;}
  else if (nPV==17){ ratio=0.317436;}
  else if (nPV==18){ ratio=0.250104;}
  else if (nPV==19){ ratio=0.199625;}
  else if (nPV==20){ ratio=0.155474;}
  else if (nPV==21){ ratio=0.128502;}
  else if (nPV==22){ ratio=0.109232;}
  else if (nPV==23){ ratio=0.0970126;}
  else if (nPV==24){ ratio=0.08588;}
  else if (nPV==25){ ratio=0.0652137;}
  else if (nPV==26){ ratio=0.0613578;}
  else if (nPV==27){ ratio=0.0511663;}
  else if (nPV==28){ ratio=0.050404;}
  else if (nPV==29){ ratio=0.0391082;}
  else if (nPV==30){ ratio=0.0319129;}
  else if (nPV==31){ ratio=0.0334386;}
  else if (nPV==32){ ratio=0.0330676;}
  else if (nPV==33){ ratio=0.0156537;}
  else if (nPV==34){ ratio=0.017986;}
  else if (nPV==35){ ratio=0.0126642;}
  else if (nPV==36){ ratio=0.00809074;}
  else {ratio=1;} 

  return ratio;

}




bool EventWeight::InitggHWeightFeHiPro(int higgsMass)
{
  std::cout << std::endl << "EventWeight::InitggHWeightFeHiPro()" << std::endl;

  TString tmpFilename = _baseDir+"/kfactors/fehipro/weight_ptH_";
  tmpFilename += higgsMass;
  tmpFilename += ".root";

  std::cout << "Open file with weights: " << tmpFilename << std::endl;
  ggHWeight_file = TFile::Open(tmpFilename);

  if (!ggHWeight_file)
  {
    std::cout << "...file not found!" << std::endl;
    return false;
  }

  TString tmpHistoname = "powheg_weight/weight_hqt_fehipro_fit_";
  tmpHistoname += higgsMass;

  std::cout << "Open histogram with weights: " << tmpHistoname << std::endl;
  ggHWeightFeHiPro_histo = (TH1F *)ggHWeight_file->Get(tmpHistoname);

  if (!ggHWeightFeHiPro_histo)
  {
    std::cout << "...histogram not found!" << std::endl;
    assert(ggHWeightFeHiPro_histo);
    return false;
  }
  std::cout << std::endl;
  return ggHWeightFeHiPro_histo;
}

float EventWeight::histIntegral(TH1F * hist, float x) {

  float norm = float(hist->Integral());
  int nBins = hist->GetNbinsX();
  int nBin = hist->FindBin(x);
  if (nBin<1) nBin = 1;
  if (nBin>nBins) nBin = nBins;
  float integral = 0;
  if (nBin==1) {
    float content = hist->GetBinWidth(1);
    float xmin = hist->GetBinLowEdge(1);
    float xmax = hist->GetBinLowEdge(2);
    integral = (x-xmin)*content/(xmax-xmin);
  }
  else {
    float xmin = hist->GetBinLowEdge(nBin);
    float xmax = hist->GetBinLowEdge(nBin+1);
    float imin = hist->Integral(1,nBin-1);
    float imax = hist->Integral(1,nBin);
    integral = imin + (x-xmin)*(imax-imin)/(xmax-xmin);
  }

  return integral/norm;

}

float EventWeight::doIsoMapping(TH1F * dataH, TH1F * mcH, float x) {

  double prob[1];
  double q[1];
  prob[0] = double(histIntegral(mcH,x));
  int nprobSum = 1;

  dataH->GetQuantiles(nprobSum, q, prob);
  
  float xcor = float(q[0]);
    
  return xcor;

}

void EventWeight::InitIsoMapping(TString dir,
				 TString dataFileName,
				 TString mcFileName,
				 std::vector<TString> histNames) {

  TString fullDataFileName = dir + "/" + dataFileName;
  TString fullMcFileName = dir + "/" + mcFileName;

  TFile * fileData = new TFile(fullDataFileName);
  TFile * fileMC   = new TFile(fullMcFileName);

  if (fileData->IsZombie()) {
    std::cout << "Cannot open file " << fullDataFileName << std::endl;
    exit(-1);
  }

  if (fileMC->IsZombie()) {
    std::cout << "Cannot open file " << fullMcFileName << std::endl;
    exit(-1);
  }

  histIsoMapData.clear();
  histIsoMapMC.clear();

  nHistosIsoMap = int(histNames.size());
  for (int iH=0; iH<nHistosIsoMap; ++iH) {
    TH1F * hData = (TH1F*)fileData->Get(histNames.at(iH));
    if (hData==NULL) {
      std::cout << "Histo " << histNames.at(iH) << " is not found in data file " << fullDataFileName << std::endl;
      exit(-1);
    }
    else {
      histIsoMapData.push_back(hData);
    }
    TH1F * hMC = (TH1F*)fileMC->Get(histNames.at(iH));
     if (hMC==NULL) {
       std::cout << "Histo " << histNames.at(iH) << " is not found in MC file " << fullMcFileName << std::endl;
       exit(-1);
    }
    else {
      histIsoMapMC.push_back(hMC);
    }

  }

}

float EventWeight::isoMapping(int iHist, float x) {

  float xcor = x;

  if (iHist>=nHistosIsoMap) {
    std::cout << "iso-mapping histogram index " << iHist << " is out of range [0," << nHistosIsoMap-1 << " : value not changed" << std::endl; 
  }
  else {
    xcor = doIsoMapping(histIsoMapData.at(iHist),histIsoMapMC.at(iHist),x);
  }

  return xcor;

}

bool EventWeight::InitggHWeightHqT(int higgsMass, float renormScale = 1.0, float factScale = 1.0)
{
  std::cout << std::endl << "EventWeight::InitggHWeightHqT()" << std::endl;

  TString tmpFilename = _baseDir+"/kfactors/Kfactors_";
  tmpFilename += higgsMass;
  tmpFilename += "_AllScales.root";

  std::cout << "Open file with weights: " << tmpFilename << std::endl;
  ggHWeight_file = TFile::Open(tmpFilename);

  if (!ggHWeight_file)
  {
    std::cout << "...file not found!" << std::endl;
    return false;
  }

  TString tmpHistoname = "kfactors/kfact_mh";
  tmpHistoname += higgsMass;
  tmpHistoname += "_ren";
  tmpHistoname += int(higgsMass*renormScale);
  tmpHistoname += "_fac";
  tmpHistoname += int(higgsMass*factScale);

  std::cout << "Open histogram with weights: " << tmpHistoname << std::endl;
  ggHWeightHqT_histo = (TH1F *)ggHWeight_file->Get(tmpHistoname);

  if (!ggHWeightHqT_histo)
  {
    std::cout << "...histogram not found!" << std::endl;
    assert(ggHWeightHqT_histo);
    return false;
  }
  std::cout << std::endl;
  return ggHWeightHqT_histo;
}

void EventWeight::InitTheoryWeights(TString theoryWeightFileName, TString theoryWeightHisto) {
	TFile * fileTh = new TFile(theoryWeightFileName);
	if (fileTh->IsZombie()) {
		std::cout << "Theory weight file " << theoryWeightFileName << " is not found in base directory " << std::endl << std::flush;
	}
	if (fileTh == NULL) std::cout << "no such file \"" << theoryWeightFileName << "\"" << std::endl;
	
	_theoryWeightHisto = (TH2F*) fileTh->Get(theoryWeightHisto);
//	_theoryWeightHistoDown = (TH2F*) fileTh->Get("m_deta_down");
	if (_theoryWeightHisto == NULL) std::cout << "no such histo " << theoryWeightHisto << std::endl;
}

void EventWeight::InitPUWeightsTruth(TString puFileName) {

  TFile * filePU = new TFile(_baseDir+"/"+puFileName);
  if (filePU->IsZombie()) {
    std::cout << "File " << puFileName << " is not found in directory " <<  _baseDir << std::endl;
  }

  _puHisto = (TH1F*)filePU->Get("theirweights");
}

bool EventWeight::InitPUWeightsS6Truth(TString puData2011, TString puMCFall11) {

// pile-up reweighting with truth information:

/*
Use the following code to create the data histogram:
	hadd RooT/PileUp_2011.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_165088-167913_7TeV_PromptReco_JSON.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_172620-173692_PromptReco_JSON.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_175832-177515_PromptReco_JSON.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_177718_178078_7TeV_PromptReco_Collisons11_JSON.pileupTruth_v2_finebin.root \
	/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_178098-180252_7TeV_PromptReco_Collisions11_JSON.pileupTruth_v2_finebin.root

Use the following code to create the MC histogram:
	TChain myChain("analysis/HTauTau2Muons")
	myChain.Add("/storage/6/zeise/events/ztautau_mumu/raw/2011-12-03_428/ztautau_mumu_DYToMuMu_Z2_powheg_0020_FA11S42v14BPUS6_v1_*.root")
	myChain.Draw("nPUTruth >> pileup(100,0,25)")
	TFile outFile("RooT/PileUp_DYToMuMu_Z2_powheg_0020_FA11S42v14BPUS6_v1.root","RECREATE")
	pileup->Write();
	outFile->Close();
*/

	TFile * puFall11MCFile = new TFile(_baseDir+"/"+puMCFall11);
	TFile * puFall11DataFile = new TFile(_baseDir+"/"+puData2011);

	if (puFall11MCFile && puFall11DataFile)
	{
		std::cout << "found" << std::endl;
		_puFall11MC = (TH1F*) puFall11MCFile->Get("pileup");
		_puFall11Data = (TH1F*) puFall11DataFile->Get("pileup");

		/*
			The applied rebinning is not really necessary as the
			reweighting function uses FindBin(). The same applies to the
			range (0,25) for the MC histogram. However, the histograms for
			both data and MC have the same binning to avoid edge/binning
			effects. A coarse binning should also limit the effect of the
			missing random seed during the pile-up creation process.
		*/
		//_puFall11Data->Rebin(10.);
		
		_puFall11MC->Scale(1./_puFall11MC->Integral("width"));
		_puFall11Data->Scale(1./_puFall11Data->Integral("width"));

		/*
		for (int idx1 = 1; idx1 < _puFall11MC->GetNbinsX(); idx1++)
		{
			std::cout << idx1 << "\t" << _puFall11Data->GetBinContent(idx1) << "\t" << _puFall11MC->GetBinContent(idx1) << "\n";
		}
		*/
		
		return true;
	}
	else
	{
		_puFall11MC = 0;
		_puFall11Data = 0;
		return false;
	}
}

bool EventWeight::InitPUWeightsS10Truth(TString puDataFileName, TString puMCFileName) {

	TFile * puMCFile = new TFile(_baseDir+"/"+puMCFileName);
	TFile * puDataFile = new TFile(_baseDir+"/"+puDataFileName);

	if (puMCFile && puDataFile)
	{
		_puHCP12MC   = (TH1F*) puMCFile->Get("pileup");
		_puHCP12Data = (TH1F*) puDataFile->Get("pileup");

		_puHCP12MC->Sumw2();
		_puHCP12Data->Sumw2();

		/*
		if (_puHCP12Data->GetNbinsX()*2 == _puHCP12MC->GetNbinsX())
			_puHCP12MC->Rebin(2);

		if (_puHCP12Data->GetNbinsX() != _puHCP12MC->GetNbinsX())
		{
			std::cout << _baseDir+"/"+puMCFileName << ": " << _puHCP12MC->GetNbinsX() << " bins" << std::endl;
			std::cout << _baseDir+"/"+puDataFileName << ": " << _puHCP12Data->GetNbinsX() << " bins" << std::endl;
		}
		*/

		_puHCP12MC->Scale(1./_puHCP12MC->Integral());
		_puHCP12Data->Scale(1./_puHCP12Data->Integral());

		/*
		for (int idx1 = 1; idx1 < _puHCP12MC->GetNbinsX(); idx1++)
			std::cout << idx1 << "\t" << _puHCP12MC->GetBinContent(idx1) << "\t" << _puHCP12Data->GetBinContent(idx1) << "\n";
		*/

		return true;
	}
	else
	{
		_puHCP12MC = 0;
		_puHCP12Data = 0;
		return false;
	}


}

bool EventWeight::InitPUWeightsS7Truth(TString puData2012, TString puMCSummer12) {
	TFile * puSummer12MCFile = new TFile(_baseDir+"/"+puMCSummer12);
	TFile * puSummer12DataFile = new TFile(_baseDir+"/"+puData2012);

	if (puSummer12MCFile && puSummer12DataFile)
	{
		_puSummer12MC = (TH1F*) puSummer12MCFile->Get("pileup");
		_puSummer12Data = (TH1F*) puSummer12DataFile->Get("pileup");

		_puSummer12MC->Sumw2();
		_puSummer12Data->Sumw2();

		/*
		if (_puSummer12Data->GetNbinsX()*2 == _puSummer12MC->GetNbinsX())
			_puSummer12MC->Rebin(2);
		if (_puSummer12Data->GetNbinsX() != _puSummer12MC->GetNbinsX())
		{
			std::cout << _baseDir+"/"+puMCSummer12 << ": " << _puSummer12MC->GetNbinsX() << " bins" << std::endl;
			std::cout << _baseDir+"/"+puData2012 << ": " << _puSummer12Data->GetNbinsX() << " bins" << std::endl;
		}
		assert(_puSummer12Data->GetNbinsX() == _puSummer12MC->GetNbinsX());
		*/
		
		_puSummer12MC->Scale(1./_puSummer12MC->Integral());
		_puSummer12Data->Scale(1./_puSummer12Data->Integral());

		/*
		for (int idx1 = 1; idx1 < _puSummer12MC->GetNbinsX(); idx1++)
		{
			std::cout << idx1 << "\t" << _puSummer12MC->GetBinContent(idx1) << "\t" << _puSummer12Data->GetBinContent(idx1) << "\n";
		}
		*/
		
		return true;
	}
	else
	{
		_puSummer12MC = 0;
		_puSummer12Data = 0;
		return false;
	}
}


float EventWeight::PUWeightTruth(int nPUI) {

  float reweight = float(_puHisto->GetBinContent(nPUI+1));

  if (_dataset==1) // Run2011A_DoubleMu
    // accounts for Normalization
    reweight *= 1.033;
  if (_dataset==2) { // Run2011B
    if (nPUI==0)
      reweight *= 0.024;
    if (nPUI==1)
      reweight *= 0.387;
    // accounts for Normalization
    reweight *= 1.034;
  }

  return reweight;

}

float EventWeight::PUWeightS6Truth(float mean)
{
	if (!_puFall11MC || !_puFall11Data) return 1.;
	
	size_t binMC = _puFall11MC->FindBin(mean);
	size_t binData = _puFall11Data->FindBin(mean);
	
	float weightMC = _puFall11MC->GetBinContent(binMC) / _puFall11MC->GetBinWidth(binMC);
	float weightData = _puFall11Data->GetBinContent(binData) / _puFall11Data->GetBinWidth(binData);
	
	if (weightMC==0.) return 1.;
	else return weightData/weightMC;
}

float EventWeight::PUWeightS7Truth(float mean)
{
	if (!_puSummer12MC || !_puSummer12Data) return 1.;

	size_t binMC = _puSummer12MC->FindBin(mean);
	size_t binData = _puSummer12Data->FindBin(mean);
	
	float weightMC = _puSummer12MC->GetBinContent(binMC) / _puSummer12MC->GetBinWidth(binMC);
	float weightData = _puSummer12Data->GetBinContent(binData) / _puSummer12Data->GetBinWidth(binData);

	if (weightMC==0.) return 1.;
	else return weightData/weightMC;
}

float EventWeight::PUWeightS10Truth(float mean)
{
	if (!_puHCP12MC || !_puHCP12Data) return 1.;

	size_t binMC = _puHCP12MC->FindBin(mean);
	size_t binData = _puHCP12Data->FindBin(mean);
	
	float weightMC = _puHCP12MC->GetBinContent(binMC) / _puHCP12MC->GetBinWidth(binMC);
	float weightData = _puHCP12Data->GetBinContent(binData) / _puHCP12Data->GetBinWidth(binData);

	if (weightMC==0.) return 1.;
	else return weightData/weightMC;
}


//FIXME: this function is just a dummy to get it compiling
void EventWeight::InitMEtWeights(TString dataFile,
				TString MCFile,
				int nZPtBins,
				float * ZPtBins,
				TString * _perpZStr,
				TString * _paralZStr,
				int nJetsBins,
				TString * _nJetsStr)
{
	std::vector<float> newZPtBins;
	std::vector<std::string> newPerpZStr;
	std::vector<std::string> newParalZStr;
	std::vector<std::string> newNJetsStr;

	for (int idx=0; idx<nZPtBins+1; ++idx)
		newZPtBins.push_back(ZPtBins[idx]);

	for (int idx=0; idx<nZPtBins; ++idx)
	{
		newPerpZStr.push_back(std::string(_perpZStr[idx]));
		newParalZStr.push_back(std::string(_paralZStr[idx]));
	}

	for (int idx=0; idx<nJetsBins; ++idx)
		newNJetsStr.push_back(std::string(_nJetsStr[idx]));

	InitMEtWeights(dataFile,
		MCFile,
		newZPtBins,
		newPerpZStr,
		newParalZStr,
		newNJetsStr);
}


void EventWeight::InitMEtWeights(TString dataFile,
				TString MCFile,
				const std::vector<float>& ZPtBins,
				const std::vector<std::string>& _perpZStr,
				const std::vector<std::string>& _paralZStr,
				const std::vector<std::string>& _nJetsStr)
{

  TFile * _fileMetData = new TFile(_baseDir+"/"+dataFile);
  TFile * _fileMetMC   = new TFile(_baseDir+"/"+MCFile);

  // checking files
  if (_fileMetData->IsZombie()) {
    std::cout << "File " << dataFile << " is not found in directory " << _baseDir << std::endl;
    std::cout << "quitting program..." << std::endl;
    exit(-1);
  }

  if (_fileMetMC->IsZombie()) {
    std::cout << "File " << MCFile << " is not found in directory " << _baseDir << std::endl;
    std::cout << "quitting program..." << std::endl;
    exit(-1);
  }

  _nZPtBins = ZPtBins.size()-1; // the -1 is on purpose!
  _nJetsBins = _nJetsStr.size();
  
  _ZPtBins = ZPtBins;


  for (int ZPtBin=0; ZPtBin<_nZPtBins; ++ZPtBin) {
    for (int jetBin=0; jetBin<_nJetsBins; ++jetBin) {

      TString binStrPerp = _perpZStr[ZPtBin]+_nJetsStr[jetBin];
      TString binStrParal = _paralZStr[ZPtBin]+_nJetsStr[jetBin];

      std::cout << binStrPerp << "  *  " << binStrParal << " ----->" << std::endl;

      TString binStrPerpH = binStrPerp + TString("H");
      TString binStrParalH = binStrParal + TString("H");

      TString binStrPerpHist = binStrPerp + TString("Hist");
      TString binStrParalHist = binStrParal + TString("Hist");      

      _metZParalData[ZPtBin][jetBin] = (TF1*)_fileMetData->Get(binStrParalH);
      _metZPerpData[ZPtBin][jetBin]  = (TF1*)_fileMetData->Get(binStrPerpH);
      _metZParalMC[ZPtBin][jetBin]   = (TF1*)_fileMetMC->Get(binStrParalH);
      _metZPerpMC[ZPtBin][jetBin]    = (TF1*)_fileMetMC->Get(binStrPerpH);

      //      _metZParalDataHist[ZPtBin][jetBin] = (TH1F*)_fileMetData->Get(binStrParalHist);
      //      _metZPerpDataHist[ZPtBin][jetBin]  = (TH1F*)_fileMetData->Get(binStrPerpHist);
      //      _metZParalMCHist[ZPtBin][jetBin]   = (TH1F*)_fileMetMC->Get(binStrParalHist);
      //      _metZPerpMCHist[ZPtBin][jetBin]    = (TH1F*)_fileMetMC->Get(binStrPerpHist);


      // checking functions
      if (_metZParalData[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrParalH
		  << " is not found in file " << dataFile << "... quitting program..." << std::endl;
	exit(-1);

      }
      if (_metZPerpData[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrPerpH
		  << " is not found in file " << dataFile << "... quitting program..." << std::endl;
	exit(-1);
	
      }

      if (_metZParalMC[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrParalH
		  << " is not found in file " << MCFile << "... quitting program..." << std::endl;
	exit(-1);

      }
      if (_metZPerpMC[ZPtBin][jetBin]==NULL) {
	std::cout << "Function with name " << binStrPerpH
		  << " is not found in file " << MCFile << "... quitting program..." << std::endl;
	exit(-1);
	
      }

      // checking histograms
//       if (_metZParalDataHist[ZPtBin][jetBin]==NULL) {
// 	std::cout << "Histogram with name " << binStrParalHist
// 		  << " is not found in file " << dataFile << "... quitting program..." << std::endl;
// 	exit(-1);

//       }
//       if (_metZPerpDataHist[ZPtBin][jetBin]==NULL) {
// 	std::cout << "Histogram with name " << binStrPerpHist
// 		  << " is not found in file " << dataFile << "... quitting program..." << std::endl;
// 	exit(-1);
	
//       }

//       if (_metZParalMCHist[ZPtBin][jetBin]==NULL) {
// 	std::cout << "Histogram with name " << binStrParalHist
// 		  << " is not found in file " << MCFile << "... quitting program..." << std::endl;
// 	exit(-1);

//       }
//       if (_metZPerpMCHist[ZPtBin][jetBin]==NULL) {
// 	std::cout << "Histogram with name " << binStrPerpHist
// 		  << " is not found in file " << MCFile << "... quitting program..." << std::endl;
// 	exit(-1);
	
//       }

      // Met Paral Data
      
      double xminD,xmaxD;

      bool isGauss = true;
      //      if (jetBin>1)
      //	isGauss = false;

      _metZParalData[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZParalData[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZParalData[ZPtBin][jetBin] = float(xmaxD);

      TF1 * func = getFuncRecoil(_metZParalData[ZPtBin][jetBin],true,isGauss);
      _rmsLeftMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZParalData[ZPtBin][jetBin],false,isGauss);
      _rmsRightMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));
      delete func;
      
      _meanMetZParalData[ZPtBin][jetBin] = _metZParalData[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZParalData[ZPtBin][jetBin] = TMath::Sqrt(_metZParalData[ZPtBin][jetBin]->CentralMoment(2,_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]));

      float integral =  _metZParalData[ZPtBin][jetBin]->Integral(_xminMetZParalData[ZPtBin][jetBin],_xmaxMetZParalData[ZPtBin][jetBin]);
      std::cout << "   Data Paral ----> " << std::endl;
      std::cout << "   X0 = " << _meanMetZParalData[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZParalData[ZPtBin][jetBin] << std::endl; 
      std::cout << "   Integral [" << _xminMetZParalData[ZPtBin][jetBin] << "," << _xmaxMetZParalData[ZPtBin][jetBin] << "] = " << integral << std::endl;

      // Met Perp Data

      _metZPerpData[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZPerpData[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZPerpData[ZPtBin][jetBin] = float(xmaxD);


      func = getFuncRecoil(_metZPerpData[ZPtBin][jetBin],true,isGauss);
      _rmsLeftMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZPerpData[ZPtBin][jetBin],false,isGauss);
      _rmsRightMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));
      delete func;
      
      _meanMetZPerpData[ZPtBin][jetBin] = _metZPerpData[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZPerpData[ZPtBin][jetBin] = TMath::Sqrt(_metZPerpData[ZPtBin][jetBin]->CentralMoment(2,_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]));

      integral =  _metZPerpData[ZPtBin][jetBin]->Integral(_xminMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpData[ZPtBin][jetBin]);
      std::cout << "   Data Perp -----> " << std::endl;
      std::cout << "   X0 = " << _meanMetZPerpData[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZPerpData[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZPerpData[ZPtBin][jetBin] << "," << _xmaxMetZPerpData[ZPtBin][jetBin] << "] = " << integral << std::endl;

     
      // Met Paral MC
      
      _metZParalMC[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZParalMC[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZParalMC[ZPtBin][jetBin] = float(xmaxD);

      func = getFuncRecoil(_metZParalMC[ZPtBin][jetBin],true,isGauss);
      _rmsLeftMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZParalMC[ZPtBin][jetBin],false,isGauss);
      _rmsRightMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));
      delete func;

      _meanMetZParalMC[ZPtBin][jetBin] = _metZParalMC[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZParalMC[ZPtBin][jetBin] = TMath::Sqrt(_metZParalMC[ZPtBin][jetBin]->CentralMoment(2,_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]));

      integral =  _metZParalMC[ZPtBin][jetBin]->Integral(_xminMetZParalMC[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]);
      std::cout << "   MC Paral ------> " << std::endl;
      std::cout << "   X0 = " << _meanMetZParalMC[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZParalMC[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZParalMC[ZPtBin][jetBin] << "," << _xmaxMetZParalMC[ZPtBin][jetBin] << "] = " << integral << std::endl;
     

      // Met Perp MC

      _metZPerpMC[ZPtBin][jetBin]->GetRange(xminD,xmaxD);
      _xminMetZPerpMC[ZPtBin][jetBin] = float(xminD);
      _xmaxMetZPerpMC[ZPtBin][jetBin] = float(xmaxD);

      func = getFuncRecoil(_metZPerpMC[ZPtBin][jetBin],true,isGauss);
      _rmsLeftMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));
      delete func;

      func = getFuncRecoil(_metZPerpMC[ZPtBin][jetBin],false,isGauss);
      _rmsRightMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(func->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));
      delete func;

      _meanMetZPerpMC[ZPtBin][jetBin] = _metZPerpMC[ZPtBin][jetBin]->GetParameter(1);
      _rmsMetZPerpMC[ZPtBin][jetBin] = TMath::Sqrt(_metZPerpMC[ZPtBin][jetBin]->CentralMoment(2,_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]));

      integral =  _metZPerpMC[ZPtBin][jetBin]->Integral(_xminMetZPerpMC[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]);
      std::cout << "   MC Perp -------> " << std::endl;
      std::cout << "   X0 = " << _meanMetZPerpMC[ZPtBin][jetBin] << "  ;  RMS = " <<  _rmsMetZPerpMC[ZPtBin][jetBin] << std::endl;
      std::cout << "   Integral [" << _xminMetZPerpMC[ZPtBin][jetBin] << "," << _xmaxMetZPerpMC[ZPtBin][jetBin] << "] = " << integral << std::endl;

      _xminMetZParal[ZPtBin][jetBin] = TMath::Max(_xminMetZParalData[ZPtBin][jetBin],_xminMetZParalMC[ZPtBin][jetBin]);
      _xmaxMetZParal[ZPtBin][jetBin] = TMath::Min(_xmaxMetZParalData[ZPtBin][jetBin],_xmaxMetZParalMC[ZPtBin][jetBin]);

      _xminMetZPerp[ZPtBin][jetBin] = TMath::Max(_xminMetZPerpData[ZPtBin][jetBin],_xminMetZPerpMC[ZPtBin][jetBin]);
      _xmaxMetZPerp[ZPtBin][jetBin] = TMath::Min(_xmaxMetZPerpData[ZPtBin][jetBin],_xmaxMetZPerpMC[ZPtBin][jetBin]);
      
      std::cout << "   Perp  : min = " <<  _xminMetZPerp[ZPtBin][jetBin] << "    xmax = " <<  _xmaxMetZPerp[ZPtBin][jetBin] << std::endl;
      std::cout << "   Paral : min = " <<  _xminMetZParal[ZPtBin][jetBin] << "    xmax = " <<  _xmaxMetZParal[ZPtBin][jetBin] << std::endl;
      std::cout << std::endl;

      // specifying ranges ---->
//       float range = 55;
//       if (jetBin>0)
// 	range = 50;

//       _xminMetZPerp[ZPtBin][jetBin] = -range;
//       _xmaxMetZPerp[ZPtBin][jetBin] = range;

//       _xminMetZParal[ZPtBin][jetBin] = -range;
//       _xmaxMetZParal[ZPtBin][jetBin] = range;


    }
  }



  //  std::cout << "Met functions downloaded...." << std::endl;
  // exit(-1);

}

void EventWeight::InitDiMuEtaWeight(TString etaWeightsFile, float etaMax) {

  TFile * file = new TFile(_baseDir+"/"+etaWeightsFile);
  if (file->IsZombie()) {
    std::cout << "File " << etaWeightsFile << " does not exist in base directory... quitting program..." << std::endl;
    exit(-1);
  }

  _diLepEtaWeightFunc = (TF1*)file->Get("DiLepEtaRatio");

  if (_diLepEtaWeightFunc==NULL) {
    std::cout << "Function DiLepEtaRatio is not found in " << etaWeightsFile << " file... quitting program..." << std::endl;
    exit(-1);
  }

  _etaMax = etaMax;

}

void EventWeight::InitRhoWeights(TString rhoWeightsFile, int nPVMax) {

  _rhoFuncMC.clear();
  _rhoFuncData.clear();

  TString fileName = _baseDir+"/"+rhoWeightsFile;

  TFile * file = new TFile(fileName);

  if (file->IsZombie()) {
    std::cout << "File " << fileName << " does not exist... quitting program..." << std::endl;
    exit(-1);
  }
  else {
    std::cout <<  "File " << fileName << " successfully opened" << std::endl;
  }

  std::vector<TString> rhoFuncDataNames; rhoFuncDataNames.clear();
  std::vector<TString> rhoFuncMCNames; rhoFuncMCNames.clear();


  rhoFuncDataNames.push_back("rhoNeutralData_1");
  rhoFuncDataNames.push_back("rhoNeutralData_2");
  rhoFuncDataNames.push_back("rhoNeutralData_3");
  rhoFuncDataNames.push_back("rhoNeutralData_4");
  rhoFuncDataNames.push_back("rhoNeutralData_5");
  rhoFuncDataNames.push_back("rhoNeutralData_6");
  rhoFuncDataNames.push_back("rhoNeutralData_7");
  rhoFuncDataNames.push_back("rhoNeutralData_8");
  rhoFuncDataNames.push_back("rhoNeutralData_9");
  rhoFuncDataNames.push_back("rhoNeutralData_10");
  rhoFuncDataNames.push_back("rhoNeutralData_11");
  rhoFuncDataNames.push_back("rhoNeutralData_12");
  rhoFuncDataNames.push_back("rhoNeutralData_13");
  rhoFuncDataNames.push_back("rhoNeutralData_14");
  rhoFuncDataNames.push_back("rhoNeutralData_15");
  rhoFuncDataNames.push_back("rhoNeutralData_16");
  rhoFuncDataNames.push_back("rhoNeutralData_17");
  rhoFuncDataNames.push_back("rhoNeutralData_18");
  rhoFuncDataNames.push_back("rhoNeutralData_19");
  rhoFuncDataNames.push_back("rhoNeutralData_20");
  rhoFuncDataNames.push_back("rhoNeutralData_21");
  rhoFuncDataNames.push_back("rhoNeutralData_22");
  rhoFuncDataNames.push_back("rhoNeutralData_23");
  rhoFuncDataNames.push_back("rhoNeutralData_24");
  rhoFuncDataNames.push_back("rhoNeutralData_25");

  rhoFuncMCNames.push_back("rhoNeutralMC_1");
  rhoFuncMCNames.push_back("rhoNeutralMC_2");
  rhoFuncMCNames.push_back("rhoNeutralMC_3");
  rhoFuncMCNames.push_back("rhoNeutralMC_4");
  rhoFuncMCNames.push_back("rhoNeutralMC_5");
  rhoFuncMCNames.push_back("rhoNeutralMC_6");
  rhoFuncMCNames.push_back("rhoNeutralMC_7");
  rhoFuncMCNames.push_back("rhoNeutralMC_8");
  rhoFuncMCNames.push_back("rhoNeutralMC_9");
  rhoFuncMCNames.push_back("rhoNeutralMC_10");
  rhoFuncMCNames.push_back("rhoNeutralMC_11");
  rhoFuncMCNames.push_back("rhoNeutralMC_12");
  rhoFuncMCNames.push_back("rhoNeutralMC_13");
  rhoFuncMCNames.push_back("rhoNeutralMC_14");
  rhoFuncMCNames.push_back("rhoNeutralMC_15");
  rhoFuncMCNames.push_back("rhoNeutralMC_16");
  rhoFuncMCNames.push_back("rhoNeutralMC_17");
  rhoFuncMCNames.push_back("rhoNeutralMC_18");
  rhoFuncMCNames.push_back("rhoNeutralMC_19");
  rhoFuncMCNames.push_back("rhoNeutralMC_20");
  rhoFuncMCNames.push_back("rhoNeutralMC_21");
  rhoFuncMCNames.push_back("rhoNeutralMC_22");
  rhoFuncMCNames.push_back("rhoNeutralMC_23");
  rhoFuncMCNames.push_back("rhoNeutralMC_24");
  rhoFuncMCNames.push_back("rhoNeutralMC_25");

  _nPVMax = nPVMax;

  for (int iH=0; iH<_nPVMax; ++iH) {

    TF1 * rhoFunDataCurrent = (TF1*)file->Get(rhoFuncDataNames.at(iH));
    if (rhoFunDataCurrent==NULL) {
      std::cout << "Histogram " << rhoFuncDataNames.at(iH) 
		<< "  does not exist in file RhoFunc.root ... quitting program... " << std::endl;
      exit(-1);
    }
    _rhoFuncData.push_back(rhoFunDataCurrent);

    TF1 * rhoFunMCCurrent = (TF1*)file->Get(rhoFuncMCNames.at(iH));
    if (rhoFunMCCurrent==NULL) {
      std::cout << "Histogram " << rhoFuncMCNames.at(iH) 
		<< "  does not exist in file RhoFunc.root ... quitting program... " << std::endl;
      exit(-1);
    }
    _rhoFuncMC.push_back(rhoFunMCCurrent);
    
  }

  _rhoDataBin1 = (TH1F*)file->Get("firstBinData");
  if (_rhoDataBin1==NULL) {
    std::cout << "Histogram firstBinData does not exist in file " << rhoWeightsFile << "  ... quitting program " << std::endl;
    exit(-1);
  }

  _rhoDataBin2 = (TH1F*)file->Get("secondBinData");
  if (_rhoDataBin2==NULL) {
    std::cout << "Histogram secondBinData does not exist in file " << rhoWeightsFile << "  ... quitting program " << std::endl;
    exit(-1);
  }

  _rhoMCBin1   = (TH1F*)file->Get("firstBinMC");
  if (_rhoMCBin1==NULL) {
    std::cout << "Histogram firstBinMC does not exist in file " << rhoWeightsFile << "  ... quitting program " << std::endl;
    exit(-1);
  }

  _rhoMCBin2   = (TH1F*)file->Get("secondBinMC");
  if (_rhoMCBin2==NULL) {
    std::cout << "Histogram secondBinMC does not exist in file " << rhoWeightsFile << "  ... quitting program " << std::endl;
    exit(-1);
  }  


}

void EventWeight::InitDCAWeights2012()
{
	TString fileNameData = _baseDir+"/Data_Dca_Run2012A.root";
	TString fileNameMC   = _baseDir+"/MC_Dca_Summer12.root";

	TFile * fileDCAData = new TFile(fileNameData);
	TFile * fileDCAMC   = new TFile(fileNameMC);

	assert(fileDCAData->IsOpen() && fileDCAMC->IsOpen());

	_funcDcaData2012A = (TF1*)fileDCAData->Get("_Full_func");
	_funcDcaMCSummer12 = (TF1*)fileDCAMC->Get("_Full_func");

	assert(_funcDcaData2012A && _funcDcaMCSummer12);
}


float EventWeight::dcaWeight2012(float dca) {

  float weight = 1;

  if ( dca > -4.49 && dca < 1.19  ) {
    float numerator = _funcDcaData2012A->Eval(dca);
    float denominator = _funcDcaMCSummer12->Eval(dca);
    weight = numerator / denominator;
  }

  return weight; 

}


void EventWeight::InitDCAWeights(TString dcaFuncDataFile,
				 TString dcaFuncMCFile,
				 int nDiLepMassBins,
				 float * diLepMassBins,
				 TString * massRangeTString,
				 int nReducedLikeBins,
				 float * reducedLikeBins,
				 TString * likeCutTString)
{
	std::vector<float> newDiLepMassBins;
	for (int idx=0; idx<nDiLepMassBins+1; ++idx)
		newDiLepMassBins.push_back(diLepMassBins[idx]);

	std::vector<std::string> newMassRangeTString;
	for (int idx=0; idx<nDiLepMassBins; ++idx)
		newMassRangeTString.push_back(std::string(massRangeTString[idx]));

	std::vector<float> newReducedLikeBins;
	for (int idx=0; idx<nReducedLikeBins+1; ++idx)
		newReducedLikeBins.push_back(diLepMassBins[idx]);

	std::vector<std::string> newLikeCutTString;
	for (int idx=0; idx<nReducedLikeBins; ++idx)
		newLikeCutTString.push_back(std::string(likeCutTString[idx]));

	InitDCAWeights(dcaFuncDataFile,
				 dcaFuncMCFile,
				newDiLepMassBins,
				newMassRangeTString,
				newReducedLikeBins,
				newLikeCutTString);
}

void EventWeight::InitDCAWeights(TString dcaFuncDataFile,
				 TString dcaFuncMCFile,
				 const std::vector<float>& diLepMassBins,
				 const std::vector<std::string>& massRangeTString,
				 const std::vector<float>& reducedLikeBins,
				 const std::vector<std::string>& likeCutTString)
{
  TString fileNameData = _baseDir+"/"+dcaFuncDataFile;
  TFile * fileDataDCA = new TFile(fileNameData);
  if (fileDataDCA->IsZombie()) {
    std::cout << "File " << dcaFuncDataFile << " does not exist in base directory... quitting program..." << std::endl;
    exit(-1);
  }


  TString fileNameMc = _baseDir+"/"+dcaFuncMCFile;
  TFile * fileMcDCA = new TFile(fileNameMc);
  if (fileMcDCA->IsZombie()) {
    std::cout << "File " << dcaFuncMCFile << " does not exist in base directory... quitting program..." << std::endl;
    exit(-1);
  }

  for (size_t iM=0; iM<diLepMassBins.size()-1;++iM) {
    for (size_t iRL=0; iRL<reducedLikeBins.size()-1; iRL++) {
      _totWeightDCA[iM][iRL] = 1.0;
    }
  }

  if (_dataset==2) {

    _totWeightDCA[0][0] = 1.025 ;
    _totWeightDCA[0][1] = 1.026 ;
    _totWeightDCA[0][2] = 1.005 ;
    _totWeightDCA[0][3] = 1.030 ;

  }


  _nDiLepMassBins   = diLepMassBins.size() - 1;
  _nReducedLikeBins = reducedLikeBins.size() - 1;

  _reducedLikeBins = reducedLikeBins;
  _diLepMassBins   = diLepMassBins; 


  for (int iMassRange=0; iMassRange<_nDiLepMassBins;++iMassRange) {
    for (int iRedLike=0; iRedLike<_nReducedLikeBins; ++iRedLike) {

      TString baseString = massRangeTString[iMassRange] + likeCutTString[iRedLike];

      std::cout << baseString << "   ------> " << std::endl;

      TString funcDcaCentralString = massRangeTString[iMassRange] + likeCutTString[iRedLike] + "_Central_func";
      _funcDCACentralData[iMassRange][iRedLike] = (TF1*)fileDataDCA->Get(funcDcaCentralString);
      _funcDCACentralMC[iMassRange][iRedLike] = (TF1*)fileMcDCA->Get(funcDcaCentralString);

      TString funcDcaLeftTailString = massRangeTString[iMassRange] + likeCutTString[iRedLike] + "_LeftTail_func";
      _funcDCALeftTailData[iMassRange][iRedLike] = (TF1*)fileDataDCA->Get(funcDcaLeftTailString);
      _funcDCALeftTailMC[iMassRange][iRedLike] = (TF1*)fileMcDCA->Get(funcDcaLeftTailString);

      TString funcDcaRightTailString = massRangeTString[iMassRange] + likeCutTString[iRedLike] + "_RightTail_func";
      _funcDCARightTailData[iMassRange][iRedLike]= (TF1*)fileDataDCA->Get(funcDcaRightTailString);
      _funcDCARightTailMC[iMassRange][iRedLike] = (TF1*)fileMcDCA->Get(funcDcaRightTailString);

      double xmin;
      double xmax;

      if (_funcDCACentralData[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaCentralString << " does not exist in Data " << std::endl;
	exit(-1);
      }
      if (_funcDCACentralMC[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaCentralString << " does not exist in MC " << std::endl;	
	exit(-1);
      }
      if (_funcDCALeftTailData[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaLeftTailString << " does not exist in Data " << std::endl;
	exit(-1);
      }
      if (_funcDCALeftTailMC[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaLeftTailString << " does not exist in MC " << std::endl;	
	exit(-1);
      }
      if (_funcDCARightTailData[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaRightTailString << " does not exist in Data " << std::endl;
	exit(-1);
      }
      if (_funcDCARightTailMC[iMassRange][iRedLike]==NULL) { 
	std::cout << "  TF1 " << funcDcaRightTailString << " does not exist in MC " << std::endl;	
	exit(-1);
      }

      _funcDCALeftTailData[iMassRange][iRedLike]->GetRange(xmin,xmax);
      _cutoffLeftDcaData[iMassRange][iRedLike] = float(xmax); 
      _minDca[iMassRange][iRedLike] = float(xmin);
      std::cout << "Data Left    : xmin = " << xmin << " xmax = " << xmax << std::endl;

      _funcDCACentralData[iMassRange][iRedLike]->GetRange(xmin,xmax);
      std::cout << "Data Central : xmin = " << xmin << " xmax = " << xmax << std::endl;

      _funcDCARightTailData[iMassRange][iRedLike]->GetRange(xmin,xmax);
      _cutoffRightDcaData[iMassRange][iRedLike] = float(xmin);
      _maxDca[iMassRange][iRedLike] = float(xmax);
      std::cout << "Data Right   : xmin = " << xmin << " xmax = " << xmax << std::endl;

      std::cout << "   left cutoff = " << _funcDCALeftTailData[iMassRange][iRedLike]->Eval(_cutoffLeftDcaData[iMassRange][iRedLike]) 
		<< "   vs. " << _funcDCACentralData[iMassRange][iRedLike]->Eval(_cutoffLeftDcaData[iMassRange][iRedLike]) << std::endl;
      std::cout << "   right cutoff = " << _funcDCACentralData[iMassRange][iRedLike]->Eval(_cutoffRightDcaData[iMassRange][iRedLike])
                << "   vs. " << _funcDCARightTailData[iMassRange][iRedLike]->Eval(_cutoffRightDcaData[iMassRange][iRedLike]) << std::endl;



      _funcDCALeftTailMC[iMassRange][iRedLike]->GetRange(xmin,xmax);
      _cutoffLeftDcaMC[iMassRange][iRedLike] = float(xmax); 
      float xMinMC = float(xmin);
      if (xMinMC>_minDca[iMassRange][iRedLike])
	_minDca[iMassRange][iRedLike] = xMinMC;
      std::cout << "MC Left      : xmin = " << xmin << " xmax = " << xmax << std::endl;

      _funcDCACentralMC[iMassRange][iRedLike]->GetRange(xmin,xmax);
      std::cout << "MC Central   : xmin = " << xmin << " xmax = " << xmax << std::endl;

      _funcDCARightTailMC[iMassRange][iRedLike]->GetRange(xmin,xmax);
      _cutoffRightDcaMC[iMassRange][iRedLike] = float(xmin); 
      float xMaxMC = float(xmax);
      if (xMaxMC<_maxDca[iMassRange][iRedLike])
	_maxDca[iMassRange][iRedLike] = xMaxMC;

      std::cout << "MC Right     : xmin = " << xmin << " xmax = " << xmax << std::endl;

      std::cout << "   left cutoff = " << _funcDCALeftTailMC[iMassRange][iRedLike]->Eval(_cutoffLeftDcaMC[iMassRange][iRedLike])
                << "   vs. " << _funcDCACentralMC[iMassRange][iRedLike]->Eval(_cutoffLeftDcaMC[iMassRange][iRedLike]) << std::endl;
      std::cout << "   right cutoff = " << _funcDCACentralMC[iMassRange][iRedLike]->Eval(_cutoffRightDcaMC[iMassRange][iRedLike])
                << "   vs. " << _funcDCARightTailMC[iMassRange][iRedLike]->Eval(_cutoffRightDcaMC[iMassRange][iRedLike]) << std::endl;

      std::cout << "   min = " << _minDca[iMassRange][iRedLike]
		<< "   cutoffLeft(Data) = " << _cutoffLeftDcaData[iMassRange][iRedLike]
		<< "   cutoffRight(Data) = " << _cutoffRightDcaData[iMassRange][iRedLike]
		<< "   max = " << _maxDca[iMassRange][iRedLike] << std::endl;

      std::cout << "   min = " << _minDca[iMassRange][iRedLike]
		<< "   cutoffLeft(MC) = " << _cutoffLeftDcaMC[iMassRange][iRedLike]
		<< "   cutoffRight(MC) = " << _cutoffRightDcaMC[iMassRange][iRedLike]
		<< "   max = " << _maxDca[iMassRange][iRedLike] << std::endl;

      int nSumProb = 2;
      double q[2];
      double sumProb[2];

      sumProb[0] = 1e-10;
      sumProb[1] = 0.999999;

      _funcDCALeftTailData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesDataLeftMin = float(q[0]);
      float quantilesDataLeftMax = float(q[1]);

      _funcDCACentralData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesDataCenterMin = float(q[0]);
      float quantilesDataCenterMax = float(q[1]);

      _funcDCARightTailData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesDataRightMin = float(q[0]);
      float quantilesDataRightMax = float(q[1]);

      std::cout << "Data quantiles left   = [" << quantilesDataLeftMin << "," << quantilesDataLeftMax << "]" << std::endl;
      std::cout << "Data quantiles center = [" << quantilesDataCenterMin << "," << quantilesDataCenterMax << "]" << std::endl;
      std::cout << "Data quantiles right  = [" << quantilesDataRightMin << "," << quantilesDataRightMax << "]" << std::endl;

      _funcDCALeftTailMC[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesMCLeftMin = float(q[0]);
      float quantilesMCLeftMax = float(q[1]);

      _funcDCACentralMC[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesMCCenterMin = float(q[0]);
      float quantilesMCCenterMax = float(q[1]);

      _funcDCARightTailMC[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      float quantilesMCRightMin = float(q[0]);
      float quantilesMCRightMax = float(q[1]);

      std::cout << "MC quantiles left   = [" << quantilesMCLeftMin << "," << quantilesMCLeftMax << "]" << std::endl;
      std::cout << "MC quantiles center = [" << quantilesMCCenterMin << "," << quantilesMCCenterMax << "]" << std::endl;
      std::cout << "MC quantiles right  = [" << quantilesMCRightMin << "," << quantilesMCRightMax << "]" << std::endl;


      // normalizing integral in data ---->

      _integralLeftDcaData[iMassRange][iRedLike] = _funcDCALeftTailData[iMassRange][iRedLike]->Integral(_minDca[iMassRange][iRedLike],_cutoffLeftDcaData[iMassRange][iRedLike]) ;
      _integralCenterDcaData[iMassRange][iRedLike] = _funcDCACentralData[iMassRange][iRedLike]->Integral(_cutoffLeftDcaData[iMassRange][iRedLike],_cutoffRightDcaData[iMassRange][iRedLike]);
      _integralRightDcaData[iMassRange][iRedLike] = _funcDCARightTailData[iMassRange][iRedLike]->Integral(_cutoffRightDcaData[iMassRange][iRedLike],_maxDca[iMassRange][iRedLike]);


      float integralData = _integralLeftDcaData[iMassRange][iRedLike] + _integralCenterDcaData[iMassRange][iRedLike] + _integralRightDcaData[iMassRange][iRedLike];      

      // normalizing integral in MC ---->

      _integralLeftDcaMC[iMassRange][iRedLike] = _funcDCALeftTailMC[iMassRange][iRedLike]->Integral(_minDca[iMassRange][iRedLike],_cutoffLeftDcaMC[iMassRange][iRedLike]); 
      _integralCenterDcaMC[iMassRange][iRedLike] = _funcDCACentralMC[iMassRange][iRedLike]->Integral(_cutoffLeftDcaMC[iMassRange][iRedLike],_cutoffRightDcaMC[iMassRange][iRedLike]);
      _integralRightDcaMC[iMassRange][iRedLike] = _funcDCARightTailMC[iMassRange][iRedLike]->Integral(_cutoffRightDcaMC[iMassRange][iRedLike],_maxDca[iMassRange][iRedLike]);


      float integralMC = _integralLeftDcaMC[iMassRange][iRedLike] + _integralCenterDcaMC[iMassRange][iRedLike] + _integralRightDcaMC[iMassRange][iRedLike];


      std::cout << "  Data (int) = " << integralData 
		<< "  left = " << _integralLeftDcaData[iMassRange][iRedLike] 
		<< "  central = " << _integralCenterDcaData[iMassRange][iRedLike] 
		<< "  right   = " << _integralRightDcaData[iMassRange][iRedLike] << std::endl;
      
      std::cout << "  MC   (int) = " << integralMC 
		<< "  left = " << _integralLeftDcaMC[iMassRange][iRedLike] 
		<< "  central = " << _integralCenterDcaMC[iMassRange][iRedLike] 
		<< "  right   = " << _integralRightDcaMC[iMassRange][iRedLike] << std::endl;
      std::cout << std::endl;



    }

  }

  std::cout << "DCA Weights are downloaded..." << std::endl;
  //  exit(1);
  
}

void EventWeight::InitTriggerWeights() {


  // HLT_IsoMu17 ------>

  _hltIso17_EtaBins[0] = 0;
  _hltIso17_EtaBins[1] = 0.9;
  _hltIso17_EtaBins[2] = 1.1;
  _hltIso17_EtaBins[3] = 1.5;
  _hltIso17_EtaBins[4] = 2.1;

  _hltIso17_PtBins[0] = 15;
  _hltIso17_PtBins[1] = 20;
  _hltIso17_PtBins[2] = 30;
  _hltIso17_PtBins[3] = 40;
  _hltIso17_PtBins[4] = 60;
  _hltIso17_PtBins[5] = 100;

  //  0 < eta < 0.9;  
  //  _hltIso17_Eff[0][0] = 0.680377;
  //  _hltIso17_Eff[1][0] = 0.917608;
  _hltIso17_Eff[0][0] = 0.710377;
  _hltIso17_Eff[1][0] = 0.917608;
  _hltIso17_Eff[2][0] = 0.898158;
  _hltIso17_Eff[3][0] = 0.897567;
  _hltIso17_Eff[4][0] = 0.873485;

  //  0.9 < eta < 1.1
//   _hltIso17_Eff[0][1] = 0.63601;
//   _hltIso17_Eff[1][1] = 0.873904;
  _hltIso17_Eff[0][1] = 0.65601;
  _hltIso17_Eff[1][1] = 0.893904;
  _hltIso17_Eff[2][1] = 0.870677;
  _hltIso17_Eff[3][1] = 0.869043;
  _hltIso17_Eff[4][1] = 0.844444;

  //  1.1 < eta < 1.5  
  //  _hltIso17_Eff[0][2] = 0.606527;
  _hltIso17_Eff[0][2] = 0.636527;
  _hltIso17_Eff[1][2] = 0.874383;
  _hltIso17_Eff[2][2] = 0.86709;
  _hltIso17_Eff[3][2] = 0.871662;
  _hltIso17_Eff[4][2] = 0.850152;

  // 1.5 < eta < 2.1
  //  _hltIso17_Eff[0][3] = 0.587694;
  _hltIso17_Eff[0][3] = 0.627694;
  _hltIso17_Eff[1][3] = 0.843881;
  _hltIso17_Eff[2][3] = 0.843413;
  _hltIso17_Eff[3][3] = 0.850237;
  _hltIso17_Eff[4][3] = 0.842298;

  _hltIso17_nPtBins = 5;
  _hltIso17_nEtaBins = 4;

}

float EventWeight::hlt_IsoMu17_Efficiency(float pt, float etaIn) {

  // check ranges;

  float eta = TMath::Abs(etaIn);

  if (eta>2.1)
    return 0;

  if (pt<15)
    return 0;

  if (pt>100)
    pt = 99;


  int ptBin = binNumber(pt, _hltIso17_nPtBins, _hltIso17_PtBins);
  int etaBin = binNumber(eta, _hltIso17_nEtaBins, _hltIso17_EtaBins);

  return _hltIso17_Eff[ptBin][etaBin];

}

float EventWeight::hlt_IsoMu17_Weight(float ptPlus, float etaPlus, float ptMinus, float etaMinus) {
  
  float effPlus  = hlt_IsoMu17_Efficiency(ptPlus,etaPlus);
  float effMinus = hlt_IsoMu17_Efficiency(ptMinus,etaMinus);

  float weight = 1 - (1-effPlus)*(1-effMinus);

  return weight;

}

float EventWeight::hlt_Mu13Mu8_Efficiency13(float abseta, float pt, float *efficiency_vec){

  float kAbsEtaBinBoundaries[6] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
  float kPTBinBoundaries[18] = {5.,6.,7.,8.,9., 10.,11.,12.,14.,16., 20.,24.,30.,40.,60., 100.,180.,1000000.};
  float effAbsEtaPT_center[5][17] = {
    {0, 0, 0, 0.0294118, 0, 0, 0, 0.502283, 0.951149, 0.976574, 0.969845, 0.967369, 0.968002, 0.967468, 0.970145, 0.979042, 1, },
    {0, 0, 0, 0, 0, 0, 0, 0.553333, 0.974874, 0.951515, 0.963696, 0.966033, 0.962098, 0.96257, 0.964833, 0.968421, 1, },
    {0, 0, 0, 0, 0, 0, 0, 0.520737, 0.945055, 0.959505, 0.95935, 0.959838, 0.955042, 0.960257, 0.966518, 0.984127, 1, },
    {0, 0, 0, 0, 0, 0.00980392, 0, 0.522901, 0.931624, 0.929236, 0.945445, 0.94578, 0.944493, 0.945942, 0.941795, 0.969388, 0.8, },
    {0, 0, 0.0508475, 0.0363636, 0.0140845, 0.0166667, 0.0645161, 0.43871, 0.820276, 0.872786, 0.878665, 0.874814, 0.879268, 0.887059, 0.88961, 0.875, 1, }
  };



  float effAbsEtaPT_lower[5][17] = {
    {0, 0, 0, 0.00502098, 0, 0, 0, 0.466101, 0.9366, 0.971702, 0.966538, 0.965344, 0.966898, 0.966503, 0.967099, 0.967879, 0.902347, },
    {0, 0, 0, 0, 0, 0, 0, 0.508953, 0.958151, 0.941588, 0.95746, 0.961965, 0.95977, 0.960697, 0.958356, 0.938509, 0.690789, },
    {0, 0, 0, 0, 0, 0, 0, 0.484317, 0.927468, 0.951721, 0.953503, 0.956047, 0.952708, 0.958498, 0.960829, 0.963475, 0.39661, },
    {0, 0, 0, 0, 0, 0.00167647, 0, 0.489935, 0.915271, 0.920473, 0.939578, 0.9418, 0.941976, 0.943863, 0.934311, 0.940363, 0.474545, },
    {0, 0, 0.0232054, 0.0128404, 0.00240757, 0.00284832, 0.0337788, 0.395981, 0.790192, 0.857897, 0.867316, 0.866972, 0.874175, 0.882455, 0.872957, 0.798911, 0.157299, }
  };



  float effAbsEtaPT_higher[5][17] = {
    {0.109169, 0.0883328, 0.0883328, 0.0941216, 0.0378005, 0.0349442, 0.0205676, 0.538443, 0.962767, 0.980681, 0.972849, 0.969283, 0.969071, 0.968407, 0.97293, 0.986772, 1, },
    {0.123759, 0.0617881, 0.0883328, 0.0545069, 0.0597914, 0.0639228, 0.0378005, 0.596942, 0.985728, 0.95994, 0.969103, 0.9697, 0.964301, 0.964359, 0.970398, 0.985619, 1, },
    {0.0441099, 0.030863, 0.0378005, 0.0370435, 0.0231408, 0.0246849, 0.0212774, 0.556953, 0.95893, 0.96616, 0.964528, 0.963328, 0.957269, 0.961946, 0.971454, 0.994404, 1, },
    {0.0171375, 0.0196918, 0.0186965, 0.0186965, 0.0201202, 0.0320924, 0.0158184, 0.55568, 0.945242, 0.937147, 0.950797, 0.949512, 0.94691, 0.947949, 0.948515, 0.98606, 0.966351, },
    {0.0303565, 0.0280544, 0.0980736, 0.0825432, 0.0458772, 0.0541263, 0.11283, 0.482296, 0.847134, 0.886394, 0.889206, 0.882271, 0.884187, 0.891509, 0.904441, 0.928119, 1, }
  };



  int iaeta = 0;
  int ipt = 0;
  for (iaeta = 0; iaeta < 5; ++iaeta){
    if ( (abseta >= kAbsEtaBinBoundaries[iaeta]) && (abseta < kAbsEtaBinBoundaries[iaeta+1]) ) break;
  }
  for (ipt = 0; ipt < 17; ++ipt){
    if ( (pt >= kPTBinBoundaries[ipt]) && (pt < kPTBinBoundaries[ipt+1]) ) break;
  }
  efficiency_vec[0] = effAbsEtaPT_center[iaeta][ipt];
  efficiency_vec[1] = effAbsEtaPT_lower[iaeta][ipt];
  efficiency_vec[2] = effAbsEtaPT_higher[iaeta][ipt];


  return effAbsEtaPT_center[iaeta][ipt];

}

float EventWeight::hlt_Mu13Mu8_Efficiency8(float abseta, float pt, float *efficiency_vec){

   float kAbsEtaBinBoundaries[6] = {0., 0.9, 1.2, 1.6, 2.1, 2.4};
   float kPTBinBoundaries[18] = {5.,6.,7.,8.,9., 10.,11.,12.,14.,16., 20.,24.,30.,40.,60., 100.,180.,1000000.};
   float effAbsEtaPT_center[5][17] = {
      {0, 0, 0, 0.941176, 0.916667, 0.961538, 0.932584, 0.954338, 0.951149, 0.978038, 0.972586, 0.969315, 0.969423, 0.968602, 0.971202, 0.979042, 1, },
      {0, 0, 0.15, 0.939394, 0.866667, 0.892857, 0.9375, 0.953333, 0.974874, 0.95303, 0.966172, 0.968459, 0.964066, 0.966252, 0.966637, 0.978947, 1, },
      {0, 0.0169492, 0.0208333, 0.877551, 0.936709, 0.972973, 0.918605, 0.958525, 0.952381, 0.967379, 0.962737, 0.96264, 0.95766, 0.962542, 0.96875, 0.984127, 1, },
      {0, 0, 0.0816327, 0.897959, 0.945055, 0.931373, 0.913793, 0.938931, 0.948718, 0.94879, 0.958538, 0.955662, 0.951733, 0.952699, 0.942603, 0.969388, 0.8, },
      {0, 0, 0.135593, 0.854545, 0.971831, 0.866667, 0.870968, 0.851613, 0.875576, 0.913043, 0.915066, 0.898565, 0.902923, 0.909531, 0.904762, 0.925, 1, }
   };



   float effAbsEtaPT_lower[5][17] = {
      {0, 0, 0, 0.868463, 0.855369, 0.912817, 0.894297, 0.935394, 0.9366, 0.973291, 0.969413, 0.967346, 0.968342, 0.967653, 0.968203, 0.967879, 0.902347, },
      {0, 0, 0.0692381, 0.864638, 0.773286, 0.799061, 0.880145, 0.928979, 0.958151, 0.943225, 0.960112, 0.964519, 0.961793, 0.964465, 0.960296, 0.951722, 0.690789, },
      {0, 0.00289653, 0.00355914, 0.811253, 0.895959, 0.938268, 0.87731, 0.940064, 0.935648, 0.960235, 0.957098, 0.958966, 0.955389, 0.960831, 0.96322, 0.963475, 0.39661, },
      {0, 0, 0.0537612, 0.857408, 0.909418, 0.896188, 0.87911, 0.920219, 0.933988, 0.941105, 0.953313, 0.952017, 0.94937, 0.950743, 0.935163, 0.940363, 0.474545, },
      {0, 0, 0.0898766, 0.790348, 0.935708, 0.80721, 0.813221, 0.817196, 0.848842, 0.900148, 0.905197, 0.891358, 0.898271, 0.905337, 0.888993, 0.857065, 0.157299, }
   };



   float effAbsEtaPT_higher[5][17] = {
      {0.109169, 0.0883328, 0.0883328, 0.979193, 0.956269, 0.986417, 0.95911, 0.968441, 0.962767, 0.982015, 0.975453, 0.971173, 0.970469, 0.969525, 0.973938, 0.986772, 1, },
      {0.123759, 0.0617881, 0.275289, 0.978559, 0.929594, 0.950792, 0.97144, 0.970457, 0.985728, 0.961327, 0.971395, 0.971997, 0.966212, 0.967954, 0.972059, 0.992575, 1, },
      {0.0441099, 0.0550258, 0.0673332, 0.925224, 0.963883, 0.990463, 0.948259, 0.972022, 0.965311, 0.973365, 0.967701, 0.96601, 0.959824, 0.964184, 0.973522, 0.994404, 1, },
      {0.0171375, 0.0196918, 0.119617, 0.929032, 0.968678, 0.956442, 0.940151, 0.953842, 0.960565, 0.955574, 0.96323, 0.959052, 0.953995, 0.954583, 0.949279, 0.98606, 0.966351, },
      {0.0303565, 0.0280544, 0.195942, 0.903471, 0.990059, 0.911646, 0.914541, 0.880964, 0.898476, 0.924505, 0.924041, 0.90536, 0.90739, 0.913561, 0.918628, 0.96568, 1, }
   };



   int iaeta = 0;
   int ipt = 0;
   for (iaeta = 0; iaeta < 5; ++iaeta){
     if ( (abseta >= kAbsEtaBinBoundaries[iaeta]) && (abseta < kAbsEtaBinBoundaries[iaeta+1]) ) break;
   }
   for (ipt = 0; ipt < 17; ++ipt){
     if ( (pt >= kPTBinBoundaries[ipt]) && (pt < kPTBinBoundaries[ipt+1]) ) break;
   }
   
   //   std::cout << "iaeta = " << iaeta << "  ipt = " << ipt << std::endl;

   efficiency_vec[0] = effAbsEtaPT_center[iaeta][ipt];
   efficiency_vec[1] = effAbsEtaPT_lower[iaeta][ipt];
   efficiency_vec[2] = effAbsEtaPT_higher[iaeta][ipt];

   return effAbsEtaPT_center[iaeta][ipt];

}



float EventWeight::hlt_Mu13Mu8_Efficiency13_Run2011B(float abseta, float pt, float *efficiency_vec){

   float kAbsEtaBinBoundaries[6] = {0, 0.9, 1.2, 1.6, 2.1, 2.4};
   float kPTBinBoundaries[13] = {10, 11, 12, 14, 16, 20, 24, 30, 40, 60, 100, 180, 1000};
   float effAbsEtaPT_center[5][12] = {
      {0, 0, 0.5375, 0.859375, 0.889849, 0.918782, 0.896552, 0.890209, 0.890023, 0.913122, 0.913793, 1, },
      {0, 0.02, 0.436364, 0.926966, 0.887097, 0.927256, 0.91165, 0.912055, 0.912105, 0.906824, 0.924658, 1, },
      {0.0131579, 0.0120482, 0.523585, 0.915789, 0.905218, 0.910602, 0.915565, 0.913768, 0.910808, 0.907358, 0.897059, 1, },
      {0.00826446, 0.015748, 0.526646, 0.889807, 0.877269, 0.900804, 0.90417, 0.897865, 0.893219, 0.890883, 0.906977, 1, },
      {0, 0.0806452, 0.485915, 0.769231, 0.807808, 0.82224, 0.831616, 0.830012, 0.832739, 0.86202, 0.852941, 1, }
   };



   float effAbsEtaPT_lower[5][12] = {
      {0, 0, 0.475103, 0.821467, 0.873229, 0.909991, 0.891781, 0.887603, 0.887843, 0.906675, 0.891099, 0.690789, },
      {0, 0.00341701, 0.385159, 0.901767, 0.87096, 0.918448, 0.906216, 0.909133, 0.909828, 0.898727, 0.895868, 0.62977, },
      {0.00224935, 0.00205984, 0.486705, 0.895863, 0.894614, 0.903199, 0.911387, 0.911362, 0.908938, 0.900809, 0.871039, 0.690789, },
      {0.00141341, 0.00555213, 0.496927, 0.870735, 0.866794, 0.89404, 0.900091, 0.895333, 0.891176, 0.88364, 0.879112, 0.62977, },
      {0, 0.0461176, 0.440499, 0.737204, 0.791223, 0.810611, 0.824225, 0.825501, 0.828686, 0.847685, 0.765272, 0.39661, }
   };



   float effAbsEtaPT_higher[5][12] = {
      {0.0976528, 0.0741721, 0.598855, 0.89094, 0.904648, 0.926821, 0.90114, 0.892762, 0.892165, 0.919179, 0.932484, 1, },
      {0.0394112, 0.0647021, 0.488838, 0.946657, 0.901552, 0.935225, 0.916808, 0.914894, 0.91433, 0.91436, 0.946701, 1, },
      {0.0429048, 0.0393367, 0.560226, 0.93247, 0.914905, 0.91751, 0.919567, 0.916116, 0.912644, 0.91353, 0.918719, 1, },
      {0.0271012, 0.0362409, 0.556187, 0.906538, 0.887058, 0.907195, 0.908102, 0.900343, 0.89523, 0.897741, 0.92947, 1, },
      {0.0293918, 0.131771, 0.531548, 0.798615, 0.823416, 0.833324, 0.838764, 0.834432, 0.836718, 0.875261, 0.915194, 1, }
   };



   int iaeta = 0;
   int ipt = 0;
   for (iaeta = 0; iaeta < 5; ++iaeta){
      if ( (abseta >= kAbsEtaBinBoundaries[iaeta]) && (abseta < kAbsEtaBinBoundaries[iaeta+1]) ) break;
      }
   for (ipt = 0; ipt < 12; ++ipt){
      if ( (pt >= kPTBinBoundaries[ipt]) && (pt < kPTBinBoundaries[ipt+1]) ) break;
      }
   efficiency_vec[0] = effAbsEtaPT_center[iaeta][ipt];
   efficiency_vec[1] = effAbsEtaPT_lower[iaeta][ipt];
   efficiency_vec[2] = effAbsEtaPT_higher[iaeta][ipt];
   return  effAbsEtaPT_center[iaeta][ipt];
}


float EventWeight::hlt_Mu13Mu8_Efficiency8_Run2011B(float abseta, float pt, float *efficiency_vec){

   float kAbsEtaBinBoundaries[6] = {0, 0.9, 1.2, 1.6, 2.1, 2.4};
   float kPTBinBoundaries[13] = {10, 11, 12, 14, 16, 20, 24, 30, 40, 60, 100, 180, 1000};
   float effAbsEtaPT_center[5][12] = {
      {0.888889, 0.958333, 0.9125, 0.882812, 0.920086, 0.947547, 0.937041, 0.933653, 0.935177, 0.947964, 0.948276, 1, },
      {0.978261, 0.86, 0.927273, 0.949438, 0.931452, 0.951197, 0.953398, 0.954644, 0.954562, 0.952756, 0.979452, 1, },
      {0.934211, 0.975904, 0.95283, 0.947368, 0.957401, 0.944986, 0.961546, 0.953647, 0.955894, 0.95523, 0.97549, 1, },
      {0.92562, 0.937008, 0.924765, 0.955923, 0.923077, 0.948615, 0.953224, 0.943002, 0.942207, 0.937648, 0.953488, 1, },
      {0.83871, 0.919355, 0.880282, 0.868778, 0.882883, 0.881494, 0.884547, 0.880981, 0.883853, 0.900427, 0.882353, 1, }
   };



   float effAbsEtaPT_lower[5][12] = {
      {0.760127, 0.868534, 0.868332, 0.84695, 0.905337, 0.940192, 0.933185, 0.931563, 0.933451, 0.942771, 0.929228, 0.690789, },
      {0.929813, 0.792424, 0.89319, 0.927086, 0.918072, 0.943709, 0.949268, 0.952468, 0.95287, 0.946642, 0.959771, 0.62977, },
      {0.891946, 0.944854, 0.933283, 0.930495, 0.949704, 0.938938, 0.958587, 0.951826, 0.954535, 0.950415, 0.959168, 0.690789, },
      {0.893336, 0.907258, 0.906854, 0.942258, 0.914376, 0.943485, 0.950242, 0.941047, 0.940652, 0.931908, 0.931237, 0.62977, },
      {0.777605, 0.868229, 0.84637, 0.841855, 0.86899, 0.871503, 0.878171, 0.877069, 0.880353, 0.887749, 0.798661, 0.39661, }
   };



   float effAbsEtaPT_higher[5][12] = {
      {0.960534, 0.992894, 0.944337, 0.911858, 0.932853, 0.954089, 0.940695, 0.935684, 0.936862, 0.952729, 0.962867, 1, },
      {0.996286, 0.910364, 0.952154, 0.965857, 0.942925, 0.95779, 0.95722, 0.956728, 0.956198, 0.958235, 0.990653, 1, },
      {0.962446, 0.991499, 0.967394, 0.960668, 0.964039, 0.950493, 0.964309, 0.955404, 0.957216, 0.959619, 0.986079, 1, },
      {0.949637, 0.958605, 0.939714, 0.966739, 0.931, 0.953322, 0.95604, 0.944898, 0.943725, 0.942952, 0.969493, 1, },
      {0.887052, 0.953882, 0.908129, 0.892016, 0.895549, 0.890833, 0.89064, 0.884787, 0.887265, 0.911886, 0.938, 1, }
   };



   int iaeta = 0;
   int ipt = 0;
   for (iaeta = 0; iaeta < 5; ++iaeta){
      if ( (abseta >= kAbsEtaBinBoundaries[iaeta]) && (abseta < kAbsEtaBinBoundaries[iaeta+1]) ) break;
      }
   for (ipt = 0; ipt < 12; ++ipt){
      if ( (pt >= kPTBinBoundaries[ipt]) && (pt < kPTBinBoundaries[ipt+1]) ) break;
      }
   efficiency_vec[0] = effAbsEtaPT_center[iaeta][ipt];
   efficiency_vec[1] = effAbsEtaPT_lower[iaeta][ipt];
   efficiency_vec[2] = effAbsEtaPT_higher[iaeta][ipt];
   return  effAbsEtaPT_center[iaeta][ipt];
}






float EventWeight::hlt_Mu13Mu8_Weight(float ptPlus, float etaPlus, float ptMinus, float etaMinus) {

  float absEtaPlus = TMath::Abs(etaPlus);
  float absEtaMinus = TMath::Abs(etaMinus);

  if (absEtaPlus>2.4)
    return 0;
  
  if (absEtaMinus>2.4)
    return 0;

  float  effVPlus13[3];
  float  effVPlus8[3];
  float  effVPMinus13[3];
  float  effVMinus8[3];

  //  std::cout << "absEtaPlus = " << absEtaPlus << "   ptPlus = " << ptPlus << std::endl;
  //  std::cout << "absEtaMinus = " << absEtaMinus << "   ptMinus = " << ptMinus << std::endl;

  float effPlus13 = 1;  
  float effPlus8  = 1;

  float effMinus13 = 1;
  float effMinus8  = 1;

  int dataset = 1;

  if (dataset==1) {

    effPlus13 = hlt_Mu13Mu8_Efficiency13(absEtaPlus,ptPlus,effVPlus13);
    effPlus8  = hlt_Mu13Mu8_Efficiency8(absEtaPlus,ptPlus,effVPlus8);

    effMinus13 = hlt_Mu13Mu8_Efficiency13(absEtaMinus,ptMinus,effVPMinus13);
    effMinus8  = hlt_Mu13Mu8_Efficiency8(absEtaMinus,ptMinus,effVMinus8);

  }
  else {

    effPlus13 = hlt_Mu13Mu8_Efficiency13_Run2011B(absEtaPlus,ptPlus,effVPlus13);
    effPlus8  = hlt_Mu13Mu8_Efficiency8_Run2011B(absEtaPlus,ptPlus,effVPlus8);

    effMinus13 = hlt_Mu13Mu8_Efficiency13_Run2011B(absEtaMinus,ptMinus,effVPMinus13);
    effMinus8  = hlt_Mu13Mu8_Efficiency8_Run2011B(absEtaMinus,ptMinus,effVMinus8);

  } 

  //  std::cout << "effPlus8   = " << effPlus8 << std::endl;
  //  std::cout << "effPlus13  = " << effPlus13 << std::endl;
  //  std::cout << "effMinus8  = " << effMinus8 << std::endl;
  //  std::cout << "effMinus13 = " << effMinus13 << std::endl;

  
  float weight = effPlus13*effMinus8 + effPlus8*effMinus13 - effPlus13*effMinus13;

  return weight;

}

float EventWeight::hlt_53x_TrigScale(float pt,float eta) {

  float absEta = TMath::Abs(eta);
  float scale = 1;

  if (pt<15.0&&pt>=10) {
    if (absEta<0.8) {
      scale = 0.9818;
    }
    else if (absEta>=0.8&&absEta<1.2) {
      scale = 0.9713;
    }
    else {
      scale = 0.9675;
    }
  }
  else if (pt<20&&pt>=15) {
    if (absEta<0.8) {
      scale = 0.9781;
    }
    else if (absEta>=0.8&&absEta<1.2) {
      scale = 0.9782;
    }
    else {
      scale = 0.9587;
    }

  }
  else if (pt<25&&pt>=20) {
    if (absEta<0.8) {
      scale = 0.9873;
    }
    else if (absEta>=0.8&&absEta<1.2) {
      scale = 0.9532;
    }
    else {
      scale = 0.9605;
    }
  }
  else if (pt<30&&pt>=25) {
    if (absEta<0.8) {
      scale = 0.9755;
    }
    else if (absEta>=0.8&&absEta<1.2) {
      scale = 0.9818;
    }
    else {
      scale = 0.9632;
    }
  }
  else {
    if (absEta<0.8) {
      scale = 0.9956;
    }
    else if (absEta>=0.8&&absEta<1.2) {
      scale = 0.9644;
    }
    else {
      scale = 0.9530;
    }
  }

  return scale;


}

float EventWeight::RhoWeight(float rho, int nPV, float upper) {

  float weight = 1;

  int nPVIndex = nPV - 1;

  if (nPVIndex<_nPVMax) {

    if (rho<upper) {

      float numerator = 1;
      float denominator = 1;

      if (rho<0.2) {
	numerator   = _rhoDataBin1->GetBinContent(nPV);
	denominator = _rhoMCBin1->GetBinContent(nPV);
      }
      else if (rho>=0.2&&rho<0.4) {
	numerator   = _rhoDataBin2->GetBinContent(nPV);
	denominator = _rhoMCBin2->GetBinContent(nPV);
	}
      else if (rho<upper) {
	numerator   = _rhoFuncData.at(nPVIndex)->Eval(rho);
	denominator = _rhoFuncMC.at(nPVIndex)->Eval(rho);
      }

      if ( (numerator>0) && (denominator>0) )
	weight = numerator / denominator;

      // do not apply lower cut on weight (only upper cut)
      if (weight>_weightMax) weight = _weightMax;

    }

  }

  return weight;
}

// for reweighting MC to different MC generators/showerings/...
float EventWeight::TheoryWeight(float deta, float mjj, float shift) {
  float weight = 1.0;
  if (shift == 0) // without a shift, this function should not be called, if it is, the weight is 1
    return weight;
  else if (mjj > 0.0 && mjj < 2000.0 && deta > 0.0 && deta < 9.0)
  // outside this range no weights are defined
  {
    if (_theoryWeightHisto == NULL)
      std::cout << "No weight histogram given" << std::endl;
    int thebin = _theoryWeightHisto->FindBin(deta, mjj);
    if (shift > 0.2)
      weight = _theoryWeightHisto->GetBinContent(thebin);
   	else 
      std::cout << "The interface has changed, no more down shifts!" << std::endl;
  }
  else
  {
    weight = 1.0;
  }
  return weight;
}

float EventWeight::DiMuEtaWeight(float eta) {

  float weight = 1;

  float absEta = TMath::Abs(eta);

  if (absEta<_etaMax)
    weight = _diLepEtaWeightFunc->Eval(eta);

  return weight;

}

float EventWeight::MetWeight(float U1, float U2, float Met, float Zpt, int njets) {

  float weight = 1;

  int jetBin = njets;

  if (jetBin>=_nJetsBins)
    jetBin = _nJetsBins - 1;

  int ZptBin = binNumber(Zpt, _ZPtBins);

  //  std::cout << "ZPt   = " << Zpt << "  bin = " << ZptBin << std::endl;
  //  std::cout << "NJets = " << njets << "  bin = " << jetBin << std::endl;

  TF1 * metZParalData = _metZParalData[ZptBin][jetBin];
  TF1 * metZPerpData  = _metZPerpData[ZptBin][jetBin];
  
  TF1 * metZParalMC   = _metZParalMC[ZptBin][jetBin];
  TF1 * metZPerpMC     = _metZPerpMC[ZptBin][jetBin];

  if (U1>_xminMetZParal[ZptBin][jetBin]&&U1<_xmaxMetZParal[ZptBin][jetBin]) {

    float numerator   = metZParalData->Eval(U1);
    float denominator = metZParalMC->Eval(U1);
  
    if (denominator>0&&numerator>0)
      weight *= numerator/denominator;
  }

  
  if (U2>_xminMetZPerp[ZptBin][jetBin]&&U2<_xmaxMetZPerp[ZptBin][jetBin]) {
    float numerator   = metZPerpData->Eval(U2);
    float denominator = metZPerpMC->Eval(U2);

    if (denominator>0&&numerator>0)
      weight *= numerator/denominator;

  }

  // additional weights --->

  if (_dataset==0) { // Run2011A SingleMu
    if (jetBin==0) {
      if (Met>28&&Met<=40)
	weight *= 0.957;
      if (Met>40&&Met<=52)
	weight *= 0.879;
      if (Met>52&&Met<=64)
	weight *= 0.523;
    }
    else if (jetBin==1) {
      if (Met>40&&Met<=48)
	weight *= 0.878;
      if (Met>48&&Met<=60)
	weight *= 0.617;
    }
    else {
      if (Met>40&&Met<=48)
	weight *= 0.643;
    }
  } // *** Run2011A SingleMu


  if (_dataset==1) { // Run2011A DoubleMu
    if (jetBin==0) {
      if (Met>36&&Met<=44)
	weight *= 0.889;
      if (Met>48&&Met<=60)
	weight *= 0.677;
    }
    if (jetBin==1) {
      if (Met>40&&Met<=60)
	weight *= 0.787;
    }
  } // *** Run2011A DoubleMu

  if (_dataset==2) { // Run2011B 

    if (jetBin==0) {
      if (Met>40&&Met<=56)
	weight *= 0.939;
      if (Met>68&&Met<=80)
	weight *= 2.624;
    }
    if (jetBin==1) {
      if (Met>48&&Met<=56)
	weight *= 0.794;
      if (Met>56&&Met<68)
	weight *= 1.364;
    }


  } // *** Run2011B


  // limits on weights --->
  if (weight>20) 
    weight = 20;

  if (weight<0.02) 
    weight = 0.02;

  return weight;


}


float EventWeight::DCACorrected(float dca,
				float diLepMass,
				float reducedLike) {

  float dcaCorrected = dca;

  int iMassRange = binNumber(diLepMass, _diLepMassBins);
  int iRedLike = binNumber(reducedLike, _reducedLikeBins);

  if (dca>_minDca[iMassRange][iRedLike] && dca<_maxDca[iMassRange][iRedLike]) {

    float integralMC = 0;

    if (dca<_cutoffLeftDcaMC[iMassRange][iRedLike]) {
      integralMC = _funcDCALeftTailMC[iMassRange][iRedLike]->Integral(_minDca[iMassRange][iRedLike],dca);
    }
    else if (dca>_cutoffLeftDcaMC[iMassRange][iRedLike]&&dca<_cutoffRightDcaMC[iMassRange][iRedLike]) {
      integralMC = _integralLeftDcaMC[iMassRange][iRedLike] 
	+ _funcDCACentralMC[iMassRange][iRedLike]->Integral(_cutoffLeftDcaMC[iMassRange][iRedLike],dca);
    }
    else {
      integralMC = _integralLeftDcaMC[iMassRange][iRedLike] 
	+ _integralCenterDcaMC[iMassRange][iRedLike] 
	+ _funcDCARightTailMC[iMassRange][iRedLike]->Integral(_cutoffRightDcaMC[iMassRange][iRedLike],dca);
    }

    int nSumProb = 1; 
    double q[1];
    double sumProb[1];

    float integralDataRightTail = _integralLeftDcaData[iMassRange][iRedLike] +
      _integralCenterDcaData[iMassRange][iRedLike];

    if (integralMC<_integralLeftDcaData[iMassRange][iRedLike]) {
      sumProb[0] = integralMC/_integralLeftDcaData[iMassRange][iRedLike];
      _funcDCALeftTailData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      dcaCorrected = float(q[0]);
    }
    else if (integralMC>=_integralLeftDcaData[iMassRange][iRedLike]&&integralMC<integralDataRightTail) {
      sumProb[0] = integralMC - _integralLeftDcaData[iMassRange][iRedLike];
      sumProb[0] /= _integralCenterDcaData[iMassRange][iRedLike];
      _funcDCACentralData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      dcaCorrected = float(q[0]);
    }
    else {
      sumProb[0] = integralMC - integralDataRightTail;
      sumProb[0] /= _integralRightDcaData[iMassRange][iRedLike];
      _funcDCARightTailData[iMassRange][iRedLike]->GetQuantiles(nSumProb,q,sumProb);
      dcaCorrected = float(q[0]);
    }
  }

  return dcaCorrected;


}

float EventWeight::RecoilCorrected(float & MetPx,
				   float & MetPy,
				   float genZPx, 
				   float genZPy,
				   float diLepPx,
				   float diLepPy,
				   int njets,
				   int method) {
  
  // input parameters
  // MetPx, MetPy - missing transverse momentum 
  //                ( corrected replaces uncorrected )
  // genZPx, genZPy - generated transverse momentum of Z
  // diLepPx, diLepPy - dilepton transverse momentum 
  // njets - number of jets 
  // method : 2 - corrections by width w(MC)=w(Data)/w(MC) w(Process)
  // method : 3 - corrections by sampling 
  //              ( calculations of quantiles )

  float Zpt = TMath::Sqrt(genZPx*genZPx + genZPy*genZPy);

  HTauTauUtils utils;

  float U1 = 0;
  float U2 = 0;
  float metU1 = 0;
  float metU2 = 0;

  float weight = 1;


  utils.CalculateU1U2FromMet(MetPx,
			     MetPy,
			     genZPx,
			     genZPy,
			     diLepPx,
			     diLepPy,
			     U1,
			     U2,
			     metU1,
			     metU2);
  if (Zpt>1000)
    Zpt = 999;

  if (njets>=_nJetsBins)
    njets = _nJetsBins - 1;

  int ZptBin = binNumber(Zpt, _ZPtBins);

  if (method==3) {

    TF1 * metZParalData = _metZParalData[ZptBin][njets];
    TF1 * metZPerpData  = _metZPerpData[ZptBin][njets];
  
    TF1 * metZParalMC   = _metZParalMC[ZptBin][njets];
    TF1 * metZPerpMC     = _metZPerpMC[ZptBin][njets];

    if (U1>_xminMetZParal[ZptBin][njets]&&U1<_xmaxMetZParal[ZptBin][njets]) {

      int nSumProb = 1;
      double q[1];
      double sumProb[1];

      sumProb[0] = metZParalMC->Integral(_xminMetZParalMC[ZptBin][njets],U1);

      if (sumProb[0]<0) {
	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 1e-10;
      }
       if (sumProb[0]>1) {
	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 0.9999999;
      }
     
      
      metZParalData->GetQuantiles(nSumProb,q,sumProb);

      float U1reco = float(q[0]);

      if (U1reco>_xminMetZParal[ZptBin][njets]&&U1reco<_xmaxMetZParal[ZptBin][njets]) 
	U1 = U1reco;

    }

    if (U2>_xminMetZPerp[ZptBin][njets]&&U2<_xmaxMetZPerp[ZptBin][njets]) {

      int nSumProb = 1;
      double q[1];
      double sumProb[1];
      
      sumProb[0] = metZPerpMC->Integral(_xminMetZPerpMC[ZptBin][njets],U2);
      
      if (sumProb[0]<0) {
	//	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 1e-10;
      }
       if (sumProb[0]>1) {
	 //	std::cout << "Warning ! ProbSum[0] = " << sumProb[0] << std::endl;
	sumProb[0] = 0.9999999;
      }

      metZPerpData->GetQuantiles(nSumProb,q,sumProb);

      float U2reco = float(q[0]);
 
      if (U2reco>_xminMetZPerp[ZptBin][njets]&&U2reco<_xmaxMetZPerp[ZptBin][njets]) 
	U2 = U2reco;
      

    }



  }
  else 
    U1U2CorrectionsByWidth(U1,U2,ZptBin,njets);
  
  utils. CalculateMetFromU1U2(U1,U2,genZPx,genZPy,diLepPx,diLepPy,MetPx,MetPy);


  if (method==3) {

    float MetX = TMath::Sqrt(MetPx*MetPx+MetPy*MetPy);
    
    if (_dataset==0) { // Run2011A SingleMu

      if (njets==0) {
	if (MetX>52&&MetX<=60)
	  weight *= 0.634;
      }
      
      if (njets==2) {
	if (MetX>40&&MetX<=60)
	  weight *= 0.752;
      }

    } // *** Run2011A SingleMu


    if (_dataset==1) { // Run2011A DoubleMu

      if (njets==0) {
	if (MetX>52&&MetX<=64)
	  weight *= 0.706; 
      }


    } // *** Run2011A DoubleMu

    if (_dataset==2) { // Run2011B 

      if (njets==0) {
	if (MetX>56&&MetX<=72) 
	  weight *= 1.271;
	if (MetX>72&&MetX<=80) 
	  weight *= 6.236;
      }

      if (njets==1) {
	if (MetX>56&&MetX<=72) 
	  weight *= 1.596;
      }

      if (njets==2) {
	if (MetX>68&&MetX<=72)
	  weight *= 0.1;
      }

    } // *** Run2011B 

  }  
  
  return weight;


}

float EventWeight::CorrectionsBySampling(float x, TF1 * funcMC, TF1 * funcData) {

  int nSumProb = 1;
  double q[1];
  double sumProb[1];
  
  double xD = double(x);

  double xminD = 0;
  double xmaxD = 0;

  funcMC->GetRange(xminD,xmaxD);

  float xmin = float(xminD);

  sumProb[0] = funcMC->Integral(xmin,xD);
  
  funcData->GetQuantiles(nSumProb,q,sumProb);

  float output = float(q[0]);

  return output;

}

void EventWeight::U1U2CorrectionsByWidth(float & U1, 
					 float & U2,
					 int ZptBin,
					 int njets) {

  if (njets>=_nJetsBins)
    njets = _nJetsBins - 1;
    
  // ********* U1 *************

  float width = U1 - _meanMetZParalMC[ZptBin][njets];

  if (width<0)
    width *= _rmsLeftMetZParalData[ZptBin][njets]/_rmsLeftMetZParalMC[ZptBin][njets];
  else
    width *= _rmsRightMetZParalData[ZptBin][njets]/_rmsRightMetZParalMC[ZptBin][njets];

  U1 = _meanMetZParalData[ZptBin][njets] + width;

  // ********* U2 *************

  width = U2 - _meanMetZPerpMC[ZptBin][njets];


  if (width<0)
    width *= _rmsLeftMetZPerpData[ZptBin][njets]/_rmsLeftMetZPerpMC[ZptBin][njets];
  else 
    width *= _rmsRightMetZPerpData[ZptBin][njets]/_rmsRightMetZPerpMC[ZptBin][njets];

  U2 = _meanMetZPerpData[ZptBin][njets] + width;

}

float EventWeight::DiMuPtWeight(float diMuonPtRatio) {

  float weight = 1;

  if (_dataset==0) {
    if ( diMuonPtRatio>0.84 )
      weight = 0.758;
    
  }
  if (_dataset==1) {
    if (diMuonPtRatio>0.72 && diMuonPtRatio<0.84)
      weight = 0.90;
  }

  if (_dataset==3) {
    if (diMuonPtRatio>0.6)
      weight = 0.85;
  }

  return weight;
}

float EventWeight::DCAWeight(float diLepMass, float reducedLike, float dca) {

  float weight = 1;

  //  std::cout << "DCA = " << dca << "  Min = " << _minDca << "  Max = " << _maxDca << std::endl;

  int iMassRange = binNumber(diLepMass, _diLepMassBins);
  int iRedLike = binNumber(reducedLike, _reducedLikeBins);

  if (dca>=_minDca[iMassRange][iRedLike] && dca<_maxDca[iMassRange][iRedLike]) {
    
    float numerator = 1;
    float denominator = 1;

    if (dca<_cutoffLeftDcaData[iMassRange][iRedLike]) {
      numerator = _funcDCALeftTailData[iMassRange][iRedLike]->Eval(dca);
    }
    else if (dca>=_cutoffLeftDcaData[iMassRange][iRedLike]&&dca<_cutoffRightDcaData[iMassRange][iRedLike]) {
      numerator = _funcDCACentralData[iMassRange][iRedLike]->Eval(dca);
    }
    else {
      numerator = _funcDCARightTailData[iMassRange][iRedLike]->Eval(dca);
    }
    
    if (dca<_cutoffLeftDcaMC[iMassRange][iRedLike]) {
      denominator = _funcDCALeftTailMC[iMassRange][iRedLike]->Eval(dca);
    }
    else if (dca>=_cutoffLeftDcaMC[iMassRange][iRedLike]&&dca<_cutoffRightDcaMC[iMassRange][iRedLike]) {
      denominator = _funcDCACentralMC[iMassRange][iRedLike]->Eval(dca);
    }
    else {
      denominator = _funcDCARightTailMC[iMassRange][iRedLike]->Eval(dca);
    }

    //    std::cout << "Numerator = " << numerator << "  Denominator = " << denominator << std::endl;

    weight = numerator / denominator;

    if (weight>_weightMax)
      weight = _weightMax;
    
    if (weight<_weightMin)
      weight = _weightMin;
  
    
  }

  //  if (reducedLike>0.5&&reducedLike<0.75) {
  //    if (diLepMass>70) {
  //      if (dca>0.0&&dca<0.6) {
  //	std::cout << "iDiLepMass = " << iMassRange << " RedLike = " << reducedLike << " iRedLike = " << iRedLike
  //		  << "  dca = " << dca << " weight = " << weight << std::endl;
  //      }
  //    }
  //  }

	if (_dataset!=3)
		weight *=  _totWeightDCA[iMassRange][iRedLike];

  return weight;

}

float EventWeight::DCAWeightZTauTau(float dca, float diLepMass, float reducedLike) {

  float weight = 1;

  if (diLepMass<70) {
    
    if (_dataset==0) { // Run2011A SingleMu

      if (reducedLike<=0.25) {
	if (dca>0.0&&dca<=0.2)
	  weight = 0.899;
	if (dca>0.8&&dca<=1.0) 
	  weight = 1.5;
      }

      if (reducedLike<=0.5&&reducedLike>0.25) {
	if (dca>0.0&&dca<=0.2)
	  weight = 1.099;
	if (dca>0.8&&dca<=1.0) 
	  weight = 0.6;
      }

      if (reducedLike<=0.75&&reducedLike>0.5) {
	if (dca>0.0&&dca<=0.2)
	  weight = 1.067;
	if (dca>0.8&&dca<=1.0) 
	  weight = 0.75;
      }

    } // Run2011A SingleMu


    if (_dataset==1) { // Run2011A DoubleMu

      
      if (reducedLike<=0.5&&reducedLike>0.25) {
	if (dca>-0.6&&dca<=-0.4)
	  weight = 1.074;
	if (dca>0.4&&dca<=0.6) 
	  weight = 1.068;
        if (dca>0.8&&dca<=1.0) 
          weight = 0.8;
	if (dca>1.0&&dca<=1.2)
	  weight = 0.3;
      }




    } // ******** Run2011A DoubleMu

    if (_dataset==2) { // Run2011B

      if (reducedLike<=0.25) {
	if (dca>0.0&&dca<=0.2)
	  weight = 1.071;
	if (dca>0.8&&dca<=1.0) 
	  weight = 0.7;
      }

      if (reducedLike<=0.5&&reducedLike>0.25) {
	if (dca>0.0&&dca<=0.2) 
	  weight = 1.068;
        if (dca>0.8&&dca<=1.0) 
          weight = 0.7;
      }

      if (reducedLike<=0.75&&reducedLike>0.5) {
	if (dca>0.0&&dca<=0.2)
	  weight = 1.177;
	if (dca>0.4&&dca<=0.8)
	  weight = 0.7;
      }

    } // Run2011B

  }

  return weight;

}

float EventWeight::getggHWeightHqT(float pT)
{
	if (ggHWeightHqT_histo)
		return std::max<float>(ggHWeightHqT_histo->GetBinContent(ggHWeightHqT_histo->GetBin(pT)), 0.);
	else
		return 1.0;
}


float EventWeight::getggHWeightFeHiPro(float pT)
{
	if (ggHWeightFeHiPro_histo)
		return std::max<float>(ggHWeightFeHiPro_histo->Interpolate(pT), 0.);
	else
		return 1.0;
}


void EventWeight::InitZttEmbedDataWeight(TString ZttEmbedDataFileName) {

  TString tmpFilename = _baseDir+"/"+ZttEmbedDataFileName;
  ZttEmbedDataWeightFile = TFile::Open(tmpFilename);

  if (!ZttEmbedDataWeightFile) {
    std::cout << "File " << ZttEmbedDataFileName << " is not found in directory " <<  _baseDir << std::endl;
    exit(-1);
  }
}

float EventWeight::ZttEmbedDataWeight(float value, int catIdx) {

  float weight = 1.;

  if (catIdx == 0)      // 0 jets low pT
    ZttEmbedDataWeightHisto = (TH1F*)ZttEmbedDataWeightFile->Get("correctionFactor0JetSoft");
  else if (catIdx == 1) // 0 jets high pT
    ZttEmbedDataWeightHisto = (TH1F*)ZttEmbedDataWeightFile->Get("correctionFactor0JetHard");
  else if (catIdx == 2) // 1 jet low pT
    ZttEmbedDataWeightHisto = (TH1F*)ZttEmbedDataWeightFile->Get("correctionFactor1JetSoft");
  else if (catIdx == 3) // 1 jet high pT
    ZttEmbedDataWeightHisto = (TH1F*)ZttEmbedDataWeightFile->Get("correctionFactor1JetHard");
  else if (catIdx == 5) // VBF
    ZttEmbedDataWeightHisto = (TH1F*)ZttEmbedDataWeightFile->Get("correctionFactorVbf");
  else { // no weights available for the category or no category at all
    //std::cout << "Embedded data weight set to 1" << std::endl;
    return weight;
  }

  int bin = ZttEmbedDataWeightHisto->FindBin(value);
  weight  = ZttEmbedDataWeightHisto->GetBinContent(bin);

  return weight;
}

