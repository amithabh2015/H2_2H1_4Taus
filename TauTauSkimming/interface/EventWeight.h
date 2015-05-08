#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2F.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <assert.h>

class EventWeight {
  
 public:
  EventWeight(TString baseDir, int dataset=-1);
  ~EventWeight();

  void InitTriggerWeights();

	void InitMEtWeights(TString dataFile,
		TString MCFile,
		int nZPtBins,
		float * ZPtBins,
		TString * histPerpNames,
		TString * histParalNames,
		int nJetsBins,
		TString * histNJetsNames);

	void InitMEtWeights(TString dataFile,
		TString MCFile,
		//int nZPtBins,
		//float * ZPtBins,
		const std::vector<float>& ZPtBins,
		const std::vector<std::string>& _perpZStr,
		const std::vector<std::string>& _paralZStr,
		const std::vector<std::string>& _nJetsStr);

  bool InitggHWeightFeHiPro(int higgsMass);
  bool InitggHWeightHqT(int higgsMass, float renormScale, float factScale);

  void InitDiMuEtaWeight(TString etaWeightsFile, float etaMax);

  void InitRhoWeights(TString rhoWeightsFile, int nPVMax);

  void InitTheoryWeights(TString theoryWeightsFileName, TString theoryWeightHisto);

  void InitDCAWeights(TString dcaFuncDataFile,
		      TString dcaFuncMCFile,
		      int nDiLepMassBins,
		      float * diLepMassBins,
		      TString * massRangeTString,
		      int nReducedLikeBins,
		      float * reducedLikeBins,
		      TString * likeCutTString);

  void InitDCAWeights(TString dcaFuncDataFile,
		      TString dcaFuncMCFile,
		      const std::vector<float>& diLepMassBins,
		      const std::vector<std::string>& massRangeTString,
		      const std::vector<float>& reducedLikeBins,
		      const std::vector<std::string>& likeCutTString);

  void InitDCAWeights2012();

  float dcaWeight2012(float dca);

  void InitPUWeightsTruth(TString puFileName);

  bool InitPUWeightsS6Truth(TString puDataFile,
			    TString puFall11MCFile);

  bool InitPUWeightsS7Truth(TString puDataFile,
			    TString puFall11MCFile);


  bool InitPUWeightsS10Truth(TString puDataFile,
			     TString puMCFile);

  float PUWeightTruth(int nPUI);
  float PUWeightS6Truth(float mean);
  float PUWeightS7Truth(float mean);
  float PUWeightS10Truth(float mean);

  float RhoWeight(float rho, int nPV, float upper);

  float TheoryWeight(float deta, float mjj, float shift = 0.0);

  float DiMuEtaWeight(float eta);

  float DCAWeightZTauTau(float dca, float diLepMass, float reducedLike);

  float RecoilCorrected(float & MetPx, 
			float & MetPy, 
			float genZpt, 
			float genZphi, 
			float diLepPt,
			float diLepEta,
			int njets,
			int method);

  float DCACorrected(float dca,
		     float diLepMass,
		     float redLike);

  float DiMuPtWeight(float pt);

  float DCAWeight(float diLepMass, float redLike, float dca);

  float MetWeight(float U1,
		  float U2,
		  float Met,
		  float Zpt,
		  int njets);

  float getggHWeightHqT(float pT);
  float getggHWeightFeHiPro(float pT);

  float hlt_IsoMu17_Weight(float ptPlus, float etaPlus, float ptMinus, float etaMins);
  float hlt_Mu13Mu8_Weight(float ptPlus, float etaPlus, float ptMinus, float etaMins);
  float hlt_53x_TrigScale(float pt, float eta);

  float PUWeightS3(int nPV);
  float PUWeightS4(int nPV);

  float PUWeightS7(int nPV);
  float PUWeightS6(int nPV);

  float MuEff2012(float pt, float eta);
  float histIntegral(TH1F * h, float x);
  float doIsoMapping(TH1F * dataH, TH1F * mcH, float x);

  void InitIsoMapping(TString dir,
		      TString dataFileName,
		      TString mcFileName,
		      std::vector<TString> histNames);

  float isoMapping(int iHist, float x);

  // Iso mapping histos;
  int nHistosIsoMap;
  std::vector<TH1F*> histIsoMapData;
  std::vector<TH1F*> histIsoMapMC;

  // PUI weight histogram
  TH1F * _puHisto;
  TH1F * _puFall11MC;
  TH1F * _puFall11Data;
  TH1F * _puSummer12MC;
  TH1F * _puSummer12Data;
  TH1F * _puHCP12MC;
  TH1F * _puHCP12Data;

  void InitZttEmbedDataWeight(TString ZttEmbedDataFileName);
  float ZttEmbedDataWeight(float value, int catIdx);

  int binNumber(float x, const std::vector<float> bins) const
  {
    for (size_t iB=0; iB<bins.size(); ++iB)
      if (x>=bins[iB]&&x<bins[iB+1])
	return iB;
    return 0;
  }

  int binNumber(float x, int nbins, const float * bins) {

    int binN = 0;

    for (int iB=0; iB<nbins; ++iB) {
      if (x>=bins[iB]&&x<bins[iB+1]) {
	binN = iB;
	break;
      }
    }
    
    return binN;

  }

  float effBin(float x, int nbins, float * bins, const float * eff) {

    int bin = binNumber(x, nbins, bins);

    return eff[bin];

  }

	float effBin(float x, std::vector<float>& bins, const float * eff)
	{
		return eff[binNumber(x, bins)];
	}


 private:


  float hlt_IsoMu17_Efficiency(float pt, float eta);
  float hlt_Mu13Mu8_Efficiency13(float abseta, float pt, float *efficiency_vec);
  float hlt_Mu13Mu8_Efficiency8(float abseta, float pt, float *efficiency_vec);
  float hlt_Mu13Mu8_Efficiency13_Run2011B(float abseta, float pt, float *efficiency_vec);
  float hlt_Mu13Mu8_Efficiency8_Run2011B(float abseta, float pt, float *efficiency_vec);

  void  U1U2CorrectionsByWidth(float & U1, float & U2,
			       int nZptBin,
			       int njets);

  float CorrectionsBySampling(float x, TF1 * funcMC, TF1 * funcData);


  TF1 * getFuncRecoil(TF1* initFunc, bool left, bool isGauss);


  TString _baseDir;

  float _hltIso17_EtaBins[5];
  float _hltIso17_PtBins[6];

  int _hltIso17_nPtBins;
  int _hltIso17_nEtaBins;

  float _hltIso17_Eff[5][4];

  //
  //  TFile * _fileMetData;
  //  TFile * _fileMetMC;
  //  TFile * _fileEtaWeights; 
  //  TFile * _fileRhoWeights;
  //  TFile * _fileDCAWeights;

  //  TString _perpZStr[4];
  //  TString _paralZStr[4];
  //  TString _nJetsStr[3];
  //


  // *****************************************
  // Change here number of ZPt bins if needed
  // *****************************************

  //float * _ZPtBins;
  std::vector<float> _ZPtBins;

  int _nZPtBins;
  int _nJetsBins;

  TF1 * _metZParalData[5][3];
  TF1 * _metZPerpData[5][3];
  TF1 * _metZParalMC[5][3];
  TF1 * _metZPerpMC[5][3];

  TH1F * _metZParalDataHist[5][3];
  TH1F * _metZPerpDataHist[5][3];
  TH1F * _metZParalMCHist[5][3];
  TH1F * _metZPerpMCHist[5][3];

  float _meanMetZParalData[5][3];
  float _meanMetZParalMC[5][3];
  float _meanMetZPerpData[5][3];
  float _meanMetZPerpMC[5][3];
  
  float _rmsMetZParalData[5][3];
  float _rmsLeftMetZParalData[5][3];
  float _rmsRightMetZParalData[5][3];
  
  float _rmsMetZParalMC[5][3];
  float _rmsLeftMetZParalMC[5][3];
  float _rmsRightMetZParalMC[5][3];

  float _rmsMetZPerpData[5][3];
  float _rmsLeftMetZPerpData[5][3];
  float _rmsRightMetZPerpData[5][3];
  
  float _rmsMetZPerpMC[5][3];
  float _rmsLeftMetZPerpMC[5][3];
  float _rmsRightMetZPerpMC[5][3];
  
  float _xminMetZPerp[5][3];
  float _xmaxMetZPerp[5][3];

  float _xminMetZPerpData[5][3];
  float _xmaxMetZPerpData[5][3];

  float _xminMetZPerpMC[5][3];
  float _xmaxMetZPerpMC[5][3];

  float _xminMetZParal[5][3];
  float _xmaxMetZParal[5][3];

  float _xminMetZParalData[5][3];
  float _xmaxMetZParalData[5][3];

  float _xminMetZParalMC[5][3];
  float _xmaxMetZParalMC[5][3];

  // diLepEta weights  

  TF1 * _diLepEtaWeightFunc;
  float _etaMax;

  // Rho functions

  std::vector<TF1*> _rhoFuncMC; 
  std::vector<TF1*> _rhoFuncData; 

  TH1F * _rhoDataBin1;
  TH1F * _rhoDataBin2;
  
  TH1F * _rhoMCBin1;
  TH1F * _rhoMCBin2;

  // DCA related functions
  // and parameters

  TF1 * _funcDcaData2012A;
  TF1 * _funcDcaMCSummer12;

  int _nDiLepMassBins;
  int _nReducedLikeBins;

  std::vector<float> _diLepMassBins;
  std::vector<float> _reducedLikeBins;

  TF1 * _funcDCACentralData[5][4];
  TF1 * _funcDCACentralMC[5][4];

  TF1 * _funcDCALeftTailData[5][4];
  TF1 * _funcDCALeftTailMC[5][4];

  TF1 * _funcDCARightTailData[5][4];
  TF1 * _funcDCARightTailMC[5][4];

  float _cutoffLeftDcaData[5][4];
  float _cutoffLeftDcaMC[5][4];
  float _cutoffRightDcaData[5][4];
  float _cutoffRightDcaMC[5][4];
  
  float _minDca[5][4];
  float _maxDca[5][4];

  float _integralLeftDcaData[5][4];
  float _integralCenterDcaData[5][4];
  float _integralRightDcaData[5][4];

  float _integralLeftDcaMC[5][4];
  float _integralCenterDcaMC[5][4];
  float _integralRightDcaMC[5][4];

  float _totWeightDCA[5][4];

  // # of primary vertices max
 
  int _nPVMax;

  int _dataset;

  float _weightMin;
  float _weightMax;

  
  TFile * ggHWeight_file;
  TH1F * ggHWeightHqT_histo;
  TH1F * ggHWeightFeHiPro_histo;
  TH2F * _theoryWeightHisto;

  TFile * ZttEmbedDataWeightFile;
  TH1F * ZttEmbedDataWeightHisto;
};
