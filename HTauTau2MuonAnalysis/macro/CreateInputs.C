#include "HtoHold.H"
#include "tdrstyle.C"
void CreateInputs(TString signalFileName = "Signal_mH1_4_New",
		  TString dataFileName = "Data_New",
		  TString corrCentralFileName = "Data_RegionA_corr",
		  TString corrUpFileName = "Data_RegionB_corr",
		  TString binVersion = "_v0",
		  int contaminationTreatment = -1,
		  float sigma = 1.) {
  // signalFileName : file name of the signal sample w/o root extension 
  // dataFileName : file name of the data sample w/o root extension
  // corrCentralFileName : file with correlation coefficients (central estimates)
  // corrUpFileName : file with correlation coefficients (upward variation)
  // binVersion     : binning option : _v0 _v1 _v2
  // contaminationTreatment : 0 - one-sided gaussian error due to signal contamination of 
  //                              background side-band
  // sigma : assumed sigma*BR for signal (default value 1 pb)


  TString suffix("");
  if (contaminationTreatment==0) 
    suffix = "_contam0";
  if (contaminationTreatment==1)
    suffix = "_contam1";
  if (contaminationTreatment==2)
    suffix = "_contam2";

  setTDRStyle();
  gStyle->SetOptStat(0);
  float lumi = 19710; // integrated luminosity

  // correlation coefficients

  TString dir("Correlations/");
  TFile * fileCorrCentral = new TFile(dir+corrCentralFileName+binVersion+".root");
  TFile * fileCorrCentralContamination = new TFile(dir+corrCentralFileName+binVersion+"_"+signalFileName+".root");
  TFile * fileCorrUp = new TFile(dir+corrUpFileName+binVersion+".root"); 
  TH2F * corr2Dcentral = (TH2F*)fileCorrCentral->Get("Corr_2D");
  TH2F * corr2Dcontamination = (TH2F*)fileCorrCentralContamination->Get("Corr_2D");
  TH2F * corr2Dup = (TH2F*)fileCorrUp->Get("Corr_2D");


  TFile * file = new TFile(dataFileName+".root");
  TFile * fileSignal = new TFile(signalFileName+".root");
  TH1F * inputEventsH = (TH1F*)fileSignal->Get("inputEventsH");
  float xGen = float(inputEventsH->GetSumOfWeights());
  float normSignal = sigma * lumi / xGen;

  TString histNamePos2Q("InvMassTrackPlusNeg_2QH");
  TString histNameNeg2Q("InvMassTrackPlusPos_2QH");
  TString histNamePos3Q("InvMassTrackPlusNeg_3QH");
  TString histNameNeg3Q("InvMassTrackPlusPos_3QH");

  TH1F * histOldPos2Q = (TH1F*)file->Get(histNamePos2Q);
  TH1F * histOldNeg2Q = (TH1F*)file->Get(histNameNeg2Q);
  TH1F * histOldPos3Q = (TH1F*)file->Get(histNamePos3Q);
  TH1F * histOldNeg3Q = (TH1F*)file->Get(histNameNeg3Q);
  TH1F * histOld2Q = histOldPos2Q->Clone("histOld2Q");
  TH1F * histOld3Q = histOldPos3Q->Clone("histOld3Q");
  histOld2Q->Add(histOld2Q,histOldNeg2Q);
  histOld3Q->Add(histOld3Q,histOldNeg3Q);
  TH1F * hist1DOld = histOld2Q->Clone("hist1D");
  hist1DOld->Add(hist1DOld,histOld3Q);

  int nBinsX = hist1DOld->GetNbinsX();
  float xLower = hist1DOld->GetBinLowEdge(1);
  float xUpper = hist1DOld->GetBinLowEdge(nBinsX+1);

  std::cout << std::endl;
  std::cout << " nBins = " << nBinsX
	    << "  ; lower = " << xLower
	    << "  ; upper = " << xUpper << std::endl;
  std::cout << std::endl;

  int nBins = 10;
  float xBins[20]; 

  TH1F * binsH = fileCorrCentral->Get("binsH");
  nBins = binsH->GetNbinsX();
  for (int iB=0; iB<=nBins; ++iB) 
    xBins[iB] = binsH->GetXaxis()->GetBinLowEdge(iB+1);

  TH1F * hist1D = TH1toTH1(hist1DOld,nBins,xBins,true,"_rebin"); 
  std::cout << std::endl;
  for (int iB=1; iB<=nBins; ++iB)
  hist1D->Scale(1/hist1D->GetSumOfWeights());

  std::cout << std::endl;
  for (int iB=1; iB<=nBins; ++iB) {
    std::cout << iB << " " << hist1D->GetBinContent(iB) << std::endl;
  }
  std::cout << std::endl;

  TH2F * histData2DOld = (TH2F*)file->Get("InvMassTrackPlusMuon2D_H");
  TH2F * histSignal2DOld = (TH2F*)fileSignal->Get("InvMassTrackPlusMuon2D_H");

  TH2F * histData2D   = TH2toTH2(histData2DOld, nBins, xBins, nBins, xBins, true, "_data");
  TH2F * histSignal2D = TH2toTH2(histSignal2DOld, nBins, xBins, nBins, xBins, true, "_signal");

  TH2F * histBkg2DCentral = new TH2F("histBkg2DCentral","",nBins, xBins, nBins, xBins);
  TH2F * histBkg2DDown = new TH2F("histBkg2DDown","",nBins, xBins, nBins, xBins);
  TH2F * histBkg2DUp = new TH2F("histBkg2DUp","",nBins, xBins, nBins, xBins);
  TH2F * histBkgContamination2DDown = new TH2F("histBkgContamination2DDown","",nBins,xBins,nBins,xBins);
  TH2F * histBkgContamination2DUp   = new TH2F("histBkgContamination2DUp","",nBins,xBins,nBins,xBins);
  

  for (int iB=1; iB<=nBins; ++iB) {
    for (int jB=1; jB<=nBins; ++jB) {
      float x = hist1D->GetBinContent(iB);
      float y = hist1D->GetBinContent(jB);
      float corrCentral = corr2Dcentral->GetBinContent(iB,jB);
      float corrCentralE = corr2Dcentral->GetBinError(iB,jB);
      float corrUp = corr2Dup->GetBinContent(iB,jB);
      float corrContamination = corr2Dcontamination->GetBinContent(iB,jB);

      float central_tmp = corrCentral*x*y;
      float central_err = corrCentralE*x*y;
      
      float central = central_tmp;
      float downContamination = corrContamination*x*y;
      float upContamination = central_tmp;

      if (contaminationTreatment==1) {
	central = 0.5*(corrContamination+corrCentral)*x*y;
	upContamination = central_tmp;
      }
      if (contaminationTreatment==2) {
	upContamination = 2*central_tmp - downContamination;
      }

      float up = corrUp*x*y;
      float down = 2*central - up;

      if (down<1e-10) down = 1e-10; // enforce positive bin content
      if (up<1e-10) up = 1e-10; // enforce positive bin content
      if (downContamination<1e-10) downContamination = 1e-10; // enforce positive bin content
      if (upContamination<1e-10) upContamination = 1e-10; // enforce positive bin content
      if (central<1e-10) central = 1e-10;
      std::cout << "[" << iB << "," << jB << "]   ";
      printf(" up = %4.2f   down = %4.2f\n",up/central, down/central);
      histBkg2DDown->SetBinContent(iB,jB,down);
      histBkg2DCentral->SetBinContent(iB,jB,central);
      histBkg2DCentral->SetBinError(iB,jB,central_err);
      histBkg2DUp->SetBinContent(iB,jB,up);
      histBkgContamination2DDown->SetBinContent(iB,jB,downContamination);
      histBkgContamination2DUp->SetBinContent(iB,jB,upContamination);
    }
  }

  //  TString inputsFileName = + signalFileName+"_inputs"+binVersion+suffix+".root";
  TString inputsFileName = signalFileName+"_inputs.root";
  TFile * fileInputs = new TFile("inputs/"+inputsFileName,"recreate");
  fileInputs->cd("");
  
  int nBins1D =  nBins * (nBins+1) / 2 ;
  TH1F * data = new TH1F("data_obs","",nBins1D,0.,float(nBins1D));
  TH1F * bkgd = new TH1F("bkgd","",nBins1D,0.,float(nBins1D));
  TH1F * bkgdUp = new TH1F("bkgd_QCDShapeUp","",nBins1D,0.,float(nBins1D));
  TH1F * bkgdDown = new TH1F("bkgd_QCDShapeDown","",nBins1D,0.,float(nBins1D));
  TH1F * bkgdContaminationUp = new TH1F("bkgd_SigContamUp","",nBins1D,0.,float(nBins1D));
  TH1F * bkgdContaminationDown = new TH1F("bkgd_SigContamDown","",nBins1D,0.,float(nBins1D));
  //  std::vector<TH1F*> bkgdUpBBB;
  //  std::vector<TH1F*> bkgdDownBBB;

  std::cout << std::endl;
  std::cout << "NBins1D = " << nBins1D << std::endl;
  std::cout << std::endl;

  TH1F * signal = new TH1F("sig","",nBins1D,0.,float(nBins1D));

  float yieldBkgd = histBkg2DCentral->GetSum();
  float yieldBkgdUp = histBkg2DUp->GetSum();
  float yieldBkgdDown = histBkg2DDown->GetSum();
  float yieldBkgdContaminationUp = histBkgContamination2DUp->GetSum();
  float yieldBkgdContaminationDown = histBkgContamination2DDown->GetSum();

  float yieldData = histData2D->GetSum();
  float estimatedBkgd = yieldData;
  float normBkgd = estimatedBkgd / yieldBkgd;

  int iBin = 1;

  for (int iB=1; iB<=nBins; ++iB) {
    for (int jB=iB; jB<=nBins; ++jB) {
      float dataX = histData2D->GetBinContent(iB,jB);
      float bkgdX = histBkg2DCentral->GetBinContent(iB,jB);
      float bkgdE = histBkg2DCentral->GetBinError(iB,jB);
      float bkgdUpX = histBkg2DUp->GetBinContent(iB,jB);
      float bkgdDownX = histBkg2DDown->GetBinContent(iB,jB);
      float signalX = histSignal2D->GetBinContent(iB,jB);
      float signalE2 = histSignal2D->GetBinError(iB,jB)*histSignal2D->GetBinError(iB,jB);
      float bkgdContaminationUpX = histBkgContamination2DUp->GetBinContent(iB,jB);
      float bkgdContaminationDownX = histBkgContamination2DDown->GetBinContent(iB,jB);
      if (iB!=jB) {
	dataX += histData2D->GetBinContent(jB,iB);
	bkgdX += histBkg2DCentral->GetBinContent(jB,iB);
	bkgdE += histBkg2DCentral->GetBinError(jB,iB);
	bkgdUpX += histBkg2DUp->GetBinContent(jB,iB);
	bkgdDownX += histBkg2DDown->GetBinContent(jB,iB);
	signalX += histSignal2D->GetBinContent(jB,iB);
	signalE2 += histSignal2D->GetBinError(jB,iB)*histSignal2D->GetBinError(jB,iB);
	bkgdContaminationUpX += histBkgContamination2DUp->GetBinContent(jB,iB);
	bkgdContaminationDownX += histBkgContamination2DDown->GetBinContent(jB,iB);
      }
      float signalE = TMath::Sqrt(signalE2);
      data->SetBinContent(iBin,dataX);
      bkgd->SetBinContent(iBin,normBkgd*bkgdX);
      bkgd->SetBinError(iBin,normBkgd*bkgdE);
      bkgdUp->SetBinContent(iBin,normBkgd*bkgdUpX);
      bkgdDown->SetBinContent(iBin,normBkgd*bkgdDownX);
      bkgdContaminationUp->SetBinContent(iBin,normBkgd*bkgdContaminationUpX);
      bkgdContaminationDown->SetBinContent(iBin,normBkgd*bkgdContaminationDownX);
      signal->SetBinContent(iBin,normSignal*signalX);
      signal->SetBinError(iBin,normSignal*signalE);
      if (bkgd->GetBinContent(iBin)<1e-4)
	bkgd->SetBinContent(iBin,1e-4);
      if (bkgdUp->GetBinContent(iBin)<1e-4)
	bkgdUp->SetBinContent(iBin,1e-4);
      if (bkgdDown->GetBinContent(iBin)<1e-4)
	bkgdDown->SetBinContent(iBin,1e-4);
      if (bkgdContaminationUp->GetBinContent(iBin)<1e-4)
	bkgdContaminationUp->SetBinContent(iBin,1e-4);
      if (bkgdContaminationDown->GetBinContent(iBin)<1e-4)
	bkgdContaminationDown->SetBinContent(iBin,1e-4);
      printf("bin %2i   up = %4.2f   down = %4.2f\n",iBin,bkgdUp->GetBinContent(iBin)/bkgd->GetBinContent(iBin),bkgdDown->GetBinContent(iBin)/bkgd->GetBinContent(iBin));
      iBin++;
    }
  }


  float yieldUp   = bkgdUp->GetSumOfWeights();
  float yieldDown = bkgdDown->GetSumOfWeights();
  bkgdUp->Scale( estimatedBkgd / yieldUp);
  bkgdDown->Scale( estimatedBkgd / yieldDown);

  float yieldContaminationUp   = bkgdContaminationUp->GetSumOfWeights();
  float yieldContaminationDown = bkgdContaminationDown->GetSumOfWeights();
  bkgdContaminationUp->Scale( estimatedBkgd / yieldContaminationUp);
  bkgdContaminationDown->Scale( estimatedBkgd / yieldContaminationDown);


  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Data     = " << data->GetSum() << std::endl;
  std::cout << "Bkgd     = " << bkgd->GetSum() << std::endl;
  std::cout << "BkgdUp   = " << bkgdUp->GetSum() << std::endl;
  std::cout << "BkgdDown = " << bkgdDown->GetSum() << std::endl;
  std::cout << "BkgdUp (Contamination)   = " << bkgdContaminationUp->GetSum() << std::endl;
  std::cout << "BkgdDown (Contamination) = " << bkgdContaminationDown->GetSum() << std::endl;
  std::cout << "Signal   = " << signal->GetSum() << std::endl;
  float significance = signal->GetSum() / TMath::Sqrt(data->GetSum());
  std::cout << "Signif   = " << significance << std::endl;

  int obsData = int(data->GetSum()+0.01);
  float sigEvents = signal->GetSum();
  float bkgEvents = bkgd->GetSum(); 

  ostringstream str;
  str << "inputs/" << signalFileName << "_inputs.txt";
  string nn = str.str();
  const char * p = nn.c_str();

  std::ofstream textFile(p);
  textFile << "imax 1   number of channels" << std::endl;
  textFile << "jmax *   number of backgrounds" << std::endl;
  textFile << "kmax *   number of nuisance parameters" << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "observation " << obsData << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "shapes * * " << inputsFileName << "  $PROCESS    $PROCESS_$SYSTEMATIC " << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "bin     H2to2H1  H2to2H1  " << std::endl;
  textFile << "process     sig     bkgd  " << std::endl;      
  textFile << "process       0        1  " << std::endl; 
  textFile << "rate    " << sigEvents << "      " << bkgEvents << std::endl;
  textFile << "-----------------" << std::endl;
  textFile << "lumi            lnN    1.026    -  " << std::endl;
  textFile << "sig_accept_4tau lnN    1.033    -  " << std::endl;
  textFile << "CMS_eff_m       lnN    1.040    -  " << std::endl;
  textFile << "CMS_eff_trk     lnN    1.100    -  " << std::endl;
  textFile << "CMS_BkgNorm     lnU      -     2.0 " << std::endl;  
  textFile << "QCDShape      shape      -     1.0 " << std::endl;
  if (contaminationTreatment>=0)
    textFile << "SigContam     shape      -     1.0 " << std::endl;

  std::cout << std::endl;
  std::cout << std::endl;
  for (int iB=1; iB<=nBins1D; ++iB) {
    char number[2];
    if (iB<10) 
      sprintf(number,"%1i",iB);
    else 
      sprintf(number,"%2i",iB);
    TString nameUp = TString("bkgd_BkgBin") + TString(number) + TString("Up");
    TString nameDown = "bkgd_BkgBin" + TString(number) + TString("Down");
    TString nameSigUp = TString("sig_SigBin") + TString(number) + TString("Up");
    TString nameSigDown = "sig_SigBin" + TString(number) + TString("Down");
    TH1F * hUp = bkgd->Clone(nameUp);
    TH1F * hDown = bkgd->Clone(nameDown);
    TH1F * hSigUp = sig->Clone(nameSigUp);
    TH1F * hSigDown = sig->Clone(nameSigDown);
    float xUp = bkgd->GetBinContent(iB)+bkgd->GetBinError(iB);
    float xDown = bkgd->GetBinContent(iB)-bkgd->GetBinError(iB);
    float xSigUp = sig->GetBinContent(iB)+sig->GetBinError(iB);
    float xSigDown = sig->GetBinContent(iB)-sig->GetBinError(iB);
    if (xUp<1e-4) xUp = 1e-4;
    if (xDown<1e-4) xDown = 1e-4;
    if (xSigUp<1e-6) xSigUp = 1e-6;
    if (xSigDown<1e-6) xSigDown = 1e-6;
    hUp->SetBinContent(iB,xUp);
    hDown->SetBinContent(iB,xDown);
    hSigUp->SetBinContent(iB,xSigUp);
    hSigDown->SetBinContent(iB,xSigDown);
    float yieldUp = hUp->GetSum();
    float yieldDown = hDown->GetSum();
    hUp->Scale( estimatedBkgd / yieldUp);
    hDown->Scale( estimatedBkgd / yieldDown);
    if (iB<10) {
      textFile << "BkgBin" << iB << "       shape      -     1.0" << std::endl;
      textFile << "SigBin" << iB << "       shape     1.0     - " << std::endl;
    }
    else {
      textFile << "BkgBin" << iB << "      shape      -     1.0" << std::endl;
      textFile << "SigBin" << iB << "      shape     1.0     - " << std::endl;
    }
  }

  textFile << std::endl;

  fileInputs->Write();
  fileInputs->Close();

}
