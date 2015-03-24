#include "HttStylesNew.cc"
void PlotInputs(TString fileName = "Signal_mH1_8_New_inputs",
		TString fileName4 = "Signal_mH1_4_New_inputs",
		bool logY = false,
		bool blindData = false,
		float signalScale = 5.) {

  // RooT macro plotting unrolled 2D distributions used in the
  // statistical analysis
  // Input parameters
  // 1. fileName (TString) : file name with input histograms (w/o root extension) 
  // 2. logY (bool) : use log scale for Y axis
  // 3. blindData (bool) : blind data in the signal region (last 4 bins of the unrolled 2D distribution)
  // 4. signalScale (float) : signal strength

  int nSysErr = 4;
  float sysErr[4] = {2.6, 3.3, 4.0, 10.0};

  float totErrSig2 = 0;
  for (int iS=1; iS<4; ++iS)
    totErrSig2 += sysErr[iS]*sysErr[iS];

  float totErrSig = 0.01*TMath::Sqrt(totErrSig2);

  SetStyle();
  //  gStyle->SetOptStat(0);
  //  gStyle->SetErrorX(0.5);
  
  int bB = 0;
  if (blindData)
    bB = 4;
  int blindedBins[4] = {7,8,9,10};

  TFile * file = new TFile(fileName+".root");
  TFile * file4 = new TFile(fileName4+".root");

  TH1F * data = (TH1F*)file->Get("data_obs");
  TH1F * sig = (TH1F*)file->Get("sig");
  TH1F * sig4 = (TH1F*)file4->Get("sig");
  TH1F * bkgd = (TH1F*)file->Get("bkgd");
  TH1F * bkgdUp = (TH1F*)file->Get("bkgd_QCDShapeUp");
  TH1F * bkgdDown = (TH1F*)file->Get("bkgd_QCDShapeDown");

  TH1F * blindH = data->Clone("blindH");
  blindH->SetMarkerSize(0);
  blindH->SetFillColor(7);
  blindH->SetFillStyle(3013);
  blindH->SetLineWidth(1);
  blindH->SetLineColor(0);
  for (int iB=1; iB<=blindH->GetNbinsX(); ++iB)
    blindH->SetBinContent(iB,0);
  for (int iB=0; iB<bB; ++iB) {
    blindH->SetBinContent(blindedBins[iB],20.);
    data->SetBinContent(blindedBins[iB],-10.);
  }

  data->SetLineColor(1);
  data->SetLineWidth(3);
  data->SetMarkerSize(1.4);
  bkgd->SetLineColor(4);
  bkgd->SetLineWidth(2);
  bkgdUp->SetLineColor(4);
  bkgdUp->SetLineStyle(2);
  bkgdDown->SetLineColor(4);
  bkgdDown->SetLineStyle(3);
  sig->SetLineColor(2);
  sig->SetLineWidth(2);
  sig4->SetLineColor(kGreen+1);
  sig4->SetLineWidth(2);
  data->GetXaxis()->SetTitle("");
  data->GetYaxis()->SetTitle("Events");
  bkgd->GetXaxis()->SetTitle("");
  bkgd->GetYaxis()->SetTitle("Events");
  data->GetYaxis()->SetRangeUser(0,1.1*data->GetMaximum());
  bkgd->GetYaxis()->SetRangeUser(0,1.1*bkgd->GetMaximum());
  TString binLabel[15] = {"(1,1)","(1,2)","(1,3)","(1,4)","(1,5)",
			  "(2,2)","(2,3)","(2,4)","(2,5)",
			  "(3,3)","(3,4)","(3,5)",
			  "(4,4)","(4,5)",
			  "(5,5)"};  

  int nBins = bkgd->GetNbinsX();
  float xMin = bkgd->GetBinLowEdge(1);
  float xMax = bkgd->GetBinLowEdge(nBins+1);
  if (nBins<15) {
    binLabel[0] = "(1,1)";
    binLabel[1] = "(1,2)";
    binLabel[2] = "(1,3)";
    binLabel[3] = "(1,4)";
    binLabel[4] = "(2,2)";
    binLabel[5] = "(2,3)";
    binLabel[6] = "(2,4)";
    binLabel[7] = "(3,3)";
    binLabel[8] = "(3,4)";
    binLabel[9] = "(4,4)";
  }

  for (int iB=1;iB<=nBins;++iB) {
    data->GetXaxis()->SetBinLabel(iB,"");
    bkgd->GetXaxis()->SetBinLabel(iB,"");
  }

  data->GetXaxis()->SetLabelSize(0.08);
  bkgd->GetXaxis()->SetLabelSize(0.08);

  data->GetYaxis()->SetRangeUser(0,300);
  bkgd->GetYaxis()->SetRangeUser(0,300);

  

  if (logY) {
    data->GetYaxis()->SetRangeUser(2,2000);
    bkgd->GetYaxis()->SetRangeUser(2,2000);
  }


  TH1F * errorBand = bkgd->Clone("errorBand");
  errorBand->SetMarkerSize(0);
  errorBand->SetFillColor(4);
  errorBand->SetFillStyle(3013);
  errorBand->SetLineWidth(2);
  
  float chi2 = 0;
  
  std::cout << std::endl;
  std::cout << std::endl;

  int nBinsIncluded = 0;

  float totData = data->GetSum();
  float relErr = 1/TMath::Sqrt(totData);

  for (int iB=1; iB<=nBins; ++iB) {
    float eSys =  bkgdUp->GetBinContent(iB)-bkgd->GetBinContent(iB);
    float eSys2 = eSys*eSys;
    for (int iSys=1; iSys<=nBins; ++iSys) {
      char binChar[5];
      if (iSys<10)
	sprintf(binChar,"%1i",iSys);
      else
	sprintf(binChar,"%2i",iSys);
      TString binStr(binChar);
      TH1F * bkgBin = file->Get("bkgd_BkgBin"+binStr+"Up");
      float err = bkgBin->GetBinContent(iB) - bkgd->GetBinContent(iB);
      eSys2 += err*err;
    }
    float eNor = relErr*bkgd->GetBinContent(iB);
    eSys2 += eNor*eNor;
    float eSys = TMath::Sqrt(eSys2);
    errorBand->SetBinError(iB,eSys);
    bool include = true;
    if (blindData) {
      for (int iJ=0; iJ<bB; ++iJ) {
	if (iB==blindedBins[iJ]) {
	  include = false;
	  break;
	}
      }
    }
    if (include) {
      float dataX = data->GetBinContent(iB);
      float bkgdX = bkgd->GetBinContent(iB);
      float dataE = data->GetBinError(iB);
      float dataE2 = dataX; 
      if (dataX < 1.0) dataE2 = 1.0; 
      float err2 = dataE2 + eSys*eSys;
      float eTot = TMath::Sqrt(err2);
      
      float chi2X = (dataX-bkgdX)*(dataX-bkgdX)/err2;
      chi2 += chi2X;

      float pull = (dataX-bkgdX)/eTot;
      
      printf("bin %2i -> Data = %4i ; B = %8.4f ; Error = %7.4f ; Pull = %4.1f\n",iB,int(dataX+0.1),bkgdX,eTot,pull);
      nBinsIncluded++;
    }

  }

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "chi2/ndof = " << chi2 << "/" << nBinsIncluded << "  Prob = " << TMath::Prob(chi2,float(nBinsIncluded)) << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;

  // ratio plots
  TH1F * dataRatioH = data->Clone("dataRatioH");
  TH1F * bkgRatioErrH = errorBand->Clone("bkgRatioH");
  TH1F * unitH = bkgd->Clone("unitH");

  for (int iB=1; iB<=nBins; ++iB) {
    float datX = data->GetBinContent(iB);
    float datE = data->GetBinError(iB);
    float bkgX = errorBand->GetBinContent(iB);
    float bkgE = errorBand->GetBinError(iB);
    float datRatioX = datX/bkgX;
    float datRatioE = datE/bkgX;
    float bkgErr = bkgE / bkgX;
    dataRatioH->SetBinContent(iB,datRatioX);
    dataRatioH->SetBinError(iB,datRatioE);
    bkgRatioErrH->SetBinContent(iB,1);
    bkgRatioErrH->SetBinError(iB,bkgErr);
    unitH->SetBinContent(iB,1);
    unitH->SetBinError(iB,0);
    bkgd->SetBinError(iB,0);
    float xSig = sig->GetBinContent(iB);
    float eSig = sig->GetBinError(iB);
    float esysSig = xSig*totErrSig;
    float xSig4 = sig4->GetBinContent(iB);
    float eSig4 = sig4->GetBinError(iB);
    sig->SetBinContent(iB,xSig*signalScale);
    sig4->SetBinContent(iB,xSig4*signalScale);
    char str[100];
    sprintf(str,"|  %4i  |  %8.2f+/-%5.2f  |  %6.4f+/-%6.4f  |",datX,bkgX,bkgE,xSig,esysSig);
    std::cout << binLabel[iB-1] << str << std::endl;
  }

  unitH->SetLineColor(4);
  unitH->SetLineWidth(2);

  bkgRatioErrH->SetFillStyle(3013);
  bkgRatioErrH->SetFillColor(4);
  bkgRatioErrH->SetMarkerStyle(21);
  bkgRatioErrH->SetMarkerSize(0);
  bkgRatioErrH->SetLineColor(4);
  bkgRatioErrH->SetLineWidth(2);

  dataRatioH->GetYaxis()->SetRangeUser(0.4,1.6);
  dataRatioH->GetYaxis()->SetNdivisions(505);
  dataRatioH->GetXaxis()->SetLabelFont(42);
  dataRatioH->GetXaxis()->SetLabelOffset(0.04);
  dataRatioH->GetXaxis()->SetLabelSize(0.14);
  dataRatioH->GetXaxis()->SetTitleSize(0.13);
  dataRatioH->GetXaxis()->SetTitleOffset(1.15);
  dataRatioH->GetXaxis()->SetTitle("Bins");
  dataRatioH->GetYaxis()->SetTitle("Data/Bkg");
  dataRatioH->GetXaxis()->CenterTitle();
  dataRatioH->GetYaxis()->CenterTitle();
  dataRatioH->GetYaxis()->SetLabelFont(42);
  dataRatioH->GetYaxis()->SetLabelOffset(0.015);
  dataRatioH->GetYaxis()->SetLabelSize(0.11);
  dataRatioH->GetYaxis()->SetTitleSize(0.14);
  dataRatioH->GetYaxis()->SetTitleOffset(0.55);
  
  for (int iB=1;iB<=nBins;++iB) {
    dataRatioH->GetXaxis()->SetBinLabel(iB,binLabel[iB-1]);
  }

  bkgd->SetTitle("CMS Preliminary,  L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
  bkgd->GetXaxis()->SetLabelOffset(0.02);
  bkgd->GetYaxis()->SetTitle("Events");
  bkgd->GetYaxis()->SetTitleSize(0.07);
  bkgd->GetYaxis()->SetTitleOffset(1.1);
  bkgd->SetMarkerSize(0);
  bkgd->SetMarkerStyle(0);

  TCanvas * canv = new TCanvas("canv","",700,700);
  canv->SetLeftMargin(0.45);
  TPad * upper =  new TPad("upper", "pad",0,0.26,1,1);
  //upper->SetBottomMargin(0);                                                                                                                              
  upper->SetLeftMargin(0.2);
  upper->Draw();
  upper->cd();

  //  data->Draw("e1");
  //  blindH->Draw("same");
  //  bkgd->Draw("same");
  //  errorBand->Draw("e2same");
  bkgd->Draw("h");
  sig->Draw("same");
  sig4->Draw("same");
  errorBand->Draw("e2same");
  data->Draw("e1same");
  TLegend * leg = new TLegend(0.55,0.68,0.95,0.92);
  leg->SetFillColor(0);
  //  leg->SetHeader("                  Prefit");
  leg->SetTextSize(0.04);
  leg->AddEntry(data,"Data","lp");
  leg->AddEntry(errorBand,"Bkg(+uncertainty)","lf");
  //  leg->AddEntry(bkgd,"Bkgd","lf");
  leg->AddEntry(sig4,"Signal, m(#phi_{1}) = 4 GeV","l");
  leg->AddEntry(sig, "Signal, m(#phi_{1}) = 8 GeV","l");
  //  leg->AddEntry(blindH,"blinded bins","f");
  upper->SetLogy(logY);
  upper->RedrawAxis();
  leg->Draw();
  upper->SetGridx(true);
  upper->SetGridy(true);
  upper->Update();
  canv->cd();

  TPad * lower = new TPad("lower", "pad",0,0.0,1,0.34);
  lower->Draw();
  //lower->SetTopMargin(0);
  lower->SetBottomMargin(0.35);
  lower->SetLeftMargin(0.2);
  lower->cd();

  dataRatioH->Draw("e1");
  bkgRatioErrH->Draw("e2same");
  unitH->Draw("hsame");

  TLine * line = new TLine(0,1,10,1);
  line->SetLineColor(4);
  line->SetLineWidth(2);
  line->Draw();
  dataRatioH->Draw("e1same");

  lower->Modified();
  lower->RedrawAxis();
  lower->SetGridx(true);
  lower->SetGridy(true);
  canv->cd();
  canv->Modified();
  canv->cd();
  canv->SetSelected(canv);


  if (logY) {
    canv->Print("Prefit_log.png");
    canv->Print("Prefit_log.pdf","Portrait pdf");
  }
  else {
    canv->Print("Prefit.png");
    canv->Print("Prefit.pdf","Portrait pdf");
  }

}
