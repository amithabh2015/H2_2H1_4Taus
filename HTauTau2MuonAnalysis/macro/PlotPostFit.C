#include "HttStylesNew.cc"
void PlotPostFit(TString fileName = "Signal_mH1_8_New_inputs",
		 TString fileName4 = "Signal_mH1_4_New_inputs",
		 bool logY = false,
		 bool blindData = true,
		 float signalScale = 5.) {

  // RooT macro plotting unrolled 2D distributions used in the
  // statistical analysis
  // Input parameters
  // 1. fileName (TString) : file name with input histograms (w/o root extension) 
  // 2. logY (bool) : use log scale for Y axis
  // 3. blindData (bool) : blind data in the signal region (last 4 bins of the unrolled 2D distribution)
  // 4. signalScale (float) : signal strength

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

  sig->Scale(signalScale);
  sig4->Scale(signalScale);

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
  data->SetLineWidth(2);
  data->SetMarkerSize(1.2);
  bkgd->SetLineColor(4);
  bkgd->SetLineWidth(2);
  bkgdUp->SetLineColor(4);
  bkgdUp->SetLineStyle(2);
  bkgdDown->SetLineColor(4);
  bkgdDown->SetLineStyle(3);
  sig->SetLineColor(2);
  sig->SetLineWidth(2);
  sig4->SetLineColor(3);
  sig4->SetLineWidth(2);
  data->GetXaxis()->SetTitle("Bins");
  data->GetYaxis()->SetTitle("Events");
  bkgd->GetXaxis()->SetTitle("Bins");
  bkgd->GetYaxis()->SetTitle("Events");
  data->GetYaxis()->SetRangeUser(0,1.1*data->GetMaximum());
  bkgd->GetYaxis()->SetRangeUser(0,1.1*bkgd->GetMaximum());
  TString binLabel[15] = {"(1,1)","(1,2)","(1,3)","(1,4)","(1,5)",
			  "(2,2)","(2,3)","(2,4)","(2,5)",
			  "(3,3)","(3,4)","(3,5)",
			  "(4,4)","(4,5)",
			  "(5,5)"};  

  int nBins = bkgd->GetNbinsX();
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
    data->GetXaxis()->SetBinLabel(iB,binLabel[iB-1]);
    bkgd->GetXaxis()->SetBinLabel(iB,binLabel[iB-1]);
  }

  data->GetXaxis()->SetLabelSize(0.08);
  bkgd->GetXaxis()->SetLabelSize(0.08);

  if (logY) {
    data->GetYaxis()->SetRangeUser(1,20*data->GetMaximum());
    bkgd->GetYaxis()->SetRangeUser(1,20*bkgd->GetMaximum());
  }

  TCanvas * canv = new TCanvas("canv","",700,500);

  TH1F * errorBand = bkgd->Clone("errorBand");
  errorBand->SetMarkerSize(0);
  errorBand->SetFillColor(4);
  errorBand->SetFillStyle(3013);
  errorBand->SetLineWidth(2);


  
  float chi2 = 0;
  
  std::cout << std::endl;
  std::cout << std::endl;

  int nBinsIncluded = 0;

  float pullBins[15] = {2.7374e-01,
			-4.0241e-01,
			-1.6041e-01,
			2.3417e-01,
			-1.4423e-01,
			2.3732e-01,
			1.5813e-01,
			-2.9314e-01,
			5.6583e-01,
			2.3727e-01,
			0,0,0,0,0};

  float errBins[15] = {9.41e-01,
		       9.55e-01,
		       9.65e-01,
		       9.32e-01,
		       9.42e-01,
		       9.48e-01,
		       9.44e-01,
		       9.56e-01,
		       9.19e-01,
		       9.18e-01,
		       0,0,0,0,0};

  float shapeCentral = 2.7196e-01;
  float shapeErr = 4.82e-01;

  for (int iB=1; iB<=nBins; ++iB) {
    float eSys =  bkgdUp->GetBinContent(iB)-bkgd->GetBinContent(iB);
    float content = bkgd->GetBinContent(iB) + shapeCentral*eSys;
    float eSys2 = eSys*eSys*shapeErr*shapeErr;
    for (int iSys=1; iSys<=nBins; ++iSys) {
      char binChar[5];
      if (iSys<10)
	sprintf(binChar,"%1i",iSys);
      else
	sprintf(binChar,"%2i",iSys);
      TString binStr(binChar);
      TH1F * bkgBin = file->Get("bkgd_BkgBin"+binStr+"Up");
      float err = bkgBin->GetBinContent(iB) - bkgd->GetBinContent(iB);
      content += err*pullBins[iB-1];
      eSys2 += errBins[iB-1]*errBins[iB-1]*err*err;
    }
    errorBand->SetBinContent(iB,content);
    bkgd->SetBinContent(iB,content);
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

  //  data->Draw("e1");
  //  blindH->Draw("same");
  //  bkgd->Draw("same");
  //  errorBand->Draw("e2same");
  bkgd->Draw();
  sig->Draw("same");
  //  sig4->Draw("same");
  errorBand->Draw("e2same");
  data->Draw("e1same");
  TLegend * leg = new TLegend(0.55,0.65,0.95,0.92);
  leg->SetFillColor(0);
  leg->SetHeader("                  Postfit");
  leg->SetTextSize(0.045);
  leg->AddEntry(data,"Data","lp");
  leg->AddEntry(errorBand,"Bkgd(+uncertainty)","lf");
  //  leg->AddEntry(bkgd,"Bkgd","lf");
  //  leg->AddEntry(sig4,"Signal m(H_{1})=4GeV","l");
  leg->AddEntry(sig,"Signal m(#phi_{1}) = 8 GeV","l");
  //  leg->AddEntry(blindH,"blinded bins","f");
  leg->Draw();
  TPad * pad = canv->GetPad(0);
  //  pad->RedrawAxis();

  TLatex * cms = new TLatex(0.25,0.94,"CMS Preliminary   L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");

  cms->SetNDC();
  cms->SetTextSize(0.05);
  cms->Draw();

  canv->Update();
  if (logY)
    canv->SetLogy();
  canv->SetGridx();
  canv->SetGridy();
  canv->Print("Postfit.png");
  canv->Print("Postfit.pdf","Portrait pdf");

}
