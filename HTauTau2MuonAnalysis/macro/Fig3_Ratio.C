#include "HtoH.h"
#include "HttStylesNew.cc"
#include "CMS_lumi.C"
#include <string>
#include <vector>
#include <iostream>
#include "Rtypes.h"
#include "TH1.h"
#include "TH2.h"
#include <TROOT.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
// .x Fig3_Ratio.C("invMassTrackSoftestMuonNtrk1H","invMassTrackSoftestMuonNtrk23H","m_{#mu1, softest trk} (GeV)","1/N #times dN/dm (GeV^{-1})","Data_MCtest2",0,10,0.5,1,0)
// .x Fig3_Ratio.C("ptTrackSoftestMuonNtrk1H","ptTrackSoftestMuonNtrk23H","softest track p_{T} (GeV)","1/N #times dN/dp_{T} (GeV^{-1})","Data_MCtest2",0,50,0.5,1,0)
// .x Fig3_Ratio.C("etaTrackSoftestMuonNtrk1H","etaTrackSoftestMuonNtrk23H","softest track #eta","1/N #times dN/d#eta","Data_MCtest2",-2.5,2.5,0.5,1,0)
// .x Fig3_Ratio.C("dRTrackSoftestMuonNtrk1H","dRTrackSoftestMuonNtrk23H","#DeltaR(#mu, softest trk)","1/N #times dN/dR","Data_MCtest2",0,1,0.3,1,0)

// .x Fig3_Ratio.C("invMassTrackHardestMuonNtrk1H","invMassTrackHardestMuonNtrk23H","m_{#mu1, hardest trk} (GeV)","1/N #times dN/dm (GeV^{-1})","Data_MCtest2",0,10,0.5,1,0)
// .x Fig3_Ratio.C("ptTrackHardestMuonNtrk1H","ptTrackHardestMuonNtrk23H","hardest track p_{T} (GeV)","1/N #times dN/dp_{T} (GeV^{-1})","Data_MCtest2",0,50,0.5,1,0)
// .x Fig3_Ratio.C("etaTrackHardestMuonNtrk1H","etaTrackHardestMuonNtrk23H","softest track #eta","1/N #times dN/d#eta","Data_MCtest2",-2.5,2.5,0.5,1,0)
// .x Fig3_Ratio.C("dRTrackHardestMuonNtrk1H","dRTrackHardestMuonNtrk23H","#DeltaR(#mu, softest trk)","1/N #times dN/dR","Data_MCtest2",0,1,0.3,1,0)

// .x Fig3_Ratio.C("invMassTrackRandomMuonNtrk1H","invMassTrackRandomMuonNtrk23H","m_{#mu1, random trk} (GeV)","1/N #times dN/dm (GeV^{-1})","Data_MCtest2",0,10,0.5,1,0)
// .x Fig3_Ratio.C("ptTrackRandomMuonNtrk1H","ptTrackRandomMuonNtrk23H","random track p_{T} (GeV)","1/N #times dN/dp_{T} (GeV^{-1})","Data_MCtest2",0,50,0.5,1,0)
// .x Fig3_Ratio.C("etaTrackRandomMuonNtrk1H","etaTrackRandomMuonNtrk23H","random track #eta","1/N #times dN/d#eta","Data_MCtest2",-2.5,2.5,0.5,1,0)
// .x Fig3_Ratio.C("dRTrackRandomMuonNtrk1H","dRTrackRandomMuonNtrk23H","#DeltaR(#mu1, random trk)","1/N #times dN/dR","Data_MCtest2",0,1,0.3,1,0)

void Fig3_Ratio(TString histoName1 = "InvMassSoftestNtrk1H",
		TString histoName2 = "InvMassSoftestNtrk23H",
		TString xTitle = "m_{#mu1, softest track} (GeV)",
		TString yTitle = "1/N #times dN/dm (GeV^{-1})",
		TString fileName = "Data2015D",
		float xMin = 0,
		float xMax = 10,
		float yMax = 0.5,
		bool isNorm = true) {


  SetStyle();
  gStyle->SetOptStat(0000);
  // gStyle->SetErrorX(0.5);

  TFile * file = new TFile(fileName+".root");

  TH1F * hist1Old = (TH1F*)file->Get(histoName1);
  TH1F * hist2Old = (TH1F*)file->Get(histoName2);

  int nBins = hist1Old->GetNbinsX();
  float xmin = hist1Old->GetBinLowEdge(1);
  float xmax = hist1Old->GetBinLowEdge(nBins+1);

  std::cout << std::endl;
  std::cout << "Histogram " << histoName1 << " : " << "nbins = " << nBins
	    << " , min = " << xmin
	    << " , max = " << xmax << std::endl;
  std::cout << std::endl;

  float bins[201];
  int nBinsNew = nBins;
  std::cout << "New number of bins : ";
  std::cin >> nBinsNew;

  if (nBins % nBinsNew >0) { 
    std::cout << "new number of bins = " << nBinsNew 
	      << "  not multiple of " << nBins << std::endl;
    return;
  }
  float binWidth = (xmax-xmin)/float(nBinsNew);
  for (int iB=0; iB<=nBinsNew; ++iB)
    bins[iB] = xmin + float(iB)*binWidth;

  TH1F * hist1 = TH1toTH1(hist1Old,nBinsNew,bins,true,"_1_new");
  TH1F * hist2 = TH1toTH1(hist2Old,nBinsNew,bins,true,"_2_new");

  int nBins1 = hist1->GetNbinsX();
  float xLower = hist1->GetXaxis()->GetBinLowEdge(1);
  float xUpper = hist1->GetXaxis()->GetBinLowEdge(nBins1+1);

  float norm1x = hist1->GetSum();
  float norm2x = hist2->GetSum();

  std::cout << "ratio of normalizations >>> " << norm1x/norm2x << std::endl;

  float norm1 = hist1->GetEntries();
  float norm2 = hist2->GetEntries();

  std::cout << "ratio of entries        >>> " << norm1/norm2 << std::endl;
  std::cout << std::endl;


  if (isNorm) {
    for (int iB=1; iB<=nBins1; ++iB) {
      hist1->SetBinContent(iB,hist1->GetBinContent(iB)/norm1x); 
      hist2->SetBinContent(iB,hist2->GetBinContent(iB)/norm2x);
      hist1->SetBinError(iB,hist1->GetBinError(iB)/norm1x); 
      hist2->SetBinError(iB,hist2->GetBinError(iB)/norm2x);
    }
  }

  hist1->GetYaxis()->SetRangeUser(0,yMax);

  hist1->SetLineColor(2);
  hist1->SetMarkerColor(2);
  hist1->SetMarkerStyle(21);
  hist1->SetMarkerSize(1.7);
  hist1->SetLineWidth(3);

  hist2->SetLineColor(4);
  hist2->SetMarkerColor(4);
  hist2->SetMarkerStyle(20);
  hist2->SetMarkerSize(1.7);
  hist2->SetLineWidth(3);

  float mean1 = hist1->GetMean();
  float rms1 = hist1->GetRMS();

  float mean2 = hist2->GetMean();
  float rms2 = hist2->GetRMS();

  //  hist1->SetTitle("CMS Simulation");
  //  hist1->SetTitleSize(1);
  //  hist1->SetTitleOffset(0);
 
  hist1->GetXaxis()->SetTitle(xTitle);
  hist1->GetYaxis()->SetTitle(yTitle);
  hist1->GetXaxis()->SetTitleSize(0.06);
  hist1->GetYaxis()->SetTitleSize(0.06);
  hist1->GetXaxis()->SetTitleOffset(1.0);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist1->GetXaxis()->SetRangeUser(xMin,xMax);
  hist1->GetYaxis()->SetNdivisions(505);
  hist1->GetXaxis()->SetNdivisions(505);
  hist1->GetXaxis()->SetLabelSize(0);
  std::cout << "Label offset = " << hist1->GetYaxis()->GetLabelOffset() << std::endl;
  
  TH1F * ratioH = (TH1F*)hist1->Clone("ratioH");
  ratioH->SetMarkerColor(1);
  ratioH->SetMarkerStyle(20);
  ratioH->SetMarkerSize(1.5);
  ratioH->SetLineColor(1);
  ratioH->GetYaxis()->SetRangeUser(0.3,1.7);
  ratioH->GetYaxis()->SetNdivisions(505);
  ratioH->GetXaxis()->SetLabelFont(42);
  ratioH->GetXaxis()->SetLabelOffset(0.04);
  ratioH->GetXaxis()->SetLabelSize(0.14);
  ratioH->GetXaxis()->SetTitleSize(0.13);
  ratioH->GetXaxis()->SetTitleOffset(1.2);
  ratioH->GetYaxis()->SetTitle("ratio");
  ratioH->GetYaxis()->SetLabelFont(42);
  ratioH->GetYaxis()->SetLabelOffset(0.015);
  ratioH->GetYaxis()->SetLabelSize(0.11);
  ratioH->GetYaxis()->SetTitleSize(0.14);
  ratioH->GetYaxis()->SetTitleOffset(0.5);
  ratioH->GetXaxis()->SetTickLength(0.07);
  ratioH->GetYaxis()->SetTickLength(0.04);
  ratioH->GetYaxis()->SetLabelOffset(0.01);

  for (int iB=1; iB<=nBinsNew; ++iB) {
    float x1 = hist1->GetBinContent(iB);
    float x2 = hist2->GetBinContent(iB);
    if (x1>0&&x2>0) {
      float e1 = hist1->GetBinError(iB);
      float e2 = hist2->GetBinError(iB);
      float r1 = e1/x1;
      float r2 = e2/x2;
      float r = TMath::Sqrt(r1*r1+r2*r2);
      float ratio = x1/x2;
      float eratio = r * ratio;
      ratioH->SetBinContent(iB,ratio);
      ratioH->SetBinError(iB,eratio);
    }
    else {
      ratioH->SetBinContent(iB,1000);
    }
  }

  TCanvas * c1 = new TCanvas("c1","",700,800);
  // ------------>Primitives in pad: upper
  TPad *upper = new TPad("upper", "pad",0,0.30,1,1);
  upper->Draw();
  upper->cd();
  upper->SetFillColor(0);
  upper->SetBorderMode(0);
  upper->SetBorderSize(10);
  upper->SetTickx(1);
  upper->SetTicky(1);
  upper->SetLeftMargin(0.17);
  upper->SetRightMargin(0.05);
  upper->SetBottomMargin(0.02);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  upper->SetFrameFillStyle(0);
  upper->SetFrameLineStyle(0);
  upper->SetFrameLineWidth(2);
  upper->SetFrameBorderMode(0);
  upper->SetFrameBorderSize(10);
  
  hist1->Draw("e1");
  hist2->Draw("e1same");
  

  TLegend * leg = NULL;
  if (histoName1.Contains("etaTrack")) 
    leg = new TLegend(0.20,0.65,0.52,0.92);
  else
    leg = new TLegend(0.6,0.4,0.95,0.72);

  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  leg->SetBorderSize(0);
  if (fileName.Contains("Data2015D"))
    leg->SetHeader("Data");
  else
    leg->SetHeader("QCD multijets");
  leg->AddEntry(hist1,"#it{N}_{trk,2} = 1","lp");
  leg->AddEntry(hist2,"#it{N}_{trk,2} = 2,3","lp");
  leg->Draw();

  writeExtraText = true;
  extraText   = "Simulation";
  if (fileName.Contains("Data2015D"))
    extraText   = "";
  if (fileName.Contains("Data2015D"))
    CMS_lumi(upper,4,33); 
  else
    CMS_lumi(upper,0,33); 

  upper->Draw("SAME");
  upper->RedrawAxis();

  upper->Modified();
  upper->Update();
  c1->cd();

// ------------>Primitives in pad: lower
  TPad * lower = new TPad("lower", "pad",0,0,1,0.30);
  lower->Draw();
  lower->cd();
  lower->SetFillColor(0);
  lower->SetBorderMode(0);
  lower->SetBorderSize(10);
  lower->SetGridy();
  lower->SetTickx(1);
  lower->SetTicky(1);
  lower->SetLeftMargin(0.17);
  lower->SetRightMargin(0.05);
  lower->SetTopMargin(0.026);
  lower->SetBottomMargin(0.35);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  lower->SetFrameFillStyle(0);
  lower->SetFrameLineStyle(0);
  lower->SetFrameLineWidth(2);
  lower->SetFrameBorderMode(0);
  lower->SetFrameBorderSize(10);
  
  ratioH->Draw("e1");

  lower->Modified();
  lower->RedrawAxis();
  c1->cd();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);

  if (fileName.Contains("Data2015D")) {
    if (histoName1.Contains("Softest")) {
      c1->Print("FiguresPaper/DataSoftest.png");
      c1->Print("FiguresPaper/DataSoftest.pdf","Portrait pdf");
    }
    else {
      c1->Print("FiguresPaper/DataHardest.png");
      c1->Print("FiguresPaper/DataHardest.pdf","Portrait pdf");
    }
  }
  else {
    if (histoName1.Contains("Softest")) {
      c1->Print("FiguresPaper/MCSoftest.png");
      c1->Print("FiguresPaper/MCSoftest.pdf","Portrait pdf");
    }
    else {
      c1->Print("FiguresPaper/MCHardest.png");
      c1->Print("FiguresPaper/MCHardest.pdf","Portrait pdf");
    }
  }

}
