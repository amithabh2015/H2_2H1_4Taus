#include "HttStylesNew.cc"
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
#include "TMath.h"
#include "HtoHold.H"

void correlation(TString fileName = "Data_RegionA",
		 TString fileSignalName = "Signal_mH1_8",
		 bool bkgSideBand = true,
		 TString legName = "Region A",
		 float sigma = 6,
		 int binningOption = 0,
		 bool subtractSignal = true) {

  
  float scale4 = 0.874; 
  float scale8 = 0.941;

  float mH1 = 4;

  if (fileSignalName.Contains("mH1_5"))
    mH1 = 5;
  else if (fileSignalName.Contains("mH1_6"))
    mH1 = 6;
  else if (fileSignalName.Contains("mH1_7"))
    mH1 = 7;
  else if (fileSignalName.Contains("mH1_8"))
    mH1 = 8;
  

  float scale = scale4 + (mH1-4)*(scale8-scale4)/4;


  float lumi = 19700;
 
  TString suffix("");
  if (subtractSignal) 
    suffix = "_"+fileSignalName;

  int newNBinsX = 9;
  float newBins[11] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,10.,20.};
  if (binningOption==0) {
    newNBinsX = 4;
    newBins[0] = 0;
    newBins[1] = 1;
    newBins[2] = 2;
    newBins[3] = 3;
    newBins[4] = 10;
  }
  else if (binningOption==1) {
    newNBinsX = 5;
    newBins[0] = 0;
    newBins[1] = 1;
    newBins[2] = 2;
    newBins[3] = 3;
    newBins[4] = 4;
    newBins[5] = 10;
  }
  else if (binningOption==2) {
    newNBinsX = 6;
    newBins[0] = 0;
    newBins[1] = 1;
    newBins[2] = 2;
    newBins[3] = 3;
    newBins[4] = 4;
    newBins[5] = 6;
    newBins[6] = 10;
  }
  char binOpt[10];
  sprintf(binOpt,"_v%1i",binningOption);
  TString BinVersion(binOpt);
  
  TFile * f = new TFile (fileName+".root");
  TFile * fs = new TFile (fileSignalName+".root");

  SetStyle();
  gStyle->SetOptStat(0000);

  TFile * fd = new TFile("Data.root");
  
  TH1F * h1  = (TH1F*)fd->Get("InvMassTrackPlusPos_2QH");
  TH1F * h2  = (TH1F*)fd->Get("InvMassTrackPlusPos_3QH");
  TH1F * h3  = (TH1F*)fd->Get("InvMassTrackPlusNeg_2QH");
  TH1F * h4  = (TH1F*)fd->Get("InvMassTrackPlusNeg_3QH");
  
  h4->Add(h4,h1);
  h4->Add(h4,h2);
  h4->Add(h4,h3);
  
  TH1F * mass1DH = TH1toTH1(h4,newNBinsX,newBins,true,"_rebin1D");
  mass1DH->Scale(1/mass1DH->GetSumOfWeights());
  std::cout << std::endl;
  for (int iB=1; iB<=newNBinsX; ++iB) {
    std::cout << iB << " " << mass1DH->GetBinContent(iB) << std::endl;
  }
  
  if (bkgSideBand) {
    h1  = (TH1F*)f->Get("InvMassTrackPlusPos_2Q_ControlH");
    h2  = (TH1F*)f->Get("InvMassTrackPlusPos_3Q_ControlH");
    h3  = (TH1F*)f->Get("InvMassTrackPlusNeg_2Q_ControlH");
    h4  = (TH1F*)f->Get("InvMassTrackPlusNeg_3Q_ControlH");
  }
  
  TH1F * a1 = new TH1F("add1","",100,0.,10.);
  TH1F * a2 = new TH1F("add2","",100,0.,10.);
  //   TH1F * N1D_control_i_H= new TH1F("N1D_control_i_H","",100,0.,10.);
  
  //   a1->Add(h1,h2);
  //   a2->Add(h3,h4);
  //   N1D_control_i_H->Add(a1,a2);
  
  TH1F * N1D_control_i_H = (TH1F*)f->Get("InvMassTrackPlusMuon1D_ControlH");
  
  N1D_control_i_H->Scale(1./N1D_control_i_H->GetSumOfWeights());   
  // rebinning --->
  TH1F * N1Drebin = (TH1F*)TH1toTH1(N1D_control_i_H,newNBinsX,newBins,true,"_rebin");

//    float xBin;
//    float yBin;
//    float x1D = N1D_control_i_H->GetBinContent(xBin);
//    float y1D = N1D_control_i_H->GetBinContent(yBin);

   

   TH2F * True2DH = NULL;
   TH2F * True2DsignalH = NULL;

   TH2F * bkgd2DH = new TH2F("bkgd2DH","",newNBinsX,newBins,newNBinsX,newBins);

   if (bkgSideBand) {

     True2DH = new TH2F("True2DH","",100,0.,10.,100,0.,10.);
     True2DsignalH = new TH2F("True2DsignalH","",100,0.,10.,100,0.,10.);

     TH2F * H1  = (TH2F*)f->Get("InvMassTrackPlusMuon2D_ControlPosH");
     TH2F * H2  = (TH2F*)f->Get("InvMassTrackPlusMuon2D_ControlNegH");
     TH2F * H3  = (TH2F*)f->Get("InvMassTrackPlusMuon2D_ControlBothH");   
     
     float entriesBkg = H1->GetEntries() + H2->GetEntries() + H3->GetEntries();
     TH1F * counterFinalH = (TH1F*)f->Get("counterFinalH");
     float entriesSig = counterFinalH->GetEntries();
     float ratio = entriesBkg / entriesSig;     
     std::cout << std::endl;
     std::cout << "ratio(Bkg/Sig) = " << entriesBkg << "/" << entriesSig << " = " << ratio << std::endl; 

     TH2F * A1 = new TH2F("A1","",100,0.,10.,100,0.,10.);
     A1->Add(H1,H2);
     True2DH->Add(A1,H3);

     TH1F * inputs = (TH1F*)fs->Get("inputEventsH");

     float events = inputs->GetSum();

     TH2F * SigH1  = (TH2F*)fs->Get("InvMassTrackPlusMuon2D_ControlPosH");
     TH2F * SigH2  = (TH2F*)fs->Get("InvMassTrackPlusMuon2D_ControlNegH");
     TH2F * SigH3  = (TH2F*)fs->Get("InvMassTrackPlusMuon2D_ControlBothH");   
     TH2F * Sig    = (TH2F*)fs->Get("InvMassTrackPlusMuon2D_H");

     TH2F * tempSig = new TH2F("tempSig","",100,0.,10.,100,0.,10.);

     tempSig->Add(SigH1,SigH2);
     True2DsignalH->Add(tempSig,SigH3);

     float normSig = scale*sigma*lumi/events;
     float norm1pb = lumi/events;

     float acceptance = scale*True2DsignalH->GetSum()/events;
     float acceptanceSR = Sig->GetSum()/events;

     std::cout << std::endl;
     std::cout << "Acceptance (side-band)     = " << 100*acceptance << " %" << std::endl;
     std::cout << "Acceptance (signal region) = " << 100*acceptanceSR << " %" << std::endl;

     for (int iB=1; iB<=100; ++iB) {
       for (int jB=1; jB<=100; ++jB) {
	 float xHisto = True2DsignalH->GetBinContent(iB,jB);
	 float eHisto = True2DsignalH->GetBinError(iB,jB);
	 True2DsignalH->SetBinContent(iB,jB,xHisto*normSig);
	 True2DsignalH->SetBinError(iB,jB,eHisto*normSig);
	 xHisto = Sig->GetBinContent(iB,jB);
	 Sig->SetBinContent(iB,jB,xHisto*norm1pb);
       }
     }

     std::cout << "Number of signal events = " << True2DsignalH->GetSum();
     std::cout << std::endl;

   }
   else {

     True2DH = (TH2F*)f->Get("InvMassTrackPlusMuon2D_H");
     

   }

   TH2F * True2Drebin = (TH2F*)TH2toTH2(True2DH,newNBinsX,newBins,newNBinsX,newBins,true,"_rebin");
   TH2F * True2DrebinNoNorm = (TH2F*)TH2toTH2(True2DH,newNBinsX,newBins,newNBinsX,newBins,true,"_rebinNoNorm");
   TH2F * True2DrebinSig = (TH2F*)TH2toTH2(True2DsignalH,newNBinsX,newBins,newNBinsX,newBins,true,"_rebinNoNorm");
   TH2F * Sig2D = (TH2F*)TH2toTH2(Sig,newNBinsX,newBins,newNBinsX,newBins,false,"_SigRebinned");

   True2DrebinSig->SetBinContent(4,4,True2DrebinSig->GetBinContent(4,4)-1);
   True2DrebinSig->SetBinContent(3,3,True2DrebinSig->GetBinContent(3,3)+0.5);
   True2DrebinSig->SetBinContent(3,4,True2DrebinSig->GetBinContent(3,4)-1);
   True2DrebinSig->SetBinContent(4,3,True2DrebinSig->GetBinContent(4,3)-1);
   True2DrebinSig->SetBinContent(2,4,True2DrebinSig->GetBinContent(2,4)-0.8);
   True2DrebinSig->SetBinContent(4,2,True2DrebinSig->GetBinContent(4,2)-0.8);

   True2DrebinSig->SetBinContent(3,3,True2DrebinSig->GetBinContent(3,3)+0.5);
   True2DrebinSig->SetBinContent(2,2,True2DrebinSig->GetBinContent(1,4)+0.25);
   True2DrebinSig->SetBinContent(2,2,True2DrebinSig->GetBinContent(4,1)+0.25);


   if (subtractSignal) {
     for (int iB=1; iB<=newNBinsX; ++iB) {
       for (int jB=1; jB<=newNBinsX; ++jB) {
	 float X1 = True2Drebin->GetBinContent(iB,jB);
	 float X2 = True2DrebinSig->GetBinContent(iB,jB);
	 float E  = True2Drebin->GetBinError(iB,jB);
	 float X = X1 - X2;
	 if (X<0) X = 0;
	 printf("[%1i,%1i] : %4i - %4.2f = %4i\n",iB,jB,floor(X1+0.5),X2,floor(X+0.5));
	 True2Drebin->SetBinContent(iB,jB,floor(X+0.5));
	 True2Drebin->SetBinError(iB,jB,E);
	 True2DrebinNoNorm->SetBinContent(iB,jB,floor(X+0.5));
	 True2DrebinNoNorm->SetBinError(iB,jB,E);
       }
     }
   }
   float Norm = True2Drebin->GetSumOfWeights();

   std::cout << std::endl;
   std::cout << "Entries  (True2D) = " << True2DH->GetEntries() << std::endl;
   std::cout << "Norm (Tru2Drebin) = " << Norm << std::endl;
   
   TH2F * True2DrebinTmp = (TH2F*)True2Drebin->Clone("Tmp2D");
   TH2F * True2DrebinSigTmp = (TH2F*)True2DrebinSig->Clone("SignalTmp2D");
   TH2F * Contamination2D = (TH2F*)True2DrebinSig->Clone("Contamination2D");
   TH2F * Sig2DTmp = (TH2F*)Sig2D->Clone("Sig2D");
   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=1; jB<=newNBinsX; ++jB) {
       if (iB!=jB) {
	 // bkgd 
	 float c1 = True2DrebinTmp->GetBinContent(iB,jB);
	 float c2 = True2DrebinTmp->GetBinContent(jB,iB);
	 float e1 = True2DrebinTmp->GetBinError(iB,jB);
	 float e2 = True2DrebinTmp->GetBinError(jB,iB);
	 float c = c1 + c2;
	 float e = TMath::Sqrt(e1*e1+e2*e2);
	 True2Drebin->SetBinContent(iB,jB,c);
	 True2Drebin->SetBinError(iB,jB,e);
	 True2DrebinNoNorm->SetBinContent(iB,jB,c);
	 True2DrebinNoNorm->SetBinError(iB,jB,e);
	 // signal
	 c1 = True2DrebinSigTmp->GetBinContent(iB,jB);
	 c2 = True2DrebinSigTmp->GetBinContent(jB,iB);
	 e1 = True2DrebinSigTmp->GetBinError(iB,jB);
	 e2 = True2DrebinSigTmp->GetBinError(jB,iB);
	 float cSig = c1 + c2;
	 float eSig = TMath::Sqrt(e1*e1+e2*e2);
	 True2DrebinSig->SetBinContent(iB,jB,cSig);
         True2DrebinSig->SetBinError(iB,jB,eSig);
	 c1 = Sig2DTmp->GetBinContent(iB,jB);
	 c2 = Sig2DTmp->GetBinContent(jB,iB);
	 float sig = c1 + c2;
	 Sig2D->SetBinContent(iB,jB,sig);
       }
       if (True2Drebin->GetBinContent(iB,jB)<1e-10) True2Drebin->SetBinError(iB,jB,1.0);
       True2Drebin->SetBinContent(iB,jB,True2Drebin->GetBinContent(iB,jB)/Norm);
       True2Drebin->SetBinError(iB,jB,True2Drebin->GetBinError(iB,jB)/Norm);
       float numerator = True2DrebinSig->GetBinContent(iB,jB);
       float denominator = True2DrebinNoNorm->GetBinContent(iB,jB);
       float ratio = 1;
       Sig2D->SetBinError(iB,jB,0);
       if (denominator>0) 
	 ratio = numerator / denominator;
       Contamination2D->SetBinContent(iB,jB,ratio);
       Contamination2D->SetBinError(iB,jB,0);
     }
   }

   for (int iB=1; iB<newNBinsX; ++iB) {
     for (int jB=iB+1; jB<=newNBinsX; ++jB) {
       True2DrebinNoNorm->SetBinContent(jB,iB,0);
       True2DrebinNoNorm->SetBinError(jB,iB,0);
       True2DrebinSig->SetBinContent(jB,iB,0);
       True2DrebinSig->SetBinError(jB,iB,0);
       Contamination2D->SetBinContent(jB,iB,0);
       Contamination2D->SetBinError(jB,iB,0);
       Sig2D->SetBinContent(jB,iB,0);
       Sig2D->SetBinError(jB,iB,0);
     }
   }


   TH2F * prod2DH= new TH2F("InvMass2DTrackPlusProbe_2H","",newNBinsX,newBins,newNBinsX,newBins);
   for(int ix=1; ix<=newNBinsX; ix++){
    for(int iy=1; iy<=newNBinsX; iy++){
      float x = N1Drebin->GetBinContent(ix);
      float y = N1Drebin->GetBinContent(iy);
      float prod = x * y;
      if (ix==iy)
	prod2DH->SetBinContent(ix,iy,prod);
      else
	prod2DH->SetBinContent(ix,iy,2*prod);
    }
   }


   TH2F * Corr2DH = new TH2F("Corr2DH","",newNBinsX,newBins,newNBinsX,newBins);
   
   int nBins1D = (newNBinsX+1) * newNBinsX / 2;
   TH1F * Corr1DH = new TH1F("Corr1DH","",nBins1D,0.,float(nBins1D));

   int counter = 0;
   for(int ix=1; ix<=newNBinsX; ix++){
    for(int iy=ix; iy<=newNBinsX; iy++){
      float XC = prod2DH->GetBinContent(ix,iy);
      float true2DC = True2Drebin->GetBinContent(ix,iy);
      float err2DC  = True2Drebin->GetBinError(ix,iy);
      float corrC = 1;
      float errC = 0;
      if (XC>0) 
	corrC = true2DC / XC ;

      if (true2DC>0)
	errC = corrC * err2DC / true2DC;
      else 
	errC = err2DC / XC;

      Corr2DH->SetBinContent(ix,iy,corrC);
      Corr2DH->SetBinError(ix,iy,errC);
      Corr2DH->SetBinContent(iy,ix,corrC);
      Corr2DH->SetBinError(iy,ix,errC);
      counter++;
      Corr1DH->SetBinContent(counter,corrC);
      Corr1DH->SetBinError(counter,errC);
      char strLab[20];
      sprintf(strLab,"(%1i,%1i)",ix,iy);
      Corr1DH->GetXaxis()->SetBinLabel(counter,strLab);
    }
   }   

   // 2D correlations : rounding numbers for plotting
   TH2F * Corr2DplotH = (TH2F*)Corr2DH->Clone("Corr2DplotH");
   std::cout << std::endl;
   std::cout << "Correlation coefficients -> " << std::endl;
   for(int ix=1; ix<=newNBinsX; ix++){
     for(int iy=1; iy<=newNBinsX; iy++){

       float corrC = Corr2DplotH->GetBinContent(ix,iy);
       float errC  = Corr2DplotH->GetBinError(ix,iy);

      char strLab[20];
      sprintf(strLab,"(%1i,%1i)",ix,iy);
       char strVal[20];
       sprintf(strVal,"%6.2f +/- %6.2f",corrC,errC); 
       std::cout << strLab << " : " << strVal << std::endl;

       if (corrC<10) {
	 corrC = floor(100*corrC+0.5)/100;
	 errC  = floor(100*errC+0.5)/100;
       }
       else {
	 corrC = floor(corrC+0.5);
	 errC  = floor(errC+0.5);
       }
       
       Corr2DplotH->SetBinContent(ix,iy,corrC);
       Corr2DplotH->SetBinError(ix,iy,errC);

     }
   }
   std::cout << std::endl;

   // Producing background estimates
   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {
       float x = mass1DH->GetBinContent(iB);
       float y = mass1DH->GetBinContent(jB);
       float prod = x * y * Corr2DH->GetBinContent(iB,jB);
       float error = x * y * Corr2DH->GetBinError(iB,jB);
       if (iB!=jB) {
   	 prod  = 2*prod;
   	 error = 2*error;
       }
       bkgd2DH->SetBinContent(iB,jB,prod);
       bkgd2DH->SetBinError(iB,jB,error);
     }
   }

   float normBkgd2D = bkgd2DH->GetSum();
   float dataYield  = 873;
   float scaleYield = dataYield/normBkgd2D;
   TH2F * bkgd2DnoRoundH = (TH2F*)bkgd2DH->Clone("bkgd2DnoRoundH");
   // rounding errors
   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {

       float x = scaleYield*bkgd2DH->GetBinContent(iB,jB);
       float e = scaleYield*bkgd2DH->GetBinError(iB,jB);

       bkgd2DnoRoundH->SetBinContent(iB,jB,x);
       bkgd2DnoRoundH->SetBinError(iB,jB,e);

       if (x<1) {
	 x = floor(100*x+0.5)/100; 
       	 e = floor(100*e+0.5)/100;
       }
       else if (x<100) {
       	 x = floor(10*x+0.5)/10; 
       	 e = floor(10*e+0.5)/10;
       }
       else {
       	 x = floor(x+0.5);
         e = floor(e+0.5);
       }

       bkgd2DH->SetBinContent(iB,jB,x);
       bkgd2DH->SetBinError(iB,jB,e);

     }
   }

   // *********************************************
   // canv -> unrolled distribution of correlations
   // *********************************************
   // TCanvas * canv = new TCanvas("canv","Canvas",1000,700); 
   // Corr1DH->GetXaxis()->SetLabelSize(0.07);
   // Corr1DH->SetMarkerSize(1.4);
   // Corr1DH->GetYaxis()->SetTitle("Correlation coefficient");
   // Corr1DH->GetXaxis()->SetTitle("Bin");
   // Corr1DH->GetYaxis()->SetTitleOffset(1);
   // Corr1DH->GetXaxis()->SetTitleOffset(1);
   // Corr1DH->SetLineWidth(2);
   // Corr1DH->GetYaxis()->SetRangeUser(0.,2.);
   // Corr1DH->Draw("e1");
   // canv->SetGridx();
   // canv->SetGridy();
   // canv->Update();
   // canv->Print("figures/Counts/"+fileName+"_1Dcorr"+suffix+BinVersion+".png");

   // for (int iB=1; iB<newNBinsX; ++iB) {
   //   for (int jB=iB+1; jB<=newNBinsX; ++jB) {
   //     Corr2DH->SetBinContent(iB,jB,0);
   //   }
   // }


   // ************************
   // canv1 -> 2D correlations
   // ************************

   TCanvas * canv1 = new TCanvas("canv1","Canvas",900,700);
   Corr2DplotH->GetXaxis()->SetTitle("m_{1} [GeV]");
   Corr2DplotH->GetYaxis()->SetTitle("m_{2} [GeV]");
   Corr2DplotH->GetXaxis()->SetTitleOffset(1.1);
   Corr2DplotH->GetYaxis()->SetTitleOffset(1.1);
   Corr2DplotH->GetXaxis()->SetTickLength(0.01);
   Corr2DplotH->GetYaxis()->SetTickLength(0.01);
   
   Corr2DplotH->GetXaxis()->SetNdivisions(110);
   Corr2DplotH->GetYaxis()->SetNdivisions(110);
   Corr2DplotH->Draw("texte");
   Corr2DplotH->SetMarkerSize(1.8);
   int nMax = newNBinsX - 1;
   for (int i=1; i<=nMax; ++i) {
     TLine * lineX = new TLine(0., newBins[i], 10., newBins[i]); 
     lineX->SetLineStyle(1);
     lineX->SetLineWidth(2);
     lineX->Draw();
     TLine * lineY = new TLine(newBins[i], 0., newBins[i], 10.);  
     lineY->SetLineStyle(1);
     lineY->SetLineWidth(2);
     lineY->Draw();
   }
   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=1; jB<=newNBinsX; ++jB) {
       float xCorr = Corr2DplotH->GetBinContent(iB,jB);
       float eCorr = Corr2DplotH->GetBinError(iB,jB);
       if (xCorr<1e-10) {
         char fff[20];
         if (eCorr<10) {
           sprintf(fff,"0^{+%3.1f}",eCorr);

         }
         else if {
           sprintf(fff,"0^{+%2.0f}",eCorr);
         }
         float xcoor = 0.5*(Corr2DplotH->GetXaxis()->GetBinLowEdge(iB)+Corr2DplotH->GetXaxis()->GetBinLowEdge(iB+1));
         float ycoor = 0.5*(Corr2DplotH->GetYaxis()->GetBinLowEdge(jB)+Corr2DplotH->GetYaxis()->GetBinLowEdge(jB+1));
         TLatex * binContent = new TLatex(xcoor-0.3,ycoor-0.2,TString(fff));
         binContent->SetTextSize(0.04);
         binContent->Draw();
       }
     }
   }

   TLatex * cms = new TLatex(0.16,0.94,"CMS Preliminary   L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
   
   cms->SetNDC();
   cms->SetTextSize(0.04);
   cms->Draw();

   TLegend * leg = new TLegend(0.72,0.92,0.95,1.0);
   leg->SetFillColor(0);
   leg->SetTextSize(0.07);
   leg->SetHeader(legName);
   TPad * pad1 = canv1->GetPad(0);
   pad1->RedrawAxis();
   leg->Draw();
   canv1->Update();
   canv1->Print("figures/Counts/"+fileName+"_2Dcorr"+suffix+BinVersion+".png");
   canv1->Print("figures/Counts/"+fileName+"_2Dcorr"+suffix+BinVersion+".pdf","Portrait pdf");


   // **************************************
   // canv2 -> true data counts in side-band
   // **************************************

   std::cout << std::endl;
   std::cout << "Sideband : data counts in bins = " << True2DrebinNoNorm->GetSumOfWeights() << std::endl;
   std::cout << "Sideband : total data counts   = " << True2DrebinNoNorm->GetSum() << std::endl;
   std::cout << std::endl;

   TCanvas * canv2 = new TCanvas("canv2","Canvas",900,700);
   True2DrebinNoNorm->GetXaxis()->SetTitle("m_{1} [GeV]");
   True2DrebinNoNorm->GetYaxis()->SetTitle("m_{2} [GeV]");
   True2DrebinNoNorm->GetXaxis()->SetTitleOffset(1.1);
   True2DrebinNoNorm->GetYaxis()->SetTitleOffset(1.1);
   True2DrebinNoNorm->GetXaxis()->SetTickLength(0.01);
   True2DrebinNoNorm->GetYaxis()->SetTickLength(0.01);
   
   True2DrebinNoNorm->GetXaxis()->SetNdivisions(110);
   True2DrebinNoNorm->GetYaxis()->SetNdivisions(110);
   True2DrebinNoNorm->Draw("text");
   True2DrebinNoNorm->SetMarkerSize(1.8);
   int nMax = newNBinsX - 1;
   for (int i=1; i<=nMax; ++i) {
     TLine * lineX = new TLine(0., newBins[i], 10., newBins[i]); 
     lineX->SetLineStyle(1);
     lineX->SetLineWidth(2);
     lineX->Draw();
     TLine * lineY = new TLine(newBins[i], 0., newBins[i], 10.);  
     lineY->SetLineStyle(1);
     lineY->SetLineWidth(2);
     lineY->Draw();
   }

   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {
       float xCount = True2DrebinNoNorm->GetBinContent(iB,jB);
       printf("[%1i,%1i] : %4i\n",iB,jB,int(xCount+0.01));
       if (xCount<1e-10) {
         float xcoor = 0.5*(Corr2DplotH->GetXaxis()->GetBinLowEdge(iB)+Corr2DplotH->GetXaxis()->GetBinLowEdge(iB+1));
         float ycoor = 0.5*(Corr2DplotH->GetYaxis()->GetBinLowEdge(jB)+Corr2DplotH->GetYaxis()->GetBinLowEdge(jB+1));
         TLatex * binContent = new TLatex(xcoor-0.1,ycoor-0.2,"0");
         binContent->SetTextSize(0.04);
         binContent->Draw();
       }
     }
   }

   cms->Draw();

   TPad * pad2 = canv2->GetPad(0);
   pad2->RedrawAxis();
   leg->Draw();
   canv2->Update();
   canv2->Print("figures/Counts/"+fileName+"_2D"+suffix+BinVersion+".png");
   canv2->Print("figures/Counts/"+fileName+"_2D"+suffix+BinVersion+".pdf","Portrait pdf");


   // ******************************
   // canv3 -> true counts in signal
   // ******************************

   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=1; jB<=newNBinsX; ++jB) {
       float xHist = True2DrebinSig->GetBinContent(iB,jB);
       float yHist = floor(100*xHist+0.5)/100;
       True2DrebinSig->SetBinContent(iB,jB,yHist);
     }
   }
 
   std::cout << std::endl;
   std::cout << "Signal counts in bins = " << True2DrebinSig->GetSumOfWeights() << std::endl;
   std::cout << "Total signal counts = " << True2DrebinSig->GetSum() << std::endl;
   std::cout << std::endl;

   TCanvas * canv3 = new TCanvas("canv3","Canvas",900,700);
   True2DrebinSig->GetXaxis()->SetTitle("m_{1} [GeV]");
   True2DrebinSig->GetYaxis()->SetTitle("m_{2} [GeV]");
   True2DrebinSig->GetXaxis()->SetTitleOffset(1.1);
   True2DrebinSig->GetYaxis()->SetTitleOffset(1.1);
   True2DrebinSig->GetXaxis()->SetTickLength(0.01);
   True2DrebinSig->GetYaxis()->SetTickLength(0.01);
   
   True2DrebinSig->GetXaxis()->SetNdivisions(110);
   True2DrebinSig->GetYaxis()->SetNdivisions(110);
   True2DrebinSig->Draw("text");
   True2DrebinSig->SetMarkerSize(1.8);
   int nMax = newNBinsX - 1;
   for (int i=1; i<=nMax; ++i) {
     TLine * lineX = new TLine(0., newBins[i], 10., newBins[i]); 
     lineX->SetLineStyle(1);
     lineX->SetLineWidth(2);
     lineX->Draw();
     TLine * lineY = new TLine(newBins[i], 0., newBins[i], 10.);  
     lineY->SetLineStyle(1);
     lineY->SetLineWidth(2);
     lineY->Draw();
   }

   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {
       float xCount = True2DrebinSig->GetBinContent(iB,jB);
       if (xCount<1e-10) {
         float xcoor = 0.5*(Corr2DplotH->GetXaxis()->GetBinLowEdge(iB)+Corr2DplotH->GetXaxis()->GetBinLowEdge(iB+1));
         float ycoor = 0.5*(Corr2DplotH->GetYaxis()->GetBinLowEdge(jB)+Corr2DplotH->GetYaxis()->GetBinLowEdge(jB+1));
         TLatex * binContent = new TLatex(xcoor-0.1,ycoor-0.2,"0");
         binContent->SetTextSize(0.04);
         binContent->Draw();
       }
     }
   }

   cms->Draw();

   TPad * pad3 = canv3->GetPad(0);
   pad3->RedrawAxis();
   leg->Draw();
   canv3->Update();
   canv3->Print("figures/Counts/"+fileSignalName+"_2D"+suffix+BinVersion+".png");
   canv3->Print("figures/Counts/"+fileSignalName+"_2D"+suffix+BinVersion+".pdf","Portrait pdf");


   // *******************************
   // canv4 -> background predictions
   // *******************************

   TCanvas * canv4 = new TCanvas("canv4","Canvas",900,700);
   bkgd2DH->GetXaxis()->SetTitle("m_{1} [GeV]");
   bkgd2DH->GetYaxis()->SetTitle("m_{2} [GeV]");
   bkgd2DH->GetXaxis()->SetTitleOffset(1.1);
   bkgd2DH->GetYaxis()->SetTitleOffset(1.1);
   bkgd2DH->GetXaxis()->SetTickLength(0.01);
   bkgd2DH->GetYaxis()->SetTickLength(0.01);
   
   bkgd2DH->GetXaxis()->SetNdivisions(110);
   bkgd2DH->GetYaxis()->SetNdivisions(110);
   bkgd2DH->Draw("texte");
   bkgd2DH->SetMarkerSize(1.8);
   int nMax = newNBinsX - 1;
   for (int i=1; i<=nMax; ++i) {
     TLine * lineX = new TLine(0., newBins[i], 10., newBins[i]); 
     lineX->SetLineStyle(1);
     lineX->SetLineWidth(2);
     lineX->Draw();
     TLine * lineY = new TLine(newBins[i], 0., newBins[i], 10.);  
     lineY->SetLineStyle(1);
     lineY->SetLineWidth(2);
     lineY->Draw();
   }
   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {
       float xCorr = bkgd2DH->GetBinContent(iB,jB);
       float eCorr = bkgd2DH->GetBinError(iB,jB);
       if (xCorr<1e-10) {
         char fff[20];
         if (eCorr<10) {
           sprintf(fff,"0^{+%4.2f}",eCorr);

         }
         else if {
           sprintf(fff,"0^{+%2.0f}",eCorr);
         }
         float xcoor = 0.5*(bkgd2DH->GetXaxis()->GetBinLowEdge(iB)+bkgd2DH->GetXaxis()->GetBinLowEdge(iB+1));
         float ycoor = 0.5*(bkgd2DH->GetYaxis()->GetBinLowEdge(jB)+bkgd2DH->GetYaxis()->GetBinLowEdge(jB+1));
         TLatex * binContent = new TLatex(xcoor-0.3,ycoor-0.2,TString(fff));
         binContent->SetTextSize(0.04);
         binContent->Draw();
       }
     }
   }

   std::cout << std::endl;
   std::cout << "Total background prediction = " << bkgd2DH->GetSum() << std::endl;
   std::cout << "Total background prediction (no rounding) = " << bkgd2DnoRoundH->GetSum() << std::endl;

   TLatex * cms = new TLatex(0.16,0.94,"CMS Preliminary   L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
   
   cms->SetNDC();
   cms->SetTextSize(0.04);
   cms->Draw();

   TLegend * leg = new TLegend(0.72,0.92,0.95,1.0);
   leg->SetFillColor(0);
   leg->SetTextSize(0.07);
   leg->SetHeader(legName);
   TPad * pad4 = canv4->GetPad(0);
   pad4->RedrawAxis();
   leg->Draw();
   canv4->Update();
   canv4->Print("figures/Counts/"+fileName+"_bkg2D"+suffix+BinVersion+".png");
   canv4->Print("figures/Counts/"+fileName+"_bkg2D"+suffix+BinVersion+".pdf","Portrait pdf");

   // **********************************************
   // canv5 -> true counts in signal (signal region)
   // **********************************************

   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=1; jB<=newNBinsX; ++jB) {
       float xHist = Sig2D->GetBinContent(iB,jB);
       float yHist = floor(100*xHist+0.5)/100;
       Sig2D->SetBinContent(iB,jB,yHist);
     }
   }
 
   std::cout << std::endl;
   std::cout << "Signal counts in bins = " << Sig2D->GetSumOfWeights() << std::endl;
   std::cout << "Total signal counts = " << Sig2D->GetSum() << std::endl;
   std::cout << std::endl;
   std::cout << std::endl;

   TCanvas * canv5 = new TCanvas("canv5","Canvas",900,700);
   Sig2D->GetXaxis()->SetTitle("m_{1} [GeV]");
   Sig2D->GetYaxis()->SetTitle("m_{2} [GeV]");
   Sig2D->GetXaxis()->SetTitleOffset(1.1);
   Sig2D->GetYaxis()->SetTitleOffset(1.1);
   Sig2D->GetXaxis()->SetTickLength(0.01);
   Sig2D->GetYaxis()->SetTickLength(0.01);
   
   Sig2D->GetXaxis()->SetNdivisions(110);
   Sig2D->GetYaxis()->SetNdivisions(110);
   Sig2D->Draw("text");
   Sig2D->SetMarkerSize(1.8);
   int nMax = newNBinsX - 1;
   for (int i=1; i<=nMax; ++i) {
     TLine * lineX = new TLine(0., newBins[i], 10., newBins[i]); 
     lineX->SetLineStyle(1);
     lineX->SetLineWidth(2);
     lineX->Draw();
     TLine * lineY = new TLine(newBins[i], 0., newBins[i], 10.);  
     lineY->SetLineStyle(1);
     lineY->SetLineWidth(2);
     lineY->Draw();
   }

   for (int iB=1; iB<=newNBinsX; ++iB) {
     for (int jB=iB; jB<=newNBinsX; ++jB) {
       float xCount = Sig2D->GetBinContent(iB,jB);
       if (xCount<1e-10) {
         float xcoor = 0.5*(Corr2DplotH->GetXaxis()->GetBinLowEdge(iB)+Corr2DplotH->GetXaxis()->GetBinLowEdge(iB+1));
         float ycoor = 0.5*(Corr2DplotH->GetYaxis()->GetBinLowEdge(jB)+Corr2DplotH->GetYaxis()->GetBinLowEdge(jB+1));
         TLatex * binContent = new TLatex(xcoor-0.1,ycoor-0.2,"0");
         binContent->SetTextSize(0.04);
         binContent->Draw();
       }
     }
   }

   cms->Draw();

   TPad * pad5 = canv5->GetPad(0);
   pad5->RedrawAxis();
   leg->Draw();
   canv5->Update();
   canv5->Print("figures/Counts/"+fileSignalName+"_2D_SR"+BinVersion+".png");
   canv5->Print("figures/Counts/"+fileSignalName+"_2D_SR"+BinVersion+".pdf","Portrait pdf");

   // ****************
   // saving file --->
   // ****************
   TFile * outputFile = new TFile("Correlations/"+fileName+"_corr"+BinVersion+suffix+".root","recreate");
   outputFile->cd("");
   TH1F * binsH = new TH1F("binsH","",newNBinsX,newBins);
   Corr1DH->Clone("Corr_1D");
   Corr2DH->Clone("Corr_2D");

   outputFile->Write();
   outputFile->Close();

}
   




