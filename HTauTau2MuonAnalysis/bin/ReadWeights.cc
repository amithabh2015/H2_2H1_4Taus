#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "THStack.h"


gStyle->SetOptStat(00000);
gStyle->SetOptStat(kFalse);
gPad->SetOptStat(00000);
void ReadWeights(){

TString fileName("HPt_Spectrum.root");
TFile *file = TFile::Open(fileName);

 TH1F * Original = (TH1F*) file->Get("scale_nominal_step5/PythiaPtSpectrumHisto");

 TH1F * Nominal = (TH1F*) file->Get("scale_nominal_step5/RewgtPt");
 TH1F * SU = (TH1F*) file->Get("scale_up_step5/RewgtPt");
 TH1F * SD = (TH1F*) file->Get("scale_down_step5/RewgtPt");

Original->SetLineColor(kBlack);
Original->SetLineStyle(2);
Nominal->SetLineColor(kBlue);
SU->SetLineColor(kRed);
SD->SetLineColor(kGreen);

TCanvas *c;
c = new TCanvas("c", "Plots", 900, 900);
c->Divide(2,1);
c->cd(1);
Original->SetStats(00000);
Nominal->SetStats(00000);
SU->SetStats(00000);
SD->SetStats(00000);

Original->SetTitle("Pythia H_pT spectrum");
Nominal->SetTitle("Pythia + rwgt with HqT -nominal");
SU->SetTitle("Pythia + rwgt with HqT -ScaleUp");
SD->SetTitle("Pythia + rwgt with HqT -ScaleDown");



Original->GetXaxis()->SetTitle("p_{T} (GeV)");
Original->GetYaxis()->SetTitle("a.u");
Original->Draw();
Original->SetMaximum(1.1*Original->GetMaximum());
Nominal->Draw("same");
SU->Draw("same");
SD->Draw("same");

c->cd(2);
gPad->SetLogy();

Original->Draw();
Original->SetMaximum(1.5*Original->GetMaximum());
Nominal->Draw("same");
SU->Draw("same");
SD->Draw("same");
gPad->BuildLegend();





}
