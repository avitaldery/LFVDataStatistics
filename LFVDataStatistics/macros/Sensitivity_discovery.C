#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include <iostream>
#include <sstream>
using namespace std;
using namespace RooFit ;

#include "LFVDataStatistics/Functions_ToyData.h"
#include "LFVDataStatistics/Functions_BkgEstimator.h"
#include "LFVDataStatistics/Constants.h"

void Sensitivity_discovery(TString dataFile, TString signalFile)
{
   InitExterns();
   RAND.SetSeed(11);

   int Polydegree = 2;

   //get two histos
   TFile* f = new TFile(dataFile);
   TH1D* h_ME = (TH1D*)f->Get("nom/em_LFV_ME_noS_McollHiggs_Unblind");
   TH1D* h_EM = (TH1D*)f->Get("nom/em_LFV_EM_noS_McollHiggs_Unblind");

   //get corresponding signal histo
   TFile* fsig = new TFile(signalFile);
   TH1D* h_sig =(TH1D*)fsig->Get("nom/em_LFV_ME_noS_McollHiggs_Unblind");

   int nbins = h_ME->GetXaxis()->GetNbins();
   const TArrayD* binsarray = h_ME->GetXaxis()->GetXbins();
   const Double_t * bins = binsarray->GetArray();
   double muHatB[nbins+2+Polydegree];

   //get BG estimation
//   TH1D* h_EMPlusSig = ToyData::appendDataGaussian(h_EM,10,"h_EMPlusSig");
   TH1D* h_b = BkgEstimator::meanDataEstimator(h_EM,h_ME,muHatB,h_sig,Polydegree,"h_b");

   //create blind histos to draw
   	TH1D* h_ME_blind = ToyData::getBlindHisto(h_ME,8,12,"ME blind");
   	TH1D* h_EM_blind = ToyData::getBlindHisto(h_EM,8,12,"EM blind");

    //draw the two histos and the BG estimation
   TLegend* leg = new TLegend(0.5,0.7,0.7,0.9);
   leg->SetFillColor(kWhite);
   leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);

   TH1D* h_bPlusSig = ToyData::appendDataGaussian(h_b,1.9,"h_b_signal");

   TCanvas* c1 = new TCanvas("BG estimation","BG estimation",600,600);
   h_ME_blind->GetXaxis()->SetTitle("M_{Collinear} (GeV)");
   h_ME_blind->SetLineColor(kRed);h_ME_blind->SetLineWidth(2); h_EM_blind->SetLineColor(kBlue);h_b->SetLineColor(kBlack);
   h_EM_blind->SetLineWidth(2); h_b->SetLineWidth(2); h_bPlusSig->SetLineWidth(2);
   h_ME_blind->Draw("e1");h_EM_blind->Draw("e1 sames"); h_bPlusSig->SetLineColor(kGreen+2);h_bPlusSig->Draw("sames"); h_b->Draw("e1 sames");

   leg->AddEntry(h_ME_blind,"#mue","l"); leg->AddEntry(h_EM_blind,"e#mu","l"); leg->AddEntry(h_b,"BG estimation","l");
   leg->AddEntry(h_bPlusSig,"BR(h#rightarrow#tau#mu) = 10%","l");  leg->Draw();

   //generate brazil plot
   TH1D* h_vanilla = new TH1D("vanilla","3#sigma sensitivity for discovery;#mu",150,0,5);

   //Toy MC to generate sensitivity width
   int numMC = 1000;

   for (int i=0;i<numMC;i++){

	   TH1D* h_brand = ToyData::drawFromHistoGaus(h_b,"h_brand");
	   TH1D* h_EMrand = ToyData::drawFromHisto(h_brand,"h_EMrand");
	   TH1D* h_MErand = ToyData::drawFromHisto(h_brand,"h_MErand");

	   //find mu for which Lambda_mu is 3^2
	   double mu = BkgEstimator::muCL_discovery2(h_EMrand,h_MErand,h_b,h_sig);

	   h_vanilla->Fill(mu,1./numMC);
   }

   TCanvas* c2 = new TCanvas("sensitivity","sensitivity",600,600);

   h_vanilla->Draw();

   Double_t quantile;
   Double_t quantile_1sigma;
   Double_t quantile_1sigma2;
   Double_t quantile_2sigma;
   Double_t quantile_2sigma2;
   const Double_t prob = 0.5;
   const Double_t prob2 = 0.683;
   const Double_t prob3 = 0.317;
   const Double_t prob4 = 0.954;
   const Double_t prob5 = 0.046;
   h_vanilla->GetQuantiles(1,&quantile,&prob);
   h_vanilla->GetQuantiles(1,&quantile_1sigma,&prob2);
   h_vanilla->GetQuantiles(1,&quantile_1sigma2,&prob3);
   h_vanilla->GetQuantiles(1,&quantile_2sigma,&prob4);
   h_vanilla->GetQuantiles(1,&quantile_2sigma2,&prob5);

   TH1D* h_green = new TH1D("greenquantile","greenquantile",150,0,5);
   int bin_med = h_green->GetXaxis()->FindBin(quantile);
   int bin_1sig = h_green->GetXaxis()->FindBin(quantile_1sigma);

   for (int k=bin_med;k<=bin_1sig;k++){
	   h_green->SetBinContent(k,h_vanilla->GetBinContent(k));
   }
   h_green->SetFillColor(kGreen);
   h_green->Draw("sames");

   TH1D* h_green2 = new TH1D("greenquantile2","greenquantile2",150,0,5);
   int bin_1sig2 = h_green2->GetXaxis()->FindBin(quantile_1sigma2);

   for (int k=bin_1sig2;k<=bin_med;k++){
	   h_green2->SetBinContent(k,h_vanilla->GetBinContent(k));
   }
   h_green2->SetFillColor(kGreen);
   h_green2->Draw("sames");
   TH1D* h_yellow = new TH1D("yellowquantile","yellowquantile",150,0,5);
   int bin_2sig = h_yellow->GetXaxis()->FindBin(quantile_2sigma);

   for (int k=bin_1sig;k<=bin_2sig;k++){
	   h_yellow->SetBinContent(k,h_vanilla->GetBinContent(k));
   }
   h_yellow->SetFillColor(kYellow);
   h_yellow->Draw("sames");
   TH1D* h_yellow2 = new TH1D("yellowquantile2","yellowquantile2",150,0,5);
   int bin_2sig2 = h_yellow2->GetXaxis()->FindBin(quantile_2sigma2);

   for (int k=bin_2sig2;k<=bin_1sig2;k++){
	   h_yellow2->SetBinContent(k,h_vanilla->GetBinContent(k));
   }
   h_yellow2->SetFillColor(kYellow);
   h_yellow2->Draw("sames");

   //print values
   cout<<"median value = "<<quantile<<endl;
   cout<<"1sigma values = "<<quantile_1sigma2<<", "<<quantile_1sigma<<endl;
   cout<<"2sigma values = "<<quantile_2sigma2<<", "<<quantile_2sigma<<endl;

}
