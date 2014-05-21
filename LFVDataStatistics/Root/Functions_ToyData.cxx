
#include <iostream>
#include "LFVDataStatistics/Functions_ToyData.h"
#include "LFVDataStatistics/Constants.h"
#include "TRandom.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"

namespace ToyData
{

TH1D* appendSmoothTestGaussian(TH1D* h, double signalStrength, TString name)
{
	delete gDirectory->FindObject(name);
	TH1D* hplus = (TH1D*)h->Clone(name);
	hplus->Add(TEST_GAUSSIAN,signalStrength);
	return hplus;
}

TH1D* appendDataGaussian(TH1D* h, double signalStrength, TString name)
{
	delete gDirectory->FindObject(name);
	TH1D* hplus = (TH1D*)h->Clone(name);
	hplus->Add(DATA_GAUSSIAN,signalStrength);
	return hplus;
}

TH1D* appendSignal(TH1D* h, TH1D* hsig, double signalStrength, TString name)
{
	delete gDirectory->FindObject(name);
	TH1D* hplus = (TH1D*)h->Clone(name);
	hplus->Add(hsig,signalStrength);
	return hplus;
}

TH1D* drawFromHisto(TH1D* h_source,TString name)
{
	int nbins = h_source->GetXaxis()->GetNbins();
	double min = h_source->GetXaxis()->GetXmin();
	double max = h_source->GetXaxis()->GetXmax();

	delete gDirectory->FindObject(name);
	TH1D* h_new = new TH1D(name,name,nbins,min,max);
//	 cout<<"seed5 = "<<RAND.GetSeed()<<endl;
	for (int i=1; i<=nbins; i++){
		double v = RAND.Poisson(h_source->GetBinContent(i));
		if (v<0){v = 0;}
		h_new->SetBinContent(i,v);
	}
//	cout<<"seed6 = "<<RAND.GetSeed()<<endl;

	return h_new;
}

TH1D* drawFromHistoGaus(TH1D* h_source, TString name)
{
	int nbins = h_source->GetXaxis()->GetNbins();
	double min = h_source->GetXaxis()->GetXmin();
	double max = h_source->GetXaxis()->GetXmax();
	delete gDirectory->FindObject(name);
	TH1D* h_new = new TH1D(name,name,nbins,min,max);

	for (int i=1; i<=nbins; i++){
		double n = h_source->GetBinContent(i);
		double err = h_source->GetBinError(i);
		double v = RAND.Gaus(n,err);//sqrt(0.5*n));
		if (v<0){v = 0;}
		h_new->SetBinContent(i,v);

	}

	return h_new;
}

void setBErrors(TH1D* h_1, TH1D* h_2,TH1D* h_b)
{
	for (int i=1; i<=h_b->GetXaxis()->GetNbins(); i++){
		double n1 = h_1->GetBinContent(i);
		double n2 = h_2->GetBinContent(i);
		h_b->SetBinError(i,0.5*sqrt(n1+n2));
	}
}

TH1D* getBlindHisto(TH1D* h_source,int bin1, int bin2,TString name)
{
	const TArrayD* binsarray = h_source->GetXaxis()->GetXbins();
	const Double_t * massBins = binsarray->GetArray();
	int nbins = h_source->GetXaxis()->GetNbins();
	TH1D* h_blind = new TH1D(name,name,nbins,massBins);
   	for (int j=1; j<=bin1; j++){
   		h_blind->SetBinContent(j,h_source->GetBinContent(j));
   	}
   	for (int j=bin2; j<=nbins; j++){
   		h_blind->SetBinContent(j,h_source->GetBinContent(j));
   	}
   	return h_blind;
}
}
