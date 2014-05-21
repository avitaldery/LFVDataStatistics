#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TLine.h"
#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <stdio.h>
#include "TMinuit.h"

#include "LFVDataStatistics/Functions_ToyData.h"
#include "LFVDataStatistics/Functions_BkgEstimator.h"
#include "LFVDataStatistics/Constants.h"


using namespace std;
using namespace ToyData;

namespace BkgEstimator{


TMinuit* minuit = 0;
TMinuit* minuit2 = 0;
TMinuit* minuit3 = 0;

void Newf1(int &npar, double *gin, double &f, double *par, int iflag) {
	//par[0-degree] are the polynomial coefficients

	//get histograms and polynomial degree
	HistoObject *data;
	data = (HistoObject*)minuit2->GetObjectFit();
	TH1D *h1 = data->m_h1;
	TH1D *h2 = data->m_h2;
	int degree = data->m_degree;
	int nbins = h1->GetXaxis()->GetNbins();

	// go blind
	//h1 = BkgEstimator::GetBlindHisto(h1);
	//h2 = BkgEstimator::GetBlindHisto(h2);

	f = 0;
	//add poly to h1 and calculate poisson for side bands
	for (int i=1; i<=7; i++){
		double n1 = h1->GetBinContent(i);
		double n2 = h2->GetBinContent(i);
		double Mcoll = h1->GetBinCenter(i);
		double Width = h1->GetXaxis()->GetBinWidth(i);
		double g = 0;
		for (int deg = 0; deg <= degree; deg++){
			g += par[deg] * TMath::Power(Mcoll, deg) * Width;
		}
		//double b = (n1+g+n2)/2;
		//f += -2*(TMath::Log(TMath::Poisson(b,n1+g))+TMath::Log(TMath::Poisson(b,n2)));
//		cout<<"newf1: for i = "<<i<<", Mcoll = "<<Mcoll<<",width = "<<Width<<", g = "<<g<<endl;
		f += -2*(TMath::Log(TMath::Poisson(n2,n1+g)));
	}

	for (int i=13; i<=nbins; i++){
		double n1 = h1->GetBinContent(i);
		double n2 = h2->GetBinContent(i);
		double Mcoll = h1->GetBinCenter(i);
		double Width = h1->GetXaxis()->GetBinWidth(i);
		double g = 0;
		for (int deg = 0; deg <= degree; deg++){
			g += par[deg] * TMath::Power(Mcoll, deg) * Width;
		}
		//double b = (n1+g+n2)/2;
		//f += -2*(TMath::Log(TMath::Poisson(n1+g,b))+TMath::Log(TMath::Poisson(n2,b)));
		f += -2*(TMath::Log(TMath::Poisson(n2,n1+g)));
	}
}

//void Newf1(int &npar, double *gin, double &f, double *par, int iflag) {
//	//par[i] are the polynomial coefficients
//
//	//get histograms and polynomial degree
//	HistoObject *data;
//	data = (HistoObject*)minuit2->GetObjectFit();
//	TH1D *h1 = data->m_h1;
//	TH1D *h2 = data->m_h2;
//	int degree = data->m_degree;
//
//	f = 0;
//	int nbins = h1->GetXaxis()->GetNbins();
//
//	//add poly to h1 and calculate poisson for side bands
//	for (int i=1; i<=7; i++){
//		double n1 = h1->GetBinContent(i);
//		double n2 = h2->GetBinContent(i);
//		double Mcoll = h1->GetBinCenter(i);
//		double width = h1->GetXaxis()->GetBinWidth(i);
//		double g = 0;
//		for (int deg = 0; deg <= degree; deg++){
//			g += par[deg] * TMath::Power(Mcoll, deg)*width;
//		}
//		double b = (n1+g+n2)/2;
//		f += -2*(TMath::Log(TMath::Poisson(n1+g,b))+TMath::Log(TMath::Poisson(n2,b)));
//	}
//	for (int i=13; i<=nbins; i++){
//		double n1 = h1->GetBinContent(i);
//		double n2 = h2->GetBinContent(i);
//		double Mcoll = h1->GetBinCenter(i);
//		double width = h1->GetXaxis()->GetBinWidth(i);
//		double g = 0;
//		for (int deg = 0; deg <= degree; deg++){
//			g += par[deg] * TMath::Power(Mcoll, deg)*width;
//		}
//		double b = (n1+g+n2)/2;
//		f += -2*(TMath::Log(TMath::Poisson(n1+g,b))+TMath::Log(TMath::Poisson(n2,b)));
//	}
//
//}

void PrintPolyCoefficients(TH1D* h1, TH1D* h2, int degree)
{
	InitExterns();

	HistoObject d;
	d.seth1(h1);
	d.seth2(h2);
	d.setDegree(degree);

	minuit2 = new TMinuit(degree+1);
	minuit2->SetObjectFit(&d);
	minuit2->SetFCN(Newf1);

	h1->Draw();h2->Draw("sames");

	// initialize the parameters:
	TString varname;
	double startValue = 0.1;
	double stepSize = 0.0001;
	double min = -1;
	double max = 1;
	Int_t ierflg = 0;

	//wider range for the free coefficient
	minuit2->mnparm(0,"a_0",startValue,stepSize,min,max,ierflg);

	cout<<"ierflg = "<<ierflg<<endl;

	//loop over the other coefficients
	for (int deg = 1; deg <= degree; deg++) {
		varname.Form("a_%d", deg);
		minuit2->mnparm(deg, varname, startValue, stepSize, min, max, ierflg);
	}
	minuit2->SetMaxIterations(5000);
	minuit2->Migrad();

	//print
	double value[degree+1];
	double err;
	double g = 0;
	TString printThis;
	double a[degree+1];
 	for (int i = 0; i <= degree; i++) {
//		a[i]=minuit2->GetParameter(i, value[i], err);
 		minuit2->GetParameter(i, value[i], err);
		printThis.Form("a_%d   %f", i, value[i]);
    	cout << printThis << endl;
    }
}

//
//void PrintPolyCoefficients(TH1D* h1, TH1D* h2, int degree)
//{
//	InitExterns();
//
//	HistoObject d;
//	d.seth1(h1);
//	d.seth2(h2);
//	d.setDegree(degree);
//
//	minuit3 = new TMinuit(degree + 1);
//	minuit3->SetObjectFit(&d);
//
//	minuit3->SetFCN(Newf1);
//
//	// initialize the parameters:
//	TString varname;
//	double startValue = 0;
//	double stepSize = 0.0001;
//	double min = -1;
//	double max = 1;
//	Int_t ierflg = 0;
//
//	//wider range for the free coefficient
//	minuit3->mnparm(0,"a_0",startValue,stepSize,10*min,10*max,ierflg);
//
//	cout<<"ierflg = "<<ierflg<<endl;
//
//	//loop over the other coefficients
//	for (int deg = 1; deg <= degree; deg++) {
//		varname.Form("a_%d", deg);
//		minuit3->mnparm(deg, varname, startValue, stepSize, min, max, ierflg);
//	}
//	minuit3->SetMaxIterations(500);
//	minuit3->Migrad();
//
//	//print
//	double value, err;
//	TString printThis;
//	for (int i = 0; i <= degree; i++) {
//		minuit3->GetParameter(i, value, err);
//		printThis.Form("a_%d   %f", i, value);
//		cout << printThis << endl;
//	}
//}

TH1D* AddPolonHisto(TH1D* h, double value1, double value2)
{
	const TArrayD* binsarray = h->GetXaxis()->GetXbins();
	const Double_t * massBins = binsarray->GetArray();	int nbins = h->GetXaxis()->GetNbins();
	TH1D* hplus = new TH1D("MEplus","MEplus",nbins,massBins);
	double n[nbins+2];
		for (int i=1; i<=nbins; i++){
			n[i] = h->GetBinContent(i);
			double Mcoll = h->GetBinCenter(i);
			double Width = h->GetXaxis()->GetBinWidth(i);
			n[i] += value1 * Width + value2 * Mcoll * Width;
		}
	n[0] = 0; n[nbins+1] = 0;
	hplus->SetContent(n);
	return hplus;
}

void meanDataEstimatorFunc(int &npar, double *gin, double &f, double *par, int iflag)
{
	//par[0] is muHat, par[1-bins] is B, par[bins+1 - bins+degree+1] is polynomial coefficients
	NewDataObject *data;
	data = (NewDataObject*)minuit->GetObjectFit();

	double *n1 = data->m_n1;
	double *n2 = data->m_n2;
	double *S = data->m_S;
	double *bins = data->m_bins;
	int degree = data->m_degree;
	f=0;

	int size = data->m_size;
	int nbins = size-2;

	for (int i=1; i<size-1; i++){
		if (n1[i]+n2[i]!=0){
			double Mcoll = 0.5*(bins[i]+bins[i-1]);
			double width = bins[i]-bins[i-1];
			double g = 0;
			for (int deg = 0; deg <= degree; deg++){
				g += par[nbins+1+deg] * TMath::Power(Mcoll, deg)*width;
			}
			f -= 2*(TMath::Log(TMath::Poisson(n1[i],par[i]+g))
			+TMath::Log(TMath::Poisson(n2[i],par[i]+par[0]*S[i])));
		}
	}
}

double LikelihoodTestValue(TH1D* h_1, TH1D* h_2, double* muHatB, TH1D* signalGaus, int degree)
{
	int nbins = h_1->GetXaxis()->GetNbins();
	double maxBin = h_1->GetXaxis()->GetXmax();

	const TArrayD* binsarray = h_1->GetXaxis()->GetXbins();
	const Double_t * massBins = binsarray->GetArray();
	double n1[nbins+2];
	double n2[nbins+2];
	double S[nbins+2];
	double min[nbins+2];
	double max[nbins+2];
	double startValue[nbins+2];

	double errorVector[nbins+1];

	for(int i=0; i<=nbins+1; i++){
		n1[i] = h_1->GetBinContent(i);
		n2[i] = h_2->GetBinContent(i);
		S[i] = signalGaus->GetBinContent(i);
		double tempmin = (n1[i] > n2[i] ? n2[i] : n1[i]);
		double tempmax = (n1[i] > n2[i] ? n1[i] : n2[i]);
		if (tempmin == tempmax){ tempmax = tempmax+1;}
		min[i] = tempmin;//+ 0.25*(tempmax-tempmin);;
		max[i] = tempmax;//tempmin+ 0.75*(tempmax-tempmin);
		startValue[i] = min[i];
	}

	NewDataObject d;
	d.setN1(n1,nbins+2);
	d.setN2(n2,nbins+2);
	d.setS(S,nbins+2);
	d.setSize(nbins+2);
	d.setBins(massBins,nbins+1);
	d.setDegree(degree);

	minuit = new TMinuit(nbins+1+degree+1);
	minuit->SetObjectFit(&d);
	minuit->SetFCN(meanDataEstimatorFunc);

	Int_t ierflg = 0;
	minuit->SetPrintLevel(-1);

	// initialize the parameters:

	double stepSize = 0.01;

	TString varname;
	minuit->mnparm(0,"muHat",0,0.0001,0,10,ierflg);
	for (int j=1;j<=nbins;j++){
		varname.Form("B%d",j);
		minuit->mnparm(j,varname,startValue[j],stepSize,min[j],max[j],ierflg);
	}
	for (int d = 0; d<=degree; d++){
		varname.Form("a%d",d);
		minuit->mnparm(nbins+1+d,varname,0,0.0001,-1,1,ierflg);
	}

	minuit->SetMaxIterations(500);
	minuit->Migrad();

	for (int j=0; j<=nbins+degree+1; j++){
		minuit->GetParameter(j,muHatB[j],errorVector[j]);
	}

	return minuit->fAmin;

}



TH1D* meanDataEstimator(TH1D* h_1, TH1D* h_2, double* muHatB, TH1D* signalGaus, int degree, TString name)
{
	int nbins = h_1->GetXaxis()->GetNbins();
	const TArrayD* binsarray = h_1->GetXaxis()->GetXbins();
	const Double_t * massBins = binsarray->GetArray();

	//get two options: signal in ME or signal in EM and choose best fit
	double testmuHatB1[nbins+2+degree];
	double testmuHatB2[nbins+2+degree];
	double v1 = LikelihoodTestValue(h_1,h_2,testmuHatB1,signalGaus,degree);
	double v2 = LikelihoodTestValue(h_2,h_1,testmuHatB2,signalGaus,degree);
//	cout<<"muHat1 = "<<testmuHatB1[0]<<", muHat2 = "<<testmuHatB2[0]<<endl;
	if (v1<v2){memcpy(muHatB,testmuHatB1,(nbins+2+degree)*sizeof(double));}//cout<<"chosen is 1"<<endl;}
	else {memcpy(muHatB,testmuHatB2,(nbins+2+degree)*sizeof(double));}//cout<<"chosen is 2"<<endl;}

	TH1D* h_bkg = new TH1D(name,name,nbins,massBins);

	double Bcontent[nbins+2];
	memcpy(Bcontent,muHatB,(nbins+1)* sizeof(double));
	Bcontent[0] = 0; Bcontent[nbins+1] = 0;
	h_bkg->SetContent(Bcontent);

	for (int d = 0; d<=degree; d++){
		cout<<"coefficient"<<d<<" = "<<muHatB[nbins+1+d]<<endl;
	}

	ToyData::setBErrors(h_1,h_2,h_bkg);

	return h_bkg;

}



double LikelihoodRatio_mu(TH1D* h1, TH1D* h2, TH1D* h_b, double mu)
{
	double lambda = 0;
	//find mu Hat
	double mu_Hat = muHat(h1,h2,h_b);
//	cout<<"muhat = "<< mu_Hat<<endl;
	//create B+muHatS
	TH1D* h_bmuHatS = ToyData::appendSmoothTestGaussian(h_b,mu_Hat,"h_bmuHatS");
	//create B+muS
	TH1D* h_bmuS = ToyData::appendSmoothTestGaussian(h_b,mu,"h_bmuS");

	for(int k=COLL_MASS_MinSigBin; k<=COLL_MASS_MaxSigBin; k++)
	{
		double n1 = h1->GetBinContent(k);
		double n2 = h2->GetBinContent(k);
		double b = h_b->GetBinContent(k);
		double bmuS = h_bmuS->GetBinContent(k);
		double bmuHS = h_bmuHatS->GetBinContent(k);

		double l_mu = -2.*(TMath::Log(TMath::Poisson(n1,b))+TMath::Log(TMath::Poisson(n2,bmuS)));
		double l_muHat = -2.*(TMath::Log(TMath::Poisson(n1,b))+TMath::Log(TMath::Poisson(n2,bmuHS)));

		if(TMath::Poisson(n1,b)==0){
			cout<<"LikelihoodRatio_mu: in bin:"<<k<<", b = "<<b<<", n1 = "<<n1<<endl;
		}
		if(TMath::Poisson(n2,bmuS)==0){
			cout<<"LikelihoodRatio_mu: in bin:"<<k<<", bmuS = "<<bmuS<<", n2 = "<<n2<<endl;
		}
		lambda += l_mu-l_muHat;
	}
//	cout<<"LikelihoodRatio_mu: Lambda = "<<lambda<<endl;
	return lambda;
}

double LikelihoodRatio_mu2(TH1D* h1, TH1D* h2, TH1D* h_b, TH1D* h_signal,double mu)
{
	double lambda = 0;
	//find mu Hat
	double mu_Hat = muHat2(h1,h2,h_b,h_signal);
//	cout<<"muhat = "<< mu_Hat<<endl;
	//create B+muHatS
	TH1D* h_bmuHatS = ToyData::appendSignal(h_b,h_signal,mu_Hat,"h_bmuHatS");
	//create B+muS
	TH1D* h_bmuS = ToyData::appendSignal(h_b,h_signal,mu,"h_bmuS");

	for(int k=COLL_MASS_MinSigBin; k<=COLL_MASS_MaxSigBin; k++)
	{
		double n1 = h1->GetBinContent(k);
		double n2 = h2->GetBinContent(k);
		double b = h_b->GetBinContent(k);
		double bmuS = h_bmuS->GetBinContent(k);
		double bmuHS = h_bmuHatS->GetBinContent(k);

		if (n1+n2 !=0){
			double l_mu = -2.*(TMath::Log(TMath::Poisson(n1,b))+TMath::Log(TMath::Poisson(n2,bmuS)));
			double l_muHat = -2.*(TMath::Log(TMath::Poisson(n1,b))+TMath::Log(TMath::Poisson(n2,bmuHS)));

			if(TMath::Poisson(n1,b)==0){
				cout<<"LikelihoodRatio_mu2: in bin:"<<k<<", b = "<<b<<", n1 = "<<n1<<endl;
			}
			if(TMath::Poisson(n2,bmuS)==0){
				cout<<"LikelihoodRatio_mu: in bin:"<<k<<", bmuS = "<<bmuS<<", n2 = "<<n2<<endl;
			}
			lambda += l_mu-l_muHat;
		}
	}
//	cout<<"LikelihoodRatio_mu: Lambda = "<<lambda<<endl;
	return lambda;
}


double Likelihood(TH1D* h1, TH1D* h2, TH1D* hb, TH1D* hbS)
{
	double lambda = 0;

	for(int k=COLL_MASS_MinSigBin; k<=COLL_MASS_MaxSigBin; k++)
	{
		double n1 = h1->GetBinContent(k);
		double n2 = h2->GetBinContent(k);
		double b = hb->GetBinContent(k);
		double bS = hbS->GetBinContent(k);

		lambda += -2.*(TMath::Log(TMath::Poisson(n1,b))+TMath::Log(TMath::Poisson(n2,bS)));
	}
	return lambda;
}


TH1D* Pdf_ratio_qMu_zero(TH1D* h_bg, double mu)
{
	delete gDirectory->FindObject("h1_likelihoodRatioPdf");
	TH1D* h_testStatisticPdf = new TH1D("h1_likelihoodRatioPdf","h1_likelihoodRatioPdf",
			LIKELIHOODRATIOMU_PDF_NBINS,LIKELIHOODRATIOMU_PDF_MIN,LIKELIHOODRATIOMU_PDF_MAX);

	int NUM_ITERATIONS  = 10000;

	for(int i=0; i<NUM_ITERATIONS; i++)
	{
		//generate toy histos
		delete gDirectory->FindObject("htoy1");
		delete gDirectory->FindObject("htoy2");
		delete gDirectory->FindObject("hbg_rand");

		TH1D *hbg_rand = drawFromHistoGaus(h_bg,"hbg_rand");
		TH1D *htoy1 = ToyData::drawFromHisto(hbg_rand,"htoy1");
		TH1D *htoy2 = ToyData::drawFromHisto(hbg_rand,"htoy2");

		//caluculate test statistic: -2Log[L(mu)/L(muHat)]
		double t = LikelihoodRatio_mu(htoy1,htoy2,h_bg,mu);
//		cout<<"L_mu/L_muHat = "<<t<<endl;

		h_testStatisticPdf->Fill(t,1./NUM_ITERATIONS);
	}
	if ((1-h_testStatisticPdf->Integral(1,LIKELIHOODRATIO_PDF_NBINS))>0.0001){
		cout<<"normalization problem? "<< h_testStatisticPdf->Integral(1,LIKELIHOODRATIO_PDF_NBINS)<<endl;
	}
	return h_testStatisticPdf;
}

TH1D* Pdf_ratio_qMu_signal(TH1D* h_bg, TH1D* h_bS, double mu)
{
	delete gDirectory->FindObject("h1_Pdf");
	TH1D* h_testStatisticPdf2 = new TH1D("h1_Pdf","h1_Pdf",
			LIKELIHOODRATIOMU_PDF_NBINS,LIKELIHOODRATIOMU_PDF_MIN,LIKELIHOODRATIOMU_PDF_MAX);

	int NUM_ITERATIONS  = 100000;

	for(int i=0; i<NUM_ITERATIONS; i++)
	{
		//generate toy histos
		delete gDirectory->FindObject("htoy1");
		delete gDirectory->FindObject("htoy2");
		delete gDirectory->FindObject("hbg_rand");
		delete gDirectory->FindObject("hsig_rand");

		TH1D* hsig_rand = drawFromHistoGaus(h_bS,"bsig_rand");
		TH1D *hbg_rand = drawFromHistoGaus(h_bg,"hbg_rand");
		TH1D *htoy1 = ToyData::drawFromHisto(hbg_rand,"htoy1");
		TH1D *htoy2 = ToyData::drawFromHisto(hsig_rand,"htoy2");

		//caluculate test statistic: -2Log[L(mu)/L(muHat)]
		double t = LikelihoodRatio_mu(htoy1,htoy2,h_bg,mu);
//		cout<<"L_mu/L_muHat = "<<t<<endl;

		h_testStatisticPdf2->Fill(t,1./NUM_ITERATIONS);
	}
	if ((1-h_testStatisticPdf2->Integral(1,LIKELIHOODRATIO_PDF_NBINS))>0.0001){
		cout<<"normalization problem? "<< h_testStatisticPdf2->Integral(1,LIKELIHOODRATIO_PDF_NBINS)<<endl;
	}
	return h_testStatisticPdf2;
}


double muHat(TH1D* h_1,TH1D* h_2,TH1D* h_b)
{
	Minimizer minim;

	LikelihoodClass LL;
	LL.seth1(h_1);
	LL.seth2(h_2);
	LL.setB(h_b);

	double min = 0;
	double max = 20;

	double mu_hat = minim.minimize(LL,min,max,0.00001);

//	cout<<"muHat:muHat = "<<mu_hat<<endl;

	return mu_hat;

}


double muHat2(TH1D* h_1,TH1D* h_2,TH1D* h_b, TH1D* h_sig)
{
	Minimizer minim;

	LikelihoodClass2 LL;
	LL.seth1(h_1);
	LL.seth2(h_2);
	LL.setB(h_b);
	LL.setS(h_sig);

	double min = 0;
	double max = 20;

	double mu_hat = minim.minimize(LL,min,max,0.00001);

//	cout<<"muHat2:muHat = "<<mu_hat<<endl;

	return mu_hat;

}


double pValue(TH1D* h_pdf, double value)
{
	int bin_i = h_pdf->GetXaxis()->FindBin(value);

	return h_pdf->Integral(bin_i,LIKELIHOOD_PDF_NBINS);
}


double muCL_discovery(TH1D* h1,TH1D* h2,TH1D* hb)
{
	Minimizer minim;

	VanillaClass_discovery V;
	V.seth1(h1);
	V.seth2(h2);
	V.sethb(hb);

	double min = 0;
	double max = 3;

	double mu = minim.minimize(V,min,max,0.00001);

//	cout<<"muCL_discovery: mu = "<<mu<<", abs(L_mu - 9) = "<<minim.minF<<endl;

	return mu;
}

double muCL_discovery2(TH1D* h1,TH1D* h2,TH1D* hb,TH1D* hs)
{
	Minimizer minim;

	VanillaClass_discovery2 V;
	V.seth1(h1);
	V.seth2(h2);
	V.sethb(hb);
	V.seths(hs);

	double min = 0;
	double max = 5;

	double mu = minim.minimize(V,min,max,0.00001);

//	cout<<"muCL_discovery: mu = "<<mu<<", abs(L_mu - 9) = "<<minim.minF<<endl;

	return mu;
}



double muCL_limit(TH1D* h,TH1D* hb)
{
	Minimizer minim;

	VanillaClass_limit V;
	V.sethb(hb);
	V.seth(h);

	double min = 0.1;
	double max = 0.4;

	double mu = minim.minimize(V,min,max,0.0000001);

//	cout<<"muCL_limit: mu = "<<mu<<", abs(L_mu - 3.84) = "<<minim.minF<<endl;

	return mu;
}


double VanillaClass_discovery::function(double a)
{
	double l = LikelihoodRatio_mu(m_h1,m_h2,m_hb,a);
	return abs(l-9);
}

double VanillaClass_discovery2::function(double a)
{
	double l = LikelihoodRatio_mu2(m_h1,m_h2,m_hb,m_hs,a);
	return abs(l-9);
}

double VanillaClass_limit::function(double a)
{
	TH1D* h_sig = ToyData::appendSmoothTestGaussian(m_hb,a,"h_sig");
	TH1D* h_sigRand = ToyData::drawFromHistoGaus(h_sig,"h_sigRand");
	double t = LikelihoodRatio_mu(m_h,h_sigRand,m_hb,0);
	return abs(t-pow(1.96,2));
}


double LikelihoodClass::function(double a)
{
	int nbins = m_B->GetXaxis()->GetNbins();
	double min = m_B->GetXaxis()->GetXmin();
	double max = m_B->GetXaxis()->GetXmax();

	delete gDirectory->FindObject("Bnew");

	TH1D *B_new = new TH1D("Bnew","Bnew",nbins,min,max);

	//Add smooth gaussian to the bg estimator with strength mu
	B_new = ToyData::appendSmoothTestGaussian(m_B,a,"Bnew");

	double l = Likelihood(m_h1,m_h2,m_B,B_new);

	return l;

}

double LikelihoodClass2::function(double a)
{
	int nbins = m_B->GetXaxis()->GetNbins();
	double min = m_B->GetXaxis()->GetXmin();
	double max = m_B->GetXaxis()->GetXmax();

	delete gDirectory->FindObject("Bnew");

//	TH1D *B_new = new TH1D("Bnew","Bnew",nbins,min,max);

	//Add smooth gaussian to the bg estimator with strength mu
	TH1D* B_new = ToyData::appendSignal(m_B,m_S,a,"Bnew");

	double l = Likelihood(m_h1,m_h2,m_B,B_new);

	return l;

}



double Median(const TH1D * h1)
{
   int n = h1->GetXaxis()->GetNbins();
   std::vector<double>  x(n);
   h1->GetXaxis()->GetCenter( &x[0] );
   const double * y = h1->GetArray();
   // exclude underflow/overflows from bin content array y
   return TMath::Median(n, &x[0], &y[1]);
}

}
