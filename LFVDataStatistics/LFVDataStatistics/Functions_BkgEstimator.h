#include "TH1.h"
#include "LFVDataStatistics/minim.h"
#include "TMath.h"
#include "LFVDataStatistics/Functions_ToyData.h"

namespace BkgEstimator{

TH1D* H_DEFAULT;

TH1D* meanDataEstimator(TH1D* h_1, TH1D* h_2, double* muHatB, TH1D* signalGauss,int degree, TString name);

void meanDataEstimatorFunc(int &npar, double *gin, double &f, double *par, int iflag);

double LikelihoodTestValue(TH1D* h_1, TH1D* h_2, double* muHatB, TH1D* signalGaus,int degree);


class NewDataObject : public TObject
{
public:

void setN1(double n1[],int size){m_n1 = new double[size]; for (int i=0; i<size; i++){m_n1[i]=n1[i];}}
void setN2(double n2[],int size){m_n2 = new double[size]; for (int i=0; i<size; i++){m_n2[i]=n2[i];}}
void setS(double S[],int size){m_S = new double[size]; for (int i=0; i<size; i++){m_S[i]=S[i];}}
void setSize(int size){m_size = size;}
void setBins(const Double_t* bins,int size){m_bins = new Double_t[size]; for (int i=0; i<size; i++){m_bins[i] = bins[i];}}
void setDegree(int deg){m_degree = deg;}

private:
int m_size;
double *m_n1;
double *m_n2;
double *m_S;
const Double_t *m_bins;
int m_degree;
};


class HistoObject : public TObject
{
public:

	void seth1(TH1D* h1){m_h1 = h1;}
	void seth2(TH1D* h2){m_h2 = h2;}
	void setDegree(int d){m_degree = d;}

private:
	TH1D* m_h1;
	TH1D* m_h2;
	int m_degree;
};

class LikelihoodClass : public IObjective
{
	public:
	void seth1(TH1D* h1){m_h1 = h1;}
	void seth2(TH1D* h2){m_h2 = h2;}
	void setB(TH1D* B){m_B = B;}

	virtual double function(double a);


	private:
	TH1D* m_h1;
	TH1D* m_h2;
	TH1D* m_B;
};

class LikelihoodClass2 : public IObjective
{
	public:
	void seth1(TH1D* h1){m_h1 = h1;}
	void seth2(TH1D* h2){m_h2 = h2;}
	void setB(TH1D* B){m_B = B;}
	void setS(TH1D* S){m_S = S;}

	virtual double function(double a);


	private:
	TH1D* m_h1;
	TH1D* m_h2;
	TH1D* m_B;
	TH1D* m_S;
};


class VanillaClass_discovery :public IObjective
{
	public:
	void seth1(TH1D* h1){m_h1 = h1;}
	void seth2(TH1D* h2){m_h2 = h2;}
	void sethb(TH1D* hb){m_hb = hb;}

	virtual double function(double a);

	private:
	TH1D* m_h1;
	TH1D* m_h2;
	TH1D* m_hb;
};


class VanillaClass_discovery2 :public IObjective
{
	public:
	void seth1(TH1D* h1){m_h1 = h1;}
	void seth2(TH1D* h2){m_h2 = h2;}
	void sethb(TH1D* hb){m_hb = hb;}
	void seths(TH1D* hs){m_hs = hs;}

	virtual double function(double a);

	private:
	TH1D* m_h1;
	TH1D* m_h2;
	TH1D* m_hb;
	TH1D* m_hs;
};



class VanillaClass_limit :public IObjective
{
	public:
	void sethb(TH1D* hb){m_hb = hb;}
	void seth(TH1D* h){m_h = h;}

	virtual double function(double a);

	private:
	TH1D* m_h;
	TH1D* m_hb;

};


void PrintPolyCoefficients(TH1D* h1, TH1D* h2, int degree);
double pValue(TH1D* h_pdf, double value);

double muHat(TH1D* h1,TH1D* h2);

double Likelihood(TH1D* h1, TH1D* h2, TH1D* hb, TH1D* hbS);
double muHat(TH1D* h_1,TH1D* h_2,TH1D* h_b);
double muHat2(TH1D* h_1,TH1D* h_2,TH1D* h_b, TH1D* h_signal);

TH1D* Pdf_ratio_qMu_zero(TH1D* h_bg, double mu);
TH1D* Pdf_ratio_qMu_signal(TH1D* h_bg, TH1D* h_bS, double mu);
double LikelihoodRatio_mu(TH1D* h1, TH1D* h2, TH1D* h_b, double mu);
double LikelihoodRatio_mu2(TH1D* h1, TH1D* h2, TH1D* h_b, TH1D* h_signal,double mu);
double Median(const TH1D * h1);

double muCL_discovery(TH1D* h1,TH1D* h2,TH1D* hb);
double muCL_discovery2(TH1D* h1,TH1D* h2,TH1D* hb,TH1D* hs);
double muCL_limit(TH1D* h,TH1D* hb);

void f1(double n_1,double n_2,double a,double& f);
void printPoissonMin(double n1,double n2);


}
