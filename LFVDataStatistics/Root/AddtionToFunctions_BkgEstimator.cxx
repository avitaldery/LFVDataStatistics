#include "TMinuit.h"

TMinuit* minuit = 0;
void f1(int &npar, double *gin, double &f, double *par, int iflag)
{
	DataObject* data = (DataObject*) minuit->GetObjectFit();
	f = 1 - TMath::Poisson(data->m_n1,*par)*TMath::Poisson(data->m_n2,*par);
}

void printPoissonMin(double n1,double n2)
{
	// put your data points in an object so you can access them in the function
	DataObject d;
	d.setN1(n1);
	d.setN2(n2);

	int number_of_params = 1;
	// initialize a TMinuit object
	minuit = new TMinuit(number_of_params);
	minuit->SetObjectFit(&d);
	// set the function you want to minimize
	minuit->SetFCN(f1);
	// initialize the parameters
	double startValue = n1;
	double stepSize = 0.0001;
	double min = n1;
	double max = n2;
	Int_t ierflg = 0;
	minuit->mnparm(0,"b",startValue,stepSize,min,max,ierflg);
	// next one would look like:
	// minuit->mnparm(1,"nameOfSecond",startValue[1],stepSize[1],min[1],max[1],ierflg);

	// use the MIGRAD algorithm to minimize
	minuit->SetMaxIterations(500);
	minuit->Migrad();
}

