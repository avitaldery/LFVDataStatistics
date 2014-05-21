#include "LFVDataStatistics/Constants.h"

#include "TRandom.h"
#include "TH1.h"

double COLL_MASS_MIN;
double COLL_MASS_MAX;
int COLL_MASS_NBINS;
int COLL_MASS_MinSigBin;
int COLL_MASS_MaxSigBin;
int TOY_MC_NUMOFEVENTS;
int LIKELIHOOD_PDF_NBINS;
double LIKELIHOOD_PDF_MIN;
double LIKELIHOOD_PDF_MAX;
int LIKELIHOODRATIO_PDF_NBINS;
double LIKELIHOODRATIO_PDF_MIN;
double LIKELIHOODRATIO_PDF_MAX;
int LIKELIHOODRATIOMU_PDF_NBINS;
double LIKELIHOODRATIOMU_PDF_MIN;
double LIKELIHOODRATIOMU_PDF_MAX;
int LIKELIHOOD_PDF_SEED;
double GAUSS_WIDTH;
double GAUSS_MEAN;
TH1D* TEST_GAUSSIAN;
TH1D* DATA_GAUSSIAN;
TRandom RAND;
unsigned int SEED;
TH1D* H_DEFAULT;

void InitExterns()
{
	COLL_MASS_MIN = 0;
	COLL_MASS_MAX = 500;
	COLL_MASS_NBINS = 125;
	COLL_MASS_MinSigBin = 6;
	COLL_MASS_MaxSigBin = 40;

	TOY_MC_NUMOFEVENTS=1000;
	LIKELIHOOD_PDF_NBINS=500;
	LIKELIHOOD_PDF_MIN=-200;
	LIKELIHOOD_PDF_MAX=200;
	LIKELIHOODRATIO_PDF_NBINS=1000;
	LIKELIHOODRATIO_PDF_MIN=-50;
	LIKELIHOODRATIO_PDF_MAX=50;
	LIKELIHOODRATIOMU_PDF_NBINS=550;
	LIKELIHOODRATIOMU_PDF_MIN=-5;
	LIKELIHOODRATIOMU_PDF_MAX=50;
	LIKELIHOOD_PDF_SEED=12;

	GAUSS_WIDTH = 20;
	GAUSS_MEAN = 300;

	H_DEFAULT = new TH1D("default","default",COLL_MASS_NBINS,COLL_MASS_MIN,COLL_MASS_MAX);


	TEST_GAUSSIAN = new TH1D("test_gauss","test_gauss",COLL_MASS_NBINS,COLL_MASS_MIN,COLL_MASS_MAX);
	TRandom rand;
	rand.SetSeed(12);
	for (int i=0; i<100000; i++)
	{
		TEST_GAUSSIAN->Fill(rand.Gaus(126,8.6));
		// 126, 8.63 after Delphes 20 GeV
		// 125, 7.3 for regular cutflow without detector effects
		// 125, 6.9 for 20 GeV cut without detector effects
		// 125.8, 8.8 after Delphes 10 GeV
		// 125.8, 10.7 for Delphes plus MET smearing 20%
		// 128.4, 21.5 for Delphes plus MET smearing 80%
		// 125, 4.6 for Empty Delphes
	}
	TEST_GAUSSIAN->Scale(49.7/100000);
	// 98.2 for Delphes 20 Gev (smearing 20% also)
	// 177.8 for regular cutflow without detector effects
	// 101.3 for 20 GeV cut without detector effects
	// 159.7 for Delphes 10 GeV
	// 89.9 for Delphes plush MET smearing 80%
	// 49.7 with correct c factor Delphes 20 GeV
	// 78 for Empty Delphes

//	double massBins[] = {0, 20, 40, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
//		               225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 700};
//	DATA_GAUSSIAN = new TH1D("data_gauss","data_gauss",30,massBins);
//	for (int i=0; i<100000; i++)
//	{
//		DATA_GAUSSIAN->Fill(rand.Gaus(126,8.6));
//	}
//	DATA_GAUSSIAN->Scale(50./100000);


}
