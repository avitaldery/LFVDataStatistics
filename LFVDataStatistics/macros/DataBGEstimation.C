#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>

using namespace std;

#include "LFVDataStatistics/Functions_BkgEstimator.h"
#include "LFVDataStatistics/Constants.h"


void DataBGEstimation(TString filename, TString sigFile)
{
	InitExterns();

	int Polydegree = 2;

	TFile* f = new TFile(filename);

	//get unblinded histos
	TH1D* h_ME = (TH1D*)f->Get("nom/em_LFV_ME_McollHiggs_Unblind");
	TH1D* h_EM = (TH1D*)f->Get("nom/em_LFV_EM_McollHiggs_Unblind");

	TFile* fsig = new TFile(sigFile);
	TH1D* h_sig = (TH1D*)fsig->Get("nom/em_LFV_ME_McollHiggs_Unblind");

	const TArrayD* binsarray = h_ME->GetXaxis()->GetXbins();
	const Double_t * massBins = binsarray->GetArray();
	int nbins = h_ME->GetXaxis()->GetNbins();

	double muHatB[nbins+2+Polydegree];

	//trial: add polynomial to one histo
	double n[nbins+2];
	for (int i=1; i<=nbins; i++){
		n[i] = h_ME->GetBinContent(i);
		double Mcoll = h_ME->GetBinCenter(i);
		double width = h_ME->GetXaxis()->GetBinWidth(i);
		n[i] += (0.5 + 0.0002*Mcoll)*width;
	}

	n[0] = 0; n[nbins+1] = 0;

	TH1D *h_MEplus = new TH1D("hMEplus","hMEplus",nbins, massBins);
	h_MEplus->SetContent(n);

	TH1D* h_MESignal = ToyData::appendSignal(h_ME,h_sig,5,"hMEsignal");

	//get BG estimation
	TH1D* h_B =  BkgEstimator::meanDataEstimator(h_MESignal,h_EM,muHatB,h_sig,Polydegree,"B");

	//create blind histos to draw
	TH1D* h_ME_blind = ToyData::getBlindHisto(h_ME,8,12,"ME blind");
	TH1D* h_EM_blind = ToyData::getBlindHisto(h_EM,8,12,"EM blind");

	TLegend* leg = new TLegend(0.5,0.7,0.7,0.9);
	leg->SetFillColor(kWhite); leg->SetBorderSize(1); leg->SetLineColor(0); leg->SetTextFont(42);

	h_B->GetXaxis()->SetTitle("M_{Collinear} (GeV)");
	h_B->SetLineColor(kBlack); h_B->SetLineWidth(2); h_B->Draw();
	h_ME_blind->SetLineColor(kRed); h_ME_blind->SetLineWidth(2); //h_ME_blind->Draw("e1 sames");
	h_EM_blind->SetLineColor(kBlue); h_EM_blind->SetLineWidth(2);// h_EM_blind->Draw("e1 sames");
	h_EM->SetLineColor(kBlue); h_EM->SetLineWidth(2); h_EM->Draw("e1 sames");
	h_MESignal->SetLineColor(kRed);h_MESignal->SetLineWidth(2); h_MESignal->Draw("sames");
//	h_MEplus->Draw("sames");
	leg->AddEntry(h_ME_blind,"#mue","l");
	leg->AddEntry(h_EM_blind,"e#mu","l");
	leg->AddEntry(h_B,"BG estimation","l");
	leg->Draw();

	//get best polynomial coefficients
	BkgEstimator::PrintPolyCoefficients(h_ME,h_MEplus,2);

}
