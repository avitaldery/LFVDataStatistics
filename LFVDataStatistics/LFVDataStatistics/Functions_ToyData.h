
#include "TH1D.h"
#include "TString.h"

namespace ToyData{

TH1D* appendSmoothTestGaussian(TH1D* h, double signalStrength, TString name);
TH1D* appendDataGaussian(TH1D* h, double signalStrength, TString name);
TH1D* appendSignal(TH1D* h, TH1D* hsig, double signalStrength, TString name);
TH1D* drawFromHisto(TH1D* h_source,TString name);
TH1D* drawFromHistoGaus(TH1D* h_source, TString name);
void setBErrors(TH1D* h_1, TH1D* h_2,TH1D* h_b);
TH1D* getBlindHisto(TH1D* h_source,int bin1, int bin2,TString name);

}
