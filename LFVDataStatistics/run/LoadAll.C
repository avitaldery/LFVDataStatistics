#include "TROOT.h"
#include "TSystem.h"

void LoadAll()
{

    gROOT->ProcessLine("#include <vector>");

    // Set the include path for all packages
    gSystem->AddIncludePath("-I/home/avitald/LFV/LFVDataStatistics");
    cout << gSystem->GetIncludePath() << endl;

    // Load package libraries
    gROOT->ProcessLine(".L ../cmt/LFVDataStatistics.C+");

    //Load additional macros
    gROOT->ProcessLine(".L ../macros/DataBGEstimation.C+");

    gROOT->ProcessLine(".L ../macros/Sensitivity_discovery.C+");


}
