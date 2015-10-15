#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include "TPaveText.h"
#include "../../bin/triggers.h" 

#define LIGHT_RED TColor::GetColor(0xcf, 0xa0, 0xa1)

using namespace std;

int main(int argc, char* argv[]) {

  TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_12_patch4/src/JetMETCorrections/GammaJetFilter/analysis/tuples/GJET_MC/PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25nsSample_25nsV5.root"); 
  TTree* AnalysisTree_mc = (TTree*) f1.Get("gammaJet/analysis");
  // TTree* PhotonTree_mc = (TTree*) f1.Get("gammaJet/photon");
  TTree* METTree_mc = (TTree*) f1.Get("gammaJet/PFlowAK4chs/met");
  uint64_t totalEvents_mc = AnalysisTree_mc->GetEntries();

  cout<< totalEvents_mc << endl;
  //loop

  for (uint64_t i = 0; i < totalEvents_mc; i++) {      
    if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
       
    // PhotonTree_mc->GetEntry(i);
    AnalysisTree_mc->GetEntry(i);
    METTree_mc->GetEntry(i);
    //    float ptPhot;
    //    PhotonTree_mc->SetBranchAddress("pt",&ptPhot);
    UInt_t run;
    UInt_t lumi_block;
    AnalysisTree_mc->SetBranchAddress("run", &run); 
    AnalysisTree_mc->SetBranchAddress("lumi_block", &lumi_block); 

    Float_t et; 
  
    METTree_mc->SetBranchAddress("et", &et);     

    if(et > 10000) cout<< "run " <<run << "   lumi   "<< lumi_block<< "    met.et()  "<< et<< endl;
    
  }
  
}//main





