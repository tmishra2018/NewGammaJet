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

  TFile* dataFile=TFile::Open("pu_truth_data2015_50bins.root"); 
  TFile* mcFile = TFile::Open("computed_mc_GJET_plus_QCD_pu_truth_50bins.root");

  TH1* dataHisto = static_cast<TH1*>(dataFile->Get("pileup"));
  TH1* mcHisto = static_cast<TH1*>(mcFile->Get("pileup"));

  TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/analysis/tuples/GJET_MC/PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns_ReReco.root");                                       
  TTree* AnalysisTree_mc = (TTree*) f1.Get("gammaJet/analysis");                                                        
  uint64_t totalEvents_mc = AnalysisTree_mc->GetEntries();  

  
    ///////////////////////////////////////////////////////////////////

    // Normalize
    dataHisto->Scale(1.0 / dataHisto->Integral());
    mcHisto->Scale(1.0 / mcHisto->Integral());
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    dataHisto-> SetLineColor(kBlue);
    mcHisto -> SetLineColor(kRed);
    dataHisto -> Draw();
    mcHisto -> Draw("same");
    c1-> SaveAs("Plot/pileup_truth.png");
    //   c1-> SaveAs("pileup_truth.root");
        
    // MC * data / MC = data, so the weights are data/MC:
    TH1* ratioHisto = static_cast<TH1*>(dataHisto->Clone());
    ratioHisto->Divide(mcHisto);

    int NBins = ratioHisto->GetNbinsX();

    for (int ibin = 1; ibin < NBins + 1; ++ibin) {
      std::cout << "   " << ibin - 1 << " " << ratioHisto->GetBinContent(ibin) << std::endl;
    }
        
    TCanvas *c2 = new TCanvas("c2","c2",800,800);
    ratioHisto -> Draw();
    c2-> SaveAs("Plot/ratioHisto.png");
    //    c2-> SaveAs("ratioHisto.root");
    
    TH1D* mcHisto_reweighted = new TH1D("mcHisto_reweighted", "mc Reweighted", 50, 0, 50);
    
    for (uint64_t i = 0; i < totalEvents_mc; i++) {      
      if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
      AnalysisTree_mc->GetEntry(i);
      Float_t interaction;
      AnalysisTree_mc->SetBranchAddress("ntrue_interactions", &interaction);
      //cout<< ptPhot << endl;
      
	int bin = ratioHisto ->GetXaxis()-> FindBin(interaction);
	float PUweight = ratioHisto -> GetBinContent(bin);
	
	double generator_weight;
	AnalysisTree_mc->SetBranchAddress("generator_weight", &generator_weight);   
	if(generator_weight == 0) generator_weight = 1;
	float event_weight;
	AnalysisTree_mc->SetBranchAddress("evtWeightTot", &event_weight);   
	
	double Weight = PUweight* generator_weight * event_weight;
	
	//      cout<< "nvertex "<< nvertex_mc<<endl;
	//      cout<< "PUWeight "<< PUweight<<endl;	
	mcHisto_reweighted -> Fill(interaction, Weight);     
    }
    
    mcHisto_reweighted ->Scale(1.0 / mcHisto_reweighted->Integral());
    
    TCanvas *c3 = new TCanvas("c3","c3",800,800);
    dataHisto -> SetLineColor(kBlue);
    mcHisto_reweighted -> SetLineColor(kRed);
    dataHisto -> Draw();
    mcHisto_reweighted -> Draw("same");
    c3-> SaveAs("Plot/Histo_reweighted.png");
    //    c3-> SaveAs("Histo_reweighted.root");


        
}//main





