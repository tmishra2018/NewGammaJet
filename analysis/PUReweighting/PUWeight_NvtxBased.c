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

  TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/analysis/tuples/GJET_MC/PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_25ns_ReReco.root"); 
  TTree* AnalysisTree_mc = (TTree*) f1.Get("gammaJet/analysis");
  TTree* PhotonTree_mc = (TTree*) f1.Get("gammaJet/photon");
  uint64_t totalEvents_mc = AnalysisTree_mc->GetEntries();
  
  TFile f2("/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/analysis/tuples/Data/PhotonJet_2ndLevel_SinglePhoton__Run2015D-2015-10-27.root");
  TTree* AnalysisTree_data = (TTree*) f2.Get("gammaJet/analysis");
  TTree* PhotonTree_data = (TTree*) f2.Get("gammaJet/photon");
  uint64_t totalEvents_data = AnalysisTree_data->GetEntries();
  
    ///////////////////////////////////////////////////////////////////
  
  TH1D *h_nvertex_mc = new TH1D("nvertex_mc", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_40_60= new TH1D("nvertex_mc_ptPhot_40_60", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_60_85= new TH1D("nvertex_mc_ptPhot_60_85", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_85_100= new TH1D("nvertex_mc_ptPhot_85_100", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_100_130= new TH1D("nvertex_mc_ptPhot_100_130", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_130_175= new TH1D("nvertex_mc_ptPhot_130_175", "nvertex_mc", 51, -0.5, 50.5);
  TH1D *h_nvertex_mc_ptPhot_175_5000= new TH1D("nvertex_mc_ptPhot_175_5000", "nvertex_mc", 51, -0.5, 50.5);
  
  
  cout<< totalEvents_mc << endl;
  //loop
  for (uint64_t i = 0; i < totalEvents_mc; i++) {      
    //    for (uint64_t i = 0; i < 10000; i++) {      
       
    if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
    
    PhotonTree_mc->GetEntry(i);
    AnalysisTree_mc->GetEntry(i);
    float ptPhot;
    PhotonTree_mc->SetBranchAddress("pt",&ptPhot);
    UInt_t nvertex_mc;
    AnalysisTree_mc->SetBranchAddress("nvertex", &nvertex_mc);   
    
    double generator_weight;
    AnalysisTree_mc->SetBranchAddress("generator_weight", &generator_weight);   
    if(generator_weight == 0) generator_weight = 1;
    
    if(ptPhot < 40) continue;
    
    //      cout<< "generator weight  "<<generator_weight<<endl;
    
    //      float event_weight;
    //      AnalysisTree_mc->SetBranchAddress("evtWeightTot", &event_weight);   
    //      cout<< "event weight  "<<event_weight<<endl;
    
    double Weight = generator_weight;
    
    //      cout<< "Weight Tot  "<<Weight<<endl;
    
    h_nvertex_mc -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 40 && ptPhot < 60)        h_nvertex_mc_ptPhot_40_60 -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 60 && ptPhot < 85)        h_nvertex_mc_ptPhot_60_85 -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 85 && ptPhot < 100)      h_nvertex_mc_ptPhot_85_100 -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 100 && ptPhot < 130)    h_nvertex_mc_ptPhot_100_130 -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 130 && ptPhot < 175)    h_nvertex_mc_ptPhot_130_175 -> Fill(nvertex_mc, Weight);
    if(ptPhot >= 175 && ptPhot <= 5000)  h_nvertex_mc_ptPhot_175_5000 -> Fill(nvertex_mc, Weight);
  }
  
  h_nvertex_mc      ->Scale(1.0 / h_nvertex_mc->Integral());
  h_nvertex_mc_ptPhot_40_60      ->Scale(1.0 / h_nvertex_mc_ptPhot_40_60->Integral());
  h_nvertex_mc_ptPhot_60_85      ->Scale(1.0 / h_nvertex_mc_ptPhot_60_85->Integral());
  h_nvertex_mc_ptPhot_85_100    ->Scale(1.0 / h_nvertex_mc_ptPhot_85_100->Integral());
  h_nvertex_mc_ptPhot_100_130  ->Scale(1.0 / h_nvertex_mc_ptPhot_100_130->Integral());
  h_nvertex_mc_ptPhot_130_175  ->Scale(1.0 / h_nvertex_mc_ptPhot_130_175->Integral());
  h_nvertex_mc_ptPhot_175_5000->Scale(1.0 / h_nvertex_mc_ptPhot_175_5000->Integral());
  
  
  
  //////////////////////////////////////////////////////////     
  
  TH1D *h_nvertex_data = new TH1D("nvertex_data", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_40_60= new TH1D("nvertex_data_ptPhot_40_60", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_60_85= new TH1D("nvertex_data_ptPhot_60_85", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_85_100= new TH1D("nvertex_data_ptPhot_85_100", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_100_130= new TH1D("nvertex_data_ptPhot_100_130", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_130_175= new TH1D("nvertex_data_ptPhot_130_175", "nvertex_data", 51, -0.5, 50.5);
  TH1D *h_nvertex_data_ptPhot_175_5000= new TH1D("nvertex_data_ptPhot_175_5000", "nvertex_data", 51, -0.5, 50.5);
  
    cout<< totalEvents_data << endl;
    //loop
    for (uint64_t i = 0; i < totalEvents_data; i++) {      
      // for (uint64_t i = 0; i < 10000; i++) {      
      if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
      
      AnalysisTree_data->GetEntry(i);
      PhotonTree_data->GetEntry(i);
      float ptPhot;
      PhotonTree_data->SetBranchAddress("pt",&ptPhot);
      UInt_t nvertex_data;     
      AnalysisTree_data->SetBranchAddress("nvertex", &nvertex_data);
      //cout<< ptPhot << endl;
      
      if(ptPhot < 40) continue;
      
      h_nvertex_data -> Fill(nvertex_data);     
      if(ptPhot >= 40 && ptPhot < 60)        h_nvertex_data_ptPhot_40_60 -> Fill(nvertex_data);
      if(ptPhot >= 60 && ptPhot < 85)        h_nvertex_data_ptPhot_60_85 -> Fill(nvertex_data);
      if(ptPhot >= 85 && ptPhot < 100)      h_nvertex_data_ptPhot_85_100 -> Fill(nvertex_data);
      if(ptPhot >= 100 && ptPhot < 130)    h_nvertex_data_ptPhot_100_130 -> Fill(nvertex_data);
      if(ptPhot >= 130 && ptPhot < 175)    h_nvertex_data_ptPhot_130_175 -> Fill(nvertex_data);
      if(ptPhot >= 175 && ptPhot <= 5000)  h_nvertex_data_ptPhot_175_5000 -> Fill(nvertex_data);
    }      
    
    h_nvertex_data      ->Scale(1.0 / h_nvertex_data->Integral());
    h_nvertex_data_ptPhot_40_60      ->Scale(1.0 / h_nvertex_data_ptPhot_40_60->Integral());
    h_nvertex_data_ptPhot_60_85      ->Scale(1.0 / h_nvertex_data_ptPhot_60_85->Integral());
    h_nvertex_data_ptPhot_85_100    ->Scale(1.0 / h_nvertex_data_ptPhot_85_100->Integral());
    h_nvertex_data_ptPhot_100_130  ->Scale(1.0 / h_nvertex_data_ptPhot_100_130->Integral());
    h_nvertex_data_ptPhot_130_175  ->Scale(1.0 / h_nvertex_data_ptPhot_130_175->Integral());
    h_nvertex_data_ptPhot_175_5000->Scale(1.0 / h_nvertex_data_ptPhot_175_5000->Integral());
    
    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    h_nvertex_mc_ptPhot_60_85-> SetLineColor(kBlue);
    h_nvertex_data_ptPhot_60_85 -> SetLineColor(kRed);
    h_nvertex_data_ptPhot_60_85 -> Draw();
    h_nvertex_mc_ptPhot_60_85 -> Draw("same");
    c1-> SaveAs("Nvtx_ptPhot_60_85.png");
    
    
    //calcolo rapporto
    TH1D *h_ratio = (TH1D*)h_nvertex_data->Clone("h_ratio");
    h_ratio->Divide(h_nvertex_mc);
    
    TH1D *h_ratio_ptPhot_40_60 = (TH1D*)h_nvertex_data_ptPhot_40_60->Clone("h_ratio_ptPhot_40_60");
    h_ratio_ptPhot_40_60->Divide(h_nvertex_mc_ptPhot_40_60);
    
    TH1D *h_ratio_ptPhot_60_85 = (TH1D*)h_nvertex_data_ptPhot_60_85->Clone("h_ratio_ptPhot_60_85");
    h_ratio_ptPhot_60_85->Divide(h_nvertex_mc_ptPhot_60_85);
    
    TH1D *h_ratio_ptPhot_85_100 = (TH1D*)h_nvertex_data_ptPhot_85_100->Clone("h_ratio_ptPhot_85_100");
    h_ratio_ptPhot_85_100->Divide(h_nvertex_mc_ptPhot_85_100);
    
    TH1D *h_ratio_ptPhot_100_130 = (TH1D*)h_nvertex_data_ptPhot_100_130->Clone("h_ratio_ptPhot_100_130");
    h_ratio_ptPhot_100_130->Divide(h_nvertex_mc_ptPhot_100_130);
    
    TH1D *h_ratio_ptPhot_130_175 = (TH1D*)h_nvertex_data_ptPhot_130_175->Clone("h_ratio_ptPhot_130_175");
    h_ratio_ptPhot_130_175->Divide(h_nvertex_mc_ptPhot_130_175);
    
    TH1D *h_ratio_ptPhot_175_5000 = (TH1D*)h_nvertex_data_ptPhot_175_5000->Clone("h_ratio_ptPhot_175_5000");
    h_ratio_ptPhot_175_5000->Divide(h_nvertex_mc_ptPhot_175_5000);
    
    
    /////////////////////////////////////////////////////////////////////////
    TCanvas *c2 = new TCanvas("c2","c2",800,800);
    h_ratio_ptPhot_60_85 -> Draw();
    c2-> SaveAs("Ratio_nvertex_ptPhot_60_85.png");
    
    TH1D *h_nvertex_mc_reweighted = new TH1D("nvertex_mc_reweighted", "nvertex reweighted", 51, -0.5, 50.5);
    
    for (uint64_t i = 0; i < totalEvents_mc; i++) {      
      if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
      AnalysisTree_mc->GetEntry(i);
      PhotonTree_mc->GetEntry(i);
      float ptPhot;
      PhotonTree_mc->SetBranchAddress("pt",&ptPhot);
      UInt_t nvertex_mc;
      AnalysisTree_mc->SetBranchAddress("nvertex", &nvertex_mc);
      //cout<< ptPhot << endl;
      
      if(ptPhot>=60 && ptPhot<85){
	
	//	int bin = nvertex_mc +1;
	int bin = h_ratio_ptPhot_60_85 -> FindBin(nvertex_mc);
	float PUweight = h_ratio_ptPhot_60_85 -> GetBinContent(bin);
	
	double generator_weight;
	AnalysisTree_mc->SetBranchAddress("generator_weight", &generator_weight);   
	if(generator_weight == 0) generator_weight = 1;
	float event_weight;
	AnalysisTree_mc->SetBranchAddress("evtWeightTot", &event_weight);   
	
	double Weight = PUweight;
	
	//      cout<< "nvertex "<< nvertex_mc<<endl;
	//      cout<< "PUWeight "<< PUweight<<endl;	
	h_nvertex_mc_reweighted -> Fill(nvertex_mc, Weight);     
      }
    }
    
    h_nvertex_mc_reweighted      ->Scale(1.0 / h_nvertex_mc_reweighted->Integral());
    
    TCanvas *c3 = new TCanvas("c3","c3",800,800);
    h_nvertex_mc_reweighted -> SetLineColor(kBlue);
    h_nvertex_data_ptPhot_60_85 -> SetLineColor(kRed);
    h_nvertex_mc_reweighted -> Draw();
    h_nvertex_data_ptPhot_60_85 -> Draw("same");
    c3-> SaveAs("Nvtx_reweighted_ptPhot_60_85.png");
    
    
    TFile f_new("NvertexPU_ReReco_27Oct2015.root", "recreate");          
    
    //    h_nvertex_mc->Write();
    h_nvertex_mc_ptPhot_40_60  -> Write();      
    h_nvertex_mc_ptPhot_60_85  -> Write();      
    h_nvertex_mc_ptPhot_85_100  -> Write();      
    h_nvertex_mc_ptPhot_100_130  -> Write();      
    h_nvertex_mc_ptPhot_130_175  -> Write();      
    h_nvertex_mc_ptPhot_175_5000  -> Write();      
    //    h_nvertex_data->Write();
    h_nvertex_data_ptPhot_40_60  -> Write();      
    h_nvertex_data_ptPhot_60_85  -> Write();      
    h_nvertex_data_ptPhot_85_100  -> Write();      
    h_nvertex_data_ptPhot_100_130  -> Write();      
    h_nvertex_data_ptPhot_130_175  -> Write();      
    h_nvertex_data_ptPhot_175_5000  -> Write();      
    //    h_ratio -> Write();
    h_ratio_ptPhot_40_60 -> Write();
    h_ratio_ptPhot_60_85 -> Write();
    h_ratio_ptPhot_85_100 -> Write();
    h_ratio_ptPhot_100_130 -> Write();
    h_ratio_ptPhot_130_175 -> Write();
    h_ratio_ptPhot_175_5000 -> Write();
    f_new.Close();
    
}//main





