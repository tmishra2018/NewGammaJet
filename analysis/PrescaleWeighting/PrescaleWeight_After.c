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
#include <TParameter.h>
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


  TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_12_patch4/src/JetMETCorrections/GammaJetFilter/analysis/PrescaleWeighting/MC_ptPhot_scaled.root");
  TH1D *h_mc = (TH1D*)f1.Get("h_mc");

  TFile f2("/cmshome/fpreiato/GammaJet/CMSSW_7_4_12_patch4/src/JetMETCorrections/GammaJetFilter/analysis/tuples/Data/PhotonJet_SinglePhoton_25ns_Run2015D_09Oct_NoPrescale_alphacut030_PFlowAK4chs.root");
  TTree* PhotonTree_data = (TTree*) f2.Get("photon");
  uint64_t totalEvents_data = PhotonTree_data->GetEntries();
  
  cout<< totalEvents_data << endl;  
  
  double lumi = -1;
  TParameter<double>* dLumi = static_cast<TParameter<double>*>(f2.Get("analysis/luminosity"));
  lumi = dLumi->GetVal();
  cout<< "lumi  " << lumi<< endl;  
  
  ///////////////////////////////////////////////////////////////////
  
  double arraybins[7] = {40, 60, 85,100,130,175,5000};
      
  TH1D *h_data = new TH1D("h_data", "h_data", 6, arraybins);
  
    //loop
  for (uint64_t i = 0; i < totalEvents_data; i++) {      
    //   for (uint64_t i = 0; i < 10000; i++) {      
    if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
    
    PhotonTree_data->GetEntry(i);
    float ptPhot;
    PhotonTree_data->SetBranchAddress("pt",&ptPhot);
    
    h_data -> Fill(ptPhot);        
  }      
  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->SetLogy();
  h_mc-> SetLineColor(kBlue);
  h_data -> SetLineColor(kRed);
  h_mc -> Draw();
  h_data -> Draw("same");
  c1-> SaveAs("Distributions.png");
  ///////////////////////////////////////////////////
  
  //calcolo rapporto
  
  TH1D *h_ratio = (TH1D*)h_mc->Clone("h_ratio");                                                                                                            
  h_ratio->Divide(h_data);                                                                                                                                      
  
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  c2->SetLogx();
  h_ratio -> Draw();
  c2-> SaveAs("Ratio.png");
  /////////////////////////////////////////////////////////////////////////    
  
  TH1D *h_data_reweighted = new TH1D("h_data_reweighted", "data reweighted", 6, arraybins);
  
  for (uint64_t i = 0; i < totalEvents_data; i++) {      
    //  for (uint64_t i = 0; i < 10000; i++) {      
    if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
    PhotonTree_data->GetEntry(i);
    
    float ptPhot;
    PhotonTree_data->SetBranchAddress("pt",&ptPhot);
    
    int bin;
    
    if(ptPhot >=40 && ptPhot <60)       bin =1;
    if(ptPhot >=60 && ptPhot <85)       bin =2;
    if(ptPhot >=85 && ptPhot <100)     bin =3; 
    if(ptPhot >=100 && ptPhot <130)   bin =4;
    if(ptPhot >=130 && ptPhot <175)   bin =5;
    if(ptPhot >=175 && ptPhot <=5000) bin =6;
    
    
    double Prescale = h_ratio->GetBinContent(bin);
    
    //    cout<< "prescale applicato  " << Prescale <<endl;
    
    h_data_reweighted -> Fill(ptPhot, Prescale);     
  }
  
  
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  c3->SetLogy();
  //      c3->SetLogx();
  h_mc -> SetLineColor(kBlue);
  h_data_reweighted -> SetLineColor(kRed);
  h_data_reweighted -> Draw();
  h_mc -> Draw("same");
  c3-> SaveAs("Data_reweighted.png");
  
  
  TFile f_new("Prescale_Run2015D_09Oct_alphacut030.root", "recreate");          
  
  h_mc  -> Write();      
  h_data->Write();
  h_ratio -> Write();
  h_data_reweighted->Write();
  f_new.Close();
  
}//main





