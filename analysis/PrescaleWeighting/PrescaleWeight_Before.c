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

    TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_5/src/JetMETCorrections/GammaJetFilter/analysis/tuples/GJET_MC/PhotonJet_2ndLevel_GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_v2.root"); 
    TTree* AnalysisTree_mc = (TTree*) f1.Get("gammaJet/analysis");
    TTree* PhotonTree_mc = (TTree*) f1.Get("gammaJet/photon");
    uint64_t totalEvents_mc = AnalysisTree_mc->GetEntries();
    cout<< totalEvents_mc << endl;

    TFile f2("/cmshome/fpreiato/GammaJet/CMSSW_7_4_5/src/JetMETCorrections/GammaJetFilter/analysis/tuples/Data/PhotonJet_2ndLevel_SinglePhoton_GoldenJson_v2.root ");
    TTree* AnalysisTree_data = (TTree*) f2.Get("gammaJet/analysis");
    TTree* PhotonTree_data = (TTree*) f2.Get("gammaJet/photon");
    uint64_t totalEvents_data = AnalysisTree_data->GetEntries();
    //    TParameter<double> luminosity = <TParameter<double>(f2.Get("analysis/luminosity")->GetVal());

    //    double lumi = f2.Get("gammaJet/total_luminosity")->GetVal();

    double lumi = 40.028;

    TFile f_PU("/cmshome/fpreiato/GammaJet/CMSSW_7_4_5/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting/NvertexPU.root"); 
    TH1D *h_PU;

 
   ///////////////////////////////////////////////////////////////////
    
    double arraybins[7] = {40, 60, 85,100,130,175,1000};

    TH1D *h_mc = new TH1D("h_mc", "h_mc", 6, arraybins);

    //loop
    for (uint64_t i = 0; i < totalEvents_mc; i++) {      
	  //for (uint64_t i = 0; i < 10000; i++) {      
      if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
      
      PhotonTree_mc->GetEntry(i);
      AnalysisTree_mc->GetEntry(i);

      float ptPhot;
      PhotonTree_mc->SetBranchAddress("pt",&ptPhot);

      double generator_weight;
      AnalysisTree_mc->SetBranchAddress("generator_weight", &generator_weight);   
      if(generator_weight == 0) generator_weight = 1;      
      float event_weight;
      AnalysisTree_mc->SetBranchAddress("evtWeightTot", &event_weight);         
      UInt_t nvertex_mc;                                                                                                                                        
      AnalysisTree_mc->SetBranchAddress("nvertex", &nvertex_mc);  


      //      cout<< ptPhot << endl;

      if(ptPhot<40) continue;

      if(ptPhot >= 40 && ptPhot < 60)          h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_40_60");  
      if(ptPhot >= 60 && ptPhot < 85)          h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_60_85");  
      if(ptPhot >= 85 && ptPhot < 100)        h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_85_100");  
      if(ptPhot >= 100 && ptPhot < 130)      h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_100_130");  
      if(ptPhot >= 130 && ptPhot < 175)      h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_130_175");  
      if(ptPhot >= 175 && ptPhot < 1000)    h_PU = (TH1D*)f_PU.Get("h_ratio_ptPhot_175_1000");  

      //      cout << "OKKK"<<endl;

      int bin = nvertex_mc+1;
      double  PUWeight = h_PU->GetBinContent(bin);

      //      cout << "OKKK"<<endl;

      //      cout<< "event weight  "<<event_weight<<endl;
      //      cout<< "generator weight  "<<generator_weight<<endl;
      //      cout<< "lumi  "<<lumi<<endl;
      //      cout<< "PU weight  "<<PUWeight<<endl;

      double Weight = generator_weight * event_weight * lumi * PUWeight;

      //      cout<< "Weight Tot  "<<Weight<<endl;
      


      h_mc -> Fill(ptPhot, Weight);
      
    }
    
    
    TH1D *h_data = new TH1D("h_data", "h_data", 6, arraybins);

    //loop
    for (uint64_t i = 0; i < totalEvents_data; i++) {      
      //   for (uint64_t i = 0; i < 10000; i++) {      
      if(i == 0 || i % 5000 == 0) cout<< "Events processed "<< i <<endl;
      
      AnalysisTree_data->GetEntry(i);
      PhotonTree_data->GetEntry(i);
      float ptPhot;
      PhotonTree_data->SetBranchAddress("pt",&ptPhot);
      //cout<< ptPhot << endl;
      
      if(ptPhot<40) continue;

      h_data -> Fill(ptPhot);     
      
    }      
    

    
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->SetLogx();
    h_mc-> SetLineColor(kBlue);
    h_data -> SetLineColor(kRed);
    h_data -> Draw();
    h_mc -> Draw("same");
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
      AnalysisTree_data->GetEntry(i);
      PhotonTree_data->GetEntry(i);
      
      float ptPhot;
      PhotonTree_data->SetBranchAddress("pt",&ptPhot);
     
      int bin;

      if(ptPhot<40) continue;

      if(ptPhot >=40 && ptPhot <60)       bin =1;
      if(ptPhot >=60 && ptPhot <85)       bin =2;
      if(ptPhot >=85 && ptPhot <100)     bin =3; 
      if(ptPhot >=100 && ptPhot <130)   bin =4;
      if(ptPhot >=130 && ptPhot <175)   bin =5;
      if(ptPhot >=175 && ptPhot <1000) bin =6;


      double Prescale = h_ratio->GetBinContent(bin);

      cout<< "prescale applicato  " << Prescale <<endl;
	
      h_data_reweighted -> Fill(ptPhot, Prescale);     
      }
      

      
      //    h_nvertex_mc_reweighted      ->Scale(1.0 / h_nvertex_mc_reweighted->Integral());
      
      TCanvas *c3 = new TCanvas("c3","c3",800,800);
      //    c3->SetLogy();
      c3->SetLogx();
      h_mc -> SetLineColor(kBlue);
      h_data_reweighted -> SetLineColor(kRed);
      h_data_reweighted -> Draw();
      h_mc -> Draw("same");
      c3-> SaveAs("Data_reweighted.png");


      TFile f_new("Prescale.root", "recreate");          
    
      h_mc  -> Write();      
      h_data->Write();
      h_ratio -> Write();
      h_data_reweighted->Write();
      f_new.Close();
	
}//main





