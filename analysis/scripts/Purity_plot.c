#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
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

#define LIGHT_RED TColor::GetColor(0xcf, 0xa0, 0xa1)


int main(int argc, char* argv[]) {
  

  TFile GJetFile("../tuples/GJET_Pythia/PhotonJet_GJet_plus_QCD_Pythia_2016-02-26_alphacut030_PFlowAK4chs.root");
  TFile QCDFile("../tuples/QCD_Pythia/PhotonJet_QCD_Pythia_2015-02-26_alphacut030_PFlowAK4chs.root");

  TH1D *h_GJet_ptPhoton    = (TH1D*)GJetFile.Get("analysis/ptPhoton_Binned");
  TH1D *h_QCD_ptPhoton    = (TH1D*)QCDFile.Get("analysis/ptPhoton_Binned");

  TH1D *h_GJet_plus_QCD = (TH1D*)h_GJet_ptPhoton->Clone("h_GJet_plus_QCD");
  h_GJet_plus_QCD-> Add(h_QCD_ptPhoton, +1);

  TH1D *h_purity = (TH1D*)h_GJet_ptPhoton->Clone("h_purity");
  h_purity -> Divide(h_GJet_plus_QCD);



  for(int ii = 1 ; ii<15 ; ii++){                                                                                                                 
    double Nevents_GJet = h_GJet_ptPhoton->GetBinContent(ii);                                                                                     
    double Nevents_QCD = h_QCD_ptPhoton->GetBinContent(ii);                                                                                       
    double NeventsTot = h_GJet_plus_QCD ->GetBinContent(ii);                                                                                               
    double purity = h_purity ->GetBinContent(ii);                                                                                               
    std::cout<<"NEVENTS"<<std::endl;                                                                                                              
    std::cout<<ii<<"th bin" << std::endl;                                                                                                         
    std::cout<<"GJET= "<< Nevents_GJet  << std::endl;                                                                                             
    std::cout<<"QCD= "<<  Nevents_QCD  << std::endl;                                                                                              
    std::cout<<"GJET+QCD= "<< NeventsTot  << std::endl;
    std::cout<<"purity = "<< purity*100 <<"%"<< std::endl;
  }  

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  h_purity -> SetStats(0);  
  h_purity -> SetTitle(" ");  
  h_purity -> SetXTitle("pT(#gamma) [GeV]");  
  h_purity -> SetYTitle("Purity");  
  h_purity -> GetYaxis()->SetTitleOffset(1.35);  

  h_purity -> SetMarkerStyle(20);                                                                                                              
  h_purity -> SetMarkerSize(0.5);                        
  //  h_purity -> SetMarkerColor(kBlue - 6); 

  //  h_purity -> Draw("ZAP");
  h_purity -> Draw();
  c1->SaveAs("Purity_plot.png");

}






