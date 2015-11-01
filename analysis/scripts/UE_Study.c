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


  TFile f1("/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/analysis/draw/Plot/2015-10-28/AlphaCut030/PhotonJet_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4chs.root");
  TH1F* MC_rho = (TH1F*) f1.Get("analysis/rho_passedID");
  TH1F* MC_mu = (TH1F*) f1.Get("analysis/mu");
  TH1F* MC_npvGood = (TH1F*) f1.Get("analysis/npvGood");
  TH2F* MC_rho_vs_mu = (TH2F*) f1.Get("analysis/rho_vs_mu");
  TH2F* MC_npvGood_vs_mu = (TH2F*) f1.Get("analysis/npvGood_vs_mu");
  

  TFile f2("/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/analysis/draw/Plot/2015-10-28/AlphaCut030/PhotonJet_SinglePhoton__Run2015D_2015-10-28_alphacut030_PFlowAK4chs.root");
  TH1F* Data_rho = (TH1F*) f2.Get("analysis/rho_passedID");
  TH1F* Data_mu = (TH1F*) f2.Get("analysis/mu");
  TH1F* Data_npvGood = (TH1F*) f2.Get("analysis/npvGood");
  TH2F* Data_rho_vs_mu = (TH2F*) f2.Get("analysis/rho_vs_mu");
  TH2F* Data_npvGood_vs_mu = (TH2F*) f2.Get("analysis/npvGood_vs_mu");

  double lumi = -1;
  TParameter<double>* dLumi = static_cast<TParameter<double>*>(f2.Get("analysis/luminosity"));
  lumi = dLumi->GetVal();
  cout<< "lumi  " << lumi<< endl;  

  MC_rho -> Scale(lumi);
  MC_mu -> Scale(lumi);
  MC_npvGood -> Scale(lumi);
  MC_rho_vs_mu -> Scale(lumi);
  MC_npvGood_vs_mu -> Scale(lumi);

  /////////////////////////////////////////////

  TCanvas *c1 = new TCanvas("c1","" , 800, 800);
  Data_rho ->Draw();
  MC_rho -> SetLineColor(kRed);
  MC_rho -> Draw("same"); 
  c1->SaveAs("Test.png");

  TFile * f_new = new TFile("UE_Study.root", "RECREATE");
  f_new->cd();
  f_new->mkdir("data");
  f_new->mkdir("MC");

  f_new -> cd("data");
  Data_rho                      -> Write();
  Data_mu                      -> Write();
  Data_npvGood             -> Write();
  Data_rho_vs_mu           -> Write();
  Data_npvGood_vs_mu  -> Write();
  f_new -> cd("MC");
  MC_rho                      -> Write();
  MC_mu                      -> Write();
  MC_npvGood             -> Write();
  MC_rho_vs_mu           -> Write();
  MC_npvGood_vs_mu  -> Write();

  f_new->Close();








}
