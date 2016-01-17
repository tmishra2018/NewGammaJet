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

#define LIGHT_RED TColor::GetColor(0xcf, 0xa0, 0xa1)

using namespace std;


int main(int argc, char* argv[]) {

  TFile f1("GammaJet_ppCollision2015_17-11-2015.root"); 
  TGraphErrors *Balancing = (TGraphErrors*)f1.Get("resp_PtBalchs_a30_eta00_13");
  TGraphErrors *MPF         = (TGraphErrors*)f1.Get("resp_MPFchs_a30_eta00_13");

  TGraphErrors *Balancing_Extrap = (TGraphErrors*)f1.Get("resp_PtBalchs_extrap_a30_eta00_13");
  TGraphErrors *MPF_Extrap         = (TGraphErrors*)f1.Get("resp_MPFchs_extrap_a30_eta00_13");

  TLegend *legend=new TLegend(0.9, 0.9, 0.6, 0.75);
  legend->AddEntry(Balancing,"Balancing","LP");
  legend->AddEntry(MPF,"MPF","LP");

  TPaveText* NoExtrap = new TPaveText(0.9, 0.1, 0.6, 0.2, "brNDC");
  //  NoExtrap->SetTextSize(0.08);
  NoExtrap->SetFillColor(0);
  NoExtrap->AddText("No extrapolation");

  TPaveText* SiExtrap = new TPaveText(0.9, 0.1, 0.6, 0.2, "brNDC");
  //  NoExtrap->SetTextSize(0.08);
  SiExtrap->SetFillColor(0);
  SiExtrap->AddText("With extrapolation");
  


  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1 -> SetLogx();
  Balancing ->GetXaxis()->SetMoreLogLabels();  
  Balancing ->GetXaxis()->SetNoExponent();  
  Balancing ->GetYaxis()-> SetTitle("Jet response (ratio)");  
  Balancing ->GetYaxis()-> SetTitleOffset(1.50);  
  Balancing ->GetXaxis()-> SetTitle("p_{T} (GeV)");  
  Balancing ->SetMinimum(0.9);  
  Balancing -> SetMaximum(1.05);
  Balancing -> SetLineColor(kRed);  
  Balancing -> SetMarkerColor(kRed);
  Balancing -> SetMarkerStyle(23);  
  MPF -> SetLineColor(kBlue);
  MPF -> SetMarkerColor(kBlue);
  MPF -> SetMarkerStyle(20);
  Balancing -> Draw("ZAP");
  MPF -> Draw("PSAME");
  legend->Draw();
  NoExtrap->Draw("same");
  c1-> SaveAs("NoExtrapolation.png");
 
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  c2 -> SetLogx();
  Balancing_Extrap ->GetXaxis()->SetMoreLogLabels();  
  Balancing_Extrap ->GetXaxis()->SetNoExponent();  
  Balancing_Extrap ->GetYaxis()-> SetTitle("Jet response (ratio)");  
  Balancing_Extrap ->GetYaxis()-> SetTitleOffset(1.50);  
  Balancing_Extrap ->GetXaxis()-> SetTitle("p_{T} (GeV)");  
  Balancing_Extrap ->SetMinimum(0.9);  
  Balancing_Extrap -> SetMaximum(1.05);
  Balancing_Extrap ->GetYaxis()-> SetLimits(0.85, 1.2);  
  Balancing_Extrap -> SetLineColor(kRed);  
  Balancing_Extrap -> SetMarkerColor(kRed);  
  Balancing_Extrap -> SetMarkerStyle(23);  
  MPF_Extrap -> SetLineColor(kBlue);
  MPF_Extrap -> SetMarkerColor(kBlue);
  MPF_Extrap -> SetMarkerStyle(20);
  Balancing_Extrap -> Draw("ZAP");
  MPF_Extrap -> Draw("PSAME");
  legend->Draw();
  SiExtrap->Draw("same");
  c2-> SaveAs("Extrapolation.png");




}//main





