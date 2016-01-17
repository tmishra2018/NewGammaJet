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

  TGraphErrors* gr_Bal_Par0_vs_run = new TGraphErrors(0);                                                                   
  gr_Bal_Par0_vs_run ->SetName("Bal_Par0_vs_run");  
  TGraphErrors* gr_MPF_Par0_vs_run = new TGraphErrors(0);                                                                   
  gr_MPF_Par0_vs_run ->SetName("MPF_Par0_vs_run");  

  TGraphErrors* gr_Bal_Par1_vs_run = new TGraphErrors(0);                                                                   
  gr_Bal_Par1_vs_run ->SetName("Bal_Par1_vs_run");  
  TGraphErrors* gr_MPF_Par1_vs_run = new TGraphErrors(0);                                                                   
  gr_MPF_Par1_vs_run ->SetName("MPF_Par1_vs_run");  

  std::vector<double> Bal_Par0;
  std::vector<double> Bal_Par0_Error;
  Bal_Par0.push_back( 0.988 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.987 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.992 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.989 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.982 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.991 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.984 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.983 );
  Bal_Par0_Error.push_back( 0.004 );
  Bal_Par0.push_back( 0.981 );
  Bal_Par0_Error.push_back( 0.003 );

  std::vector<double> Bal_Par1;
  std::vector<double> Bal_Par1_Error;
  Bal_Par1.push_back( 0.012 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.011 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.021 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.018 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.006 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.020 );
  Bal_Par1_Error.push_back( 0.005 );
  Bal_Par1.push_back( 0.010 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.005 );
  Bal_Par1_Error.push_back( 0.006 );
  Bal_Par1.push_back( 0.006 );
  Bal_Par1_Error.push_back( 0.005 );

  std::vector<double> MPF_Par0;
  std::vector<double> MPF_Par0_Error;
  MPF_Par0.push_back( 0.977 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.974 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.974 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.973 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.970 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.974 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.974 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.970 );
  MPF_Par0_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.967 );
  MPF_Par0_Error.push_back( 0.002 );

  std::vector<double> MPF_Par1;
  std::vector<double> MPF_Par1_Error;
  MPF_Par1.push_back( 0.018 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.019 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.015 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.015 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.018 );
  MPF_Par1_Error.push_back( 0.005 );
  MPF_Par1.push_back( 0.023 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.026 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.015 );
  MPF_Par1_Error.push_back( 0.004 );
  MPF_Par1.push_back( 0.012 );
  MPF_Par1_Error.push_back( 0.004 );

  int run_index;

  for( int iplot =0 ; iplot<9 ; iplot++){                                                                                                                                     

    run_index = iplot+1;                                                         
                                                                                                                     
    gr_Bal_Par0_vs_run -> SetPoint(iplot, run_index, Bal_Par0.at(iplot));
    gr_Bal_Par0_vs_run -> SetPointError(iplot, 0, Bal_Par0_Error.at(iplot));
    gr_MPF_Par0_vs_run -> SetPoint(iplot, run_index, MPF_Par0.at(iplot));
    gr_MPF_Par0_vs_run -> SetPointError(iplot, 0, MPF_Par0_Error.at(iplot));

    gr_Bal_Par1_vs_run -> SetPoint(iplot, run_index, Bal_Par1.at(iplot));
    gr_Bal_Par1_vs_run -> SetPointError(iplot, 0, Bal_Par1_Error.at(iplot));
    gr_MPF_Par1_vs_run -> SetPoint(iplot, run_index, MPF_Par1.at(iplot));
    gr_MPF_Par1_vs_run -> SetPointError(iplot, 0, MPF_Par1_Error.at(iplot));


  }
                                                                                                                                                                              
  TLegend* legend = new TLegend(0.6, 0.9, 0.4, 0.7);                                                                                                                      
  legend->SetTextFont(42);                                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                                   
  legend->SetFillColor(kWhite);                                                                                                                                               
  legend->SetFillStyle(0);                                                                                                                                                    
  legend->SetTextSize(0.036);                                                                 
  legend->SetHeader("|#eta| < 1.3");                                                                                                        
  legend->AddEntry(gr_Bal_Par0_vs_run, "Balancing", "PL");                                                                                                                                        
  legend->AddEntry(gr_MPF_Par0_vs_run, "MPF", "PL");     

  gStyle -> SetOptStat(kFALSE); 


  TF1* Bal_Par0_Fit = new TF1("Bal_Par0_Fit", "[0]", 0, 10);
  Bal_Par0_Fit->SetParameter(0, 0.);
  gr_Bal_Par0_vs_run->Fit(Bal_Par0_Fit, "RQ");  
  double Bal_Par0_FitValue = Bal_Par0_Fit->GetParameter(0);
  double Bal_Par0_FitError = Bal_Par0_Fit->GetParError(0);

  TF1* MPF_Par0_Fit = new TF1("MPF_Par0_Fit", "[0]", 0, 10);
  MPF_Par0_Fit->SetParameter(0, 0.);
  MPF_Par0_Fit->SetLineColor(kBlue);
  gr_MPF_Par0_vs_run->Fit(MPF_Par0_Fit, "RQ");  
  double MPF_Par0_FitValue = MPF_Par0_Fit->GetParameter(0);
  double MPF_Par0_FitError = MPF_Par0_Fit->GetParError(0);
  
  TPaveText* fitlabel = new TPaveText(0.65, 0.85, 0.85, 0.7, "brNDC");
  fitlabel->SetTextSize(0.03);
  fitlabel->SetFillColor(0);
  TString Bal_fitLabelText = TString::Format("Fit: %.3f #pm %.3f", Bal_Par0_FitValue, Bal_Par0_FitError);
  TString MPF_fitLabelText = TString::Format("Fit: %.3f #pm %.3f", MPF_Par0_FitValue, MPF_Par0_FitError);
   fitlabel->AddText(Bal_fitLabelText);
  fitlabel->AddText(MPF_fitLabelText);
 
    
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gr_Bal_Par0_vs_run-> SetTitle();  
  gr_Bal_Par0_vs_run ->GetHistogram()-> SetXTitle("run index");  
  gr_Bal_Par0_vs_run ->GetHistogram()-> SetYTitle("Par [0]");  
  gr_Bal_Par0_vs_run -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_Bal_Par0_vs_run -> GetYaxis()->SetRangeUser(0.9, 1.05);  
  gr_Bal_Par0_vs_run -> GetXaxis()->SetRangeUser(0., 10);  

  gr_Bal_Par0_vs_run -> SetMarkerStyle(20);  
  gr_Bal_Par0_vs_run -> SetMarkerColor(kRed);  
  gr_Bal_Par0_vs_run -> SetLineColor(kRed); 

  gr_MPF_Par0_vs_run -> SetMarkerStyle(20);  
  gr_MPF_Par0_vs_run -> SetMarkerColor(kBlue);  
  gr_MPF_Par0_vs_run -> SetLineColor(kBlue); 

  gr_Bal_Par0_vs_run -> Draw("ZAP");
  gr_MPF_Par0_vs_run -> Draw("PSAME");

 fitlabel->Draw("same");
 legend->Draw();

  c1->SaveAs("Par0_vs_run.png");

  ////////////////////////

  TF1* Bal_Par1_Fit = new TF1("Bal_Par1_Fit", "[0]", 0, 10);
  Bal_Par1_Fit->SetParameter(0, 0.);
  gr_Bal_Par1_vs_run->Fit(Bal_Par1_Fit, "RQ");  
  double Bal_Par1_FitValue = Bal_Par1_Fit->GetParameter(0);
  double Bal_Par1_FitError = Bal_Par1_Fit->GetParError(0);

  TF1* MPF_Par1_Fit = new TF1("MPF_Par1_Fit", "[0]", 0, 10);
  MPF_Par1_Fit->SetParameter(0, 0.);
  MPF_Par1_Fit->SetLineColor(kBlue);
  gr_MPF_Par1_vs_run->Fit(MPF_Par1_Fit, "RQ");  
  double MPF_Par1_FitValue = MPF_Par1_Fit->GetParameter(0);
  double MPF_Par1_FitError = MPF_Par1_Fit->GetParError(0);
  
  TPaveText* fitlabel_Par1 = new TPaveText(0.65, 0.85, 0.85, 0.7, "brNDC");
  fitlabel_Par1->SetTextSize(0.03);
  fitlabel_Par1->SetFillColor(0);
  TString Bal_fitLabelText_Par1 = TString::Format("Fit: %.3f #pm %.3f", Bal_Par1_FitValue, Bal_Par1_FitError);
  TString MPF_fitLabelText_Par1 = TString::Format("Fit: %.3f #pm %.3f", MPF_Par1_FitValue, MPF_Par1_FitError);
  //  MPF_fitLabelText.SetLabelColor(kBlue);
  fitlabel_Par1->AddText(Bal_fitLabelText_Par1);
  fitlabel_Par1->AddText(MPF_fitLabelText_Par1);


  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  gr_Bal_Par1_vs_run -> SetTitle();  
  gr_Bal_Par1_vs_run ->GetHistogram()-> SetXTitle("run index");  
  gr_Bal_Par1_vs_run ->GetHistogram()-> SetYTitle("Par [1]");  
  gr_Bal_Par1_vs_run -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_Bal_Par1_vs_run -> GetYaxis()->SetRangeUser(-0.01, 0.05);  
  gr_Bal_Par1_vs_run -> GetXaxis()->SetRangeUser(0., 10);  

  gr_Bal_Par1_vs_run -> SetMarkerStyle(20);  
  gr_Bal_Par1_vs_run -> SetMarkerColor(kRed);  
  gr_Bal_Par1_vs_run -> SetLineColor(kRed); 

  gr_MPF_Par1_vs_run -> SetMarkerStyle(20);  
  gr_MPF_Par1_vs_run -> SetMarkerColor(kBlue);  
  gr_MPF_Par1_vs_run -> SetLineColor(kBlue); 

  gr_Bal_Par1_vs_run -> Draw("ZAP");
  gr_MPF_Par1_vs_run -> Draw("PSAME");

  fitlabel_Par1->Draw("same");
  legend->Draw();

  c2->SaveAs("Par1_vs_run.png");


}//main





