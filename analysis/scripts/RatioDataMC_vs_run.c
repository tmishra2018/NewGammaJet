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

  TGraphErrors* gr_Bal_vs_run = new TGraphErrors(0);                                                                   
  gr_Bal_vs_run ->SetName("Bal_vs_run");  
  TGraphErrors* gr_MPF_vs_run = new TGraphErrors(0);                                                                   
  gr_MPF_vs_run ->SetName("MPF_vs_run");  

  std::vector<double> Bal_Par0;
  std::vector<double> Bal_Error;
  Bal_Par0.push_back( 0.984 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.983 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.986 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.983 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.980 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.986 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.981 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.981 );
  Bal_Error.push_back( 0.003 );
  Bal_Par0.push_back( 0.979 );
  Bal_Error.push_back( 0.003 );


  std::vector<double> MPF_Par0;
  std::vector<double> MPF_Error;
  MPF_Par0.push_back( 0.973 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.971 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.971 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.971 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.967 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.972 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.970 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.969 );
  MPF_Error.push_back( 0.002 );
  MPF_Par0.push_back( 0.966 );
  MPF_Error.push_back( 0.002 );

  int run_index;

  for( int iplot =0 ; iplot<9 ; iplot++){                                                                                                                                     

    run_index = iplot+1;                                                         
                                                                                                                     
    gr_Bal_vs_run -> SetPoint(iplot, run_index, Bal_Par0.at(iplot));
    gr_Bal_vs_run -> SetPointError(iplot, 0, Bal_Error.at(iplot));
    gr_MPF_vs_run -> SetPoint(iplot, run_index, MPF_Par0.at(iplot));
    gr_MPF_vs_run -> SetPointError(iplot, 0, MPF_Error.at(iplot) ) ;

  }

  TF1* ratioFit = new TF1("ratioFit", "[0]", 0, 10);
  ratioFit->SetParameter(0, 0.);
  ratioFit->SetLineColor(46);
  ratioFit->SetLineWidth(2);
  gr_Bal_vs_run->Fit(ratioFit, "RQ");

  double fitValue = ratioFit->GetParameter(0);
  double fitError = ratioFit->GetParError(0);
  double chiSquare = ratioFit->GetChisquare();
  double NDF = ratioFit->GetNDF();
  
 std::cout << "Balancing: " << std::endl;
 std::cout << "Fit value: " << fitValue << " #pm " << fitError << std::endl;
  std::cout << "-> ChiSquare: " << chiSquare << "   NDF: " << NDF << std::endl;

 TPaveText* fitlabel = new TPaveText(0.4, 0.2, 0.6, 0.3, "brNDC");
 fitlabel->SetTextSize(0.03);
 fitlabel->SetFillColor(0);
 TString fitLabelText = TString::Format("Fit: %.3f #pm %.3f", fitValue, fitError);
 TString fitLabelText_2 = TString::Format("#chi^{2}/NDF: %.3f / %.3f", chiSquare, NDF);
 fitlabel->AddText(fitLabelText);
 fitlabel->AddText(fitLabelText_2);

 ///////////////////////////////////////////////////////////////////

 TF1* ratioFit_MPF = new TF1("ratioFit_MPF", "[0]", 0, 10);
  ratioFit_MPF->SetParameter(0, 0.);
  ratioFit_MPF->SetLineColor(46);
  ratioFit_MPF->SetLineWidth(2);
  gr_MPF_vs_run->Fit(ratioFit_MPF, "RQ");

  double fitValue_MPF = ratioFit_MPF->GetParameter(0);
  double fitError_MPF = ratioFit_MPF->GetParError(0);
  double chiSquare_MPF = ratioFit_MPF->GetChisquare();
  double NDF_MPF = ratioFit_MPF->GetNDF();


 TPaveText* fitlabel_MPF = new TPaveText(0.4, 0.2, 0.6, 0.3, "brNDC");
 fitlabel_MPF->SetTextSize(0.03);
 fitlabel_MPF->SetFillColor(0);
 TString fitLabel_MPFText = TString::Format("Fit: %.3f #pm %.3f", fitValue_MPF, fitError_MPF);
 TString fitLabel_MPFText_2 = TString::Format("#chi^{2}/NDF: %.3f / %.3f", chiSquare_MPF, NDF_MPF);
 fitlabel_MPF->AddText(fitLabel_MPFText);
 fitlabel_MPF->AddText(fitLabel_MPFText_2);

 std::cout << "MPF: " << std::endl;
 std::cout << "Fit value: " << fitValue_MPF << " #pm " << fitError_MPF << std::endl;
 std::cout << "-> ChiSquare: " << chiSquare_MPF << "   NDF: " << NDF_MPF << std::endl;


  ///////////////////////////// drawer                                                                                                                                                                              
  TLegend* legend = new TLegend(0.9, 0.9, 0.7, 0.7);                                                                                                                      
  legend->SetTextFont(42);                                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                                   
  legend->SetFillColor(kWhite);                                                                                                                                               
  legend->SetFillStyle(0);                                                                                                                                                    
  legend->SetTextSize(0.036);                                                                 
  legend->SetHeader("|#eta| < 1.3");                                                                                                        
  //  legend->AddEntry(gr_Bal_vs_run, "Balancing", "PL");                                                                                                                                        
  //  legend->AddEntry(gr_MPF_vs_run, "MPF", "PL");     

  gStyle -> SetOptStat(kFALSE); 

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gr_Bal_vs_run -> SetTitle();  
  gr_Bal_vs_run ->GetHistogram()-> SetXTitle("run index");  
  gr_Bal_vs_run ->GetHistogram()-> SetYTitle("Data/MC");  
  gr_Bal_vs_run -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_Bal_vs_run -> GetYaxis()->SetRangeUser(0.9, 1.05);  
  gr_Bal_vs_run -> GetXaxis()->SetRangeUser(0., 10);  

  gr_Bal_vs_run -> SetMarkerStyle(20);  
  gr_Bal_vs_run -> SetMarkerColor(kRed);  
  gr_Bal_vs_run -> SetLineColor(kRed); 

  gr_Bal_vs_run -> Draw("ZAP");

  fitlabel->Draw("same");
  legend->Draw();

  c1->SaveAs("Balancing_DataMCRatio_vs_run.png");

  TCanvas *c2 = new TCanvas("c2","c2",800,800);

  gr_MPF_vs_run -> SetTitle();  
  gr_MPF_vs_run ->GetHistogram()-> SetXTitle("run index");  
  gr_MPF_vs_run ->GetHistogram()-> SetYTitle("Data/MC");  
  gr_MPF_vs_run -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_MPF_vs_run -> GetYaxis()->SetRangeUser(0.9, 1.05);  
  gr_MPF_vs_run -> GetXaxis()->SetRangeUser(0., 10);  

  gr_MPF_vs_run -> SetMarkerStyle(20);  
  gr_MPF_vs_run -> SetMarkerColor(kRed);  
  gr_MPF_vs_run -> SetLineColor(kRed); 

  gr_MPF_vs_run -> Draw("ZAP");

  fitlabel_MPF->Draw("same");
  legend->Draw();

  c2->SaveAs("MPF_DataMCRatio_vs_run.png");


}//main





