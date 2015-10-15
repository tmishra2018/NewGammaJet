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

  TFile f1("../tuples/Data/PhotonJet_TEST_PFlowAK4chs.root"); 
  TH1F *resp_balancing_eta0013_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta0013_run_256630_256678");
  TH1F *resp_balancing_eta0013_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta0013_run_256728_256730");
  TH1F *resp_balancing_eta0013_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta0013_run_256734_256844");
  TH1F *resp_balancing_eta0013_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta0013_run_256866_257600");

  TH1F *resp_balancing_eta1319_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1319_run_256630_256678");
  TH1F *resp_balancing_eta1319_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1319_run_256728_256730");
  TH1F *resp_balancing_eta1319_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1319_run_256734_256844");
  TH1F *resp_balancing_eta1319_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1319_run_256866_257600");

  TH1F *resp_balancing_eta1925_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1925_run_256630_256678");
  TH1F *resp_balancing_eta1925_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1925_run_256728_256730");
  TH1F *resp_balancing_eta1925_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1925_run_256734_256844");
  TH1F *resp_balancing_eta1925_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta1925_run_256866_257600");

  TH1F *resp_balancing_eta2530_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta2530_run_256630_256678");
  TH1F *resp_balancing_eta2530_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta2530_run_256728_256730");
  TH1F *resp_balancing_eta2530_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta2530_run_256734_256844");
  TH1F *resp_balancing_eta2530_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta2530_run_256866_257600");

  TH1F *resp_balancing_eta3032_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3032_run_256630_256678");
  TH1F *resp_balancing_eta3032_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3032_run_256728_256730");
  TH1F *resp_balancing_eta3032_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3032_run_256734_256844");
  TH1F *resp_balancing_eta3032_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3032_run_256866_257600");

  TH1F *resp_balancing_eta3252_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3252_run_256630_256678");
  TH1F *resp_balancing_eta3252_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3252_run_256728_256730");
  TH1F *resp_balancing_eta3252_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3252_run_256734_256844");
  TH1F *resp_balancing_eta3252_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_balancing_eta3252_run_256866_257600");

  ///////////////// MPF

  TH1F *resp_mpf_eta0013_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta0013_run_256630_256678");
  TH1F *resp_mpf_eta0013_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta0013_run_256728_256730");
  TH1F *resp_mpf_eta0013_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta0013_run_256734_256844");
  TH1F *resp_mpf_eta0013_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta0013_run_256866_257600");

  TH1F *resp_mpf_eta1319_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1319_run_256630_256678");
  TH1F *resp_mpf_eta1319_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1319_run_256728_256730");
  TH1F *resp_mpf_eta1319_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1319_run_256734_256844");
  TH1F *resp_mpf_eta1319_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1319_run_256866_257600");

  TH1F *resp_mpf_eta1925_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1925_run_256630_256678");
  TH1F *resp_mpf_eta1925_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1925_run_256728_256730");
  TH1F *resp_mpf_eta1925_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1925_run_256734_256844");
  TH1F *resp_mpf_eta1925_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta1925_run_256866_257600");

  TH1F *resp_mpf_eta2530_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta2530_run_256630_256678");
  TH1F *resp_mpf_eta2530_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta2530_run_256728_256730");
  TH1F *resp_mpf_eta2530_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta2530_run_256734_256844");
  TH1F *resp_mpf_eta2530_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta2530_run_256866_257600");

  TH1F *resp_mpf_eta3032_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3032_run_256630_256678");
  TH1F *resp_mpf_eta3032_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3032_run_256728_256730");
  TH1F *resp_mpf_eta3032_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3032_run_256734_256844");
  TH1F *resp_mpf_eta3032_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3032_run_256866_257600");

  TH1F *resp_mpf_eta3252_run_256630_256678 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3252_run_256630_256678");
  TH1F *resp_mpf_eta3252_run_256728_256730 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3252_run_256728_256730");
  TH1F *resp_mpf_eta3252_run_256734_256844 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3252_run_256734_256844");
  TH1F *resp_mpf_eta3252_run_256866_257600 = (TH1F*)f1.Get("analysis/run/resp_mpf_eta3252_run_256866_257600");
  ////////////////////
                                                                                                                                                                             
  std::vector<std::pair<float, float> > etaBins;                                                                                                                                
  etaBins.push_back(std::make_pair(0., 1.3) );                                                                                                                                
  etaBins.push_back(std::make_pair(1.3, 2.0) );                                                                                                                               
  etaBins.push_back(std::make_pair(2.0, 2.5) );                                                                                                                               
  etaBins.push_back(std::make_pair(2.5, 3.0) );                                                                                                                               
  etaBins.push_back(std::make_pair(3.0, 3.2) );                                                                                                                               
  etaBins.push_back(std::make_pair(3.2, 5.2) );    

  std::vector<double> etaBinsError;                                                                                                                                
  etaBinsError.push_back( 0.65 );                                                                                                                                
  etaBinsError.push_back( 0.35 );                                                                                                                               
  etaBinsError.push_back( 0.25 );                                                                                                                               
  etaBinsError.push_back( 0.25 );                                                                                                                               
  etaBinsError.push_back( 0.1 );                                                                                                                               
  etaBinsError.push_back( 1.0 );    

  std::vector<double > valueMean_run1;           
  valueMean_run1.push_back(resp_balancing_eta0013_run_256630_256678 ->GetMean() );
  valueMean_run1.push_back(resp_balancing_eta1319_run_256630_256678 ->GetMean() );
  valueMean_run1.push_back(resp_balancing_eta1925_run_256630_256678 ->GetMean() );
  valueMean_run1.push_back(resp_balancing_eta2530_run_256630_256678 ->GetMean() );
  valueMean_run1.push_back(resp_balancing_eta3032_run_256630_256678 ->GetMean() );
  valueMean_run1.push_back(resp_balancing_eta3252_run_256630_256678 ->GetMean() );
  std::vector<double > valueMeanError_run1;           
  valueMeanError_run1.push_back(resp_balancing_eta0013_run_256630_256678 ->GetMeanError() );
  valueMeanError_run1.push_back(resp_balancing_eta1319_run_256630_256678 ->GetMeanError() );
  valueMeanError_run1.push_back(resp_balancing_eta1925_run_256630_256678 ->GetMeanError() );
  valueMeanError_run1.push_back(resp_balancing_eta2530_run_256630_256678 ->GetMeanError() );
  valueMeanError_run1.push_back(resp_balancing_eta3032_run_256630_256678 ->GetMeanError() );
  valueMeanError_run1.push_back(resp_balancing_eta3252_run_256630_256678 ->GetMeanError() );
  ///////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMean_run2;           
  valueMean_run2.push_back(resp_balancing_eta0013_run_256728_256730 ->GetMean() );
  valueMean_run2.push_back(resp_balancing_eta1319_run_256728_256730 ->GetMean() );
  valueMean_run2.push_back(resp_balancing_eta1925_run_256728_256730 ->GetMean() );
  valueMean_run2.push_back(resp_balancing_eta2530_run_256728_256730 ->GetMean() );
  valueMean_run2.push_back(resp_balancing_eta3032_run_256728_256730 ->GetMean() );
  valueMean_run2.push_back(resp_balancing_eta3252_run_256728_256730 ->GetMean() );
  std::vector<double > valueMeanError_run2;           
  valueMeanError_run2.push_back(resp_balancing_eta0013_run_256728_256730 ->GetMeanError() );
  valueMeanError_run2.push_back(resp_balancing_eta1319_run_256728_256730 ->GetMeanError() );
  valueMeanError_run2.push_back(resp_balancing_eta1925_run_256728_256730 ->GetMeanError() );
  valueMeanError_run2.push_back(resp_balancing_eta2530_run_256728_256730 ->GetMeanError() );
  valueMeanError_run2.push_back(resp_balancing_eta3032_run_256728_256730 ->GetMeanError() );
  valueMeanError_run2.push_back(resp_balancing_eta3252_run_256728_256730 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMean_run3;           
  valueMean_run3.push_back(resp_balancing_eta0013_run_256734_256844 ->GetMean() );
  valueMean_run3.push_back(resp_balancing_eta1319_run_256734_256844 ->GetMean() );
  valueMean_run3.push_back(resp_balancing_eta1925_run_256734_256844 ->GetMean() );
  valueMean_run3.push_back(resp_balancing_eta2530_run_256734_256844 ->GetMean() );
  valueMean_run3.push_back(resp_balancing_eta3032_run_256734_256844 ->GetMean() );
  valueMean_run3.push_back(resp_balancing_eta3252_run_256734_256844 ->GetMean() );
  std::vector<double > valueMeanError_run3;           
  valueMeanError_run3.push_back(resp_balancing_eta0013_run_256734_256844 ->GetMeanError() );
  valueMeanError_run3.push_back(resp_balancing_eta1319_run_256734_256844 ->GetMeanError() );
  valueMeanError_run3.push_back(resp_balancing_eta1925_run_256734_256844 ->GetMeanError() );
  valueMeanError_run3.push_back(resp_balancing_eta2530_run_256734_256844 ->GetMeanError() );
  valueMeanError_run3.push_back(resp_balancing_eta3032_run_256734_256844 ->GetMeanError() );
  valueMeanError_run3.push_back(resp_balancing_eta3252_run_256734_256844 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMean_run4;           
  valueMean_run4.push_back(resp_balancing_eta0013_run_256866_257600 ->GetMean() );
  valueMean_run4.push_back(resp_balancing_eta1319_run_256866_257600 ->GetMean() );
  valueMean_run4.push_back(resp_balancing_eta1925_run_256866_257600 ->GetMean() );
  valueMean_run4.push_back(resp_balancing_eta2530_run_256866_257600 ->GetMean() );
  valueMean_run4.push_back(resp_balancing_eta3032_run_256866_257600 ->GetMean() );
  valueMean_run4.push_back(resp_balancing_eta3252_run_256866_257600 ->GetMean() );
  std::vector<double > valueMeanError_run4;           
  valueMeanError_run4.push_back(resp_balancing_eta0013_run_256866_257600 ->GetMeanError() );
  valueMeanError_run4.push_back(resp_balancing_eta1319_run_256866_257600 ->GetMeanError() );
  valueMeanError_run4.push_back(resp_balancing_eta1925_run_256866_257600 ->GetMeanError() );
  valueMeanError_run4.push_back(resp_balancing_eta2530_run_256866_257600 ->GetMeanError() );
  valueMeanError_run4.push_back(resp_balancing_eta3032_run_256866_257600 ->GetMeanError() );
  valueMeanError_run4.push_back(resp_balancing_eta3252_run_256866_257600 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////

  //////////////////////////////////// MPF

  std::vector<double > valueMeanMPF_run1;           
  valueMeanMPF_run1.push_back(resp_mpf_eta0013_run_256630_256678 ->GetMean() );
  valueMeanMPF_run1.push_back(resp_mpf_eta1319_run_256630_256678 ->GetMean() );
  valueMeanMPF_run1.push_back(resp_mpf_eta1925_run_256630_256678 ->GetMean() );
  valueMeanMPF_run1.push_back(resp_mpf_eta2530_run_256630_256678 ->GetMean() );
  valueMeanMPF_run1.push_back(resp_mpf_eta3032_run_256630_256678 ->GetMean() );
  valueMeanMPF_run1.push_back(resp_mpf_eta3252_run_256630_256678 ->GetMean() );
  std::vector<double > valueMeanMPFError_run1;           
  valueMeanMPFError_run1.push_back(resp_mpf_eta0013_run_256630_256678 ->GetMeanError() );
  valueMeanMPFError_run1.push_back(resp_mpf_eta1319_run_256630_256678 ->GetMeanError() );
  valueMeanMPFError_run1.push_back(resp_mpf_eta1925_run_256630_256678 ->GetMeanError() );
  valueMeanMPFError_run1.push_back(resp_mpf_eta2530_run_256630_256678 ->GetMeanError() );
  valueMeanMPFError_run1.push_back(resp_mpf_eta3032_run_256630_256678 ->GetMeanError() );
  valueMeanMPFError_run1.push_back(resp_mpf_eta3252_run_256630_256678 ->GetMeanError() );
  ///////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMeanMPF_run2;           
  valueMeanMPF_run2.push_back(resp_mpf_eta0013_run_256728_256730 ->GetMean() );
  valueMeanMPF_run2.push_back(resp_mpf_eta1319_run_256728_256730 ->GetMean() );
  valueMeanMPF_run2.push_back(resp_mpf_eta1925_run_256728_256730 ->GetMean() );
  valueMeanMPF_run2.push_back(resp_mpf_eta2530_run_256728_256730 ->GetMean() );
  valueMeanMPF_run2.push_back(resp_mpf_eta3032_run_256728_256730 ->GetMean() );
  valueMeanMPF_run2.push_back(resp_mpf_eta3252_run_256728_256730 ->GetMean() );
  std::vector<double > valueMeanMPFError_run2;           
  valueMeanMPFError_run2.push_back(resp_mpf_eta0013_run_256728_256730 ->GetMeanError() );
  valueMeanMPFError_run2.push_back(resp_mpf_eta1319_run_256728_256730 ->GetMeanError() );
  valueMeanMPFError_run2.push_back(resp_mpf_eta1925_run_256728_256730 ->GetMeanError() );
  valueMeanMPFError_run2.push_back(resp_mpf_eta2530_run_256728_256730 ->GetMeanError() );
  valueMeanMPFError_run2.push_back(resp_mpf_eta3032_run_256728_256730 ->GetMeanError() );
  valueMeanMPFError_run2.push_back(resp_mpf_eta3252_run_256728_256730 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMeanMPF_run3;           
  valueMeanMPF_run3.push_back(resp_mpf_eta0013_run_256734_256844 ->GetMean() );
  valueMeanMPF_run3.push_back(resp_mpf_eta1319_run_256734_256844 ->GetMean() );
  valueMeanMPF_run3.push_back(resp_mpf_eta1925_run_256734_256844 ->GetMean() );
  valueMeanMPF_run3.push_back(resp_mpf_eta2530_run_256734_256844 ->GetMean() );
  valueMeanMPF_run3.push_back(resp_mpf_eta3032_run_256734_256844 ->GetMean() );
  valueMeanMPF_run3.push_back(resp_mpf_eta3252_run_256734_256844 ->GetMean() );
  std::vector<double > valueMeanMPFError_run3;           
  valueMeanMPFError_run3.push_back(resp_mpf_eta0013_run_256734_256844 ->GetMeanError() );
  valueMeanMPFError_run3.push_back(resp_mpf_eta1319_run_256734_256844 ->GetMeanError() );
  valueMeanMPFError_run3.push_back(resp_mpf_eta1925_run_256734_256844 ->GetMeanError() );
  valueMeanMPFError_run3.push_back(resp_mpf_eta2530_run_256734_256844 ->GetMeanError() );
  valueMeanMPFError_run3.push_back(resp_mpf_eta3032_run_256734_256844 ->GetMeanError() );
  valueMeanMPFError_run3.push_back(resp_mpf_eta3252_run_256734_256844 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////
  std::vector<double > valueMeanMPF_run4;           
  valueMeanMPF_run4.push_back(resp_mpf_eta0013_run_256866_257600 ->GetMean() );
  valueMeanMPF_run4.push_back(resp_mpf_eta1319_run_256866_257600 ->GetMean() );
  valueMeanMPF_run4.push_back(resp_mpf_eta1925_run_256866_257600 ->GetMean() );
  valueMeanMPF_run4.push_back(resp_mpf_eta2530_run_256866_257600 ->GetMean() );
  valueMeanMPF_run4.push_back(resp_mpf_eta3032_run_256866_257600 ->GetMean() );
  valueMeanMPF_run4.push_back(resp_mpf_eta3252_run_256866_257600 ->GetMean() );
  std::vector<double > valueMeanMPFError_run4;           
  valueMeanMPFError_run4.push_back(resp_mpf_eta0013_run_256866_257600 ->GetMeanError() );
  valueMeanMPFError_run4.push_back(resp_mpf_eta1319_run_256866_257600 ->GetMeanError() );
  valueMeanMPFError_run4.push_back(resp_mpf_eta1925_run_256866_257600 ->GetMeanError() );
  valueMeanMPFError_run4.push_back(resp_mpf_eta2530_run_256866_257600 ->GetMeanError() );
  valueMeanMPFError_run4.push_back(resp_mpf_eta3032_run_256866_257600 ->GetMeanError() );
  valueMeanMPFError_run4.push_back(resp_mpf_eta3252_run_256866_257600 ->GetMeanError() );
  ////////////////////////////////////////////////////////////////////////////////////////////


  TGraphErrors* gr_response_vs_run1 = new TGraphErrors(0);                                                                   
  gr_response_vs_run1->SetName("response_vs_run1");  
  TGraphErrors* gr_response_vs_run2 = new TGraphErrors(0);                                                                   
  gr_response_vs_run2->SetName("response_vs_run2");  
  TGraphErrors* gr_response_vs_run3 = new TGraphErrors(0);                                                                   
  gr_response_vs_run3->SetName("response_vs_run3");  
  TGraphErrors* gr_response_vs_run4 = new TGraphErrors(0);                                                                   
  gr_response_vs_run4->SetName("response_vs_run4");  

  TGraphErrors* gr_responseMPF_vs_run1 = new TGraphErrors(0);                                                                   
  gr_responseMPF_vs_run1->SetName("responseMPF_vs_run1");  
  TGraphErrors* gr_responseMPF_vs_run2 = new TGraphErrors(0);                                                                   
  gr_responseMPF_vs_run2->SetName("responseMPF_vs_run2");  
  TGraphErrors* gr_responseMPF_vs_run3 = new TGraphErrors(0);                                                                   
  gr_responseMPF_vs_run3->SetName("responseMPF_vs_run3");  
  TGraphErrors* gr_responseMPF_vs_run4 = new TGraphErrors(0);                                                                   
  gr_responseMPF_vs_run4->SetName("responseMPF_vs_run4");  


  for( int iplot =0 ; iplot<6 ; iplot++){                                                                                                                                     
                                                                                                                                                                              
    std::pair<float, float> currentBin = etaBins[iplot];                                                                                                                      
    float etaMean = (currentBin.first + currentBin.second) / 2.;    

    gr_response_vs_run1-> SetPoint(iplot, etaMean, valueMean_run1.at(iplot));
    gr_response_vs_run1-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanError_run1.at(iplot));

    gr_response_vs_run2-> SetPoint(iplot, etaMean, valueMean_run2.at(iplot));
    gr_response_vs_run2-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanError_run2.at(iplot));

    gr_response_vs_run3-> SetPoint(iplot, etaMean, valueMean_run3.at(iplot));
    gr_response_vs_run3-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanError_run3.at(iplot));

    gr_response_vs_run4-> SetPoint(iplot, etaMean, valueMean_run4.at(iplot));
    gr_response_vs_run4-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanError_run4.at(iplot));

    /////////////////////// MPF
    gr_responseMPF_vs_run1-> SetPoint(iplot, etaMean, valueMeanMPF_run1.at(iplot));
    gr_responseMPF_vs_run1-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanMPFError_run1.at(iplot));

    gr_responseMPF_vs_run2-> SetPoint(iplot, etaMean, valueMeanMPF_run2.at(iplot));
    gr_responseMPF_vs_run2-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanMPFError_run2.at(iplot));

    gr_responseMPF_vs_run3-> SetPoint(iplot, etaMean, valueMeanMPF_run3.at(iplot));
    gr_responseMPF_vs_run3-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanMPFError_run3.at(iplot));

    gr_responseMPF_vs_run4-> SetPoint(iplot, etaMean, valueMeanMPF_run4.at(iplot));
    gr_responseMPF_vs_run4-> SetPointError(iplot, etaBinsError.at(iplot), valueMeanMPFError_run4.at(iplot));

  }

                                                                                                                                                                              
  TLegend* legend = new TLegend(0.15, 0.15, 0.35, 0.35);                                                                                                                      
  legend->SetTextFont(42);                                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                                   
  legend->SetFillColor(kWhite);                                                                                                                                               
  legend->SetFillStyle(0);                                                                                                                                                    
  legend->SetTextSize(0.036);                                                                                                                                                 
  legend->SetHeader("Run: ");                                                                                                                                  
  legend->AddEntry(gr_response_vs_run1, "256630-256677", "L");                                                                                                                                        
  legend->AddEntry(gr_response_vs_run2, "256729", "L");                                                                                                                                        
  legend->AddEntry(gr_response_vs_run3, "256734-256843", "L");                                                                                                                                        
  legend->AddEntry(gr_response_vs_run4, "256866-257599", "L");                                                                                                                                        

  gStyle -> SetOptStat(kFALSE); 

  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  gr_response_vs_run1 -> SetTitle();  
  gr_response_vs_run1 ->GetHistogram()-> SetXTitle("#eta (jet)");  
  gr_response_vs_run1 ->GetHistogram()-> SetYTitle("Balancing Response");  
  gr_response_vs_run1 -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_response_vs_run1 -> GetYaxis()->SetRangeUser(0.3, 1.2);  
  gr_response_vs_run1 -> GetXaxis()->SetRangeUser(0., 5.2);  

  gr_response_vs_run1 -> SetMarkerColor(kRed);  
  gr_response_vs_run1 -> SetLineColor(kRed); 
  gr_response_vs_run2 -> SetMarkerColor(kBlue);  
  gr_response_vs_run2 -> SetLineColor(kBlue);  
  gr_response_vs_run3 -> SetMarkerColor(kGreen);  
  gr_response_vs_run3 -> SetLineColor(kGreen);  
  gr_response_vs_run4 -> SetMarkerColor(kBlack);  
  gr_response_vs_run4 -> SetLineColor(kBlack);  


  gr_response_vs_run1-> Draw("ZAP");
  gr_response_vs_run2-> Draw("PSAME");
  gr_response_vs_run3-> Draw("PSAME");
  gr_response_vs_run4-> Draw("PSAME");

  legend->Draw();

  c2->SaveAs("Balancing_vs_run.png");


  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  gr_responseMPF_vs_run1 -> SetTitle();  
  gr_responseMPF_vs_run1 ->GetHistogram()-> SetXTitle("#eta (jet)");  
  gr_responseMPF_vs_run1 ->GetHistogram()-> SetYTitle("MPF Response");  
  gr_responseMPF_vs_run1 -> GetYaxis()-> SetTitleOffset(1.4);  
  gr_responseMPF_vs_run1 -> GetYaxis()->SetRangeUser(0.3, 1.4);  
  gr_responseMPF_vs_run1 -> GetXaxis()->SetRangeUser(0., 5.2);  

  gr_responseMPF_vs_run1 -> SetMarkerColor(kRed);  
  gr_responseMPF_vs_run1 -> SetLineColor(kRed); 
  gr_responseMPF_vs_run2 -> SetMarkerColor(kBlue);  
  gr_responseMPF_vs_run2 -> SetLineColor(kBlue);  
  gr_responseMPF_vs_run3 -> SetMarkerColor(kGreen);  
  gr_responseMPF_vs_run3 -> SetLineColor(kGreen);  
  gr_responseMPF_vs_run4 -> SetMarkerColor(kBlack);  
  gr_responseMPF_vs_run4 -> SetLineColor(kBlack);  


  gr_responseMPF_vs_run1-> Draw("ZAP");
  gr_responseMPF_vs_run2-> Draw("PSAME");
  gr_responseMPF_vs_run3-> Draw("PSAME");
  gr_responseMPF_vs_run4-> Draw("PSAME");

  legend->Draw();

  c3->SaveAs("MPF_vs_run.png");




  ////////////////////////////////////////////////////
  /*
  resp_balancing_eta0013_run_256630_256678  ->Scale(1.0 / resp_balancing_eta0013_run_256630_256678->Integral()); 
  resp_balancing_eta0013_run_256728_256730  ->Scale(1.0 / resp_balancing_eta0013_run_256728_256730->Integral()); 
  resp_balancing_eta0013_run_256734_256844  ->Scale(1.0 / resp_balancing_eta0013_run_256734_256844->Integral()); 
  resp_balancing_eta0013_run_256866_257600  ->Scale(1.0 / resp_balancing_eta0013_run_256866_257600->Integral()); 
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  resp_balancing_eta0013_run_256630_256678  -> SetLineColor(kRed);  
  resp_balancing_eta0013_run_256728_256730 -> SetLineColor(kBlue);  
  resp_balancing_eta0013_run_256734_256844  -> SetLineColor(kGreen);  
  resp_balancing_eta0013_run_256866_257600  -> SetLineColor(kBlack);  
  resp_balancing_eta0013_run_256630_256678 -> Draw();
  resp_balancing_eta0013_run_256728_256730 -> Draw("same");  
  resp_balancing_eta0013_run_256734_256844  ->  Draw("same");  
  resp_balancing_eta0013_run_256866_257600  -> Draw("same");  

  c1-> SaveAs("Test.png");
  */
  ////////////////////////////////////////////////////////



}//main





