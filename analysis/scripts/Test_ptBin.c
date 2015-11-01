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

using namespace std;


int main(int argc, char* argv[]) {
  
  TFile file1("TestFinalized.root");    
  TH1D *ptPhot = (TH1D*)file1.Get("analysis/ptPhoton_passedID");
  
  double MeanAll=  ptPhot-> GetMean();
  
  cout<< "Stamp Mean ALL distribution  "<< MeanAll << endl;

  std::vector<std::pair<float, float> > mPtBins;   
  mPtBins.push_back(std::make_pair(40., 60.));                                                                                                                         
  mPtBins.push_back(std::make_pair(60., 85.));                                                                                                                         
  mPtBins.push_back(std::make_pair(85., 100.));                                                                                                                        
  mPtBins.push_back(std::make_pair(100., 130.));                                                                                                                       
  mPtBins.push_back(std::make_pair(130., 175.));                                                                                                                       
  mPtBins.push_back(std::make_pair(175., 250.));                                                                                                                       
  mPtBins.push_back(std::make_pair(250., 300.));                                                                                                                       
  mPtBins.push_back(std::make_pair(300., 400.));                                                                                                                       
  mPtBins.push_back(std::make_pair(400., 500.));                                                                                                                       
  mPtBins.push_back(std::make_pair(500., 1200.)); 
  
  for( int i = 0 ; i< 10 ; i++){ 

    std::pair<float, float> currentBin = mPtBins[i];
 
    ptPhot ->GetXaxis()->SetRangeUser(currentBin.first, currentBin.second);
    double Mean = ptPhot->GetMean();
    cout<< "Bin " << currentBin.first<< "-"<<currentBin.second<<" -> Mean  "<< Mean << endl;

  }


  cout<<"Calcolo a mano"<< endl;

  int bin_500 = ptPhot -> FindBin(500);
  int bin_1100 = ptPhot -> FindBin(1100);
  int start = bin_500;
  int finish = bin_1100+1;

  cout<<"bin(500) "<<bin_500<<endl;
  cout<<"bin(1100) "<<bin_1100<<endl;

  double N_tot = 0 ;
  double num = 0;
  double My_Mean = 0;

  for(int i = start; i<finish; i++){
    double bin = ptPhot-> GetBin(i);
    double N_i = ptPhot-> GetBinContent(i);
    double x_i = ptPhot-> GetBinCenter(i);
    //    cout<<"bin "<<bin<<endl;
    //    cout<<"N_i "<<N_i<<endl;
    //    cout<<"x_i "<<x_i<<endl;
    num += N_i * x_i ;
    N_tot +=N_i;
  }
  
  My_Mean = num / N_tot;

  cout << "My Mean  "<<My_Mean << endl;


}






