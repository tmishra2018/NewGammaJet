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
  //  TGraphErrors *Balancing_a10 = (TGraphErrors*)f1.Get("resp_PtBalchs_a10_eta00_13");
  // TGraphErrors *Balancing_a15 = (TGraphErrors*)f1.Get("resp_PtBalchs_a15_eta00_13");
  // TGraphErrors *Balancing_a20 = (TGraphErrors*)f1.Get("resp_PtBalchs_a20_eta00_13");
  // TGraphErrors *Balancing_a30 = (TGraphErrors*)f1.Get("resp_PtBalchs_a30_eta00_13");
  
  TGraphErrors *Balancing_a10 = (TGraphErrors*)f1.Get("resp_MPFchs_a10_eta00_13");
  TGraphErrors *Balancing_a15 = (TGraphErrors*)f1.Get("resp_MPFchs_a15_eta00_13");
  TGraphErrors *Balancing_a20 = (TGraphErrors*)f1.Get("resp_MPFchs_a20_eta00_13");
  TGraphErrors *Balancing_a30 = (TGraphErrors*)f1.Get("resp_MPFchs_a30_eta00_13");


  //  TLegend *legend=new TLegend(0.9, 0.9, 0.6, 0.75);
  //  legend->AddEntry(Balancing,"Balancing","LP");
  //  legend->AddEntry(MPF,"MPF","LP");

  //  TPaveText* NoExtrap = new TPaveText(0.9, 0.1, 0.6, 0.2, "brNDC");
  //  NoExtrap->SetTextSize(0.08);
  //  NoExtrap->SetFillColor(0);
  //  NoExtrap->AddText("No extrapolation");

  //  TPaveText* SiExtrap = new TPaveText(0.9, 0.1, 0.6, 0.2, "brNDC");
  //  NoExtrap->SetTextSize(0.08);
  //  SiExtrap->SetFillColor(0);
  //  SiExtrap->AddText("With extrapolation");

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
  mPtBins.push_back(std::make_pair(500., 1500.));
  
  Int_t n;

  for(int j = 0; j<4; j++){
    if(j==0)  n=Balancing_a10->GetN(); //get ploted array dimention
    if(j==1)  n=Balancing_a15->GetN(); //get ploted array dimention
    if(j==2)  n=Balancing_a20->GetN(); //get ploted array dimention
    if(j==3)  n=Balancing_a30->GetN(); //get ploted array dimention
    cout<<"Numero di punti = "<< n<<endl;
  }  

  
  Double_t x;
  Double_t y;
  Double_t yErr;

  std::vector<double> pT_a30;
  std::vector<double> ln_pT;
  std::vector<double> R_a10;
  std::vector<double> RErr_a10;
  std::vector<double> R_a15;
  std::vector<double> RErr_a15;
  std::vector<double> R_a20;
  std::vector<double> RErr_a20;
  std::vector<double> R_a30;
  std::vector<double> RErr_a30;
  std::vector<double> kFSR;
  std::vector<double> kFSRErr;
  
  for(int j = 0; j<4; j++){
      std::cout<<std::endl;
      std::cout<<"New graphic:"<<std::endl;
      for(Int_t i=0; i<n; i++) {
	if(j==0){
	  Balancing_a10->GetPoint(i, x, y);
	  yErr = Balancing_a10->GetErrorY(i);
	  R_a10.push_back(y);    
	  RErr_a10.push_back(yErr);    

	}
	if(j==1){
	  Balancing_a15->GetPoint(i, x, y);
	  yErr = Balancing_a15->GetErrorY(i);	  
	  R_a15.push_back(y);    
	  RErr_a15.push_back(yErr);    
	}	
	if(j==2){
	  Balancing_a20->GetPoint(i, x, y);
	  yErr = Balancing_a20->GetErrorY(i);	  
	  R_a20.push_back(y);    
	  RErr_a20.push_back(yErr);    
	}
	if(j==3){
	  Balancing_a30->GetPoint(i, x, y);
	  yErr = Balancing_a30->GetErrorY(i);	  
	  pT_a30.push_back(x);    
	  R_a30.push_back(y);    
	  RErr_a30.push_back(yErr);    
	}
	cout<<i+1<<"th element of X array: "<< x<<endl;
	cout<<i+1<<"th element of Y array: "<< y<<endl;
       	std::cout<< yErr <<std::endl;
      }
  }
  

  std::cout<< R_a10.size() <<std::endl;
  
  for(Int_t i=0; i<n; i++) {

    const std::pair<float, float> bin = mPtBins.at(i);
    std::cout<< bin.first << " "<<bin.second<<std::endl;
    TString plotName = TString::Format("Ratio_vs_alpha_ptPhot_%d_%d.png", (int) bin.first, (int) bin.second);
    
    TGraphErrors* gr = new TGraphErrors(0);    

    std::cout<< "R_a10 = "<< R_a10.at(i)<<std::endl;
    std::cout<< "R_a15 = "<< R_a15.at(i)<<std::endl;
    std::cout<< "R_a20 = "<< R_a20.at(i)<<std::endl;
    std::cout<< "R_a30 = "<< R_a30.at(i)<<std::endl;

    gr->SetPoint(0, 0.10, R_a10.at(i));
    gr->SetPointError(0, 0, RErr_a10.at(i));
    gr->SetPoint(1, 0.15, R_a15.at(i));
    gr->SetPointError(1, 0, RErr_a15.at(i));
    gr->SetPoint(2, 0.20, R_a20.at(i));
    gr->SetPointError(2, 0, RErr_a20.at(i));
    gr->SetPoint(3, 0.30, R_a30.at(i));
    gr->SetPointError(3, 0, RErr_a30.at(i));
    
    TF1* ratioFit = new TF1("ratioFit", "[0]+[1]*x", 0, 0.30);
    ratioFit->SetParameter(0, 1.);
    ratioFit->SetParameter(1, 0.);
    ratioFit->SetLineColor(kBlue);
    ratioFit->SetLineWidth(2);
    gr -> Fit(ratioFit, "RQ");
    //  gr_resp_ratio->Fit(ratioFit, "RQ");
    //  std::cout << "-> ChiSquare: " << constline->GetChisquare() << "   NDF: " << constline->GetNDF() << std::endl;

    double Par0 = ratioFit->GetParameter(0);
    double Par0Err = ratioFit->GetParError(0);
    double Par1 = ratioFit->GetParameter(1);
    double Par1Err = ratioFit->GetParError(1);
    double ChiSquare = ratioFit->GetChisquare();
    double NDF = ratioFit->GetNDF();
    std::cout << "Par 0= "<<Par0<<" #pm "<<Par0Err     <<std::endl;
    std::cout << "Par 1= "<<Par1<<" #pm "<<Par1Err     <<std::endl;

    std::cout << "ChiSquare/NDF " << ChiSquare << "/" << NDF << std::endl;

    double k = (R_a30.at(i) - Par0) / 0.30 ;
    double Err_k = sqrt( (  (RErr_a30.at(i)*RErr_a30.at(i) ) + (Par0Err*Par0Err) )/0.30 );

    std::cout << "kFSR " << k <<std::endl;

    kFSR.push_back(k);
    kFSRErr.push_back(Err_k);

    
    TPaveText* fitText = new TPaveText(0.9, 0.9, 0.6, 0.7, "brNDC");                                                                                           
    //  NoExtrap->SetTextSize(0.08);                                                                                        
    fitText->SetFillColor(0);                                                                                                                                  
    TString fitLabelText_1 = TString::Format("[0]: %.3f #pm %.3f", Par0, Par0Err);
    TString fitLabelText_2 = TString::Format("[1]: %.3f #pm %.3f", Par1, Par1Err);
    TString fitLabelText_3 = TString::Format("ChiSquare/NDF: %.3f / %.3f", ChiSquare, NDF);
    fitText->AddText(fitLabelText_1);
    fitText->AddText(fitLabelText_2);
    fitText->AddText(fitLabelText_3);

    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    gr ->SetMinimum(0.94);                                                                                                                           
    gr -> SetMaximum(1.00);
    gr -> SetMarkerStyle(20);
    gr -> SetMarkerSize(1.5);
    gr -> SetMarkerColor(kRed);  
    gr -> Draw("ZAP");
    fitText->Draw("same");

    c1->SaveAs(plotName);
    c1->Destructor();
  }

    TGraphErrors* gr_kFSR = new TGraphErrors(0);

    for(Int_t i=0; i<mPtBins.size(); i++) {

      std::cout<<"kFSR plot"<<std::endl;
      std::cout<<"x = "<<pT_a30.at(i)<<std::endl;
      std::cout<<"y = "<<kFSR.at(i)<<std::endl;
 
      gr_kFSR->SetPoint(i, pT_a30.at(i), kFSR.at(i));
      gr_kFSR->SetPointError(i, 0, kFSRErr.at(i));
      
    }
   
    TF1* ratioFit = new TF1("ratioFit", "[0]+[1]*log(x)+[2]*log(x)*log(x)", 60, 1500);
    //    TF1* ratioFit = new TF1("ratioFit", "[0]+[1]*x+[2]*x*x", 0, 1500);
    //    ratioFit->SetParameter(0, 1.);
    //    ratioFit->SetParameter(1, 0.);
    ratioFit->SetLineColor(kBlue);
    ratioFit->SetLineWidth(2);
    //  gr_resp_ratio->Fit(ratioFit, "RQ");
    //  std::cout << "-> ChiSquare: " << constline->GetChisquare() << "   NDF: " << constline->GetNDF() << std::endl;
    
    TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
    c1->SetLogx();
    gr_kFSR ->GetXaxis()->SetMoreLogLabels();                                                                                                                 
    gr_kFSR ->GetXaxis()->SetNoExponent(); 
    gr_kFSR ->SetMinimum(-0.15);                                                                                                                           
    gr_kFSR -> SetMaximum(0.08);
    gr_kFSR -> SetMarkerStyle(20);
    gr_kFSR -> SetMarkerSize(1.5);
    gr_kFSR -> SetMarkerColor(kRed);  
    gr_kFSR -> Draw("ZAP");
    gr_kFSR -> Fit(ratioFit, "RQ");
    std::cout << "-> ChiSquare: " << ratioFit->GetChisquare() << "   NDF: " << ratioFit->GetNDF() << std::endl;
    
     c1->SaveAs("kFSR_vs_pt_NoFirstPoint.png");
     //  c1->SaveAs("kFSR_vs_pt.png");
    
    



}//main





