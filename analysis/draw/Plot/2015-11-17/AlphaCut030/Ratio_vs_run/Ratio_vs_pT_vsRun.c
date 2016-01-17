#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
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
  
  TFile file1("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run256630-257613_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run1 = (TGraphErrors*)file1.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run1 = (TGraphErrors*)file1.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file2("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run257614-257969_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run2 = (TGraphErrors*)file2.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run2 = (TGraphErrors*)file2.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file3("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run258129-258177_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run3 = (TGraphErrors*)file3.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run3 = (TGraphErrors*)file3.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file4("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run258211-258448_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run4 = (TGraphErrors*)file4.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run4 = (TGraphErrors*)file4.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file5("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run258655-258713_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run5 = (TGraphErrors*)file5.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run5 = (TGraphErrors*)file5.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file6("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run258714-259685_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run6 = (TGraphErrors*)file6.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run6 = (TGraphErrors*)file6.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file7("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run259686-259891_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run7 = (TGraphErrors*)file7.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run7 = (TGraphErrors*)file7.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file8("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run260373-260532_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run8 = (TGraphErrors*)file8.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run8 = (TGraphErrors*)file8.Get("resp_mpf_eta0013_ratio_vs_pt");

  TFile file9("PhotonJetGraphs_SinglePhoton__Run2015D_2015-11-17_alphacut030_run260533-260627_vs_GJet_Pt15To6000_ReReco_Finalized_alphacut030_PFlowAK4_LUMI.root");
  TGraphErrors* Bal_ratio_vs_run9 = (TGraphErrors*)file9.Get("resp_balancing_eta0013_ratio_vs_pt");
  TGraphErrors* MPF_ratio_vs_run9 = (TGraphErrors*)file9.Get("resp_mpf_eta0013_ratio_vs_pt");

  ///////////// grafici dei valori del fit
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
  std::vector<double> Bal_Par1;
  std::vector<double> Bal_Par1_Error;

  std::vector<double> MPF_Par0;
  std::vector<double> MPF_Par0_Error;
  std::vector<double> MPF_Par1;
  std::vector<double> MPF_Par1_Error;
  ///////////////////////////////////// faccio il fit

  bool pTTuned = true;

  TF1* Bal_ratioFit;

  if(pTTuned){
    Bal_ratioFit = new TF1("Bal_ratioFit", "[0]+[1]*log(x/150.)", 0, 700);  
  }else{
    Bal_ratioFit = new TF1("Bal_ratioFit", "[0]+[1]*log(x/200.)", 0, 700);
  } 
  
  TCanvas* c1 = new TCanvas("c1","",800,800);
  //  c1->SetLogx();
  Bal_ratio_vs_run1->SetTitle("");
  Bal_ratio_vs_run1->GetXaxis()->SetTitle("pT(#gamma) [GeV]");
  Bal_ratio_vs_run1->GetYaxis()->SetTitle("Data/MC");
  Bal_ratio_vs_run1->SetMarkerColor(1);
  Bal_ratio_vs_run1->SetLineColor(1);
  Bal_ratio_vs_run1->GetXaxis()->SetRangeUser(0,700);

  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(1);
  TFitResultPtr fitres_run1= Bal_ratio_vs_run1->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 1:"<<endl;
  fitres_run1->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run1->Draw("AP");

  Bal_ratio_vs_run2->SetMarkerColor(2);
  Bal_ratio_vs_run2->SetLineColor(2);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(2);
  TFitResultPtr fitres_run2= Bal_ratio_vs_run2->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 2:"<<endl;
  fitres_run2->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run2->Draw("PSAME");
 
  Bal_ratio_vs_run3->SetMarkerColor(3);
  Bal_ratio_vs_run3->SetLineColor(3);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(3);
  TFitResultPtr fitres_run3= Bal_ratio_vs_run3->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 3:"<<endl;
  fitres_run3->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run3->Draw("PSAME");

  Bal_ratio_vs_run4->SetMarkerColor(4);
  Bal_ratio_vs_run4->SetLineColor(4);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(4);
  TFitResultPtr fitres_run4= Bal_ratio_vs_run4->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 4:"<<endl;
  fitres_run4->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run4->Draw("PSAME");

  Bal_ratio_vs_run5->SetMarkerColor(6);
  Bal_ratio_vs_run5->SetLineColor(6);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(6);
  TFitResultPtr fitres_run5= Bal_ratio_vs_run5->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 5:"<<endl;
  fitres_run5->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run5->Draw("PSAME");

  Bal_ratio_vs_run6->SetMarkerColor(7);
  Bal_ratio_vs_run6->SetLineColor(7);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(7);
  TFitResultPtr fitres_run6= Bal_ratio_vs_run6->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 6:"<<endl;
  fitres_run6->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run6->Draw("PSAME");

  Bal_ratio_vs_run7->SetMarkerColor(8);
  Bal_ratio_vs_run7->SetLineColor(8);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(8);
  TFitResultPtr fitres_run7= Bal_ratio_vs_run7->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 7:"<<endl;
  fitres_run7->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run7->Draw("PSAME");

  Bal_ratio_vs_run8->SetMarkerColor(9);
  Bal_ratio_vs_run8->SetLineColor(9);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(9);
  TFitResultPtr fitres_run8= Bal_ratio_vs_run8->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 8:"<<endl;
  fitres_run8->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run8->Draw("PSAME");

  Bal_ratio_vs_run9->SetMarkerColor(46);
  Bal_ratio_vs_run9->SetLineColor(46);
  Bal_ratioFit->SetParameter(0, 0.);
  Bal_ratioFit->SetParameter(1, 0.);
  Bal_ratioFit->SetLineColor(46);
  TFitResultPtr fitres_run9= Bal_ratio_vs_run9->Fit(Bal_ratioFit, "RQS");
  Bal_Par0.push_back( Bal_ratioFit->GetParameter(0) );
  Bal_Par0_Error.push_back( Bal_ratioFit->GetParError(0) );
  Bal_Par1.push_back( Bal_ratioFit->GetParameter(1) );
  Bal_Par1_Error.push_back( Bal_ratioFit->GetParError(1) );
  cout<<"Run Range 9:"<<endl;
  fitres_run9->PrintCovMatrix(std::cout);
  Bal_ratio_vs_run9->Draw("PSAME");

  if(pTTuned){  
    c1->SaveAs("BALANCING_Fit_pTRef150GeV.png");
  }else{
    c1->SaveAs("BALANCING_Fit_pTRef200GeV.png");
  }
  ////////////////////////////////////////////////////

  TF1* MPF_ratioFit;

  if(pTTuned){
    MPF_ratioFit = new TF1("MPF_ratioFit", "[0]+[1]*log(x/170.)", 0, 700);
  }else{
    MPF_ratioFit = new TF1("MPF_ratioFit", "[0]+[1]*log(x/200.)", 0, 700);
  }
  
  TCanvas* c2 = new TCanvas("c2","",800,800);
  MPF_ratio_vs_run1->SetTitle("");
  MPF_ratio_vs_run1->GetXaxis()->SetTitle("pT(#gamma) [GeV]");
  MPF_ratio_vs_run1->GetYaxis()->SetTitle("Data/MC");
  MPF_ratio_vs_run1->SetMarkerColor(1);
  MPF_ratio_vs_run1->SetLineColor(1);
  MPF_ratio_vs_run1->GetXaxis()->SetRangeUser(0,700);

  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(1);
  TFitResultPtr MPF_fitres_run1= MPF_ratio_vs_run1->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 1:"<<endl;
  MPF_fitres_run1->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run1->Draw("AP");

  MPF_ratio_vs_run2->SetMarkerColor(2);
  MPF_ratio_vs_run2->SetLineColor(2);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(2);
  TFitResultPtr MPF_fitres_run2= MPF_ratio_vs_run2->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 2:"<<endl;
  MPF_fitres_run2->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run2->Draw("PSAME");

  MPF_ratio_vs_run3->SetMarkerColor(3);
  MPF_ratio_vs_run3->SetLineColor(3);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(3);
  TFitResultPtr MPF_fitres_run3= MPF_ratio_vs_run3->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 3:"<<endl;
  MPF_fitres_run3->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run3->Draw("PSAME");

  MPF_ratio_vs_run4->SetMarkerColor(4);
  MPF_ratio_vs_run4->SetLineColor(4);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(4);
  TFitResultPtr MPF_fitres_run4= MPF_ratio_vs_run4->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 4:"<<endl;
  MPF_fitres_run4->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run4->Draw("PSAME");

  MPF_ratio_vs_run5->SetMarkerColor(6);
  MPF_ratio_vs_run5->SetLineColor(6);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(6);
  TFitResultPtr MPF_fitres_run5= MPF_ratio_vs_run5->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 5:"<<endl;
  MPF_fitres_run5->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run5->Draw("PSAME");

  MPF_ratio_vs_run6->SetMarkerColor(7);
  MPF_ratio_vs_run6->SetLineColor(7);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(7);
  TFitResultPtr MPF_fitres_run6= MPF_ratio_vs_run6->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 6:"<<endl;
  MPF_fitres_run6->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run6->Draw("PSAME");

  MPF_ratio_vs_run7->SetMarkerColor(8);
  MPF_ratio_vs_run7->SetLineColor(8);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(8);
  TFitResultPtr MPF_fitres_run7= MPF_ratio_vs_run7->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 7:"<<endl;
  MPF_fitres_run7->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run7->Draw("PSAME");

  MPF_ratio_vs_run8->SetMarkerColor(9);
  MPF_ratio_vs_run8->SetLineColor(9);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(9);
  TFitResultPtr MPF_fitres_run8= MPF_ratio_vs_run8->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 8:"<<endl;
  MPF_fitres_run8->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run8->Draw("PSAME");

  MPF_ratio_vs_run9->SetMarkerColor(46);
  MPF_ratio_vs_run9->SetLineColor(46);
  MPF_ratioFit->SetParameter(0, 0.);
  MPF_ratioFit->SetParameter(1, 0.);
  MPF_ratioFit->SetLineColor(46);
  TFitResultPtr MPF_fitres_run9= MPF_ratio_vs_run9->Fit(MPF_ratioFit, "RQS");
  MPF_Par0.push_back( MPF_ratioFit->GetParameter(0) );
  MPF_Par0_Error.push_back( MPF_ratioFit->GetParError(0) );
  MPF_Par1.push_back( MPF_ratioFit->GetParameter(1) );
  MPF_Par1_Error.push_back( MPF_ratioFit->GetParError(1) );
  cout<<"Run Range 9:"<<endl;
  MPF_fitres_run9->PrintCovMatrix(std::cout);
  MPF_ratio_vs_run9->Draw("PSAME");
 
  if(pTTuned){
    c2->SaveAs("MPF_Fit_pTRef170GeV.png");
  }else{ 
    c2->SaveAs("MPF_Fit_pTRef200GeV.png");
  }
  ////////////////////////////////////////////////////
  int run_index;

  for( int iplot =0 ; iplot<9 ; iplot++){                                                                                                                                     

    run_index = iplot+1;                                                         
                                                                                                                     
    gr_Bal_Par0_vs_run -> SetPoint(iplot, run_index, Bal_Par0.at(iplot));
    gr_Bal_Par0_vs_run -> SetPointError(iplot, 0, Bal_Par0_Error.at(iplot));
    gr_Bal_Par1_vs_run -> SetPoint(iplot, run_index, Bal_Par1.at(iplot));
    gr_Bal_Par1_vs_run -> SetPointError(iplot, 0, Bal_Par1_Error.at(iplot));

    gr_MPF_Par0_vs_run -> SetPoint(iplot, run_index, MPF_Par0.at(iplot));
    gr_MPF_Par0_vs_run -> SetPointError(iplot, 0, MPF_Par0_Error.at(iplot));
    gr_MPF_Par1_vs_run -> SetPoint(iplot, run_index, MPF_Par1.at(iplot));
    gr_MPF_Par1_vs_run -> SetPointError(iplot, 0, MPF_Par1_Error.at(iplot));
  }

  ///////////////////////////////////////////

  TLegend* legend = new TLegend(0.2, 0.9, 0.4, 0.7);                                                                                                                      
  legend->SetTextFont(42);                                                                                                                                                    
  legend->SetBorderSize(0);                                                                                                                                                   
  legend->SetFillColor(kWhite);                                                                                                                                               
  legend->SetFillStyle(0);                                                                                                                                                    
  legend->SetTextSize(0.036);                                                                 
  legend->SetHeader("|#eta| < 1.3");                                                                                                        
  legend->AddEntry(gr_Bal_Par0_vs_run, "Balancing", "PL");                                                                                                                                        
  legend->AddEntry(gr_MPF_Par0_vs_run, "MPF", "PL");     

  gStyle -> SetOptStat(kFALSE); 

  bool pol1 = true;

  TF1* Bal_Par0_Fit;
  
  if(pol1){
    Bal_Par0_Fit = new TF1("Bal_Par0_Fit", "pol1", 0, 10);
  }else{
    Bal_Par0_Fit = new TF1("Bal_Par0_Fit", "[0]", 0, 10);
  }  
  Bal_Par0_Fit->SetParameter(0, 1.);
  if(pol1)  Bal_Par0_Fit->SetParameter(1, 0.);
  gr_Bal_Par0_vs_run->Fit(Bal_Par0_Fit, "RQ");  
  double Bal_Par0_FitValue = Bal_Par0_Fit->GetParameter(0);
  double Bal_Par0_FitError = Bal_Par0_Fit->GetParError(0);
  double Bal_Par0_FitValue_1 = Bal_Par0_Fit->GetParameter(1);
  double Bal_Par0_FitError_1 = Bal_Par0_Fit->GetParError(1);
  double Bal_Par0_ChiSquare = Bal_Par0_Fit->GetChisquare();
  double Bal_Par0_NDF = Bal_Par0_Fit->GetNDF();

  TF1* MPF_Par0_Fit;

  if(pol1){
    MPF_Par0_Fit = new TF1("MPF_Par0_Fit", "pol1", 0, 10);
  }else{
    MPF_Par0_Fit = new TF1("MPF_Par0_Fit", "[0]", 0, 10);
  } 
  MPF_Par0_Fit->SetParameter(0, 1.);
  if(pol1)  MPF_Par0_Fit->SetParameter(1, 0.);
  MPF_Par0_Fit->SetLineColor(kBlue);
  gr_MPF_Par0_vs_run->Fit(MPF_Par0_Fit, "RQ");  
  double MPF_Par0_FitValue = MPF_Par0_Fit->GetParameter(0);
  double MPF_Par0_FitError = MPF_Par0_Fit->GetParError(0);
  double MPF_Par0_FitValue_1 = MPF_Par0_Fit->GetParameter(1);
  double MPF_Par0_FitError_1 = MPF_Par0_Fit->GetParError(1);
  double MPF_Par0_ChiSquare = MPF_Par0_Fit->GetChisquare();
  double MPF_Par0_NDF = MPF_Par0_Fit->GetNDF();
  
  TF1* Bal_Par1_Fit;

  if(pol1){
    Bal_Par1_Fit = new TF1("Bal_Par1_Fit", "pol1", 0, 10);
  }else{
    Bal_Par1_Fit = new TF1("Bal_Par1_Fit", "[0]", 0, 10);
  }
  Bal_Par1_Fit->SetParameter(0, 0.);
  if(pol1)  Bal_Par1_Fit->SetParameter(1, 0.);
  gr_Bal_Par1_vs_run->Fit(Bal_Par1_Fit, "RQ");  
  double Bal_Par1_FitValue = Bal_Par1_Fit->GetParameter(0);
  double Bal_Par1_FitError = Bal_Par1_Fit->GetParError(0);
  double Bal_Par1_FitValue_1 = Bal_Par1_Fit->GetParameter(1);
  double Bal_Par1_FitError_1 = Bal_Par1_Fit->GetParError(1);
  double Bal_Par1_ChiSquare = Bal_Par1_Fit->GetChisquare();
  double Bal_Par1_NDF = Bal_Par1_Fit->GetNDF();

  TF1* MPF_Par1_Fit;

  if(pol1){
    MPF_Par1_Fit = new TF1("MPF_Par1_Fit", "pol1", 0, 10); 
  }else{
    MPF_Par1_Fit = new TF1("MPF_Par1_Fit", "[0]", 0, 10);
  }
  MPF_Par1_Fit->SetParameter(0, 0.);
  if(pol1)  MPF_Par1_Fit->SetParameter(1, 0.);
  MPF_Par1_Fit->SetLineColor(kBlue);
  gr_MPF_Par1_vs_run->Fit(MPF_Par1_Fit, "RQ");  
  double MPF_Par1_FitValue = MPF_Par1_Fit->GetParameter(0);
  double MPF_Par1_FitError = MPF_Par1_Fit->GetParError(0);
  double MPF_Par1_FitValue_1 = MPF_Par1_Fit->GetParameter(1);
  double MPF_Par1_FitError_1 = MPF_Par1_Fit->GetParError(1);
  double MPF_Par1_ChiSquare = MPF_Par1_Fit->GetChisquare();
  double MPF_Par1_NDF = MPF_Par1_Fit->GetNDF();


  TPaveText* fitlabel = new TPaveText(0.60, 0.85, 0.80, 0.7, "brNDC");
  fitlabel->SetTextSize(0.03);
  fitlabel->SetFillColor(0);
  TString Bal_fitLabelText = TString::Format("Par[0]: %.5f #pm %.5f", Bal_Par0_FitValue, Bal_Par0_FitError);
  TString Bal_fitLabelText_1;
  if(pol1) Bal_fitLabelText_1 = TString::Format("Par[1]: %.5f #pm %.5f", Bal_Par0_FitValue_1, Bal_Par0_FitError_1);
  TString Bal_chiLabelText = TString::Format("ChiSquare/NDF: %.1f / %.1f", Bal_Par0_ChiSquare, Bal_Par0_NDF);
  TString MPF_fitLabelText = TString::Format("Par[0]: %.5f #pm %.5f", MPF_Par0_FitValue, MPF_Par0_FitError);
  TString MPF_fitLabelText_1;
  if(pol1) MPF_fitLabelText_1 = TString::Format("Par[1]: %.5f #pm %.5f", MPF_Par0_FitValue_1, MPF_Par0_FitError_1);
  TString MPF_chiLabelText = TString::Format("ChiSquare/NDF: %.1f / %.1f", MPF_Par0_ChiSquare, MPF_Par0_NDF);
  fitlabel->AddText(Bal_fitLabelText);
  if(pol1)  fitlabel->AddText(Bal_fitLabelText_1);
  fitlabel->AddText(Bal_chiLabelText);
  fitlabel->AddText(MPF_fitLabelText);
  if(pol1)  fitlabel->AddText(MPF_fitLabelText_1);
  fitlabel->AddText(MPF_chiLabelText);
  
  TPaveText* fitlabel_Par1 = new TPaveText(0.60, 0.85, 0.80, 0.7, "brNDC");
  fitlabel_Par1->SetTextSize(0.03);
  fitlabel_Par1->SetFillColor(0);
  TString Bal_fitLabelText_Par1 = TString::Format("Par[0] %.5f #pm %.5f", Bal_Par1_FitValue, Bal_Par1_FitError);
  TString Bal_fitLabelText_Par1_1;
  if(pol1) Bal_fitLabelText_Par1_1 = TString::Format("Par[1] %.5f #pm %.5f", Bal_Par1_FitValue_1, Bal_Par1_FitError_1);
  TString Bal_chiLabelText_Par1 = TString::Format("ChiSquare/NDF: %.1f / %.1f", Bal_Par1_ChiSquare, Bal_Par1_NDF);
  TString MPF_fitLabelText_Par1 = TString::Format("Par[0] %.5f #pm %.5f", MPF_Par1_FitValue, MPF_Par1_FitError);
  TString MPF_fitLabelText_Par1_1;
  if(pol1) MPF_fitLabelText_Par1_1 = TString::Format("Par[1] %.5f #pm %.5f", MPF_Par1_FitValue_1, MPF_Par1_FitError_1);
  TString MPF_chiLabelText_Par1 = TString::Format("ChiSquare/NDF: %.1f / %.1f", MPF_Par1_ChiSquare, MPF_Par1_NDF);
  fitlabel_Par1->AddText(Bal_fitLabelText_Par1);
  if(pol1)  fitlabel_Par1->AddText(Bal_fitLabelText_Par1_1);
  fitlabel_Par1->AddText(Bal_chiLabelText_Par1);
  fitlabel_Par1->AddText(MPF_fitLabelText_Par1);
  if(pol1)  fitlabel_Par1->AddText(MPF_fitLabelText_Par1_1);
  fitlabel_Par1->AddText(MPF_chiLabelText_Par1);

  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  gr_Bal_Par0_vs_run -> SetTitle();  
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
 
 if(pTTuned){
   if(pol1){
     c3->SaveAs("Par0_vs_run_pTTuned_pol1.png");
   }else{
     c3->SaveAs("Par0_vs_run_pTTuned.png");
   }
 }else{
   if(pol1){
   c3->SaveAs("Par0_vs_run_pTRef200_pol1.png");
   }else{
     c3->SaveAs("Par0_vs_run_pTRef200.png");
   }
 }
 
 TCanvas *c4 = new TCanvas("c4","c4",800,800);
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

 if(pTTuned){
   if(pol1){
     c4->SaveAs("Par1_vs_run_pTTuned_pol1.png");
   }else{
     c4->SaveAs("Par1_vs_run_pTTuned.png");
   }
 }else{
   if(pol1){
     c4->SaveAs("Par1_vs_run_pTRef200_pol1.png");
   }else{
     c4->SaveAs("Par1_vs_run_pTRef200.png");
   }
 }

}






