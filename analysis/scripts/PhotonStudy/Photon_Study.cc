#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
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

  std::string data_dataset(argv[1]);
  TString dataFileName;
  dataFileName = TString::Format("%s", data_dataset.c_str() );
  TFile* dataFile = TFile::Open(dataFileName);
  if (dataFile) {
    std::cout << "Opened data file" <<std::endl;
  }

  TTree *PhotonStudy; 
  dataFile->GetObject("gammaJet/photonStudy",PhotonStudy);

  uint64_t totalEvents = PhotonStudy->GetEntries();

  //  std::cout<< totalEvents << std::endl;

   std::vector<double> *NEvents = 0;
   TBranch *bNEvents = 0;
   PhotonStudy->SetBranchAddress("NEvents",&NEvents,&bNEvents);
   std::vector<double> *Area = 0;
   TBranch *bArea = 0;
   PhotonStudy->SetBranchAddress("Area",&Area,&bArea);
   std::vector<double> *Rho = 0;
   TBranch *bRho = 0;
   PhotonStudy->SetBranchAddress("rho",&Rho,&bRho);

   std::vector<double> *NPV = 0;
   TBranch *bNPV = 0;
   PhotonStudy->SetBranchAddress("nPV",&NPV,&bNPV);

   std::vector<double> *NPVGood = 0;
   TBranch *bNPVGood = 0;
   PhotonStudy->SetBranchAddress("nPV_Good",&NPVGood,&bNPVGood);

   std::vector<double> *R_min = 0;
   TBranch *bR_min = 0;
   PhotonStudy->SetBranchAddress("R_min",&R_min,&bR_min);
   std::vector<double> *R_max = 0;
   TBranch *bR_max = 0;
   PhotonStudy->SetBranchAddress("R_max",&R_max,&bR_max);
   std::vector<double> *SumE_pfCand = 0;
   TBranch *bSumE_pfCand = 0;
   PhotonStudy->SetBranchAddress("SumE_pfCandidate",&SumE_pfCand,&bSumE_pfCand);
   std::vector<double> *K_L1FastJet = 0;
   TBranch *bK_L1FastJet = 0;
   PhotonStudy->SetBranchAddress("K_L1FastJet",&K_L1FastJet,&bK_L1FastJet);
   std::vector<double> *K_L1RC = 0;
   TBranch *bK_L1RC = 0;
   PhotonStudy->SetBranchAddress("K_L1RC",&K_L1RC,&bK_L1RC);
   
   TGraphErrors* gr_1 = new TGraphErrors(0);
   gr_1->SetName("SumE_pfCand_vs_dR");
   TGraphErrors* gr_2 = new TGraphErrors(0);
   gr_2->SetName("K_L1FastJet_vs_dR");
   TGraphErrors* gr_3 = new TGraphErrors(0);
   gr_3->SetName("K_L1RC_vs_dR");

   std::cout<< "NEvents | " <<"Area(A) | "<<"rho | "<<"nPV | "<<"nPV(good) | "<<"R_min | "<<"R_max | "<<"SumE(pfCand) | "<<"SumE(pfCand)/N*A | "<<"K(L1FastJet) | "<<"K(L1FastJet)/N*A | "<<"K(L1RC) | "<<"K(L1RC)/N*A"<< std::endl;

   for (Int_t i = 0; i < totalEvents; i++) {

     Long64_t tentry = PhotonStudy->LoadTree(i);
     bNEvents->GetEntry(tentry);
     bArea->GetEntry(tentry);
     bRho->GetEntry(tentry);
     bNPV->GetEntry(tentry);
     bNPVGood->GetEntry(tentry);
     bR_min->GetEntry(tentry);
     bR_max->GetEntry(tentry);
     bSumE_pfCand->GetEntry(tentry);
     bK_L1FastJet->GetEntry(tentry);
     bK_L1RC->GetEntry(tentry);
     
     //     std::cout<< NEvents->size() << std::endl;

     for (UInt_t j = 0; j < NEvents->size(); ++j) {
       
       //       std::cout<< NEvents->at(j) << std::endl;
       //       std::cout<< Area->at(j) << std::endl;
       //       std::cout<< R_min->at(j) << std::endl;
       //       std::cout<< R_max->at(j) << std::endl;
       //       std::cout<< SumE_pfCand->at(j) << std::endl;
       //       std::cout<< K_L1FastJet->at(j) << std::endl;
       //       std::cout<< K_L1RC->at(j) << std::endl;

       double N_ = NEvents->at(j);
       double A_ = Area->at(j);
       double rho_ = Rho->at(j) / N_; 
       double nPV_ = NPV->at(j) / N_;
       double nPVGood_ = NPVGood->at(j) / N_;
       double R_min_ = R_min->at(j);
       double R_max_ = R_max->at(j);
       double SumE_ = SumE_pfCand->at(j);
       double KL1FJ_ = K_L1FastJet->at(j);
       double KL1RC_ = K_L1RC->at(j);

       double x = R_min_ +0.05 ;
       //       std::cout<<x << std::endl;
       double Err_x = 0.05;

       double y1 = SumE_ / (A_ * N_);
       double y2 = KL1FJ_ / (A_ * N_);
       double y3 = KL1RC_ / (A_ * N_);

       //       std::cout<<y1 <<" " <<y2 << " " << y3<< std::endl;

       std::cout<<std::endl;
       std::cout<< N_ <<" | "<<A_ <<" | "<<rho_<<" | "<<nPV_<<" | "<<nPVGood_<< " | "<<R_min_<<" | "<<R_max_<<" | "<<SumE_<<" | "<<y1<<" | "<<KL1FJ_<<" | "<<y2<<" | "<<KL1RC_<<" | "<<y3<<std::endl;
       std::cout<<std::endl;

       double Err_y1 = y1 / sqrt(N_) ;
       double Err_y2 = y2 / sqrt(N_) ;
       double Err_y3 = y3 / sqrt(N_) ;

       gr_1->SetPoint(j, x, y1);
       gr_1->SetPointError(j, Err_x, Err_y1 );
       gr_2->SetPoint(j, x, y2);
       gr_2->SetPointError(j, Err_x, Err_y2);
       gr_3->SetPoint(j, x, y3);
       gr_3->SetPointError(j, Err_x, Err_y3);

     }
   }


   TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
   gr_1->GetXaxis()->SetTitle("R_{cone}");
   gr_1->GetYaxis()->SetTitle("Energy density");
   gr_1->GetYaxis()->SetRangeUser(0, 21);
   gr_1->SetMarkerStyle(20);
   gr_1->SetMarkerSize(1.3);
   gr_1->SetMarkerColor(kBlue);
   gr_1 -> Draw("ZAP");
   gr_2->SetMarkerStyle(20);
   gr_2->SetMarkerSize(1.3);
   gr_2->SetMarkerColor(kRed);
   gr_2 -> Draw("Psame");
   gr_3->SetMarkerStyle(20);
   gr_3->SetMarkerSize(1.3);
   gr_3->SetMarkerColor(kBlack);
   gr_3 -> Draw("Psame");

   TLegend *leg = new TLegend(0.6, 0.1, 0.9, 0.3);
   //   leg->SetHeader("Energy density from:");
   leg -> AddEntry(gr_1,"sum of pfCand","P");
   leg -> AddEntry(gr_2,"L1FastJet","P");
   leg -> AddEntry(gr_3,"L1RC","P");
   leg->Draw();

   c1->SaveAs(dataFileName+".png");

}

