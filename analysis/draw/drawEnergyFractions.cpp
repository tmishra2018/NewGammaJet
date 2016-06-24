#include <stdlib.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"

#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>
#include <TStyle.h>
#include <THStack.h>
#include <TPad.h>
#include <TAttFill.h>

int main(int argc, char* argv[]) {

  if (argc != 7) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')]" << std::endl;
    exit(23);
  }

  std::string data_dataset(argv[1]);
  std::string mc_photonjet(argv[2]);
  std::string mc_QCD(argv[3]);
  std::string recoType(argv[4]);
  std::string jetAlgo(argv[5]);
  std::string norm(argv[6]);
  if (norm != "LUMI" && norm != "SHAPE") {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  std::string postFix = recoType + jetAlgo;
  if(recoType == "PFlow") postFix += "chs";

  TString dataFileName;
  dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  TFile* dataFile = TFile::Open(dataFileName);
  
  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
  }
  
  TString mc1FileName;
  mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;


  TFile * inputfile;
  TString type;
  // TFile * inputfile = dataFile;
  //  TString type("data");

  for(int ii = 0; ii<2; ii++){

    inputfile = (ii == 0) ? mcPhotonJetFile : dataFile; 
    type = (ii == 0) ? "mc" : "data"; 


  PtBinning ptBinning;
  std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();
  size_t ptBinningSize = ptBinning.size();
  std::pair<float, float> currentBin;

  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();
  TString histoName;
  TString histoName_tmp;
  TH1F *h_tmp;
  TString stackName;
  
  //  Float_t meanTruncFraction = 0.99;
  //  Float_t rmsTruncFraction = 0.99;  
  Float_t dataResponse = 0.;
  Float_t dataResponseErr = 0.;
  //  Float_t dataRMS = 0.;
  //  Float_t dataRMSErr = 0.;

  Float_t ptphot_bins[ptBinningSize+1];
  TH1F *hCHFrac[etaBinningSize];
  TH1F *hNHFrac[etaBinningSize];
  TH1F *hElFrac[etaBinningSize];
  TH1F *hPhFrac[etaBinningSize];
  TH1F *hMuFrac[etaBinningSize];
  THStack *hsum[etaBinningSize];

  TH1F *hCHFrac_eta0013;
  TH1F *hNHFrac_eta0013;
  TH1F *hElFrac_eta0013;
  TH1F *hPhFrac_eta0013;
  TH1F *hMuFrac_eta0013;
  THStack *hsum_eta0013;


   for(size_t j = 0; j<ptBinningSize; j++) {
     currentBin = ptBinning.getBinValue(j);
     ptphot_bins[j] = currentBin.first;
     if(j == ptBinningSize-1) ptphot_bins[j+1] = currentBin.second;
   }

   for (size_t i = 0; i < etaBinningSize; i++) {
     for(int frac=0; frac<5; frac++) {
       if(frac==0){
	 histoName = TString::Format("ChHadronFraction_%s", etaBinning.getBinName(i).c_str());
	 hCHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==1){
	 histoName = TString::Format("NHadronFraction_%s", etaBinning.getBinName(i).c_str());
	 hNHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==2) {
	 histoName = TString::Format("CEmFraction_%s", etaBinning.getBinName(i).c_str());
	 hElFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==3) {
	 histoName = TString::Format("NEmFraction_%s", etaBinning.getBinName(i).c_str());
	 hPhFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else {
	 histoName = TString::Format("MuFraction_%s", etaBinning.getBinName(i).c_str());
	 hMuFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       }

       for(size_t j = 0; j<ptBinningSize; j++) {
	 histoName_tmp=histoName;
	 currentBin = ptBinning.getBinValue(j);
	 histoName_tmp.Append(TString::Format("_ptPhot_%i",int(currentBin.first)));
	 histoName_tmp.Append(TString::Format("_%i", int(currentBin.second)));
	 //std::cout << "Getting histo: "<< histoName_tmp << std::endl;
	 h_tmp=(TH1F*)inputfile->Get("analysis/ecomposition/"+histoName_tmp);
	 //     fitTools::getTruncatedMeanAndRMS(h_tmp, dataResponse, dataResponseErr, dataRMS, dataRMSErr, meanTruncFraction, rmsTruncFraction);
	 dataResponse = h_tmp->GetMean();
	 dataResponseErr = h_tmp->GetMeanError();
	 // dataRMS = h_tmp->GetRMS();
	 // dataRMSErr = h_tmp->GetRMSError();

	 if(frac==0) {
	   hCHFrac[i]->SetBinContent(j+1,dataResponse);
	   hCHFrac[i]->SetBinError(j+1,dataResponseErr);
	 } else if(frac==1) {
	   hNHFrac[i]->SetBinContent(j+1,dataResponse);
	   hNHFrac[i]->SetBinError(j+1,dataResponseErr);
	 } else if(frac==2) {
	   hElFrac[i]->SetBinContent(j+1,dataResponse);
	   hElFrac[i]->SetBinError(j+1,dataResponseErr);
	 } else if(frac==3){
	   hPhFrac[i]->SetBinContent(j+1,dataResponse);
	   hPhFrac[i]->SetBinError(j+1,dataResponseErr);
	 } else {
	   hMuFrac[i]->SetBinContent(j+1,dataResponse);
	   hMuFrac[i]->SetBinError(j+1,dataResponseErr);
	 }

       }
     }//loop on frac 

     //stack them
     stackName = TString::Format("JetEnergyComposition_%s", etaBinning.getBinName(i).c_str());
     hsum[i] = new THStack(stackName, stackName);

     if(type == "mc"){
       hMuFrac[i]->SetFillColor(7);
       hMuFrac[i]->SetFillStyle(1001);
       hElFrac[i]->SetFillColor(4);          
       hElFrac[i]->SetFillStyle(1001);
       hNHFrac[i]->SetFillColor(8);          
       hNHFrac[i]->SetFillStyle(1001);
       hPhFrac[i]->SetFillColor(6);     
       hPhFrac[i]->SetFillStyle(1001);  
       hCHFrac[i]->SetFillColor(2);     
       hCHFrac[i]->SetFillStyle(1001); 
     }else{
       hMuFrac[i]->SetMarkerStyle(29);
       //       hMuFrac[i]->SetMarkerColor(7);
       hElFrac[i]->SetMarkerStyle(22);
       //  hElFrac[i]->SetMarkerColor(4);
       hNHFrac[i]->SetMarkerStyle(21);
       // hNHFrac[i]->SetMarkerColor(8);
       hPhFrac[i]->SetMarkerStyle(23);
       // hPhFrac[i]->SetMarkerColor(6);
      hCHFrac[i]->SetMarkerStyle(20);
      // hCHFrac[i]->SetMarkerColor(2);
     }

     hsum[i]->Add(hMuFrac[i]);     
     hsum[i]->Add(hElFrac[i]);
     hsum[i]->Add(hNHFrac[i]);
     hsum[i]->Add(hPhFrac[i]);
     hsum[i]->Add(hCHFrac[i]);

     TCanvas *c = new TCanvas("c","c",600,600);    
     gPad->SetLogx();
     if(type == "mc") {
       hsum[i]->Draw("hist");
     }else{
       hsum[i]->Draw();
     }
     c->SaveAs(type+"_"+stackName+".png");
     c->Destructor();
	 
   }// eta bins
   //special case

     for(int frac=0; frac<5; frac++) {
       if(frac==0){
	 histoName = TString::Format("ChHadronFraction_eta0013");
	 hCHFrac_eta0013 = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==1){
	 histoName = TString::Format("NHadronFraction_eta0013");
	 hNHFrac_eta0013 = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==2) {
	 histoName = TString::Format("CEmFraction_eta0013");
	 hElFrac_eta0013 = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else if(frac==3) {
	 histoName = TString::Format("NEmFraction_eta0013");
	 hPhFrac_eta0013 = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       } else {
	 histoName = TString::Format("MuFraction_eta0013");
	 hMuFrac_eta0013 = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
       }

       for(size_t j = 0; j<ptBinningSize; j++) {
	 histoName_tmp=histoName;
	 currentBin = ptBinning.getBinValue(j);
	 histoName_tmp.Append(TString::Format("_ptPhot_%i",int(currentBin.first)));
	 histoName_tmp.Append(TString::Format("_%i", int(currentBin.second)));
	 h_tmp=(TH1F*)inputfile->Get("analysis/ecomposition/"+histoName_tmp);
	 //     fitTools::getTruncatedMeanAndRMS(h_tmp, dataResponse, dataResponseErr, dataRMS, dataRMSErr, meanTruncFraction, rmsTruncFraction);
	 dataResponse = h_tmp->GetMean();
	 dataResponseErr = h_tmp->GetMeanError();
	 // dataRMS = h_tmp->GetRMS();
	 // dataRMSErr = h_tmp->GetRMSError();

	 if(frac==0) {
	   hCHFrac_eta0013->SetBinContent(j+1,dataResponse);
	   hCHFrac_eta0013->SetBinError(j+1,dataResponseErr);
	 } else if(frac==1) {
	   hNHFrac_eta0013->SetBinContent(j+1,dataResponse);
	   hNHFrac_eta0013->SetBinError(j+1,dataResponseErr);
	 } else if(frac==2) {
	   hElFrac_eta0013->SetBinContent(j+1,dataResponse);
	   hElFrac_eta0013->SetBinError(j+1,dataResponseErr);
	 } else if(frac==3){
	   hPhFrac_eta0013->SetBinContent(j+1,dataResponse);
	   hPhFrac_eta0013->SetBinError(j+1,dataResponseErr);
	 } else {
	   hMuFrac_eta0013->SetBinContent(j+1,dataResponse);
	   hMuFrac_eta0013->SetBinError(j+1,dataResponseErr);
	 }
       }
     }//loop on frac 

     //stack them
     stackName = TString::Format("JetEnergyComposition_eta0013");
     hsum_eta0013 = new THStack(stackName, stackName);

     if(type == "mc"){
       hMuFrac_eta0013->SetFillColor(7);
       hMuFrac_eta0013->SetFillStyle(1001);
       hElFrac_eta0013->SetFillColor(4);          
       hElFrac_eta0013->SetFillStyle(1001);
       hNHFrac_eta0013->SetFillColor(8);          
       hNHFrac_eta0013->SetFillStyle(1001);
       hPhFrac_eta0013->SetFillColor(6);     
       hPhFrac_eta0013->SetFillStyle(1001);  
       hCHFrac_eta0013->SetFillColor(2);     
       hCHFrac_eta0013->SetFillStyle(1001); 
     }else{
       hMuFrac_eta0013->SetMarkerStyle(29);
       //       hMuFrac_eta0013->SetMarkerColor(7);
       hElFrac_eta0013->SetMarkerStyle(22);
       //  hElFrac_eta0013->SetMarkerColor(4);
       hNHFrac_eta0013->SetMarkerStyle(21);
       // hNHFrac_eta0013->SetMarkerColor(8);
       hPhFrac_eta0013->SetMarkerStyle(23);
       // hPhFrac_eta0013->SetMarkerColor(6);
      hCHFrac_eta0013->SetMarkerStyle(20);
      // hCHFrac_eta0013->SetMarkerColor(2);
     }

     hsum_eta0013->Add(hMuFrac_eta0013);     
     hsum_eta0013->Add(hElFrac_eta0013);
     hsum_eta0013->Add(hNHFrac_eta0013);
     hsum_eta0013->Add(hPhFrac_eta0013);
     hsum_eta0013->Add(hCHFrac_eta0013);

     TCanvas *c = new TCanvas("c","c",600,600);    
     gPad->SetLogx();
     if(type == "mc") {
       hsum_eta0013->Draw("hist");
     }else{
       hsum_eta0013->Draw();
     }
     c->SaveAs(type+"_"+stackName+".png");
     c->Destructor();
	 
   TFile out(type+"_EComposition.root","recreate");
   out.cd();
   for (size_t i = 0; i < etaBinningSize; i++) {
     hCHFrac[i]->Write();
     hNHFrac[i]->Write();
     hElFrac[i]->Write();
     hPhFrac[i]->Write();
     hMuFrac[i]->Write();
     hsum[i]->Write();
     hCHFrac_eta0013->Write();
     hNHFrac_eta0013->Write();
     hElFrac_eta0013->Write();
     hPhFrac_eta0013->Write();
     hMuFrac_eta0013->Write();
     hsum_eta0013->Write();
   }
   out.Close();
  }

}


