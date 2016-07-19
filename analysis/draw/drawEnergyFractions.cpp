#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>

#include <TColor.h>
#include <TStyle.h>
#include <THStack.h>
#include <TPad.h>
#include <TAttFill.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"
#include "etaBinning.h"
#include "ptBinning.h"

int main(int argc, char* argv[]) {

  if (argc != 5) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  std::string data_dataset(argv[1]);
  std::string mc_photonjet(argv[2]);
  std::string recoType(argv[3]);
  std::string jetAlgo(argv[4]);
  recoType = (recoType == "pf") ? "PFlow" : "PUPPI";
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
  if (mcPhotonJetFile) {
    std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;
  }
  
  TString outputDir = TString::Format("PhotonJetPlots_%s_vs_%s_%s_EnergyFractions",data_dataset.c_str(), mc_photonjet.c_str(), postFix.c_str());
  mkdir(outputDir, 0777);
  
  TFile * inputfile;
  TString type;  
  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size() + 1;
  
  for(int ii = 0; ii<2; ii++){
    
    inputfile = (ii == 0) ? mcPhotonJetFile : dataFile; 
    type = (ii == 0) ? "mc" : "data"; 
    
    PtBinning ptBinning;
    std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();
    size_t ptBinningSize = ptBinning.size();
    std::pair<float, float> currentBin;
    
    TString histoName;
    TString histoName_tmp;
    TH1F *h_tmp;
    TString stackName;

    Float_t dataResponse = 0.;
    Float_t dataResponseErr = 0.;    
    //  Float_t meanTruncFraction = 0.99;
    //  Float_t rmsTruncFraction = 0.99;  
    //  Float_t dataRMS = 0.;
    //  Float_t dataRMSErr = 0.;
    
    Float_t ptphot_bins[ptBinningSize+1];
    TH1F *hCHFrac[etaBinningSize];
    TH1F *hNHFrac[etaBinningSize];
    TH1F *hElFrac[etaBinningSize];
    TH1F *hPhFrac[etaBinningSize];
    TH1F *hMuFrac[etaBinningSize];
    THStack *hsum[etaBinningSize];
    
    for(size_t j = 0; j<ptBinningSize; j++) {
      currentBin = ptBinning.getBinValue(j);
      ptphot_bins[j] = currentBin.first;
      if(j == ptBinningSize-1) ptphot_bins[j+1] = currentBin.second;
    }
    
    for (size_t i = 0; i < etaBinningSize; i++) {
      for(int frac=0; frac<5; frac++) {
	if(i == etaBinningSize-1){
	  if(frac==0){
	    histoName = TString::Format("ChHadronFraction_eta0013");
	    hCHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
	  } else if(frac==1){
	    histoName = TString::Format("NHadronFraction_eta0013");
	    hNHFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
	  } else if(frac==2) {
	    histoName = TString::Format("CEmFraction_eta0013");
	    hElFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
	  } else if(frac==3) {
	    histoName = TString::Format("NEmFraction_eta0013");
	    hPhFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
	  } else {
	    histoName = TString::Format("MuFraction_eta0013");
	    hMuFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
	  }
	}else{
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
      if(i == etaBinningSize-1){
	stackName = TString::Format("JetEnergyComposition_eta0013");
      }else{
	stackName = TString::Format("JetEnergyComposition_%s", etaBinning.getBinName(i).c_str());
      }
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
	hElFrac[i]->SetMarkerStyle(22);
	hNHFrac[i]->SetMarkerStyle(21);
	hPhFrac[i]->SetMarkerStyle(23);
	hCHFrac[i]->SetMarkerStyle(20);
      }

      hsum[i]->Add(hCHFrac[i]);      
      hsum[i]->Add(hPhFrac[i]);
      hsum[i]->Add(hNHFrac[i]);
      hsum[i]->Add(hElFrac[i]);
      hsum[i]->Add(hMuFrac[i]);     
      
      TCanvas *c = new TCanvas("c","c",600,600);    
      gPad->SetLogx();
      if(type == "mc") {
	hsum[i]->Draw("hist");
      }else{
	hsum[i]->Draw();
      }
      c->SaveAs(outputDir+"/"+type+"_"+stackName+".png");
      c->Destructor();
	 
    }// eta bins
 
    TFile out(outputDir+"/"+type+"_EComposition.root","recreate");
    out.cd();
    for (size_t i = 0; i < etaBinningSize; i++) {
      hCHFrac[i]->Write();
      hNHFrac[i]->Write();
      hElFrac[i]->Write();
      hPhFrac[i]->Write();
      hMuFrac[i]->Write();
      hsum[i]->Write();
    }
    out.Close();
  }

  TString dataName;
  dataName = TString::Format(outputDir+"/data_EComposition.root");
  TFile* data = TFile::Open(dataName);
  if (data) {
    std::cout << "Opened data file '" << dataName << "'." << std::endl;
  }
  
  TString mcName;
  mcName = TString::Format(outputDir+"/mc_EComposition.root");  
  TFile* mcFile = TFile::Open(mcName);
  if(mcFile){
    std::cout << "Opened mc file '" << mcName << "'." << std::endl;
  }

  TString histoName;
  THStack *h_data;
  THStack *h_mc;
  TString histoNameEnFrac;
  TH1F *h_dataEnFrac;
  TH1F *h_mcEnFrac;
  
  for (size_t i = 0; i < etaBinningSize; i++) {

    if(i == etaBinningSize-1){
      histoName = TString::Format("JetEnergyComposition_eta0013");
    }else{
      histoName = TString::Format("JetEnergyComposition_%s", etaBinning.getBinName(i).c_str());
    }
    h_data=(THStack*)data->Get(histoName);
    h_mc=(THStack*)mcFile->Get(histoName);
    
    TCanvas *c = new TCanvas("c","c",800,800);    

    // Data / MC comparison
    TPad* pad_hi = new TPad("pad_hi", "", 0., 0.3333333, 0.99, 0.99);
    pad_hi->SetLogx();
    pad_hi->SetLeftMargin(0.15);
    pad_hi->SetBottomMargin(0.015);
    pad_hi->Draw();
     
    // Data - MC
    TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.3333333);
    pad_lo->SetGridy();
    pad_lo->SetLogx();
    pad_lo->SetLeftMargin(0.15);
    pad_lo->SetTopMargin(1.);
    pad_lo->SetBottomMargin(0.3);
    pad_lo->Draw();
     
    pad_hi -> cd();  
    h_mc->GetYaxis()->SetTitle("PF Energy Fraction");
    h_mc->GetYaxis()-> SetNdivisions(511);
    h_mc-> GetXaxis()->SetLabelColor(kWhite);
    h_mc->Draw("hist");
    h_data->Draw("same");
    pad_lo->cd();
     
    for(int frac=0; frac<5; frac++) {
      if(i == etaBinningSize-1){
	if(frac==0){
	  histoNameEnFrac = TString::Format("ChHadronFraction_eta0013");
	} else if(frac==1){
	  histoNameEnFrac = TString::Format("NHadronFraction_eta0013");
	} else if(frac==2) {
	  histoNameEnFrac = TString::Format("CEmFraction_eta0013");
	} else if(frac==3) {
	  histoNameEnFrac = TString::Format("NEmFraction_eta0013");
	} else {
	  histoNameEnFrac = TString::Format("MuFraction_eta0013");
	}
      }else{
	if(frac==0){
	  histoNameEnFrac = TString::Format("ChHadronFraction_%s", etaBinning.getBinName(i).c_str());
	} else if(frac==1){
	  histoNameEnFrac = TString::Format("NHadronFraction_%s", etaBinning.getBinName(i).c_str());
	} else if(frac==2) {
	  histoNameEnFrac = TString::Format("CEmFraction_%s", etaBinning.getBinName(i).c_str());
	} else if(frac==3) {
	  histoNameEnFrac = TString::Format("NEmFraction_%s", etaBinning.getBinName(i).c_str());
	} else {
	  histoNameEnFrac = TString::Format("MuFraction_%s", etaBinning.getBinName(i).c_str());
	}
      }
	h_dataEnFrac = (TH1F*)data->Get(histoNameEnFrac);
	h_mcEnFrac = (TH1F*)mcFile->Get(histoNameEnFrac);
       
	TH1D *h_diff = (TH1D*)h_dataEnFrac->Clone("h_diff");
	h_diff->Add(h_mcEnFrac, -1);
	h_diff->Scale(100);
	if(frac==0){
	  h_diff->GetXaxis()->SetLabelSize(0.085);
	  h_diff->GetYaxis()->SetLabelSize(0.07);
	  h_diff->GetYaxis()->SetTitleOffset(0.50);
	  h_diff->GetXaxis()->SetTitleSize(0.09);
	  h_diff->GetYaxis()->SetTitleSize(0.07);	 	 
	  h_diff -> SetTitle("");
	  h_diff-> GetXaxis()->SetMoreLogLabels();
	  h_diff-> GetXaxis()->SetNoExponent();
	  h_diff-> GetYaxis()->SetNdivisions(508);
	  h_diff -> SetYTitle("Data - MC (%)");
	  h_diff -> SetXTitle("p_{T} [GeV]");
	  h_diff -> SetStats(kFALSE);
	  h_diff -> SetMarkerColor(2);
	  h_diff -> GetYaxis()->SetRangeUser(-10., 10.);
	  h_diff -> Draw();
	}else{
	  if(frac==1) h_diff->SetMarkerColor(8);
	  if(frac==2) h_diff->SetMarkerColor(4);
	  if(frac==3) h_diff->SetMarkerColor(6);
	  if(frac==4) h_diff->SetMarkerColor(7);
	  h_diff -> Draw("same");
	}
      }
     
      c->SaveAs(outputDir+"/dataMC_"+histoName+".png");
      c->Destructor();
     
  }// eta bins


}


