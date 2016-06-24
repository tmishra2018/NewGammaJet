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

  TString dataFileName;
  dataFileName = TString::Format("data_EComposition.root");
  TFile* dataFile = TFile::Open(dataFileName);
  
  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
  }
  
  TString mc1FileName;
  mc1FileName = TString::Format("mc_EComposition.root");
  
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;



  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();
  TString histoName;
  TH1F *h_data;
  TH1F *h_mc;
  
   for (size_t i = 0; i < etaBinningSize; i++) {

     histoName = TString::Format("JetEnergyComposition_%s", etaBinning.getBinName(i).c_str());
     h_data=(TH1F*)dataFile->Get(histoName);
     h_mc=(TH1F*)mcPhotonJetFile->Get(histoName);
     
     TCanvas *c = new TCanvas("c","c",600,600);    
     gPad->SetLogx();
     h_mc->Draw("hist");
     //     h_data -> SetMarkerColor(kBlack);
     h_data->Draw("same");
     c->SaveAs("dataMC_"+histoName+".png");
     c->Destructor();
     
   }// eta bins
  


}


