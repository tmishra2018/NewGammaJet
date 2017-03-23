#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>

using namespace std;

//void generate_mc_pileup(const std::string& mc) {
int main(int argc, char* argv[]) {
  
  string* mc = new string(argv[1]);

  //USAGE: ./generate_mc_pileup.exe [MC sample]
  std::cout<<mc->c_str()<<std::endl;

  TChain* chain = new TChain("puvariable","puvariable");

  // Load MC files
   std::cout << "Loading files..." << std::endl;
  TString inputFile = TString::Format("%s.list", mc->c_str());
  std::ifstream file(inputFile.Data());

  do {
    std::string f;
    std::getline(file, f);
    if (! file.good())
      break;
    
    chain->AddFile(f.c_str());
  } while (true);
std::cout << "Done." << std::endl;

  int entries = chain->GetEntries();
  std::cout << "Entries: " << entries << std::endl;
  chain->SetBranchStatus("*", 0);
  chain->SetBranchStatus("TrueInteractionall", 1);
  chain->SetBranchStatus("evtWeightTot", 1);
  chain->SetBranchStatus("Generatorweight", 1);

  // Connect branches
  float n_trueInteractions;
  chain->SetBranchAddress("TrueInteractionall", &n_trueInteractions);
 
  float evtWeightTot;
  chain->SetBranchAddress("evtWeightTot", &evtWeightTot);
  double generatorWeight;
  chain->SetBranchAddress("Generatorweight", &generatorWeight);
  if (generatorWeight == 0)
    generatorWeight = 1.;

  float eventWeight = 0 ;  
  
  
  // PU histogram
  TH1F* pu = new TH1F("pileup", "MC Pileup truth", 100, 0, 100);
  
  for (int i = 0; i < entries; i++) {
    chain->GetEntry(i);

    if (i % 50000 == 0) {
      std::cout << "Iteration " << i << " over " << entries << "; " << (float) i / entries * 100 << "%" << std::endl;
    }
    eventWeight=  evtWeightTot * generatorWeight;
    pu->Fill(n_trueInteractions, eventWeight);
  }

  // Normalize histogram
  double scale = 1. / pu->Integral();
  pu->Scale(scale);

  TString filename = TString::Format("computed_mc_%s_pu_truth_100bins.root", mc->c_str());
  TFile * output = TFile::Open(filename.Data(), "recreate");
  pu->Write();
  output->Close();
  delete output;

  delete chain;
}
