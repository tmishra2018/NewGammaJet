#ifndef __parserootPileUpfromJson_C__
#define __parserootPileUpfromJson_C__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>

using namespace std;

map<int, map<int, double> > m_PU_root;

double getAvgPUfromlatest(int run, int ls) {

  return m_PU_root[run][ls];
}

//int parsePileUpJSON2(string filename="/cmshome/fpreiato/GammaJet/CMSSW_7_4_14/src/JetMETCorrections/GammaJetFilter/bin/nTrueInteractions_data_latest2015.txt") {
int parserootPileUpfromJson(string filename="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/data/lumtree.root") {

  cout << "Opening " << filename << "...";

  TFile* file = TFile::Open(filename.c_str());
  cout << "ok" << endl;
  
  Int_t run = 0 ;
  Int_t ls   = 0 ;
  Float_t mu = 0 ;
  TTree* tree = (TTree*) file->Get("t");
  tree->SetBranchAddress("run",&run,NULL);
  tree->SetBranchAddress("ls",&ls,NULL);
  tree->SetBranchAddress("mpx",&mu,NULL);
  int Nentries =  tree->GetEntries();
  
  for(int i = 0; i<Nentries ; i++  ){
    tree->GetEntry(i);
    //cout<<"test arguments for run nb : "<<run << " lumi section : "<<ls<<" PU : "<<mu<<endl;
    m_PU_root[run][ls] = mu;
  
  
  }
    
  file->Close();
  
 

  return 0;
}

#endif //__parserootPileUpfromJson_C__
