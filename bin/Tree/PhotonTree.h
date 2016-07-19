//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 28 15:43:50 2012 by ROOT version 5.32/00
// from TTree photon/photon tree
// found on file: output_mc.root
//////////////////////////////////////////////////////////

#pragma once

#include "BaseTree.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PhotonTree: public BaseTree {
  public :

    // Declaration of leaf types
    Bool_t          has_pixel_seed;
    Float_t         hadTowOverEm;
    Float_t         sigmaIetaIeta;
    Float_t         r9;
    Float_t         rho;
    Bool_t          hasMatchedPromptElectron;
    Float_t         chargedHadronsIsolation;
    Float_t         neutralHadronsIsolation;
    Float_t         photonIsolation; 
    Float_t         trkSumPtHollowConeDR03; 
    Float_t         ecalPFClusterIso; 
    Float_t         hcalPFClusterIso; 
    Float_t         SC_pt; 
    Float_t         SC_eta; 
    Float_t         SC_phi; 
    Float_t         SC_e; 

    // List of branches
    TBranch        *b_jet_area;

    virtual void     Init(TTree *tree);
};

void PhotonTree::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (! tree)
    return;

  BaseTree::Init(tree);

  fChain->SetBranchAddress("has_pixel_seed", &has_pixel_seed, NULL);
  fChain->SetBranchAddress("hadTowOverEm", &hadTowOverEm, NULL);
  fChain->SetBranchAddress("sigmaIetaIeta", &sigmaIetaIeta, NULL);
  fChain->SetBranchAddress("r9", &r9, NULL);
  fChain->SetBranchAddress("rho", &rho, NULL);
  fChain->SetBranchAddress("hasMatchedPromptElectron", &hasMatchedPromptElectron, NULL);
  fChain->SetBranchAddress("chargedHadronsIsolation", &chargedHadronsIsolation, NULL);
  fChain->SetBranchAddress("neutralHadronsIsolation", &neutralHadronsIsolation, NULL);
  fChain->SetBranchAddress("photonIsolation", &photonIsolation, NULL);
  fChain->SetBranchAddress("trkSumPtHollowConeDR03", &trkSumPtHollowConeDR03, NULL);
  fChain->SetBranchAddress("ecalPFClusterIso", &ecalPFClusterIso, NULL);
  fChain->SetBranchAddress("hcalPFClusterIso", &hcalPFClusterIso, NULL);
  fChain->SetBranchAddress("SC_pt", &SC_pt, NULL);
  fChain->SetBranchAddress("SC_eta", &SC_eta, NULL);
  fChain->SetBranchAddress("SC_phi", &SC_phi, NULL);
  fChain->SetBranchAddress("SC_e", &SC_e, NULL);

  BaseTree::InitCache();
}
