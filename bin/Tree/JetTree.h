//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 28 15:43:50 2012 by ROOT version 5.32/00
// from TTree photon/photon tree
// found on file: output_mc.root
//////////////////////////////////////////////////////////

#ifndef BaseTree_h
#define BaseTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "BaseTree.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTree: public BaseTree {
  public :

    // Declaration of leaf types
    Float_t         jet_area;
    Float_t         qg_tag_likelihood;
    Float_t         jet_Energy;
    Float_t         jet_CHEnF;
    Float_t         jet_NHEnF;
    Float_t         jet_CEmEnF;
    Float_t         jet_NEmEnF;
    Float_t         jet_MuEnF;
    Int_t           jet_CHMult;
    Int_t           jet_NHMult;
    Int_t           jet_PhMult;
    Int_t           jet_ElMult;
    Int_t           jet_MuonMult;
    Int_t           jet_ChargedMult;
    Int_t           jet_NeutralMult;

    // List of branches
    TBranch        *b_jet_area;
    TBranch        *b_jet_Energy;
    TBranch        *b_jet_CHEnF;
    TBranch        *b_jet_NHEnF; 
    TBranch        *b_jet_CEmEnF;
    TBranch        *b_jet_NEmEnF;
    TBranch        *b_jet_MuEnF;  
    TBranch        *b_jet_CHMult; 
    TBranch        *b_jet_NHMult;  
    TBranch        *b_jet_PhMult;  
    TBranch        *b_jet_ElMult;  
    TBranch        *b_jet_MuonMult; 
    TBranch        *b_jet_ChargedMult; 
    TBranch        *b_jet_NeutralMult; 

    JetTree();

    virtual void     Init(TTree *tree);
    void             DisableUnrelatedBranches();
};

#endif

JetTree::JetTree():
  BaseTree::BaseTree()
{
}

void JetTree::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (! tree)
    return;

  BaseTree::Init(tree);

  fChain->SetBranchAddress("jet_area", &jet_area, &b_jet_area);
  fChain->SetBranchAddress("qg_tag_likelihood", &qg_tag_likelihood, NULL);
  fChain->SetBranchAddress("jet_CHEnF", &jet_CHEnF, NULL);
  fChain->SetBranchAddress("jet_NHEnF", &jet_NHEnF, NULL);
  fChain->SetBranchAddress("jet_CEmEnF", &jet_CEmEnF, NULL);
  fChain->SetBranchAddress("jet_NEmEnF", &jet_NEmEnF, NULL);
  fChain->SetBranchAddress("jet_MuEnF", &jet_MuEnF, NULL);
  fChain->SetBranchAddress("jet_CHMult", &jet_CHMult, NULL);
  fChain->SetBranchAddress("jet_NHMult", &jet_NHMult, NULL);
  fChain->SetBranchAddress("jet_PhMult", &jet_PhMult, NULL);
  fChain->SetBranchAddress("jet_ElMult", &jet_ElMult, NULL);
  fChain->SetBranchAddress("jet_MuonMult", &jet_MuonMult, NULL);
  fChain->SetBranchAddress("jet_ChargedMult", &jet_ChargedMult, NULL);
  fChain->SetBranchAddress("jet_NeutralMult", &jet_NeutralMult, NULL);

  InitCache();
}

void JetTree::DisableUnrelatedBranches()
{
  fChain->SetBranchStatus("jet_area", 0);
  fChain->SetBranchStatus("qg_tag_likelihood", 0);
  fChain->SetBranchStatus("jet_CHEnF", 0);
  fChain->SetBranchStatus("jet_NHEnF", 0);
  fChain->SetBranchStatus("jet_CEmEnF", 0);
  fChain->SetBranchStatus("jet_NEmEnF", 0);
  fChain->SetBranchStatus("jet_MuEnF", 0);
  fChain->SetBranchStatus("jet_CHMult", 0);
  fChain->SetBranchStatus("jet_NHMult", 0);
  fChain->SetBranchStatus("jet_PhMult", 0);
  fChain->SetBranchStatus("jet_ElMult", 0);
  fChain->SetBranchStatus("jet_MuonMult", 0);
  fChain->SetBranchStatus("jet_ChargedMult", 0);
  fChain->SetBranchStatus("jet_NeutralMult", 0);

}
