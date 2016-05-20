//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jul 18 16:06:33 2015 by ROOT version 6.02/05
// from TTree misc/misc tree
// found on file: output_Data.root
//////////////////////////////////////////////////////////

#pragma once

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MiscTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        rho;

   // List of branches
   TBranch        *b_rho;   //!

   MiscTree();
   virtual ~MiscTree();
   virtual Int_t    GetEntry(Long64_t entry);

   virtual void     Init(TTree *tree);

};


MiscTree::MiscTree() : fChain(0) 
{
}

MiscTree::~MiscTree()
{
   if (!fChain) return;
   // federico ->*** Break *** segmentation violation
   //   delete fChain->GetCurrentFile();
}

Int_t MiscTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

void MiscTree::Init(TTree *tree)
{
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("rho", &rho, &b_rho);
 }
