#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTree.h>
#include <TParameter.h>
#include <TH2D.h>

#include <fstream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <chrono>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <DataFormats/Math/interface/deltaPhi.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

#include <PhysicsTools/FWLite/interface/TFileService.h>

#include <FWCore/FWLite/interface/AutoLibraryLoader.h>
#include <FWCore/Framework/interface/Event.h>

#include <DataFormats/Common/interface/Handle.h>
#include <DataFormats/FWLite/interface/Event.h>
#include <DataFormats/FWLite/interface/ChainEvent.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Photon.h>

#include <SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h>

#include "tclap/CmdLine.h"

#include "gammaJetFinalizer.h"
#include "PUReweighter.h"
#include "parsePileUpJSON2.h"

#include <boost/regex.hpp>

#define RESET_COLOR "\033[m"
#define MAKE_RED "\033[31m"
#define MAKE_BLUE "\033[34m"

#define ADD_TREES false

#define DELTAPHI_CUT (2.8)

#define TRIGGER_OK                    0
#define TRIGGER_NOT_FOUND            -1
#define TRIGGER_FOUND_BUT_PT_OUT     -2

boost::shared_ptr<PUReweighter> reweighter;

TFile* PUFile;

bool EXIT = false;

GammaJetFinalizer::GammaJetFinalizer():
  mRandomGenerator(0) {
  mPUWeight = 1.;
  mDoMCComparison = false; // do ptPhoton cut < 200GeV
  mNoPUReweighting = false; // do not PUReweighting
}

GammaJetFinalizer::~GammaJetFinalizer() {

}

std::string GammaJetFinalizer::buildPostfix() {
  std::string algo = mJetAlgo == AK4 ? "AK4" : "AK8";
  //  std::string type = mJetType == PF ? "PFlow" : "Calo";
 std::string type = "PFlow" ;

  std::string postfix = type + algo;

  if (mUseCHS)
    postfix += "chs";

  return postfix;
}

void GammaJetFinalizer::loadFiles(TChain& chain) {
  for (std::vector<std::string>::const_iterator it = mInputFiles.begin(); it != mInputFiles.end(); ++it) {
    chain.Add(it->c_str());
  }
}

void GammaJetFinalizer::cloneTree(TTree* from, TTree*& to) {
  to = from->CloneTree(0);
  from->CopyAddresses(to);
}

void GammaJetFinalizer::runAnalysis() {

  typedef std::chrono::high_resolution_clock clock;


  
  if (mIsMC) {
    // PU Reweighting 
    static std::string cmsswBase = getenv("CMSSW_BASE");
    // PU Reweighting                                                                                                                                                   
    static std::string puPrefix = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting", cmsswBase.c_str()).Data();                          
    static std::string puMC = TString::Format("%s/computed_mc_GJET_Pythia_pu_truth_100bins.root", puPrefix.c_str()).Data(); //GJET Flat -- pythia8                       
    static std::string puData = TString::Format("%s/pu_truth_data2016_100bins.root", puPrefix.c_str()).Data();                                                           
    reweighter = boost::shared_ptr<PUReweighter>(new PUReweighter(puData, puMC));                                                                                           
    // Trigger
    std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_mc.xml", cmsswBase.c_str()).Data();
    std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
    mMCTriggers      = new MCTriggers( TriggerFile.c_str() ) ;
  } else {
    // PU
    parsePileUpJSON2();
    //Trigger
    static std::string cmsswBase = getenv("CMSSW_BASE");
    std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_data.xml", cmsswBase.c_str()).Data();
    std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
    mTriggers      = new Triggers( TriggerFile.c_str() ) ;
  }
  
  // Initialization
  mExtrapBinning.initialize(mPtBinning, (mJetType == PF) ? "PFlow" : "Calo");
  mNewExtrapBinning.initialize(mAlphaCut);

  if (mIsMC) {
    std::cout << "Parsing triggers_mc.xml ..." << std::endl;
    if (! mMCTriggers->parse()) {
      std::cerr << "Failed to parse triggers_mc.xml..." << std::endl;
      return;
    }
    std::cout << "done." << std::endl;
  } else {
    std::cout << "Parsing triggers.xml ..." << std::endl;
    if (! mTriggers->parse()) {
      std::cerr << "Failed to parse triggers.xml..." << std::endl;
      return;
    }
    std::cout << "done." << std::endl;
  }
  std::cout << "triggers mapping:" << std::endl;
  if (mIsMC)
    mMCTriggers->print();
  else
   mTriggers->print();

  std::cout << "Opening files ..." << std::endl;

  const std::string postFix = buildPostfix();

  // Set max TTree size
  TTree::SetMaxTreeSize(429496729600LL);

  TChain analysisChain("gammaJet/analysis");
  TChain photonChain("gammaJet/photon");
  TChain genPhotonChain("gammaJet/photon_gen");
  TChain muonsChain("gammaJet/muons");
  TChain electronsChain("gammaJet/electrons");

  TString treeName = TString::Format("gammaJet/%s/first_jet", postFix.c_str());
  TChain firstJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/first_jet_raw", postFix.c_str());
  TChain firstRawJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/first_jet_gen", postFix.c_str());
  TChain firstGenJetChain(treeName);

  treeName = TString::Format("gammaJet/%s/second_jet", postFix.c_str());
  TChain secondJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/second_jet_raw", postFix.c_str());
  TChain secondRawJetChain(treeName);
  treeName = TString::Format("gammaJet/%s/second_jet_gen", postFix.c_str());
  TChain secondGenJetChain(treeName);

  treeName = TString::Format("gammaJet/%s/met", postFix.c_str());
  TChain metChain(treeName);
  treeName = TString::Format("gammaJet/%s/met_raw", postFix.c_str());
  TChain rawMetChain(treeName);
  treeName = TString::Format("gammaJet/%s/met_gen", postFix.c_str());
  TChain genMetChain(treeName);

  treeName = TString::Format("gammaJet/%s/misc", postFix.c_str());
  TChain miscChain(treeName);

  loadFiles(analysisChain);
  loadFiles(photonChain);
  if (mIsMC)
    loadFiles(genPhotonChain);
  loadFiles(muonsChain);
  loadFiles(electronsChain);

  loadFiles(firstJetChain);
  loadFiles(firstRawJetChain);
  if (mIsMC)
    loadFiles(firstGenJetChain);

  loadFiles(secondJetChain);
  if (mIsMC)
    loadFiles(secondGenJetChain);
  loadFiles(secondRawJetChain);

  loadFiles(metChain);
  if (mIsMC)
    loadFiles(genMetChain);
  loadFiles(rawMetChain);

  loadFiles(miscChain);

  analysis.Init(&analysisChain);
  photon.Init(&photonChain);
  muons.Init(&muonsChain);
  electrons.Init(&electronsChain);

  firstJet.Init(&firstJetChain);
  firstRawJet.Init(&firstRawJetChain);

#if !ADD_TREES
  firstJet.DisableUnrelatedBranches();
  firstRawJet.DisableUnrelatedBranches();
#endif

  secondJet.Init(&secondJetChain);
  secondRawJet.Init(&secondRawJetChain);

#if !ADD_TREES
  secondJet.DisableUnrelatedBranches();
  secondRawJet.DisableUnrelatedBranches();
#endif

  MET.Init(&metChain);
  rawMET.Init(&rawMetChain);

  if (mIsMC) {
    genPhoton.Init(&genPhotonChain);
    genMET.Init(&genMetChain);
    secondGenJet.Init(&secondGenJetChain);
    firstGenJet.Init(&firstGenJetChain);
  }

  misc.Init(&miscChain);

  std::cout << "done." << std::endl;

  std::cout << std::endl << "##########" << std::endl;
  std::cout << "# " << MAKE_BLUE << "Running on " << MAKE_RED << ((mIsMC) ? "MC" : "DATA") << RESET_COLOR << std::endl;
  std::cout << "##########" << std::endl << std::endl;

  // Output file
  std::string outputFile = TString::Format("PhotonJet_%s_%s.root", mDatasetName.c_str(), postFix.c_str()).Data();
  fwlite::TFileService fs(outputFile);

#if ADD_TREES
  TTree* photonTree = NULL;
  cloneTree(photon.fChain, photonTree);

  TTree* genPhotonTree = NULL;
  if (mIsMC)
    cloneTree(genPhoton.fChain, genPhotonTree);

  TTree* firstJetTree = NULL;
  cloneTree(firstJet.fChain, firstJetTree);

  TTree* firstGenJetTree = NULL;
  if (mIsMC)
    cloneTree(firstGenJet.fChain, firstGenJetTree);

  TTree* firstRawJetTree = NULL;
  cloneTree(firstRawJet.fChain, firstRawJetTree);

  TTree* secondJetTree = NULL;
  cloneTree(secondJet.fChain, secondJetTree);

  TTree* secondGenJetTree = NULL;
  if (mIsMC)
    cloneTree(secondGenJet.fChain, secondGenJetTree);

  TTree* secondRawJetTree = NULL;
  cloneTree(secondRawJet.fChain, secondRawJetTree);

  TTree* metTree = NULL;
  cloneTree(MET.fChain, metTree);

  TTree* rawMetTree = NULL;
  cloneTree(rawMET.fChain, rawMetTree);

  TTree* genMetTree = NULL;
  if (mIsMC)
    cloneTree(genMET.fChain, genMetTree);

  TTree* muonsTree = NULL;
  cloneTree(muons.fChain, muonsTree);

  TTree* electronsTree = NULL;
  cloneTree(electrons.fChain, electronsTree);

  TTree* analysisTree = NULL;
  cloneTree(analysis.fChain, analysisTree);
  analysisTree->SetName("misc"); 

  TTree *miscTree = NULL;
  cloneTree(misc.fChain, miscTree);
  miscTree->SetName("rho");
#endif

  std::cout << "Processing..." << std::endl;

  // Automatically call Sumw2 when creating an histogram
  TH1::SetDefaultSumw2(true);

  // Init some analysis variables
  TFileDirectory analysisDir = fs.mkdir("analysis");
  TH1F* h_nvertex = analysisDir.make<TH1F>("nvertex", "nvertex", 51, 0., 50.);
  TH1F* h_nvertex_reweighted = analysisDir.make<TH1F>("nvertex_reweighted", "nvertex_reweighted", 51, 0., 50.);
  TH1F* h_ntrue_interactions = analysisDir.make<TH1F>("ntrue_interactions", "ntrue_interactions", 76, 0., 75.);
  TH1F* h_ntrue_interactions_reweighted = analysisDir.make<TH1F>("ntrue_interactions_reweighted", "ntrue_interactions_reweighted", 76, 0., 75.);

  TH1F* h_mPUWeight = analysisDir.make<TH1F>("mPUWeight", "mPUWeight", 50, 0., 5.);
  //  TH1F* h_analysis_event_weight = analysisDir.make<TH1F>("analysis_event_weight", "analysis_event_weight", 50, 0., 5.);
  TH1F* h_generatorWeight = analysisDir.make<TH1F>("generatorWeight", "generatorWeight", 50, 0., 5.);
  TH1F* h_analysis_evtWeightTot = analysisDir.make<TH1F>("analysis_evtWeightTot", "analysis_evtWeightTot", 50, 0., 10.);
  TH1F* h_event_weight_used = analysisDir.make<TH1F>("event_weight_used", "event_weight_used", 20, 0., 2.);

  TH1F* h_ptPhoton = analysisDir.make<TH1F>("ptPhoton", "ptPhoton", 100, 0., 2000.);
  TH1F* h_EtaPhoton = analysisDir.make<TH1F>("EtaPhoton", "EtaPhoton", 60, -5, 5.);
  TH1F* h_PhiPhoton = analysisDir.make<TH1F>("PhiPhoton", "PhiPhoton", 60, -3.5, 3.5);
  TH1F* h_ptFirstJet = analysisDir.make<TH1F>("ptFirstJet", "ptFirstJet", 200, 0., 2000.);
  TH1F* h_EtaFirstJet = analysisDir.make<TH1F>("EtaFirstJet", "EtaFirstJet", 60, -5, 5.);
  TH1F* h_PhiFirstJet = analysisDir.make<TH1F>("PhiFirstJet", "PhiFirstJet", 60, -3.5, 3.5);
  TH1F* h_ptSecondJet = analysisDir.make<TH1F>("ptSecondJet", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet = analysisDir.make<TH1F>("EtaSecondJet", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet = analysisDir.make<TH1F>("PhiSecondJet", "PhiSecondJet", 60, -3.5, 3.5);
  TH1F* h_MET = analysisDir.make<TH1F>("MET", "MET", 150, 0., 300.);
  TH1F* h_alpha = analysisDir.make<TH1F>("alpha", "alpha", 100, 0., 2.);

  TH1F* h_deltaPhi_NoCut = analysisDir.make<TH1F>("deltaPhi_NoCut", "deltaPhi before cut", 60, M_PI / 2, M_PI);

  TH1F* h_deltaPhi = analysisDir.make<TH1F>("deltaPhi", "deltaPhi", 60, M_PI / 2, M_PI);
  TH1F* h_deltaPhi_2ndJet = analysisDir.make<TH1F>("deltaPhi_2ndjet", "deltaPhi of 2nd jet", 60, M_PI / 2., M_PI);
  TH1F* h_deltaPhi_Photon_MET = analysisDir.make<TH1F>("deltaPhi_Photon_MEt", "deltaPhi MET", 60, M_PI / 2., M_PI);

  std::vector<TH1F*> h_ptPhotonBinned = buildPtVector<TH1F>(analysisDir, "ptPhoton", 100, -1, -1);
  std::vector<TH1F*> h_responseEtaBinned_passedID = buildEtaVector<TH1F>(analysisDir, "responseEtaBinned_passedID", 150, 0, 2);

  TH1F* h_rho = analysisDir.make<TH1F>("rho", "rho", 100, 0, 50);
  TH1F* h_hadTowOverEm = analysisDir.make<TH1F>("hadTowOverEm", "hadTowOverEm", 100, 0, 0.05);
  TH1F* h_sigmaIetaIeta = analysisDir.make<TH1F>("sigmaIetaIeta", "sigmaIetaIeta", 100, 0, 0.011);
  TH1F* h_chargedHadronsIsolation = analysisDir.make<TH1F>("chargedHadronsIsolation", "chargedHadronsIsolation", 100, 0, 2);
  TH1F* h_neutralHadronsIsolation = analysisDir.make<TH1F>("neutralHadronsIsolation", "neutralHadronsIsolation", 100, 0, 100);
  TH1F* h_photonIsolation = analysisDir.make<TH1F>("photonIsolation", "photonIsolation", 100, 0, 15);

  TH1F* h_deltaPhi_passedID = analysisDir.make<TH1F>("deltaPhi_passedID", "deltaPhi", 40, M_PI / 2, M_PI);

  double ptBins[] = {0, 40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 1500};
  int  binnum = sizeof(ptBins)/sizeof(double) -1;
  TH1F* h_ptPhoton_passedID_Binned = new TH1F("ptPhoton_passedID_Binned","ptPhoton", binnum, ptBins);
  TH1F* h_ptPhoton_Binned = new TH1F("ptPhoton_Binned","ptPhoton", binnum, ptBins);

  TH1F* h_ptPhoton_passedID = analysisDir.make<TH1F>("ptPhoton_passedID", "ptPhoton", 100, 0., 2000.);
  TH1F* h_EtaPhoton_passedID = analysisDir.make<TH1F>("EtaPhoton_passedID", "EtaPhoton", 60, -5, 5.);
  TH1F* h_PhiPhoton_passedID = analysisDir.make<TH1F>("PhiPhoton_passedID", "PhiPhoton", 60, -3.5, 3.5);
  TH1F* h_ptFirstJet_passedID = analysisDir.make<TH1F>("ptFirstJet_passedID", "ptFirstJet", 200, 0., 2000.);
  TH1F* h_EtaFirstJet_passedID = analysisDir.make<TH1F>("EtaFirstJet_passedID", "EtaFirstJet", 60, -5, 5.);
  TH1F* h_PhiFirstJet_passedID = analysisDir.make<TH1F>("PhiFirstJet_passedID", "PhiFirstJet", 60, -3.5, 3.5);
  TH1F* h_ptSecondJet_passedID = analysisDir.make<TH1F>("ptSecondJet_passedID", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet_passedID = analysisDir.make<TH1F>("EtaSecondJet_passedID", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet_passedID = analysisDir.make<TH1F>("PhiSecondJet_passedID", "PhiSecondJet", 60, -3.5, 3.5);
  TH2F* h_PtEtaSecondJet_passedID = analysisDir.make<TH2F>("PtEtaSecondJet_passedID", "Pt vs Eta SecondJet", 60, -5, 5, 400, 0, 4000);

  TH1F* h_ptSecondJet_2ndJetOK = analysisDir.make<TH1F>("ptSecondJet_2ndJetOK", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet_2ndJetOK = analysisDir.make<TH1F>("EtaSecondJet_2ndJetOK", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet_2ndJetOK = analysisDir.make<TH1F>("PhiSecondJet_2ndJetOK", "PhiSecondJet", 60, -3.5, 3.5);

  TH1F* h_MET_passedID = analysisDir.make<TH1F>("MET_passedID", "MET", 75, 0., 600.);
  TH1F* h_rawMET_passedID = analysisDir.make<TH1F>("rawMET_passedID", "raw MET", 75, 0., 300.);
  TH1F* h_alpha_passedID = analysisDir.make<TH1F>("alpha_passedID", "alpha", 100, 0., 2.);
  TH1F* h_METResolution_passedID = analysisDir.make<TH1F>("METResolution_passedID", "MET", 100, 0., 600.);
  TH1F* h_MET_perp_passedID = analysisDir.make<TH1F>("MET_perp_passedID", "MET", 200, -600., 600.);
  TH1F* h_MET_par_passedID = analysisDir.make<TH1F>("MET_par_passedID", "MET", 200, -600., 600.);

  TH1F* h_mu = analysisDir.make<TH1F>("mu", "mu", 50, 0, 50);
  TH1F* h_npvGood = analysisDir.make<TH1F>("npvGood", "npvGood", 50, 0, 50);
  TH2F* h_rho_vs_mu = analysisDir.make<TH2F>("rho_vs_mu", "Rho vs mu", 50, 0, 50, 100, 0, 50);
  TH2F* h_npvGood_vs_mu = analysisDir.make<TH2F>("npvGood_vs_mu", "npv_good vs mu", 50, 0, 50, 50, 0, 50);

////resolution plots for mc only
////  if (mIsMC) {
////photon eergy resolution
//  TH1F* h_phPt_resolution = analysisDir.make<TH1F>("phPt_resolution","phPt_resolution", 300, -15, 15);
//  TH1F* h_phPx_resolution = analysisDir.make<TH1F>("phPx_resolution", "phPx_resolution", 300, -15, 15);
//  TH1F* h_phPy_resolution = analysisDir.make<TH1F>("phPy_resolution", "phPy_resolution", 300, -15, 15);
//  TH1F* h_phPt_regression_resolution = analysisDir.make<TH1F>("phPt_regression_resolution", "phPt_regression_resolution", 300, -15, 15);
//  TH1F* h_phPx_regression_resolution = analysisDir.make<TH1F>("phPx_regression_resolution", "phPx_regression_resolution", 300, -15, 15);
//  TH1F* h_phPy_regression_resolution = analysisDir.make<TH1F>("phPy_regression_resolution", "phPy_regression_resolution", 300, -15, 15);
////MET resolution
//  TH1F* h_MET_resolution = analysisDir.make<TH1F>("MET_resolution", "MET_resolution", 100, 0., 60);
//  TH1F* h_MET_par_resolution = analysisDir.make<TH1F>("MET_par_resolution", "MET_par_resolution", 100, -300, 300);
//  TH1F* h_MET_perp_resolution = analysisDir.make<TH1F>("MET_perp_resolution", "MET_perp_resolution", 100, 0, 300);
//  TH1F* h_MET_footprint_resolution = analysisDir.make<TH1F>("MET_footprint_resolution", "MET_footprint_resolution", 100, 0., 600);
//  TH1F* h_MET_par_footprint_resolution = analysisDir.make<TH1F>("MET_par_footprint_resolution", "MET_par_footprint_resolution", 100, -300., 300);
//  TH1F* h_MET_perp_footprint_resolution = analysisDir.make<TH1F>("MET_perp_footprint_resolution", "MET_perp_footprint_resolution", 100, -300, 300);
////  }

//jet composition - viola
  TH1F* h_JetCHEn_passedID = analysisDir.make<TH1F>("CHEnergy_passedID", "CHEnergy", 40, 0., 500);
  TH1F* h_JetNHEn_passedID = analysisDir.make<TH1F>("NHEnergy_passedID", "NHEnergy", 40, 0., 500);
  TH1F* h_JetCEEn_passedID = analysisDir.make<TH1F>("CEEnergy_passedID", "CEEnergy", 40, 0., 500);
  TH1F* h_JetNEEn_passedID = analysisDir.make<TH1F>("NEEnergy_passedID", "NEEnergy", 40, 0., 500);
  TH1F* h_JetPhEn_passedID = analysisDir.make<TH1F>("PhEnergy_passedID", "PhEnergy", 40, 0., 500);
  TH1F* h_JetElEn_passedID = analysisDir.make<TH1F>("ElEnergy_passedID", "ElEnergy", 40, 0., 500);
  TH1F* h_JetMuEn_passedID = analysisDir.make<TH1F>("MuEnergy_passedID", "MuEnergy", 40, 0., 500);
  // federico --- jet multiplicity
  //  TH1F* h_JetPhMult_passedID = analysisDir.make<TH1F>("PhMult_passedID", "PhMult", 50, 0., 50);
  //  TH1F* h_JetNHMult_passedID = analysisDir.make<TH1F>("NHMult_passedID", "NHMult", 50, 0., 50);
  //  TH1F* h_JetElMult_passedID = analysisDir.make<TH1F>("ElMult_passedID", "ElMult", 50, 0., 50);
  //  TH1F* h_JetCHMult_passedID = analysisDir.make<TH1F>("CHMult_passedID", "CHMult", 50, 0., 50);
//jet composition - histos vectors
  TFileDirectory ecompositionDir = analysisDir.mkdir("ecomposition");
  std::vector<std::vector<TH1F*> > ChHadronEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronEnergy", 40, 0., 500.);
  std::vector<std::vector<TH1F*> > NHadronEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronEnergy", 40, 0., 500.);
  std::vector<std::vector<TH1F*> > ElEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "ElEnergy", 40, 0., 500.);
  std::vector<std::vector<TH1F*> > PhEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "PhEnergy", 40, 0., 500.);
  std::vector<std::vector<TH1F*> > MuEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "MuEnergy", 40, 0., 500.);
  std::vector<std::vector<TH1F*> > TotJetEnergy = buildEtaPtVector<TH1F>(ecompositionDir, "TotJetEnergy", 80, 0., 1000.);
//jet composition fractions - histos vectors
  std::vector<std::vector<TH1F*> > ChHadronFraction = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadronFraction = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > ElFraction = buildEtaPtVector<TH1F>(ecompositionDir, "ElFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > PhFraction = buildEtaPtVector<TH1F>(ecompositionDir, "PhFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > MuFraction = buildEtaPtVector<TH1F>(ecompositionDir, "MuFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > LeptFraction = buildEtaPtVector<TH1F>(ecompositionDir, "LeptFraction", 40, 0., 1.);
  //
  std::vector<std::vector<TH1F*> > ChHadronFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronFraction_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadronFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronFraction_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > ElFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "ElFraction_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > PhFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "PhFraction_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > MuFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "MuFraction_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > LeptFraction_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "LeptFraction_mpf", 40, 0., 1.);
//
  std::vector<std::vector<TH1F*> > ChHadronFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadronFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > ElFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "ElFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > PhFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "PhFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > MuFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "MuFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > LeptFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "LeptFractionRaw", 40, 0., 1.);
//
  std::vector<std::vector<TH1F*> > ChHadronFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronFractionRaw_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadronFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronFractionRaw_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > ElFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "ElFractionRaw_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > PhFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "PhFractionRaw_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > MuFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "MuFractionRaw_mpf", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > LeptFractionRaw_mpf = buildEtaPtVector<TH1F>(ecompositionDir, "LeptFractionRaw_mpf", 40, 0., 1.);
//
  std::vector<std::vector<TH1F*> > ChHadron_realFraction = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadron_realFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadron_realFraction = buildEtaPtVector<TH1F>(ecompositionDir, "NHadron_realFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > El_realFraction = buildEtaPtVector<TH1F>(ecompositionDir, "El_realFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > Ph_realFraction = buildEtaPtVector<TH1F>(ecompositionDir, "Ph_realFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > Mu_realFraction = buildEtaPtVector<TH1F>(ecompositionDir, "Mu_realFraction", 40, 0., 1.);
//
  std::vector<std::vector<TH1F*> > ChHadron_realFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadron_realFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadron_realFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "NHadron_realFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > El_realFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "El_realFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > Ph_realFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "Ph_realFractionRaw", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > Mu_realFractionRaw = buildEtaPtVector<TH1F>(ecompositionDir, "Mu_realFractionRaw", 40, 0., 1.);
//jet multiplicities
  std::vector<std::vector<TH1F*> > ChHadronMult = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronMult", 20, 0, 20);
  std::vector<std::vector<TH1F*> > NHadronMult = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronMult", 20, 0, 20);
  std::vector<std::vector<TH1F*> > ElMult = buildEtaPtVector<TH1F>(ecompositionDir, "ElMult", 20, 0, 20);
  std::vector<std::vector<TH1F*> > PhMult = buildEtaPtVector<TH1F>(ecompositionDir, "PhMult", 20, 0, 20);
//check Nvtx vs ptphoton
  std::vector<std::vector<TH1F*> > Nvertices = buildEtaPtVector<TH1F>(ecompositionDir, "Nvertices", 50, 0., 50.);
//
  std::vector<TH1F*> h_ptPhotonBinned_passedID = buildPtVector<TH1F>(analysisDir, "ptPhoton_passedID", 100, -1, -1);
//
  TH1F* h_rho_passedID = analysisDir.make<TH1F>("rho_passedID", "rho", 100, 0, 50);
  TH1F* h_hadTowOverEm_passedID = analysisDir.make<TH1F>("hadTowOverEm_passedID", "hadTowOverEm", 100, 0, 0.05);
  TH1F* h_sigmaIetaIeta_passedID = analysisDir.make<TH1F>("sigmaIetaIeta_passedID", "sigmaIetaIeta", 100, 0, 0.011);
  TH1F* h_chargedHadronsIsolation_passedID = analysisDir.make<TH1F>("chargedHadronsIsolation_passedID", "chargedHadronsIsolation", 100, 0, 2.0);
  TH1F* h_neutralHadronsIsolation_passedID = analysisDir.make<TH1F>("neutralHadronsIsolation_passedID", "neutralHadronsIsolation", 100, 0, 100);
  TH1F* h_photonIsolation_passedID = analysisDir.make<TH1F>("photonIsolation_passedID", "photonIsolation", 100, 0, 15);

  TH2D* h_METvsfirstJet = analysisDir.make<TH2D>("METvsfirstJet", "MET vs firstJet", 150, 0., 300., 150, 0., 500.);
  TH2D* h_firstJetvsSecondJet = analysisDir.make<TH2D>("firstJetvsSecondJet", "firstJet vs secondJet", 60, 5., 100., 60, 5., 100.);

  // Balancing
  TFileDirectory balancingDir = analysisDir.mkdir("balancing");
  std::vector<std::vector<TH1F*> > responseBalancing = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseBalancingRaw = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_raw", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseBalancingGen;
  std::vector<std::vector<TH1F*> > responseBalancingRawGen;
  TH1F* responsePhotGamma_Integrated = analysisDir.make<TH1F>("responsePhotGamma", "responsePhotGamma", 150, 0., 2.);
  TH1F* responseBalancingGen_Integrated = analysisDir.make<TH1F>("responseBalancingGen", "responseBalancingGen", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responsePhotGamma;

  if (mIsMC) {
    responseBalancingGen = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_gen", 150, 0., 2.);
    responseBalancingRawGen = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_raw_gen", 150, 0., 2.);
    responsePhotGamma = buildEtaPtVector<TH1F>(balancingDir, "resp_photGamma", 150, 0., 2.);
  }

  // special case
  std::vector<TH1F*> responseBalancingEta013       = buildPtVector<TH1F>(balancingDir, "resp_balancing", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseBalancingRawEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_raw", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseBalancingGenEta013;
    std::vector<TH1F*> responseBalancingRawGenEta013;
    if (mIsMC) {
      responseBalancingGenEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_gen", "eta0013", 150, 0., 2.);
      responseBalancingRawGenEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_raw_gen", "eta0013", 150, 0., 2.);
    }
    // viola
    //  std::vector<TH1F*> responseBalancingEta024 = buildPtVector<TH1F>(balancingDir, "resp_balancing", "eta024", 150, 0., 2.);

  // MPF
  TFileDirectory mpfDir = analysisDir.mkdir("mpf");
  std::vector<std::vector<TH1F*> > responseMPF = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseMPFRaw = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_raw", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseMPFGen;
  if (mIsMC) {
    responseMPFGen = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_gen", 150, 0., 5.);
  }

  std::vector<TH1F*> responseMPFEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFRawEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_raw", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFGenEta013;
  if (mIsMC) {
    responseMPFGenEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_gen", "eta0013", 150, 0., 2.);
  }
  // viola
  //  std::vector<TH1F*> responseMPFEta024 = buildPtVector<TH1F>(mpfDir, "resp_mpf", "eta024", 150, 0., 2.);

  TFileDirectory trueDir = analysisDir.mkdir("trueresp");
  std::vector<std::vector<TH1F*> > responseTrue;
  std::vector<std::vector<TH1F*> > responsePLI;
  if (mIsMC) {
  responseTrue = buildEtaPtVector<TH1F>(trueDir, "true_resp", 150, 0., 2.);
  responsePLI = buildEtaPtVector<TH1F>(trueDir, "pli", 150, 0., 2.);
}
 
 // vs number of vertices
  TFileDirectory vertexDir = analysisDir.mkdir("vertex");
  std::vector<std::vector<TH1F*>> vertex_responseBalancing = buildEtaVertexVector<TH1F>(vertexDir, "resp_balancing", 150, 0., 2.);
  std::vector<std::vector<TH1F*>> vertex_responseBalancingRaw = buildEtaVertexVector<TH1F>(vertexDir, "resp_balancing_raw", 150, 0., 2.);
  std::vector<TH1F*> vertex_responseBalancingEta013 = buildVertexVector<TH1F>(vertexDir, "resp_balancing", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> vertex_responseBalancingRawEta013 = buildVertexVector<TH1F>(vertexDir, "resp_balancing_raw", "eta0013", 150, 0., 2.);
 
  std::vector<TH1F*> vertex_responseBalancingEta030 = buildVertexVector<TH1F>(vertexDir, "resp_balancing", "eta030", 150, 0., 2.);
  std::vector<TH1F*> vertex_responseBalancingEta050 = buildVertexVector<TH1F>(vertexDir, "resp_balancing", "eta050", 150, 0., 2.);

  std::vector<std::vector<TH1F*>> vertex_responseMPF = buildEtaVertexVector<TH1F>(vertexDir, "resp_mpf", 150, 0., 2.);
  std::vector<std::vector<TH1F*>> vertex_responseMPFRaw = buildEtaVertexVector<TH1F>(vertexDir, "resp_mpf_raw", 150, 0., 2.);
  std::vector<TH1F*> vertex_responseMPFEta013 = buildVertexVector<TH1F>(vertexDir, "resp_mpf", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> vertex_responseMPFRawEta013 = buildVertexVector<TH1F>(vertexDir, "resp_mpf_raw", "eta0013", 150, 0., 2.);
  /*
  // vs run number -- time dependence study
  TFileDirectory runDir;
  std::vector<std::vector<TH1F*>> run_responseBalancing ;
  std::vector<TH1F*> run_responseBalancingEta013 ;
  std::vector<std::vector<TH1F*>> run_responseMPF ;
  std::vector<TH1F*> run_responseMPFEta013 ;

  if(!mIsMC){
  runDir = analysisDir.mkdir("run");
  run_responseBalancing = buildEtaRunVector<TH1F>(runDir, "resp_balancing", 150, 0., 2.);
  run_responseBalancingEta013 = buildRunVector<TH1F>(runDir, "resp_balancing", "eta0013", 150, 0., 2.);
  run_responseMPF = buildEtaRunVector<TH1F>(runDir, "resp_mpf", 150, 0., 2.);
  run_responseMPFEta013 = buildRunVector<TH1F>(runDir, "resp_mpf", "eta0013", 150, 0., 2.);
  }
  */
  // Extrapolation
  int extrapolationBins = 50;
  double extrapolationMin = 0.;
  double extrapolationMax = 2.;
  TFileDirectory extrapDir = analysisDir.mkdir("extrapolation");
  ExtrapolationVectors<TH1F>::type extrap_responseBalancing = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingRaw = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingRawEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);

  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGen;
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingRawGen;
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGenPhot;
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGenGamma;
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingPhotGamma;

  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenEta013;
  std::vector<std::vector<TH1F*> > extrap_responseBalancingRawGenEta013;
  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenPhotEta013;
  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenGammaEta013;
  std::vector<std::vector<TH1F*> > extrap_responseBalancingPhotGammaEta013;

  if (mIsMC) {
    extrap_responseBalancingGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingRawGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_raw_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenPhot = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_phot", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenGamma = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_gamma", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingPhotGamma = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_phot_gamma", extrapolationBins, extrapolationMin, extrapolationMax);

    extrap_responseBalancingGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingRawGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_raw_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenPhotEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_phot", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenGammaEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_gamma", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingPhotGammaEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_phot_gamma", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  }
  ExtrapolationVectors<TH1F>::type extrap_responseMPF = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseMPFRaw = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseMPFEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseMPFRawEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);

  ExtrapolationVectors<TH1F>::type extrap_responseMPFGen;
  std::vector<std::vector<TH1F*> > extrap_responseMPFGenEta013;
  if (mIsMC) {
    extrap_responseMPFGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseMPFGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  }
  
  // New extrapolation
  TFileDirectory newExtrapDir = analysisDir.mkdir("new_extrapolation");
  std::vector<std::shared_ptr<GaussianProfile>> new_extrap_responseBalancing = buildNewExtrapolationEtaVector(newExtrapDir, "extrap_resp_balancing", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::shared_ptr<GaussianProfile>> new_extrap_responseBalancingRaw = buildNewExtrapolationEtaVector(newExtrapDir, "extrap_resp_balancing_raw", extrapolationBins, extrapolationMin, extrapolationMax);  
  std::shared_ptr<GaussianProfile> new_extrap_responseBalancingEta013 = buildNewExtrapolationVector(newExtrapDir, "extrap_resp_balancing", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::shared_ptr<GaussianProfile> new_extrap_responseBalancingRawEta013 = buildNewExtrapolationVector(newExtrapDir, "extrap_resp_balancing_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);

  std::vector<std::shared_ptr<GaussianProfile>> new_extrap_responseMPF = buildNewExtrapolationEtaVector(newExtrapDir, "extrap_resp_mpf", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::shared_ptr<GaussianProfile>> new_extrap_responseMPFRaw = buildNewExtrapolationEtaVector(newExtrapDir, "extrap_resp_mpf_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  std::shared_ptr<GaussianProfile> new_extrap_responseMPFEta013 = buildNewExtrapolationVector(newExtrapDir, "extrap_resp_mpf", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::shared_ptr<GaussianProfile> new_extrap_responseMPFRawEta013 = buildNewExtrapolationVector(newExtrapDir, "extrap_resp_mpf_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);

  // Viola
  //  std::vector<TH1F*> ptFirstJetEta024 = buildPtVector<TH1F>(analysisDir, "ptFirstJet", "eta024", 500, 5., 1005.);

  // Luminosity
  if (! mIsMC) {
    // For data, there's only one file, so open it in order to read the luminosity
    TFile* f = TFile::Open(mInputFiles[0].c_str());
    analysisDir.make<TParameter<double>>("luminosity", static_cast<TParameter<double>*>(f->Get("gammaJet/total_luminosity"))->GetVal());
    f->Close();
    delete f;
  }

  // Store alpha cut
  analysisDir.make<TParameter<double>>("alpha_cut", mAlphaCut);

  uint64_t totalEvents = photonChain.GetEntries();
  uint64_t passedEvents = 0;
  uint64_t passedEventsFromTriggers = 0;
  uint64_t rejectedEventsFromTriggers = 0;
  uint64_t rejectedEventsTriggerNotFound = 0;
  uint64_t rejectedEventsPtOut = 0;

  uint64_t passedPhotonJetCut = 0;
  uint64_t passedDeltaPhiCut = 0;
  uint64_t passedPixelSeedVetoCut = 0;
  uint64_t passedMuonsCut = 0;
  uint64_t passedElectronsCut = 0;
  uint64_t passedJetPtCut = 0;
  uint64_t passedAlphaCut = 0;

  uint64_t from = 0;
  uint64_t to = totalEvents;

  clock::time_point start = clock::now();


  // Loop -- from = 0, to = totalEvents
  for (uint64_t i = from; i < to; i++) {     

    //skip test bug
    //    if( i < 1318285) continue;
    
    if ( (i - from) < 10 || (i - from) % 50000 == 0) { //50000
      clock::time_point end = clock::now();
      double elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      start = end;
      std::cout << "Processing event #" << (i - from + 1) << " of " << (to - from) << " (" << (float) (i - from) / (to - from) * 100 << "%) - " << elapsedTime << " s" << std::endl;
    }
    
    //bug in crab outputs -- skip events with bugs on 50/25 ns
    if( mIsMC ){ // bug in GJET Pythia
      if ( i == 117167 || i == 662324 || i == 1318285) continue;
    }
    // No events skipped for GJet Madgraph
    
    if (EXIT) {
      break;
    }


    analysis.GetEntry(i);
    photon.GetEntry(i);
    if (mIsMC)
      genPhoton.GetEntry(i);
    muons.GetEntry(i);
    electrons.GetEntry(i);

    firstJet.GetEntry(i);
    firstRawJet.GetEntry(i);
    if (mIsMC)
      firstGenJet.GetEntry(i);

    secondJet.GetEntry(i);
    secondRawJet.GetEntry(i);
    if (mIsMC)
      secondGenJet.GetEntry(i);

    MET.GetEntry(i);
    if (mIsMC)
      genMET.GetEntry(i);
    rawMET.GetEntry(i);

    misc.GetEntry(i);

    // if you want analyze a run range
    //    if( analysis.run <260533 || analysis.run >260627 ) continue;

    // skip all events -- usefull to check the crab output
    //   if ( photon.is_present || firstJet.is_present || !photon.is_present || !firstJet.is_present)  continue;
   
    if (! photon.is_present || ! firstJet.is_present)
      continue;

    passedPhotonJetCut++;
    
    //        std::cout<<" passedPhotonJetCut  " << std::endl;
    

    // federico -- trigger throw away photon with pT < 40
    //    if(photon.pt <40) continue;

    int checkTriggerResult = 0;
    std::string passedTrigger;
    float triggerWeight = 1.;
    if(mVerbose) std::cout<<" finding trigger " << std::endl;
    if ((checkTriggerResult = checkTrigger(passedTrigger, triggerWeight)) != TRIGGER_OK) {
      switch (checkTriggerResult) {
      case TRIGGER_NOT_FOUND:
	if (mVerbose) {
	  std::cout << MAKE_RED << "[Run #" << analysis.run << ", pT: " << photon.pt << "] Event does not pass required trigger. List of passed triggers: " << RESET_COLOR << std::endl;
	  size_t size = analysis.trigger_names->size();
	  for (size_t i = 0; i < size; i++) {
	    if (analysis.trigger_results->at(i)) {
	      std::cout << "\t" << analysis.trigger_names->at(i) << std::endl;
	    }
	  }
	}
	rejectedEventsTriggerNotFound++;
	break;
      case TRIGGER_FOUND_BUT_PT_OUT:
	if (mVerbose) {
	  std::cout << MAKE_RED << "[Run #" << analysis.run << ", pT: " << photon.pt << "] Event does pass required trigger, but pT is out of range. List of passed triggers: " << RESET_COLOR << std::endl;
	  size_t size = analysis.trigger_names->size();
	  for (size_t i = 0; i < size; i++) {
	    if (analysis.trigger_results->at(i)) {
	      std::cout << "\t" << analysis.trigger_names->at(i) <<  std::endl;
	    }
	  }
	}
	rejectedEventsPtOut++;
	break;
      }
      rejectedEventsFromTriggers++;
      continue;
    }
    passedEventsFromTriggers++;
    
    if(mVerbose)    std::cout<<" passedEventFromTriggers  " << std::endl;
    
    if (mIsMC) {
      
      // PU reweighting -- official recipe
      computePUWeight(); // 2016 analysis
           
      //federico -- N vertex based
      // computePUWeight_NVtxBased(photon.pt, analysis.nvertex);      
    } else { // IsData
      // Note: if trigger bins == pT bins the prescale is not important
      // to calculate the response (they are normalized to shape)

      //federico -- compute hand-made prescale
      //      bool siprescale = true;
      //      if(siprescale){
      //       	computeTriggerWeight( photon.pt, triggerWeight);
      //      }else{
      //    	triggerWeight = 1;
      //     }
      //  get prescale from xml file
      //  triggerWeight = 1. / triggerWeight; // in the xml file there were the inverse of prescale
      //  get trigger prescale from ntupla
      triggerWeight = triggerWeight;
      //  std::cout<< triggerWeight << std::endl;
    }
    
    // Weights 
    double generatorWeight = (mIsMC) ? analysis.generator_weight : 1.;
    if (generatorWeight == 0.)
      generatorWeight = 1.;
    
    double evtWeightSum = (mIsMC) ? analysis.evtWeightTot : 1.;
    
    //debugging - print weights
    //    if( mIsMC){
    //      std::cout << "mPUWeight     "<< mPUWeight << std::endl; 
    //      std::cout << "generatorWeight   "<< generatorWeight << std::endl; 
    //      std::cout << "analysis.evtWeightTot   "<<evtWeightSum << std::endl; 
    //      std::cout << "analysis.event_weight    " << analysis.event_weight << std::endl;
    //    }
    //  if (! mIsMC)    std::cout<< "triggerWeight    "<< triggerWeight << std::endl;

    // old     
    // double eventWeight = (mIsMC) ? mPUWeight * analysis.event_weight * generatorWeight * analysis.evtWeightTot : triggerWeight;
    // new
    double eventWeight = (mIsMC) ? mPUWeight * generatorWeight * evtWeightSum : triggerWeight;


#if ADD_TREES
    if (mUncutTrees) {
      photonTree->Fill();
      if (mIsMC)
        genPhotonTree->Fill();
      firstJetTree->Fill();
      if (mIsMC)
        firstGenJetTree->Fill();
      firstRawJetTree->Fill();
      secondJetTree->Fill();
      if (mIsMC)
        secondGenJetTree->Fill();
      secondRawJetTree->Fill();
      metTree->Fill();
      rawMetTree->Fill();
      if (mIsMC)
        genMetTree->Fill();
      electronsTree->Fill();
      muonsTree->Fill();
      analysisTree->Fill();
      miscTree->Fill();
    }
#endif

    // Event selection
    // The photon is "Good" from previous step (Filter)
    // From previous step: fabs(deltaPhi(photon, firstJet)) > PI/2
    double deltaPhi = fabs(reco::deltaPhi(photon.phi, firstJet.phi));
    h_deltaPhi_NoCut -> Fill(deltaPhi, eventWeight);

    bool isBack2Back = (deltaPhi >= DELTAPHI_CUT); // 2.8
    if (! isBack2Back) {
      continue;
    }
    
    passedDeltaPhiCut++;
    //    std::cout<<"passedDeltaPhiCut"<<std::endl;
    
    // Pixel seed veto
    if (photon.has_pixel_seed)
      continue;
    
    passedPixelSeedVetoCut++;  
    //    std::cout<<"passedPixelSeedVetoCut"<<std::endl;
    
    // No muons -- changed
    //    if (muons.n != 0)
    //  continue;

    // No Loose muons 
    if (muons.nLooseMuon != 0)
      continue;

    passedMuonsCut++;
    //    std::cout<<"passedMuonsCut"<<std::endl;
    
    // Electron veto. No electron close to the photon
    bool keepEvent = true;
    for (int j = 0; j < electrons.n; j++) {
      double deltaR = fabs(reco::deltaR(photon.eta, photon.phi, electrons.eta[j], electrons.phi[j]));
      if (deltaR < 0.13) {
        keepEvent = false;
        break;
      }
    }
    
    if (! keepEvent)
      continue;

    passedElectronsCut++;    
    //    std::cout<<"passedElectronsCut"<<std::endl;
    
    if (firstJet.pt < 15)
      continue;

    passedJetPtCut++;
    //    std::cout<<"passedJetPtCut"<<std::endl;
    
    //bool secondJetOK = !secondJet.is_present || (secondJet.pt < mAlphaCut * photon.pt);    
    bool secondJetOK = !secondJet.is_present || (secondJet.pt < 10 || secondJet.pt < mAlphaCut * photon.pt);
    
    if (mDoMCComparison) { 
      // Lowest unprescaled trigger @ this pT
      if (photon.pt < 165.)
	continue;
    }
    
    if (secondJetOK)    
      passedAlphaCut++;

    //federico -- without this cut the extrapolation is always the same to different alpha cut
    //    if(!secondJetOK) continue;
    //      passedAlphaCut++;
   
    //    std::cout << " secondJetOK "<< std::endl; 

    double mu;
    if(mIsMC){
      mu = analysis.ntrue_interactions ;
    } else {
      mu = getAvgPU( analysis.run, analysis.lumi_block );
      //      std::cout<< analysis.run << "  " << analysis.lumi_block << "  " << mu<<endl;
    }
    
 
    h_mPUWeight                   ->Fill(mPUWeight);
    //   h_analysis_event_weight  ->Fill(analysis.event_weight);
    h_generatorWeight           ->Fill(generatorWeight);
    h_analysis_evtWeightTot  ->Fill(evtWeightSum);
    h_event_weight_used       ->Fill(eventWeight);
    
    h_nvertex->Fill(analysis.nvertex);
    h_ntrue_interactions->Fill(mu);
    h_nvertex_reweighted->Fill(analysis.nvertex, eventWeight);
    h_ntrue_interactions_reweighted->Fill(mu, eventWeight);
    
    h_ptPhoton               ->Fill(photon.pt, eventWeight);
    h_ptPhoton_Binned        ->Fill(photon.pt, eventWeight);
    h_EtaPhoton             ->Fill(photon.eta, eventWeight);
    h_PhiPhoton             ->Fill(photon.phi, eventWeight);
    h_ptFirstJet              ->Fill(firstJet.pt, eventWeight);
    h_EtaFirstJet            ->Fill(firstJet.eta, eventWeight);
    h_PhiFirstJet            ->Fill(firstJet.phi, eventWeight);
    
    h_deltaPhi                 ->Fill(deltaPhi, eventWeight); //first jet - photon
    
    h_ptSecondJet          ->Fill(secondJet.pt, eventWeight);
    h_EtaSecondJet        ->Fill(secondJet.eta, eventWeight);
    h_PhiSecondJet        ->Fill(secondJet.phi, eventWeight);
    double deltaPhi_2ndJet = fabs(reco::deltaPhi(secondJet.phi, photon.phi));
    h_deltaPhi_2ndJet     ->Fill(deltaPhi_2ndJet, eventWeight); //2nd jet - photon
    
    h_alpha                     ->Fill(secondJet.pt / photon.pt, eventWeight);
    h_MET                      ->Fill(MET.pt, eventWeight);
    h_rho                          ->Fill(photon.rho, eventWeight);
    h_hadTowOverEm      ->Fill(photon.hadTowOverEm, eventWeight);
    h_sigmaIetaIeta          ->Fill(photon.sigmaIetaIeta, eventWeight);
    h_chargedHadronsIsolation    ->Fill(photon.chargedHadronsIsolation, eventWeight);
    h_neutralHadronsIsolation     ->Fill(photon.neutralHadronsIsolation, eventWeight);
    h_photonIsolation                  ->Fill(photon.photonIsolation, eventWeight);
    
    // Dump to Tree
    /*photonToTree(photon);
      firstJetToTree(firstJet);
      if (secondJet.is_present) {
      secondJetToTree(*secondJet);
      }*/

    // Compute values
    // MPF
    deltaPhi_Photon_MET = reco::deltaPhi(photon.phi, MET.phi);
    respMPF = 1. + MET.et * photon.pt * cos(deltaPhi_Photon_MET) / (photon.pt * photon.pt);

    h_deltaPhi_Photon_MET     ->Fill(deltaPhi_Photon_MET, eventWeight); // MET - photon

    //    std::cout << " deltaPhi PhotonMET "<< deltaPhi_Photon_MET << std::endl; 
    //    std::cout << " respMPF "<< respMPF << std::endl; 

    if ( mIsMC){
    deltaPhi_Photon_MET_gen = reco::deltaPhi(genPhoton.phi, genMET.phi);
    respMPFGen = 1. + genMET.et * genPhoton.pt * cos(deltaPhi_Photon_MET_gen) / (genPhoton.pt * genPhoton.pt);
    //    std::cout << " deltaPhi PhotonMET   GEN "<< deltaPhi_Photon_MET_gen << std::endl; 
    //    std::cout << " respMPFGEN "<< respMPFGen << std::endl; 
    }

    deltaPhi_Photon_MET_raw = reco::deltaPhi(photon.phi, rawMET.phi);
    respMPFRaw = 1. + rawMET.et * photon.pt * cos(deltaPhi_Photon_MET_raw) / (photon.pt * photon.pt);
    //    std::cout << " deltaPhi PhotonMET  RAW "<< deltaPhi_Photon_MET_raw << std::endl; 
    //    std::cout << " respMPFRaw "<< respMPFRaw << std::endl; 

    // Balancing
    respBalancing = firstJet.pt / photon.pt;
    if ( mIsMC)    respBalancingGen = firstJet.pt / firstGenJet.pt;
    respBalancingRaw = firstRawJet.pt / photon.pt;
    if ( mIsMC)    respBalancingRawGen = firstRawJet.pt / firstGenJet.pt;

    // For DATA/MC comparison
    if ( mIsMC)    respGenPhoton = firstGenJet.pt / photon.pt;
    if ( mIsMC)    respGenGamma = firstGenJet.pt / genPhoton.pt;
    if ( mIsMC)    respPhotonGamma = photon.pt / genPhoton.pt;


    int ptBin = mPtBinning.getPtBin(photon.pt);
    if (ptBin < 0) {
      //std::cout << "Photon pt " << photon.pt() << " is not covered by our pt binning. Dumping event." << std::endl;
      continue;
    }

    h_ptPhotonBinned[ptBin]->Fill(photon.pt, eventWeight);

    if ( mIsMC)   ptBinGen = mPtBinning.getPtBin(genPhoton.pt);
    int etaBin = mEtaBinning.getBin(firstJet.eta);
    if ( mIsMC)   etaBinGen = mEtaBinning.getBin(firstGenJet.eta);
    int vertexBin = mVertexBinning.getVertexBin(analysis.nvertex);

    //    int runBin = -1;
    //    if(!mIsMC)  runBin = mRunBinning.getRunBin(analysis.run);

    float jetcalcen=0;
    float jetcalcenraw=0;


    if (secondJet.is_present) { //extrapolation is second jet is present
      do {

	//	if( fabs(firstJet.eta) < 0.783){
	//	  Counter++;
	//	  std::cout<< " Extrap counter "<<Counter << std::endl;
	//	}
	
        int extrapBin = mExtrapBinning.getBin(photon.pt, secondJet.pt, ptBin);
        int rawExtrapBin = mExtrapBinning.getBin(photon.pt, secondRawJet.pt, ptBin);

	//     float alpha = secondJet.pt / photon.pt;
	//	const std::pair<float, float> ExtrapBin = mExtrapBinning.getBinValue(extrapBin);
	//	std::cout<< std::endl;
	//	std::cout<< "alpha  " << alpha << std::endl;
	//	std::cout<< "pTPhot  " << photon.pt << std::endl;
	//	std::cout<< "extrapBin.getBinValue  " << ExtrapBin.first << " "<<ExtrapBin.second << std::endl;

        r_RecoPhot = firstJet.pt / photon.pt;
	if ( mIsMC)       r_RecoGen  = firstJet.pt / firstGenJet.pt;
	if ( mIsMC)       r_GenPhot  = firstGenJet.pt / photon.pt;
	if ( mIsMC)       r_GenGamma  = firstGenJet.pt / genPhoton.pt;
	if ( mIsMC)       r_PhotGamma  = photon.pt / genPhoton.pt;

 
	do {
          if (extrapBin < 0) {
            //std::cout << "No bin found for extrapolation: " << secondJet.pt / photon.pt << std::endl;
            break;
          }

          // Special case 
	  if (fabs(firstJet.eta) < 1.305) {
            extrap_responseBalancingEta013[ptBin][extrapBin]->Fill(r_RecoPhot, eventWeight);
            extrap_responseMPFEta013[ptBin][extrapBin]->Fill(respMPF, eventWeight);

            if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {
              extrap_responseBalancingGenEta013[ptBinGen][extrapBin]->Fill(r_RecoGen, eventWeight);
              extrap_responseBalancingGenPhotEta013[ptBinGen][extrapBin]->Fill(r_GenPhot, eventWeight);
              extrap_responseBalancingGenGammaEta013[ptBinGen][extrapBin]->Fill(r_GenGamma, eventWeight);
              extrap_responseBalancingPhotGammaEta013[ptBinGen][extrapBin]->Fill(r_PhotGamma, eventWeight);
              extrap_responseMPFGenEta013[ptBinGen][extrapBin]->Fill(respMPFGen, eventWeight);
            }
          }
	  
          if (etaBin < 0)
            break;

          extrap_responseBalancing[etaBin][ptBin][extrapBin]->Fill(r_RecoPhot, eventWeight);
          extrap_responseMPF[etaBin][ptBin][extrapBin]->Fill(respMPF, eventWeight);

          if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {
            extrap_responseBalancingGen[etaBinGen][ptBinGen][extrapBin]->Fill(r_RecoGen, eventWeight);
            extrap_responseBalancingGenPhot[etaBinGen][ptBinGen][extrapBin]->Fill(r_GenPhot, eventWeight);
            extrap_responseBalancingGenGamma[etaBinGen][ptBinGen][extrapBin]->Fill(r_GenGamma, eventWeight);
            extrap_responseBalancingPhotGamma[etaBinGen][ptBinGen][extrapBin]->Fill(r_PhotGamma, eventWeight);
            extrap_responseMPFGen[etaBinGen][ptBinGen][extrapBin]->Fill(respMPFGen, eventWeight);
          }
        } while (false);

        do {
          if (rawExtrapBin < 0) {
            //std::cout << "No bin found for extrapolation: " << secondJet.pt / photon.pt << std::endl;
            break;
          }

          float r_RecoPhotRaw = firstRawJet.pt / photon.pt;
          float r_RecoGenRaw  = firstRawJet.pt / firstGenJet.pt;

          // Special case	  
          if (fabs(firstJet.eta) < 1.305) {

            extrap_responseBalancingRawEta013[ptBin][rawExtrapBin]->Fill(r_RecoPhotRaw, eventWeight);
            extrap_responseMPFRawEta013[ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);
	    
            if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {
              extrap_responseBalancingRawGenEta013[ptBinGen][rawExtrapBin]->Fill(r_RecoGenRaw, eventWeight);
            }
	  }

          if (etaBin < 0)
            break;
	  
          extrap_responseBalancingRaw[etaBin][ptBin][rawExtrapBin]->Fill(r_RecoPhotRaw, eventWeight);
          extrap_responseMPFRaw[etaBin][ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);

          if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {
            extrap_responseBalancingRawGen[etaBinGen][ptBinGen][rawExtrapBin]->Fill(r_RecoGenRaw, eventWeight);
          }
        } while (false);
	
      } while (false);
      
      // New extrapolation -- done but not used later
      do {
	// Federico --> No cuts
        // Cut on photon pt. The first two bins are too low stats for beeing usefull
	//        if (photon.pt < 165)
	//          break;

        float r_RecoPhot        = firstJet.pt / photon.pt;
        float r_RecoPhotRaw = firstRawJet.pt / photon.pt;
        float alpha                 = secondJet.pt / photon.pt;
        float raw_alpha         = secondRawJet.pt / photon.pt;
	
        // Special case
	if (fabs(firstJet.eta) < 1.305) {
          new_extrap_responseBalancingEta013->fill(alpha, r_RecoPhot, eventWeight);
          new_extrap_responseBalancingRawEta013->fill(raw_alpha, r_RecoPhotRaw, eventWeight);
          new_extrap_responseMPFEta013->fill(alpha, respMPF, eventWeight);
          new_extrap_responseMPFRawEta013->fill(raw_alpha, respMPFRaw, eventWeight);
	}
	
        if (etaBin < 0)
          break;

	//	if( fabs(firstJet.eta) < 0.783 && alpha <0.30){
	//	  NewCounter++;
	//	  std::cout<< " NEW Extrap counter "<<NewCounter << std::endl;
	//	  std::cout<< " alpha="<< alpha << std::endl;
	//	}

        new_extrap_responseBalancing[etaBin]         ->fill(alpha, r_RecoPhot, eventWeight);
        new_extrap_responseBalancingRaw[etaBin]  ->fill(raw_alpha, r_RecoPhotRaw, eventWeight);
        new_extrap_responseMPF[etaBin]                 ->fill(alpha, respMPF, eventWeight);
        new_extrap_responseMPFRaw[etaBin]          ->fill(raw_alpha, respMPFRaw, eventWeight);

      } while (false);
    } // if(secondJet.is_present)
    
    float vpar=0.;
    //float vparRaw=0.;
    //float vparGen=0.;

    //    std::cout<<"If secondJet is Ok"<<std::endl;
       
    if (secondJetOK) { // ! is_present || pT < 10 || pT < 0.2*pT(pho)
      do {
        h_deltaPhi_passedID          ->Fill(deltaPhi, eventWeight);
        h_ptPhoton_passedID        ->Fill(photon.pt, eventWeight);
        h_ptPhoton_passedID_Binned        ->Fill(photon.pt, eventWeight);
	h_EtaPhoton_passedID      ->Fill(photon.eta, eventWeight);
	h_PhiPhoton_passedID      ->Fill(photon.phi, eventWeight);
        h_ptFirstJet_passedID       ->Fill(firstJet.pt, eventWeight);
	h_EtaFirstJet_passedID     ->Fill(firstJet.eta, eventWeight);
	h_PhiFirstJet_passedID     ->Fill(firstJet.phi, eventWeight);
        h_ptSecondJet_passedID       ->Fill(secondJet.pt, eventWeight);
	h_EtaSecondJet_passedID     ->Fill(secondJet.eta, eventWeight);
	h_PtEtaSecondJet_passedID  ->Fill(secondJet.eta, secondJet.pt, eventWeight);
	h_PhiSecondJet_passedID     ->Fill(secondJet.phi, eventWeight);
        h_firstJetvsSecondJet            ->Fill(firstJet.pt, secondJet.pt, eventWeight);

	if(secondJet.is_present ) {
        h_ptSecondJet_2ndJetOK       ->Fill(secondJet.pt, eventWeight);
	h_EtaSecondJet_2ndJetOK     ->Fill(secondJet.eta, eventWeight);
	h_PhiSecondJet_2ndJetOK     ->Fill(secondJet.phi, eventWeight);
	}

        h_alpha_passedID            ->Fill(secondJet.pt / photon.pt, eventWeight);
        h_MET_passedID              ->Fill(MET.et, eventWeight);
        h_rawMET_passedID        ->Fill(rawMET.et, eventWeight);


//jet energy composition
        h_JetCHEn_passedID              ->Fill(firstJet.jet_CHEn, eventWeight);
        h_JetNHEn_passedID              ->Fill(firstJet.jet_NHEn, eventWeight);
        h_JetCEEn_passedID              ->Fill(firstJet.jet_CEEn, eventWeight);
        h_JetNEEn_passedID              ->Fill(firstJet.jet_NEEn, eventWeight);
        h_JetPhEn_passedID               ->Fill(firstJet.jet_PhEn, eventWeight);
        h_JetElEn_passedID                ->Fill(firstJet.jet_ElEn, eventWeight);
        h_JetMuEn_passedID              ->Fill(firstJet.jet_MuEn, eventWeight);
//jet multiplicity
//        h_JetPhMult_passedID              ->Fill(firstJet.jet_PhMult, eventWeight);
//        h_JetNHMult_passedID              ->Fill(firstJet.jet_NHMult, eventWeight);
//        h_JetElMult_passedID              ->Fill(firstJet.jet_ElMult, eventWeight);
//        h_JetCHMult_passedID              ->Fill(firstJet.jet_CHMult, eventWeight);

        h_ptPhotonBinned_passedID[ptBin]->Fill(photon.pt, eventWeight);
	h_responseEtaBinned_passedID[etaBin]->Fill(respBalancing, eventWeight);

        h_METvsfirstJet                   ->Fill(MET.et, firstJet.pt, eventWeight);

        h_rho_passedID                         ->Fill(photon.rho, eventWeight);
        h_hadTowOverEm_passedID     ->Fill(photon.hadTowOverEm, eventWeight);
        h_sigmaIetaIeta_passedID         ->Fill(photon.sigmaIetaIeta, eventWeight);
        h_chargedHadronsIsolation_passedID     ->Fill(photon.chargedHadronsIsolation, eventWeight);
        h_neutralHadronsIsolation_passedID      ->Fill(photon.neutralHadronsIsolation, eventWeight);
        h_photonIsolation_passedID                    ->Fill(photon.photonIsolation, eventWeight);

	//	cout<< mu << " " << analysis.nvertexGood << " " << photon.rho << " " << eventWeight << endl;
	h_mu -> Fill(mu, eventWeight);
	h_npvGood -> Fill(analysis.nvertexGood, eventWeight);
	h_rho_vs_mu -> Fill(mu, photon.rho, eventWeight);
	h_npvGood_vs_mu -> Fill(mu, analysis.nvertexGood, eventWeight);

        // MET component
        vpar=(MET.px*photon.px + MET.py*photon.py)/photon.pt;
        h_METResolution_passedID             ->Fill(sqrt(pow(MET.px-genMET.px,2)+pow(MET.py-genMET.py,2)), eventWeight);
        h_MET_par_passedID                       ->Fill((MET.px*photon.px + MET.py*photon.py)/photon.pt, eventWeight);
        h_MET_perp_passedID                     ->Fill(MET.pt*(1.-pow(vpar/MET.pt,2)), eventWeight);

////fill resolution histos (for mc only)
//       if (mIsMC && genPhoton.pt!=0. && genMET.pt!=0.) {
//        h_METResolution_passedID->Fill(sqrt(pow(MET.px-genMET.px,2)+pow(MET.py-genMET.py,2)), eventWeight);
//         if(photon.regressionEnergy!=0.){
//         h_phPt_resolution->Fill( sqrt(pow((photon.px*photon.originalEnergy/photon.regressionEnergy)-genPhoton.px,2)+pow((photon.py*photon.originalEnergy/photon.regressionEnergy)-genPhoton.py,2)) , eventWeight); 
//         h_phPx_resolution->Fill((photon.px*photon.originalEnergy/photon.regressionEnergy)-genPhoton.px, eventWeight);
//         h_phPy_resolution->Fill((photon.py*photon.originalEnergy/photon.regressionEnergy)-genPhoton.py, eventWeight);
//          }
//         h_phPt_regression_resolution->Fill(sqrt(pow(photon.px-genPhoton.px,2)+pow(photon.py-genPhoton.py,2)) , eventWeight);
//         h_phPx_regression_resolution->Fill(photon.px-genPhoton.px, eventWeight);
//         h_phPy_regression_resolution->Fill(photon.py-genPhoton.py, eventWeight);
//         vparRaw=(rawMET.px*photon.px + rawMET.py*photon.py)/photon.pt;
//         vparGen=(genMET.px*photon.px + genMET.py*photon.py)/photon.pt;
//         h_MET_resolution->Fill(sqrt(pow(rawMET.px-genMET.px,2)+pow(rawMET.py-genMET.py,2)), eventWeight);
//         h_MET_par_resolution->Fill(vpar -  vparGen, eventWeight);
//         h_MET_perp_resolution->Fill(sqrt(pow((MET.px-(vpar*photon.px/photon.pt)),2)+pow((MET.py-(vpar*photon.py/photon.pt)),2))  , eventWeight);
//         h_MET_perp_resolution->Fill( rawMET.pt*(1.-pow(vparRaw/rawMET.pt,2))-genMET.pt*(1.-pow(vparGen/genMET.pt,2)), eventWeight);
//         h_MET_footprint_resolution->Fill(sqrt(pow(photon.footprintMExCorr-genMET.px,2) + pow(photon.footprintMEyCorr-genMET.py,2) ), eventWeight);
//         h_MET_par_footprint_resolution->Fill( vpar-(genMET.px*genPhoton.px + genMET.py*genPhoton.py)/genPhoton.pt  , eventWeight);
//         h_MET_perp_footprint_resolution->Fill( MET.pt*(1.-pow(vpar/MET.pt,2))-genMET.pt*(1.-pow(vparGen/genMET.pt,2)), eventWeight);
//      }
	
        // Special case
        if (fabs(firstJet.eta) <1.305) {
	  responseBalancingEta013[ptBin]->Fill(respBalancing, eventWeight);
	  responseBalancingRawEta013[ptBin]->Fill(respBalancingRaw, eventWeight);
	  responseMPFEta013[ptBin]->Fill(respMPF, eventWeight);
	  responseMPFRawEta013[ptBin]->Fill(respMPFRaw, eventWeight);

	  //	  if(!mIsMC){
	    //	    run_responseBalancingEta013[runBin]->Fill(respBalancing, eventWeight);
	    //	    run_responseMPFEta013[runBin]->Fill(respMPF, eventWeight);
	  //	  }

          if (vertexBin >= 0) {
            vertex_responseBalancingEta013[vertexBin]->Fill(respBalancing, eventWeight);
	    vertex_responseBalancingRawEta013[vertexBin]->Fill(respBalancingRaw, eventWeight);
	    vertex_responseMPFEta013[vertexBin]->Fill(respMPF, eventWeight);
	    vertex_responseMPFRawEta013[vertexBin]->Fill(respMPF, eventWeight);
	    //	    vertex_DeltapT[vertexBin]->Fill(firstJet.pt-(photon.pt*fabs(cos(deltaPhi))),eventWeight);
          }
	  if (mIsMC && ptBinGen >= 0) {
	    responseBalancingGenEta013[ptBinGen]->Fill(respBalancingGen, eventWeight);
            responseBalancingRawGenEta013[ptBinGen]->Fill(respBalancingRawGen, eventWeight);
            responseMPFGenEta013[ptBinGen]->Fill(respMPFGen, eventWeight);
	  }
        }
	//federico-> others eta bins
        if (fabs(firstJet.eta) >=1.305 && fabs(firstJet.eta) < 3.0 ) {
	  if (vertexBin >= 0) {
	    vertex_responseBalancingEta030[vertexBin]->Fill(respBalancing, eventWeight);
	  }
	}
        if (fabs(firstJet.eta) >=3.0 && fabs(firstJet.eta) < 5.0 ) {
	  if (vertexBin >= 0) {
	    vertex_responseBalancingEta050[vertexBin]->Fill(respBalancing, eventWeight);
	  }
	}
	/*
        if (fabs(firstJet.eta) < 2.4 && (fabs(firstJet.eta) < 1.4442 || fabs(firstJet.eta) > 1.5560)){ 
          // Viola
          ptFirstJetEta024[ptBin]->Fill(firstJet.pt, eventWeight);
          responseBalancingEta024[ptBin]->Fill(respBalancing, eventWeight);
          responseMPFEta024[ptBin]->Fill(respMPF, eventWeight);
	  }*/
	
        if (etaBin < 0) {
          //std::cout << "Jet eta " << firstJet.eta() << " is not covered by our eta binning. Dumping event." << std::endl;
          break;
        }
	
	if (mIsMC) {
	  if(firstGenJet.pt>0.) {
	    //       std::cout<< "responsetrue = "<< firstJet.pt / firstGenJet.pt<< " and weight "<< eventWeight<< std::endl;
	    responseTrue[etaBin][ptBin]->Fill(firstJet.pt / firstGenJet.pt, eventWeight);
	  }
	  if(photon.pt>0.) {
	    //       std::cout << "responsePLI = "<< firstGenJet.pt / photon.pt << " and weight "<< eventWeight<<  std::endl;
	    responsePLI[etaBin][ptBin]->Fill(firstGenJet.pt / photon.pt, eventWeight);
	  }
	}
	//fill N vertices as a function of eta/pT
	Nvertices[etaBin][ptBin]->Fill(analysis.nvertex, eventWeight);
	
	//fill jet energy composition histo vectors
       ChHadronEnergy[etaBin][ptBin]->Fill(firstJet.jet_CHEn, eventWeight);
       NHadronEnergy[etaBin][ptBin]->Fill(firstJet.jet_NHEn, eventWeight);
       ElEnergy[etaBin][ptBin]->Fill(firstJet.jet_ElEn, eventWeight);
       PhEnergy[etaBin][ptBin]->Fill(firstJet.jet_PhEn, eventWeight);
       MuEnergy[etaBin][ptBin]->Fill(firstJet.jet_MuEn, eventWeight);
       TotJetEnergy[etaBin][ptBin]->Fill(firstJet.e, eventWeight);
       //fill jet energy composition fractions histo vectors
       jetcalcen=firstJet.jet_CHEn+firstJet.jet_NHEn+firstJet.jet_ElEn+firstJet.jet_PhEn+firstJet.jet_MuEn;
       jetcalcenraw=firstRawJet.jet_CHEn+firstRawJet.jet_NHEn+firstRawJet.jet_ElEn+firstRawJet.jet_PhEn+firstRawJet.jet_MuEn;
       if(firstJet.e > 0.) {
	 ChHadron_realFraction[etaBin][ptBin]->Fill(firstJet.jet_CHEn/firstJet.e, eventWeight);
	 NHadron_realFraction[etaBin][ptBin]->Fill(firstJet.jet_NHEn/firstJet.e, eventWeight);
	 El_realFraction[etaBin][ptBin]->Fill(firstJet.jet_ElEn/firstJet.e, eventWeight);
	 Ph_realFraction[etaBin][ptBin]->Fill(firstJet.jet_PhEn/firstJet.e, eventWeight);
	 Mu_realFraction[etaBin][ptBin]->Fill(firstJet.jet_MuEn/firstJet.e, eventWeight);
       }
       if(jetcalcen > 0.) {
        ChHadronFraction[etaBin][ptBin]->Fill(firstJet.jet_CHEn/jetcalcen, eventWeight);
        NHadronFraction[etaBin][ptBin]->Fill(firstJet.jet_NHEn/jetcalcen, eventWeight);
        ElFraction[etaBin][ptBin]->Fill(firstJet.jet_ElEn/jetcalcen, eventWeight);
        PhFraction[etaBin][ptBin]->Fill(firstJet.jet_PhEn/jetcalcen, eventWeight);
        MuFraction[etaBin][ptBin]->Fill(firstJet.jet_MuEn/jetcalcen, eventWeight);
        LeptFraction[etaBin][ptBin]->Fill((firstJet.jet_MuEn+firstJet.jet_ElEn)/jetcalcen, eventWeight);
	//
        ChHadronFraction_mpf[etaBin][ptBin]->Fill(firstJet.jet_CHEn*respMPF/jetcalcen, eventWeight);
        NHadronFraction_mpf[etaBin][ptBin]->Fill(firstJet.jet_NHEn*respMPF/jetcalcen, eventWeight);
        ElFraction_mpf[etaBin][ptBin]->Fill(firstJet.jet_ElEn*respMPF/jetcalcen, eventWeight);
        PhFraction_mpf[etaBin][ptBin]->Fill(firstJet.jet_PhEn*respMPF/jetcalcen, eventWeight);
        MuFraction_mpf[etaBin][ptBin]->Fill(firstJet.jet_MuEn*respMPF/jetcalcen, eventWeight);
        LeptFraction_mpf[etaBin][ptBin]->Fill((firstJet.jet_MuEn+firstJet.jet_ElEn)*respMPF/jetcalcen, eventWeight);
        }
        if(jetcalcenraw > 0.) {
        ChHadronFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_CHEn/jetcalcenraw, eventWeight);
        NHadronFractionRaw[etaBin][ptBin]->Fill(firstJet.jet_NHEn/jetcalcenraw, eventWeight);
        ElFractionRaw[etaBin][ptBin]->Fill(firstJet.jet_ElEn/jetcalcenraw, eventWeight);
        PhFractionRaw[etaBin][ptBin]->Fill(firstJet.jet_PhEn/jetcalcenraw, eventWeight);
        MuFractionRaw[etaBin][ptBin]->Fill(firstJet.jet_MuEn/jetcalcenraw, eventWeight);
        LeptFractionRaw[etaBin][ptBin]->Fill((firstJet.jet_MuEn+firstJet.jet_ElEn)/jetcalcenraw, eventWeight);
	//
        ChHadronFractionRaw_mpf[etaBin][ptBin]->Fill(firstRawJet.jet_CHEn*respMPFRaw/jetcalcenraw, eventWeight);
        NHadronFractionRaw_mpf[etaBin][ptBin]->Fill(firstJet.jet_NHEn*respMPFRaw/jetcalcenraw, eventWeight);
        ElFractionRaw_mpf[etaBin][ptBin]->Fill(firstJet.jet_ElEn*respMPFRaw/jetcalcenraw, eventWeight);
        PhFractionRaw_mpf[etaBin][ptBin]->Fill(firstJet.jet_PhEn*respMPFRaw/jetcalcenraw, eventWeight);
        MuFractionRaw_mpf[etaBin][ptBin]->Fill(firstJet.jet_MuEn*respMPFRaw/jetcalcenraw, eventWeight);
        LeptFractionRaw_mpf[etaBin][ptBin]->Fill((firstJet.jet_MuEn+firstJet.jet_ElEn)*respMPFRaw/jetcalcenraw, eventWeight);
        }
        if(firstRawJet.e > 0.) {
        ChHadron_realFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_CHEn/firstRawJet.e, eventWeight);
        NHadron_realFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_NHEn/firstRawJet.e, eventWeight);
        El_realFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_ElEn/firstRawJet.e, eventWeight);
        Ph_realFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_PhEn/firstRawJet.e, eventWeight);
	Mu_realFractionRaw[etaBin][ptBin]->Fill(firstRawJet.jet_MuEn/firstRawJet.e, eventWeight);
        }
        //fill jet multiplicities histo vectors
        ChHadronMult[etaBin][ptBin]->Fill(firstJet.jet_CHMult, eventWeight);
        NHadronMult[etaBin][ptBin]->Fill(firstJet.jet_NHMult, eventWeight);
        ElMult[etaBin][ptBin]->Fill(firstJet.jet_ElMult, eventWeight);
        PhMult[etaBin][ptBin]->Fill(firstJet.jet_PhMult, eventWeight);

        responseBalancing[etaBin][ptBin]->Fill(respBalancing, eventWeight);
        responseBalancingRaw[etaBin][ptBin]->Fill(respBalancingRaw, eventWeight);

        responseMPF[etaBin][ptBin]->Fill(respMPF, eventWeight);
        responseMPFRaw[etaBin][ptBin]->Fill(respMPFRaw, eventWeight);

	//	if(!mIsMC){
	//	  run_responseBalancing[etaBin][runBin]->Fill(respBalancing, eventWeight);
	//	  run_responseMPF[etaBin][runBin]->Fill(respMPF, eventWeight);
	//	}

        if (vertexBin >= 0) {
          vertex_responseBalancing[etaBin][vertexBin]->Fill(respBalancing, eventWeight);
          vertex_responseBalancingRaw[etaBin][vertexBin]->Fill(respBalancingRaw, eventWeight);	  
          vertex_responseMPF[etaBin][vertexBin]->Fill(respMPF, eventWeight);
          vertex_responseMPFRaw[etaBin][vertexBin]->Fill(respMPF, eventWeight);
        }

        // Gen values
        if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {

          responseBalancingGen[etaBinGen][ptBinGen]->Fill(respBalancingGen, eventWeight);
          responseBalancingRawGen[etaBinGen][ptBinGen]->Fill(respBalancingRawGen, eventWeight);	  
          responseMPFGen[etaBinGen][ptBinGen]->Fill(respMPFGen, eventWeight);

	  responsePhotGamma[etaBinGen][ptBinGen]  ->Fill(respPhotonGamma, eventWeight);
	  responsePhotGamma_Integrated ->Fill(respPhotonGamma, eventWeight);
          responseBalancingGen_Integrated->Fill(respBalancingGen, eventWeight);

        }
      } while (false);

#if ADD_TREES
      if (! mUncutTrees) {
        photonTree->Fill();
        if (mIsMC)
          genPhotonTree->Fill();
        firstJetTree->Fill();
        if (mIsMC)
          firstGenJetTree->Fill();
        firstRawJetTree->Fill();
        secondJetTree->Fill();
        if (mIsMC)
          secondGenJetTree->Fill();
        secondRawJetTree->Fill();
        metTree->Fill();
        rawMetTree->Fill();
        if (mIsMC)
          genMetTree->Fill();
        electronsTree->Fill();
        muonsTree->Fill();
        analysisTree->Fill();
        miscTree->Fill();
      }
#endif

      passedEvents++;
      //  std::cout<<"passedEvents"<<std::endl;

    }// if secondJetOK --> end of "passedID" histo



  }

  analysisDir.cd();
  h_ptPhoton_passedID_Binned->Write();

  std::cout << std::endl;
  std::cout << "Absolute efficiency : related to initial number of event =  " << to-from << std::endl;
  std::cout << "Efficiency for photon/jet cut: " << MAKE_RED << (double) passedPhotonJetCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency for trigger selection: " << MAKE_RED << (double) passedEventsFromTriggers / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for Δφ cut: " << MAKE_RED << (double) passedDeltaPhiCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pixel seed veto cut: " << MAKE_RED << (double) passedPixelSeedVetoCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for muons cut: " << MAKE_RED << (double) passedMuonsCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for electrons cut: " << MAKE_RED << (double) passedElectronsCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pT(j1) cut: " << MAKE_RED << (double) passedJetPtCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for α cut: " << MAKE_RED << (double) passedAlphaCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency: " << MAKE_RED << (double) passedEvents / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << std::endl;
  std::cout << "Relative efficiency --- Initial events : " << to-from << std::endl;
  std::cout << "Efficiency for photon/jet cut: " << MAKE_RED << (double) passedPhotonJetCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency for trigger selection: " << MAKE_RED << (double) passedEventsFromTriggers / passedPhotonJetCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for Δφ cut: " << MAKE_RED << (double) passedDeltaPhiCut / passedEventsFromTriggers * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pixel seed veto cut: " << MAKE_RED << (double) passedPixelSeedVetoCut / passedDeltaPhiCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for muons cut: " << MAKE_RED << (double) passedMuonsCut / passedPixelSeedVetoCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for electrons cut: " << MAKE_RED << (double) passedElectronsCut / passedMuonsCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pT(j1) cut: " << MAKE_RED << (double) passedJetPtCut / passedElectronsCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for α cut: " << MAKE_RED << (double) passedAlphaCut / passedJetPtCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << std::endl;
  std::cout<< "Entries degli isto senza passedID -->    " << passedJetPtCut << std::endl;
  std::cout<< "Entries degli isto con passedID -->    " << passedEvents << std::endl;

  std::cout << std::endl;
  std::cout << "Rejected events because trigger was not found: " << MAKE_RED << (double) rejectedEventsTriggerNotFound / (rejectedEventsFromTriggers) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Rejected events because pT was out of range: " << MAKE_RED << (double) rejectedEventsPtOut / (rejectedEventsFromTriggers) * 100 << "%" << RESET_COLOR << std::endl;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

  bool appendText = (xMin >= 0 && xMax >= 0);
  std::vector<T*> vector;
  size_t ptBinningSize = mPtBinning.size();
  for (size_t j = 0; j < ptBinningSize; j++) {

    const std::pair<float, float> bin = mPtBinning.getBinValue(j);
    std::stringstream ss;
    if (appendText)
      ss << branchName << "_ptPhot_" << (int) bin.first << "_" << (int) bin.second;
    else
      ss << branchName << "_" << (int) bin.first << "_" << (int) bin.second;

    if (!appendText) {
      xMin = bin.first;
    }

    if (!appendText) {
      xMax = bin.second;
    }

    T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
    vector.push_back(object);
  }

  return vector;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {
  return buildPtVector<T>(dir, branchName + "_" + etaName, nBins, xMin, xMax);
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

  bool appendText = (xMin >= 0 && xMax >= 0);
  std::vector<T*> vector;
  size_t etaBinningSize = mEtaBinning.size();
  for (size_t j = 0; j < etaBinningSize; j++) {

    const std::pair<float, float> bin = mEtaBinning.getBinValue(j);
    std::stringstream ss;
    if (appendText)
      ss << branchName << "_" << mEtaBinning.getBinName(j);
    else
      ss << branchName << "_" << mEtaBinning.getBinName(j);

    if (!appendText) {
      xMin = bin.first;
    }

    if (!appendText) {
      xMax = bin.second;
    }

    T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
    vector.push_back(object);
  }

  return vector;
}

template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildEtaPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::vector<T*> > etaBinning;

  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);
    etaBinning.push_back(buildPtVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
  }

  return etaBinning;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildVertexVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {

  std::vector<T*> vector;
  size_t vertexBinningSize = mVertexBinning.size();
  for (size_t j = 0; j < vertexBinningSize; j++) {

    const std::pair<int, int> bin = mVertexBinning.getBinValue(j);
    std::stringstream ss;
    ss << branchName << "_" << etaName << "_nvertex_" << bin.first << "_" << bin.second;

    T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
    vector.push_back(object);
  }

  return vector;
}

template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildEtaVertexVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::vector<T*> > etaBinning;

  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);
    etaBinning.push_back(buildVertexVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
  }

  return etaBinning;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildRunVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {

  std::vector<T*> vector;
  size_t runBinningSize = mRunBinning.size();
  for (size_t j = 0; j < runBinningSize; j++) {

    const std::pair<int, int> bin = mRunBinning.getBinValue(j);
    std::stringstream ss;
    ss << branchName << "_" << etaName << "_run_" << bin.first << "_" << bin.second;

    T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
    vector.push_back(object);
  }

  return vector;
}

template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildEtaRunVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::vector<T*> > etaBinning;

  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);
    etaBinning.push_back(buildRunVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
  }

  return etaBinning;
}

template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildExtrapolationVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {

  std::vector<std::vector<T*> > vector;
  size_t ptBinningSize = mPtBinning.size();
  for (size_t j = 0; j < ptBinningSize; j++) {

    const std::pair<float, float> bin = mPtBinning.getBinValue(j);
    std::stringstream ss;
    ss << branchName << "_" << etaName;

    TString subDirectoryName = TString::Format("extrap_ptPhot_%d_%d", (int) bin.first, (int) bin.second);
    TFileDirectory subDir = dir.mkdir(subDirectoryName.Data());

    std::vector<T*> subvector;
    size_t extrapBinningSize = mExtrapBinning.size();
    for (size_t p = 0; p < extrapBinningSize; p++) {
      TString name = TString::Format("%s_%d", ss.str().c_str(), (int) p);

      T* object = subDir.make<T>(name, name, nBins, xMin, xMax);
      subvector.push_back(object);
    }

    vector.push_back(subvector);
  }

  return vector;
}

template<typename T>
std::vector<std::vector<std::vector<T*> > > GammaJetFinalizer::buildExtrapolationEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::vector<std::vector<T*> > > etaBinning;

  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);

    std::vector<std::vector<T*> > vector = buildExtrapolationVector<T>(dir, branchName, etaName, nBins, xMin, xMax);
    etaBinning.push_back(vector);
  }

  return etaBinning;
}

std::shared_ptr<GaussianProfile> GammaJetFinalizer::buildNewExtrapolationVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {

  std::stringstream ss;
  ss << branchName << "_" << etaName;

  std::shared_ptr<GaussianProfile> object(new GaussianProfile(ss.str(), mNewExtrapBinning.size(), 0, mNewExtrapBinning.size() * mNewExtrapBinning.getBinWidth(), nBins, xMin, xMax));
  object->setPrefix("alpha");
  object->initialize(dir);

  return object;
}

std::vector<std::shared_ptr<GaussianProfile>> GammaJetFinalizer::buildNewExtrapolationEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

  size_t etaBinningSize = mEtaBinning.size();
  std::vector<std::shared_ptr<GaussianProfile>> etaBinning;

  for (size_t i = 0; i < etaBinningSize; i++) {
    const std::string etaName = mEtaBinning.getBinName(i);
    

    std::shared_ptr<GaussianProfile> object = buildNewExtrapolationVector(dir, branchName, etaName, nBins, xMin, xMax);
    etaBinning.push_back(object);
  }
  return etaBinning;
}

void GammaJetFinalizer::cleanTriggerName(std::string& trigger) {
  boost::replace_first(trigger, "_.*", "");
  boost::replace_first(trigger, ".*", "");
}

//// PU Reweighting 2015
//void GammaJetFinalizer::computePUWeight(const std::string& passedTrigger, int run_period) {
void GammaJetFinalizer::computePUWeight() {

  //  std::cout<< "Using PU Reweighting CODE"<<std::endl;
  
  if (mNoPUReweighting)
    return;
  
  mPUWeight = reweighter->weight(analysis.ntrue_interactions);
  //  std::cout<<analysis.ntrue_interactions<<std::endl;  
  //  std::cout<<mPUWeight<<std::endl;

}//end compute PUReweight


void GammaJetFinalizer::computePUWeight_NVtxBased(double ptPhot, int nvertex) {

  //  std::cout << "My PU reweighting   "<< nvertex<<std::endl;

  if (mNoPUReweighting)
    return;

  TH1D *h_ratio=0;
  if(ptPhot >= 40 && ptPhot < 60)             h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_40_60");  
  if(ptPhot >= 60 && ptPhot < 85)             h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_60_85");  
  if(ptPhot >= 85 && ptPhot < 100)           h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_85_100");  
  if(ptPhot >= 100 && ptPhot < 130)         h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_100_130");  
  if(ptPhot >= 130 && ptPhot < 175)         h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_130_175");  
  if(ptPhot >= 175 )                                     h_ratio = (TH1D*)PUFile->Get("h_ratio_ptPhot_175_Inf");  
  
  int bin = h_ratio->FindBin(nvertex);
  mPUWeight = h_ratio->GetBinContent(bin);  
  if(mVerbose) std::cout<< "Nvtx  "<<nvertex<< std::endl;
  if(mVerbose) std::cout<< "bin "<<bin<< std::endl;
  if(mVerbose) std::cout<< "PU Weight  "<<mPUWeight<< std::endl;
}

void GammaJetFinalizer::checkInputFiles() {
  for (std::vector<std::string>::iterator it = mInputFiles.begin(); it != mInputFiles.end();) {
    TFile* f = TFile::Open(it->c_str());
    if (! f) {
      std::cerr << "Error: can't open '" << it->c_str() << "'. Removed from input files." << std::endl;
      it = mInputFiles.erase(it);
      continue;
    }

    TTree* analysis = static_cast<TTree*>(f->Get("gammaJet/analysis"));
    if (! analysis || analysis->GetEntry(0) == 0) {
      std::cerr << "Error: Trees inside '" << it->c_str() << "' were empty. Removed from input files." << std::endl;
      it = mInputFiles.erase(it);

      f->Close();
      delete f;

      continue;
    }

    f->Close();
    delete f;

    ++it;
  }
}

int GammaJetFinalizer::checkTrigger(std::string& passedTrigger, float& weight) {

  if (! mIsMC) {
    const PathVector& mandatoryTriggers = mTriggers->getTriggers(analysis.run);
    
    // - With the photon p_t, find the trigger it should pass
    // - Then, look on trigger list if it pass it or not (only for data)

    //    std::cout << "photon.pt  " << photon.pt << std::endl;
    
    const PathData* mandatoryTrigger = nullptr;
    for (auto& path: mandatoryTriggers) {
      if (path.second.range.in(photon.pt)) {
        mandatoryTrigger = &path;
      }
    }
   
    if (!mandatoryTrigger)
      return TRIGGER_FOUND_BUT_PT_OUT;

    //    weight = mandatoryTrigger->second.weight; // from file
    //    std::cout << "weight  " << weight << std::endl;
    
    // Photon had to pass mandatoryTrigger.first
    size_t size = analysis.trigger_names->size();
 
    for (int i = size - 1; i >= 0; i--) {
      bool passed = analysis.trigger_results->at(i);
      //      if (!passed) std::cout << "Trigger NOT passed" <<std::endl;
      if (! passed)
        continue;
      
      if (boost::regex_match(analysis.trigger_names->at(i), mandatoryTrigger->first)) {
	//	std::cout << "Triggers  matching"<<std::endl;
        passedTrigger = mandatoryTrigger->first.str();
	//	std::cout << "Trigger name   " << analysis.trigger_names->at(i) << std::endl;
	//	std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	weight = analysis.trigger_prescale->at(i); // prescale from ntupla
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
    }
  } else { // IsMC

    const std::map<Range<float>, std::vector<MCTrigger>>& triggers = mMCTriggers->getTriggers();

    //    std::cout<<photon.pt <<std::endl;
    
    const std::vector<MCTrigger>* mandatoryTrigger = nullptr;
    for (auto& path: triggers) {
      if (path.first.in(photon.pt)) {
	mandatoryTrigger = &path.second;
      }
    }
    
    if (!mandatoryTrigger)
      return TRIGGER_FOUND_BUT_PT_OUT;;
    
    /*
    if (mandatoryTrigger->size() > 1) {
      double random = mRandomGenerator.Rndm();

      double weight_low = 0;
      double weight_high = 0;
      for (int32_t i = 0; i < (int32_t) mandatoryTrigger->size(); i++) {
        if (i - 1 >= 0)
          weight_low += (*mandatoryTrigger)[i - 1].weight;
        weight_high += (*mandatoryTrigger)[i].weight;

        if (random > weight_low && random <= weight_high) {
          passedTrigger = (*mandatoryTrigger)[i].name.str();
          return TRIGGER_OK;
        } //if random
      } // for i

      throw new std::exception(); // This should NEVER happened
      }*/

    return TRIGGER_OK;

    // added require "trigger passed" also for MC
    size_t size = analysis.trigger_names->size();
    
    for (int i = size - 1; i >= 0; i--) {
      bool passed = analysis.trigger_results->at(i);
      if (!passed) std::cout << "Trigger NOT passed" <<std::endl;
      if (! passed)
	continue;
      
      //      if (boost::regex_match(analysis.trigger_names->at(i), mandatoryTrigger->first)) {
      if (boost::regex_match(analysis.trigger_names->at(i), mandatoryTrigger->at(0).name )) {
	std::cout << "Triggers  matching"<<std::endl;
	//passedTrigger = mandatoryTrigger->first.str();
	passedTrigger = mandatoryTrigger->at(0).name.str();
	std::cout << "Trigger name   " << analysis.trigger_names->at(i) << std::endl;
	//	std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	weight = analysis.trigger_prescale->at(i); // prescale from ntupla
	weight = 1;
	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
    }
    
    //    passedTrigger = mandatoryTrigger->at(0).name.str();
    //    return TRIGGER_OK;
  }
  // }else {
  
  //// Method 1
  //size_t size = analysis.trigger_names->size();
  //for (int i = size - 1; i >= 0; i--) {
  //bool passed = analysis.trigger_results->at(i);
  //if (! passed) continue;
  
  //for (const PathData& mandatoryTrigger: mandatoryTriggers) {
  //if (boost::regex_match(analysis.trigger_names->at(i), mandatoryTrigger.first)) {
  //passedTrigger = mandatoryTrigger.first.str();
  //weight = mandatoryTrigger.second.weight;
  //return TRIGGER_OK;
  //}
  //}
  //}
  //}
  return TRIGGER_NOT_FOUND;
}

std::vector<std::string> readInputFiles(const std::string& list) {
  std::ifstream f(list.c_str());
  std::string line;
  std::vector<std::string> files;
  while (std::getline(f, line)) {
    boost::algorithm::trim(line);
    if (line.length() == 0 || line[0] == '#')
      continue;

    files.push_back(line);
  }

  if (files.size() == 0) {
    throw new TCLAP::ArgException("No input files found in " + list);
  }

  return files;
}

void handleCtrlC(int s){
  EXIT = true;
}

int main(int argc, char** argv) {
  struct sigaction sigIntHandler;

  sigIntHandler.sa_handler = handleCtrlC;
  sigemptyset(&sigIntHandler.sa_mask);
  sigIntHandler.sa_flags = 0;

  sigaction(SIGINT, &sigIntHandler, NULL);

  try {
    TCLAP::CmdLine cmd("Step 3 of Gamma+Jet analysis", ' ', "0.1");

    TCLAP::ValueArg<std::string> datasetArg("d", "dataset", "Dataset name", true, "", "string", cmd);

    TCLAP::MultiArg<std::string> inputArg("i", "in", "Input file", true, "string");
    TCLAP::ValueArg<std::string> inputListArg("", "input-list", "Text file containing input files", true, "input.list", "string");
    cmd.xorAdd(inputArg, inputListArg);

    // Jet type
    std::vector<std::string> jetTypes;
    jetTypes.push_back("pf");
    //    jetTypes.push_back("calo");
    TCLAP::ValuesConstraint<std::string> allowedJetTypes(jetTypes);
    TCLAP::ValueArg<std::string> typeArg("", "type", "jet type", false, "pf", &allowedJetTypes, cmd);
    std::vector<std::string> algoTypes;
    algoTypes.push_back("ak4");
    algoTypes.push_back("ak8");
    TCLAP::ValuesConstraint<std::string> allowedAlgoTypes(algoTypes);
    TCLAP::ValueArg<std::string> algoArg("", "algo", "jet algo", true, "ak4", &allowedAlgoTypes, cmd);
    TCLAP::ValueArg<float> chsArg("", "chs", "Use CHS branches", false, true, "bool", cmd);
    // alpha cut
    TCLAP::ValueArg<float> alphaCutArg("", "alpha", "P_t^{second jet} / p_t^{photon} cut (default: 0.3)", false, 0.3, "float", cmd);

    TCLAP::SwitchArg mcArg("", "mc", "MC?", cmd);
    TCLAP::SwitchArg mcComparisonArg("", "mc-comp", "Cut photon pt to avoid trigger prescale issues", cmd);
    TCLAP::SwitchArg verboseArg("v", "verbose", "Enable verbose mode", cmd);
    TCLAP::SwitchArg uncutTreesArg("", "uncut-trees", "Fill trees before second jet cut", cmd);

    cmd.parse(argc, argv);

    //std::cout << "Initializing..." << std::endl;
    //gSystem->Load("libFWCoreFWLite.so");
    //AutoLibraryLoader::enable();
    //std::cout << "done." << std::endl;

    std::vector<std::string> files;
    if (inputArg.isSet()) {
      files = inputArg.getValue();
    } else {
      files = readInputFiles(inputListArg.getValue());
    }

    GammaJetFinalizer finalizer;
    finalizer.setInputFiles(files);
    finalizer.setDatasetName(datasetArg.getValue());
    finalizer.setJetAlgo(typeArg.getValue(), algoArg.getValue());
    finalizer.setMC(mcArg.getValue());
    finalizer.setMCComparison(mcComparisonArg.getValue());
    finalizer.setAlphaCut(alphaCutArg.getValue());
    finalizer.setCHS(chsArg.getValue());
    finalizer.setVerbose(verboseArg.getValue());
    finalizer.setUncutTrees(uncutTreesArg.getValue());
    //    if (totalJobsArg.isSet() && currentJobArg.isSet()) {
    //   finalizer.setBatchJob(currentJobArg.getValue(), totalJobsArg.getValue());
    //  }

    finalizer.runAnalysis();

  } catch (TCLAP::ArgException &e) {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    return 1;
  }
}
