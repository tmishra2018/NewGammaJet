#include <TFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TSystem.h>
#include <TTree.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TH2D.h>
#include <TLorentzVector.h>

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

GammaJetFinalizer::GammaJetFinalizer() {
  mPUWeight = 1.;
  mNoPUReweighting = false; // do not PUReweighting
}

GammaJetFinalizer::~GammaJetFinalizer() {
  
}

std::string GammaJetFinalizer::buildPostfix() {
  std::string algo = mJetAlgo == AK4 ? "AK4" : "AK8";
  std::string type = mJetType == PF ? "PFlow" : "PUPPI";
  
  std::string postfix = type + algo;
  
  if (mJetType == PF && mUseCHS) postfix += "chs";
  
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
    static std::string puPrefix = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/analysis/PUReweighting", cmsswBase.c_str()).Data();                          
    static std::string puMC = TString::Format("%s/computed_mc_GJET_Pythia_pu_truth_100bins.root", puPrefix.c_str()).Data();    
    static std::string puData = TString::Format("%s/pu_truth_data2016_100bins.root", puPrefix.c_str()).Data();                                                           
    reweighter = boost::shared_ptr<PUReweighter>(new PUReweighter(puData, puMC));                                                                                           
    // Trigger
    std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_mc.xml", cmsswBase.c_str()).Data();
    std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
    mMCTriggers      = new MCTriggers( TriggerFile.c_str() ) ;
  } else {
    // Get PU
    parsePileUpJSON2();
    //Trigger
    static std::string cmsswBase = getenv("CMSSW_BASE");
    std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_data.xml", cmsswBase.c_str()).Data();
    std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
    mTriggers      = new Triggers( TriggerFile.c_str() ) ;
  }
  
  // Initialization
  mExtrapBinning.initialize(mPtBinning, (mJetType == PF) ? "PFlow" : "PUPPI");

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



   TChain FullinfoChain("rootTupleTree/tree");


  loadFiles(FullinfoChain);
  
  fullinfo.Init(&FullinfoChain);
  
 

  std::cout << "done." << std::endl;

  std::cout << std::endl << "##########" << std::endl;
  std::cout << "# " << MAKE_BLUE << "Running on " << MAKE_RED << ((mIsMC) ? "MC" : "DATA") << RESET_COLOR << std::endl;
  std::cout << "##########" << std::endl << std::endl;

  // Output file
  std::string outputFile = TString::Format("PhotonJet_%s_%s.root", mDatasetName.c_str(), postFix.c_str()).Data();
  fwlite::TFileService fs(outputFile);

#if ADD_TREES
  
  TTree* fullinfoTree = NULL;
  cloneTree(fullinfo.fChain, fullinfoTree);

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
  TH1F* h_generatorWeight = analysisDir.make<TH1F>("generatorWeight", "generatorWeight", 50, 0., 1.);
  TH1F* h_analysis_evtWeightTot = analysisDir.make<TH1F>("analysis_evtWeightTot", "analysis_evtWeightTot", 50, 0., 10.);
  TH1F* h_event_weight_used = analysisDir.make<TH1F>("event_weight_used", "event_weight_used", 150, 0., 150.);

  TH1F* h_ptPhoton_NoCut = analysisDir.make<TH1F>("ptPhoton_NoCut", "ptPhoton NoCut", 150, 0., 3000.);

  TH1F* h_ptPhoton = analysisDir.make<TH1F>("ptPhoton", "ptPhoton", 600, 0., 3000.);
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

  TH1F* h_rho = analysisDir.make<TH1F>("rho", "rho", 100, 0, 50);
  TH1F* h_hadTowOverEm = analysisDir.make<TH1F>("hadTowOverEm", "hadTowOverEm", 100, 0, 0.05);
  TH1F* h_sigmaIetaIeta = analysisDir.make<TH1F>("sigmaIetaIeta", "sigmaIetaIeta", 100, 0, 0.011);
  TH1F* h_chargedHadronsIsolation = analysisDir.make<TH1F>("chargedHadronsIsolation", "chargedHadronsIsolation", 100, 0, 2);
  TH1F* h_neutralHadronsIsolation = analysisDir.make<TH1F>("neutralHadronsIsolation", "neutralHadronsIsolation", 100, 0, 100);
  TH1F* h_photonIsolation = analysisDir.make<TH1F>("photonIsolation", "photonIsolation", 100, 0, 15);
  
  TH2F* h_deltaPhi_vs_alpha = analysisDir.make<TH2F>("deltaPhi_vs_alpha", "deltaPhi_vs_alpha", 50, 0, 4, 50, 0, 0.5);

  TH1F* h_deltaPhi_passedID = analysisDir.make<TH1F>("deltaPhi_passedID", "deltaPhi", 40, M_PI / 2, M_PI);

  double ptBins[] = {40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 3000};
  int  binnum = sizeof(ptBins)/sizeof(double) -1;
  double etaBins[] = {0, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
  int  binnumEta = sizeof(etaBins)/sizeof(double) -1;

  TH1F* h_ptPhoton_Binned = new TH1F("ptPhoton_Binned","ptPhoton", binnum, ptBins); 
  TH1F* h_ptPhoton_passedID_Binned = new TH1F("ptPhoton_passedID_Binned","ptPhoton_passedID", binnum, ptBins);  

  TH1F* h_ptPhoton_passedID = analysisDir.make<TH1F>("ptPhoton_passedID", "ptPhoton", 600, 0., 3000.);
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

  TH1F* h_mu = analysisDir.make<TH1F>("mu", "mu", 50, 0, 50);
  TH1F* h_npvGood = analysisDir.make<TH1F>("npvGood", "npvGood", 50, 0, 50);
  TH2F* h_rho_vs_mu = analysisDir.make<TH2F>("rho_vs_mu", "Rho vs mu", 50, 0, 50, 100, 0, 50);
  TH2F* h_npvGood_vs_mu = analysisDir.make<TH2F>("npvGood_vs_mu", "npv_good vs mu", 50, 0, 50, 50, 0, 50);

  //jet composition - histos vectors
  TFileDirectory ecompositionDir = analysisDir.mkdir("ecomposition");
  //jet composition fractions - histos vectors
  std::vector<std::vector<TH1F*> > ChHadronFraction = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NHadronFraction = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > CEmFraction = buildEtaPtVector<TH1F>(ecompositionDir, "CEmFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > NEmFraction = buildEtaPtVector<TH1F>(ecompositionDir, "NEmFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > MuFraction = buildEtaPtVector<TH1F>(ecompositionDir, "MuFraction", 40, 0., 1.);
  std::vector<std::vector<TH1F*> > LeptFraction = buildEtaPtVector<TH1F>(ecompositionDir, "LeptFraction", 40, 0., 1.);
  std::vector<TH1F*> ChHadronFractionEta013    = buildPtVector<TH1F>(ecompositionDir, "ChHadronFraction", "eta0013", 40, 0., 1.);
  std::vector<TH1F*> NHadronFractionEta013    = buildPtVector<TH1F>(ecompositionDir, "NHadronFraction", "eta0013", 40, 0., 1.);
  std::vector<TH1F*> CEmFractionEta013    = buildPtVector<TH1F>(ecompositionDir, "CEmFraction", "eta0013", 40, 0., 1.);
  std::vector<TH1F*> NEmFractionEta013    = buildPtVector<TH1F>(ecompositionDir, "NEmFraction", "eta0013", 40, 0., 1.);
  std::vector<TH1F*> MuFractionEta013    = buildPtVector<TH1F>(ecompositionDir, "MuFraction", "eta0013", 40, 0., 1.);
  std::vector<TH1F*> LeptFractionEta013 = buildPtVector<TH1F>(ecompositionDir, "LeptFraction", "eta0013", 40, 0., 1.);
  //jet multiplicities
  std::vector<std::vector<TH1F*> > ChHadronMult = buildEtaPtVector<TH1F>(ecompositionDir, "ChHadronMult", 50, 0, 50);
  std::vector<std::vector<TH1F*> > NHadronMult = buildEtaPtVector<TH1F>(ecompositionDir, "NHadronMult", 50, 0, 50);
  std::vector<std::vector<TH1F*> > ElMult = buildEtaPtVector<TH1F>(ecompositionDir, "ElMult", 50, 0, 50);
  std::vector<std::vector<TH1F*> > PhMult = buildEtaPtVector<TH1F>(ecompositionDir, "PhMult", 50, 0, 50);
  std::vector<std::vector<TH1F*> > MuonMult = buildEtaPtVector<TH1F>(ecompositionDir, "MuonMult", 50, 0, 50);

  TH1F* h_rho_passedID = analysisDir.make<TH1F>("rho_passedID", "rho", 100, 0, 50);
  TH1F* h_hadTowOverEm_passedID = analysisDir.make<TH1F>("hadTowOverEm_passedID", "hadTowOverEm", 100, 0, 0.05);
  TH1F* h_sigmaIetaIeta_passedID = analysisDir.make<TH1F>("sigmaIetaIeta_passedID", "sigmaIetaIeta", 100, 0, 0.011);
  TH1F* h_chargedHadronsIsolation_passedID = analysisDir.make<TH1F>("chargedHadronsIsolation_passedID", "chargedHadronsIsolation", 100, 0, 2.0);
  TH1F* h_neutralHadronsIsolation_passedID = analysisDir.make<TH1F>("neutralHadronsIsolation_passedID", "neutralHadronsIsolation", 100, 0, 100);
  TH1F* h_photonIsolation_passedID = analysisDir.make<TH1F>("photonIsolation_passedID", "photonIsolation", 100, 0, 15);

  TH2F* h_METvsfirstJet = analysisDir.make<TH2F>("METvsfirstJet", "MET vs firstJet", 150, 0., 300., 150, 0., 500.);
  TH2F* h_firstJetvsSecondJet = analysisDir.make<TH2F>("firstJetvsSecondJet", "firstJet vs secondJet", 60, 5., 100., 60, 5., 100.);

  // Balancing
  TFileDirectory balancingDir = analysisDir.mkdir("balancing");
  std::vector<std::vector<TH1F*> > responseBalancing = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseBalancingRaw = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_raw", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseBalancingGen;
  std::vector<TH1F*> responseBalancingEta013       = buildPtVector<TH1F>(balancingDir, "resp_balancing", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseBalancingRawEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_raw", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseBalancingGenEta013;
  if (mIsMC) {
    responseBalancingGen = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_gen", 150, 0., 2.);
    responseBalancingGenEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_gen", "eta0013", 150, 0., 2.);
  }

  std::vector<std::vector<TH1F*> > responseBalancingGenPhot; // GenJet / photon
  std::vector<std::vector<TH1F*> > responseBalancingPhotGamma; // photon / GenPhoton
  std::vector<TH1F*> responseBalancingGenPhotEta013;
  std::vector<TH1F*> responseBalancingPhotGammaEta013;
  if (mIsMC) {
    responseBalancingGenPhot = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_gen_phot", 150, 0., 2.);
    responseBalancingPhotGamma = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing_photGamma", 150, 0., 2.);
    responseBalancingGenPhotEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_gen_phot", "eta0013", 150, 0., 2.);
    responseBalancingPhotGammaEta013 = buildPtVector<TH1F>(balancingDir, "resp_balancing_photGamma", "eta0013", 150, 0., 2.);
  }
  
  // MPF
  TFileDirectory mpfDir = analysisDir.mkdir("mpf");
  std::vector<std::vector<TH1F*> > responseMPF = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf", 150, 0., 2.);
   std::vector<std::vector<TH1F*> > responseMPFRaw = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_raw", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseMPFGen;
  std::vector<TH1F*> responseMPFEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFRawEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_raw", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFGenEta013;
  if (mIsMC) {
    responseMPFGen = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_gen", 150, 0., 2.);
    responseMPFGenEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_gen", "eta0013", 150, 0., 2.);
  }
  
  // vs number of vertices
  TFileDirectory vertexDir = analysisDir.mkdir("vertex");
  //check Nvtx vs ptphoton
  std::vector<std::vector<TH1F*> > Nvertices = buildEtaPtVector<TH1F>(vertexDir, "Nvertices", 50, 0., 50.);
  std::vector<std::vector<TH1F*>> vertex_responseBalancing = buildEtaVertexVector<TH1F>(vertexDir, "resp_balancing", 150, 0., 2.);
  std::vector<std::vector<TH1F*>> vertex_responseMPF = buildEtaVertexVector<TH1F>(vertexDir, "resp_mpf", 150, 0., 2.);

  // TProfile
  // response vs Pt (all eta)
  TProfile* Profile_Bal_vs_Pt = balancingDir.make<TProfile>("Bal_vs_Pt", "Balancing vs Pt", binnum, ptBins, 0, 2);
  TProfile* Profile_MPF_vs_Pt = mpfDir.make<TProfile>("MPF_vs_Pt", "MPF vs Pt", binnum, ptBins, 0, 2);
  // response vs Eta (all pT)
  TProfile* Profile_Bal_vs_Eta = balancingDir.make<TProfile>("Bal_vs_Eta", "Balancing vs Eta", binnumEta, etaBins, 0, 2);
  TProfile* Profile_MPF_vs_Eta = mpfDir.make<TProfile>("MPF_vs_Eta", "MPF vs Eta", binnumEta, etaBins, 0, 2);
  // response vs nVtx (all pT and eta)
  TProfile* Profile_Bal_vs_Nvtx = balancingDir.make<TProfile>("Bal_vs_Nvtx", "Balancing vs Nvtx", 10, 0, 40, 0, 2);
  TProfile* Profile_MPF_vs_Nvtx = mpfDir.make<TProfile>("MPF_vs_Nvtx", "MPF vs Nvtx", 10, 0, 40, 0, 2);

  // Phot SC vs Pho
  // |eta| < 1.3 --- for the others? --- vector of TProfile!
  TProfile* Profile_photon_SCPt_vs_Pt = analysisDir.make<TProfile>("PhotSCPt_vs_Pt", "SCPt vs Pt", binnum, ptBins, 0, 2000); 
  TH2F* h_photon_SCPt_vs_Pt = analysisDir.make<TH2F>("PhotonSCPt_vs_Pt", "SCPt vs Pt", binnum, ptBins, 200, 0, 2000);  
 
  //////////////////////////////////////////////////// Extrapolation
  int extrapolationBins = 50;
  double extrapolationMin = 0.;
  double extrapolationMax = 2.;
  TFileDirectory extrapDir = analysisDir.mkdir("extrapolation");
  ExtrapolationVectors<TH1F>::type extrap_responseBalancing = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingRaw = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGen;

  std::vector<std::vector<TH1F*> > extrap_responseBalancingEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingRawEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenEta013;  
  if (mIsMC) {
    extrap_responseBalancingGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  }

  ExtrapolationVectors<TH1F>::type extrap_responseMPF = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseMPFRaw = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseMPFGen;

  std::vector<std::vector<TH1F*> > extrap_responseMPFEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseMPFRawEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseMPFGenEta013;
  if (mIsMC) {
    extrap_responseMPFGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_mpf_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseMPFGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_mpf_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  }    

  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGenPhot;
  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenPhotEta013;
  if (mIsMC) {
    extrap_responseBalancingGenPhot = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_phot", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenPhotEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen_phot", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  }
  
  // Luminosity
   Double_t totalluminosity = 0;
  if (! mIsMC) {
    // For data, there's only one file, so open it in order to read the luminosity
    TFile* f = TFile::Open(mInputFiles[0].c_str());        
    analysisDir.make<TParameter<double>>("luminosity", static_cast<TParameter<double>*>(f->Get("totallumi"))->GetVal());
    f->Close();
    delete f;
  }
  
  // Store alpha cut
  analysisDir.make<TParameter<double>>("alpha_cut", mAlphaCut);
    
  uint64_t totalEvents = FullinfoChain.GetEntries();
  uint64_t passedEvents = 0;
  uint64_t passedEventsFromTriggers = 0;
  uint64_t rejectedEventsFromTriggers = 0;
  uint64_t rejectedEventsTriggerNotFound = 0;
  uint64_t rejectedEventsPtOut = 0;
  //  uint64_t passedEventsFromPhotonRequests = 0; 
  uint64_t passedPhotonJetCut = 0;
  uint64_t passedDeltaPhiCut = 0;
  uint64_t passedPixelSeedVetoCut = 0;
  uint64_t passedMuonsCut = 0;
  uint64_t passedElectronsCut = 0;
  uint64_t passedJetPtCut = 0;
  uint64_t passedAlphaCut = 0;
    
  uint64_t from = 0;
  uint64_t to = totalEvents;
    cout<<"total number of entries  "<<totalEvents<<endl;
  clock::time_point start = clock::now();
    
  // Loop -- from = 0, to = totalEvents
  for (uint64_t i = from; i < to; i++) {     
      
    //test: skip bug events in MC
    // if( i < 3194920 ) continue;
      
    if ( (i - from) < 10 || (i - from) % 50000 == 0) { //50000
      clock::time_point end = clock::now();
      double elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      start = end;
      std::cout << "Processing event #" << (i - from + 1) << " of " << (to - from) << " (" << (float) (i - from) / (to - from) * 100 << "%) - " << elapsedTime << " s" << std::endl;
    }
    
    //bug in crab outputs -- skip events with bugs
    if( mIsMC ){ // bug in GJET Pythia
      if ( i == 475910 || i == 1215426 || i == 2634046 || i == 3016299 || i == 3194920) continue;
    }
    
    if (EXIT) {
      break;
    }
      
     fullinfo.GetEntry(i);
      

      
    // if you want analyze a run range
    //    if( !mIsMC && analysis.run>274315 ) continue;
      
    //skip all events -- usefull to check the crab output
    

    if(mVerbose) std::cout<< std::endl;          
    if(mVerbose) std::cout<<"Event: "<< i << std::endl;          
     
    passedPhotonJetCut++; 
    if(mVerbose)        std::cout<<" passedPhotonJetPresence  " << std::endl;    

  //  if (!(fullinfo.passHLT_Photon30)  && !(fullinfo.passHLT_Photon50) && !(fullinfo.passHLT_Photon75) && !(fullinfo.passHLT_Photon90) && !(fullinfo.passHLT_Photon120) && !(fullinfo.passHLT_Photon165) && !mIsMC) 
  //  {            
 //     rejectedEventsFromTriggers++;
 //     continue;
 //   }
    
    
     int checkTriggerResult = 0;
    std::string passedTrigger;
    float triggerWeight = 1.;
    if(mVerbose) std::cout<<"Finding trigger... " << std::endl;
    if ((checkTriggerResult = checkTriggerfulltree(passedTrigger, fullinfo.passHLT_Photon30, fullinfo.passHLT_Photon50, fullinfo.passHLT_Photon75, fullinfo.passHLT_Photon90, fullinfo.passHLT_Photon120, fullinfo.passHLT_Photon165, triggerWeight)) != TRIGGER_OK) {
      switch (checkTriggerResult) {
      case TRIGGER_NOT_FOUND:
	if (mVerbose) {
	  std::cout << MAKE_RED << "[Run #" << fullinfo.run << ", pT: " << fullinfo.Pt_photon << "] Event does not pass required trigger. List of passed triggers: " << RESET_COLOR << std::endl;
	//  size_t size = analysis.trigger_names->size();
	//  for (size_t i = 0; i < size; i++) {
	 //   if (analysis.trigger_results->at(i)) {
	//      std::cout << "\t" << analysis.trigger_names->at(i) << std::endl;
	///    }
	///  }
	}
	rejectedEventsTriggerNotFound++;
	break;
      case TRIGGER_FOUND_BUT_PT_OUT:
	if (mVerbose) {
	  std::cout << MAKE_RED << "[Run #" << fullinfo.run << ", pT: " << fullinfo.Pt_photon << "] Event does pass required trigger, but pT is out of range. List of passed triggers: " << RESET_COLOR << std::endl;
	 // size_t size = analysis.trigger_names->size();
	 // for (size_t i = 0; i < size; i++) {
	  //  if (analysis.trigger_results->at(i)) {
	   //   std::cout << "\t" << analysis.trigger_names->at(i) <<  std::endl;
	  //  }
	 // }
	}
	rejectedEventsPtOut++;
	break;
      }
      rejectedEventsFromTriggers++;
      continue;
}
    

    passedEventsFromTriggers++;
        
    if(mVerbose)    std::cout<<" passedEventFromTriggers  " << std::endl;
    
    if (mIsMC) 
    {	
      // PU reweighting -- 2016 official recipe
      computePUWeight();
	
      // N vertex based (if official recipe not available)
      // computePUWeight_NVtxBased(photon.pt, analysis.nvertex);      
    } else { // IsData
      // Note: if trigger bins == pT bins the prescale is not important
      // to calculate the response (they are normalized to shape)
	
      //  OLD: get prescale from xml file
      //  in the xml file there were the inverse of prescale
      //  triggerWeight = 1. / triggerWeight; 

      //  NEW: get trigger prescale from ntupla
      triggerWeight = triggerWeight;
      if(mVerbose)  std::cout<< triggerWeight << std::endl;
    }
      
    // Weights 
    double generatorWeight = (mIsMC) ? fullinfo.weight : 1.;
    if (generatorWeight == 0.)
      generatorWeight = 1.;
    double evtWeightSum = (mIsMC) ? fullinfo.evtWeightTot : 1.;

    if (evtWeightSum == 0.)
      evtWeightSum = 1.;
    double eventWeight = (mIsMC) ? mPUWeight * generatorWeight * evtWeightSum : triggerWeight;
      
    if(mVerbose){
      if( mIsMC){
	std::cout << "generatorWeight   "<< generatorWeight << std::endl; 
	std::cout << "evtWeightTot   "<<evtWeightSum << std::endl; 
	std::cout << "mPUWeight     "<< mPUWeight << std::endl; 
      }else{
	std::cout<< "triggerWeight    "<< triggerWeight << std::endl;
      }
      std::cout<< "Final used weight   "<< eventWeight << std::endl;
    }

    h_ptPhoton_NoCut -> Fill(fullinfo.Pt_photon, eventWeight);

    double mu;
    if(mIsMC){
      mu = fullinfo.trueInteraction ;
    } else {
      mu = getAvgPU( fullinfo.run, fullinfo.lumi );
    }
    
    // Event selection
    // The photon is "Good" from previous step (Filter)
    // From previous step: fabs(deltaPhi(photon, firstJet)) > PI/2

    double deltaPhi = fabs(fullinfo.deltaPHIgj);
    h_deltaPhi_NoCut -> Fill(deltaPhi, eventWeight);      
   /* bool isBack2Back = (deltaPhi >= DELTAPHI_CUT);
    if (! isBack2Back) {
      continue;
    }*/
      
    passedDeltaPhiCut++;
    if(mVerbose) std::cout<<"passedDeltaPhiCut"<<std::endl;
      
    // Pixel seed veto
   // if (photon.has_pixel_seed)
    //  continue;
      
    passedPixelSeedVetoCut++;  
    if(mVerbose) std::cout<<"passedPixelSeedVetoCut"<<std::endl;
    
    // No Loose muons 
   // if (muons.nLooseMuon != 0)
   //   continue;

    passedMuonsCut++;
    if(mVerbose) std::cout<<"passedMuonsCut"<<std::endl;
    
    // Electron veto. No electron close to the photon
    bool keepEvent = true;
  //  for (int j = 0; j < electrons.n; j++) {
  //    double deltaR = fabs(reco::deltaR(photon.eta, photon.phi, electrons.eta[j], electrons.phi[j]));
  //    if (deltaR < 0.13) {
  //      keepEvent = false;
  //      break;
  //    }
  //  }
    
    if (! keepEvent)
      continue;

    passedElectronsCut++;    
    if(mVerbose) std::cout<<"passedElectronsCut"<<std::endl;
    
    if (fullinfo.pTAK4_j1 < 15)
      continue;

    passedJetPtCut++;
    if(mVerbose) std::cout<<"passedJetPtCut"<<std::endl;
    
    bool secondJetOK = /* fullinfo.pTAK4_j2!=0 || (fullinfo.pTAK4_j2 < 10 || */ fullinfo.alpha < mAlphaCut ;
    
    //federico -- without this cut the extrapolation is always the same to different alpha cut
    //    if( !secondJetOK) continue;
    //    passedAlphaCut++;
    
    if (secondJetOK)    
      passedAlphaCut++;    
    if(mVerbose) std::cout << "secondJetOK "<< std::endl; 
    
#if ADD_TREES
    if (mUncutTrees) {
      
      fullinfoTree->Fill();

    }
#endif
    
    h_mPUWeight                   ->Fill(mPUWeight);
    h_generatorWeight           ->Fill(generatorWeight);
    h_analysis_evtWeightTot  ->Fill(evtWeightSum);
    h_event_weight_used       ->Fill(eventWeight);
    
    h_nvertex->Fill(fullinfo.nVtx);
    h_ntrue_interactions->Fill(mu);
    h_nvertex_reweighted->Fill(fullinfo.nVtx, eventWeight);
    h_ntrue_interactions_reweighted->Fill(mu, eventWeight);
    
    h_ptPhoton               ->Fill(fullinfo.Pt_photon, eventWeight);
    h_ptPhoton_Binned        ->Fill(fullinfo.Pt_photon, eventWeight);
    h_EtaPhoton             ->Fill(fullinfo.Eta_photon, eventWeight);
    h_PhiPhoton             ->Fill(fullinfo.Phi_photon, eventWeight);
    h_ptFirstJet              ->Fill(fullinfo.pTAK4_j1, eventWeight);
    h_EtaFirstJet            ->Fill(fullinfo.etaAK4_j1, eventWeight);
    h_PhiFirstJet            ->Fill(fullinfo.phiAK4_j1, eventWeight);    
    h_deltaPhi                 ->Fill(deltaPhi, eventWeight); //first jet - photon
    
    h_ptSecondJet          ->Fill(fullinfo.pTAK4_j2, eventWeight);
    h_EtaSecondJet        ->Fill(fullinfo.etaAK4_j2, eventWeight);
    h_PhiSecondJet        ->Fill(fullinfo.phiAK4_j2, eventWeight);
    double deltaPhi_2ndJet = fabs(fullinfo.Phi_photon - fullinfo.phiAK4_j2);
    h_deltaPhi_2ndJet     ->Fill(deltaPhi_2ndJet, eventWeight); //2nd jet - photon
    
    h_alpha                     ->Fill(fullinfo.alpha, eventWeight);
    h_MET                      ->Fill(fullinfo.MET, eventWeight);
    h_rho                          ->Fill(fullinfo.rho, eventWeight);
    h_hadTowOverEm      ->Fill(fullinfo.hadTowOverEm, eventWeight);
    h_sigmaIetaIeta          ->Fill(fullinfo.sigmaietaieta_photon, eventWeight);
    h_chargedHadronsIsolation    ->Fill(fullinfo.CHiso_photon, eventWeight);
    h_neutralHadronsIsolation     ->Fill(fullinfo.NHiso_photon, eventWeight);
    h_photonIsolation                  ->Fill(fullinfo.Photoniso_photon, eventWeight);
    h_deltaPhi_vs_alpha    ->Fill(deltaPhi,fullinfo.alpha,eventWeight);
   

    if(mVerbose){
    std::cout<<"photon phi = "<< fullinfo.Phi_photon << std::endl;
    std::cout<<"MET phi = "<< MET.phi << std::endl;
    std::cout<<"deltaPhi PhotMET = "<< deltaPhi_Photon_MET << std::endl;
    std::cout<<"photon pT = "<< fullinfo.Pt_photon << std::endl;
    std::cout<<"MET = "<< fullinfo.MET << std::endl;
    std::cout<<"resp MPF = "<< fullinfo.RMPF << std::endl;
    }

   // deltaPhi_Photon_MET_raw = reco::deltaPhi(photon.phi, rawMET.phi);
    respMPFRaw = fullinfo.RMPFRAW;
    
    if ( mIsMC){
    //  deltaPhi_Photon_MET_gen = reco::deltaPhi(genPhoton.phi, genMET.phi);
      respMPFGen = fullinfo.RMPF;
    } // true MPF response

    // Balancing
    respBalancing = fullinfo.Rbalancing;
  //  respBalancingRaw = firstRawJet.pt / photon.pt;
    if( mIsMC && fullinfo.pTAK4_j1GEN!=0)    respBalancingGen = fullinfo.pTAK4_j1 / fullinfo.pTAK4_j1GEN; // true balancing response 
    if( mIsMC )    respGenPhot = fullinfo.pTAK4_j1GEN / fullinfo.Pt_photon; // used to constrain extrapolation fits // no more
    if( mIsMC )    respPhotGamma = fullinfo.Pt_photon / fullinfo.Pt_photonGEN; // to check photon response

    int ptBin = mPtBinning.getPtBin(fullinfo.Pt_photon);
    if (ptBin < 0) {
      if(mVerbose) std::cout << "Photon pt " << fullinfo.Pt_photon << " is not covered by our pt binning. Dumping event." << std::endl;
      continue;
    }
    if ( mIsMC)   ptBinGen = mPtBinning.getPtBin(fullinfo.Pt_photonGEN);
    
    int etaBin = mEtaBinning.getBin(fullinfo.etaAK4_j1);
    if (etaBin < 0) {
      if(mVerbose) std::cout << "Jet Bin " << fullinfo.etaAK4_j1 << " is not covered by our eta binning. Dumping event." << std::endl;
      continue;
    } 
    if ( mIsMC)   etaBinGen = mEtaBinning.getBin(fullinfo.etaAK4_j1GEN);
    
    int vertexBin = mVertexBinning.getVertexBin(fullinfo.nVtx);

    if (secondJetOK) { // ! is_present || pT < 10 || pT < 0.3*pT(pho)
      if(mVerbose) std::cout << "Filling histograms passedID"<< std::endl; 
      do {
        h_ptPhoton_passedID                 -> Fill(fullinfo.Pt_photon, eventWeight);
        h_ptPhoton_passedID_Binned    -> Fill(fullinfo.Pt_photon, eventWeight);
	h_EtaPhoton_passedID               -> Fill(fullinfo.Eta_photon, eventWeight);
	h_PhiPhoton_passedID               -> Fill(fullinfo.Phi_photon, eventWeight);
        h_rho_passedID                          -> Fill(fullinfo.rho, eventWeight);
        h_hadTowOverEm_passedID      -> Fill(fullinfo.hadTowOverEm, eventWeight);
        h_sigmaIetaIeta_passedID          -> Fill(fullinfo.sigmaietaieta_photon, eventWeight);
        h_chargedHadronsIsolation_passedID     -> Fill(fullinfo.CHiso_photon, eventWeight);
        h_neutralHadronsIsolation_passedID      -> Fill(fullinfo.NHiso_photon, eventWeight);
        h_photonIsolation_passedID                   -> Fill(fullinfo.Photoniso_photon, eventWeight);

        h_ptFirstJet_passedID       -> Fill(fullinfo.pTAK4_j1, eventWeight);
	h_EtaFirstJet_passedID     -> Fill(fullinfo.etaAK4_j1, eventWeight);
	h_PhiFirstJet_passedID     -> Fill(fullinfo.phiAK4_j1, eventWeight);

        h_deltaPhi_passedID          ->Fill(deltaPhi, eventWeight);

        h_ptSecondJet_passedID       ->Fill(fullinfo.pTAK4_j2, eventWeight);
	h_EtaSecondJet_passedID     ->Fill(fullinfo.etaAK4_j2, eventWeight);
	h_PhiSecondJet_passedID     ->Fill(fullinfo.phiAK4_j2, eventWeight);
	h_PtEtaSecondJet_passedID  ->Fill(fullinfo.etaAK4_j2, fullinfo.pTAK4_j2, eventWeight);	
	if(fullinfo.pTAK4_j2 > 0) {
	  h_ptSecondJet_2ndJetOK       ->Fill(fullinfo.pTAK4_j2, eventWeight);
	  h_EtaSecondJet_2ndJetOK     ->Fill(fullinfo.etaAK4_j2, eventWeight);
	  h_PhiSecondJet_2ndJetOK     ->Fill(fullinfo.phiAK4_j2, eventWeight);
	}
	
        h_alpha_passedID            ->Fill(fullinfo.alpha, eventWeight);
        h_MET_passedID              ->Fill(fullinfo.MET, eventWeight);
        h_rawMET_passedID        ->Fill(fullinfo.METRAW, eventWeight);
        h_METvsfirstJet                   ->Fill(fullinfo.MET, fullinfo.pTAK4_j1, eventWeight);
        h_firstJetvsSecondJet            ->Fill(fullinfo.pTAK4_j1, fullinfo.pTAK4_j2, eventWeight);      	
	h_mu -> Fill(mu, eventWeight);
	h_npvGood -> Fill(fullinfo.nVtx, eventWeight);
	h_rho_vs_mu -> Fill(mu, fullinfo.rho, eventWeight);
	h_npvGood_vs_mu -> Fill(mu,fullinfo.nVtx, eventWeight);
	
	Profile_Bal_vs_Pt     -> Fill(fullinfo.Pt_photon, fullinfo.Rbalancing, eventWeight);
	Profile_MPF_vs_Pt   -> Fill(fullinfo.Pt_photon, fullinfo.RMPF, eventWeight);
	Profile_Bal_vs_Eta   -> Fill(fabs(fullinfo.etaAK4_j1), fullinfo.Rbalancing, eventWeight);
	Profile_MPF_vs_Eta -> Fill(fabs(fullinfo.etaAK4_j1), fullinfo.RMPF, eventWeight);
	Profile_Bal_vs_Nvtx -> Fill(fullinfo.nVtx, fullinfo.Rbalancing, eventWeight);
	Profile_MPF_vs_Nvtx -> Fill(fullinfo.nVtx, fullinfo.RMPF, eventWeight);	

// 	// to be changed with SuperCluster pT
 	if (fabs(fullinfo.etaAK4_j1) <1.305) { //only the special case now
 	  Profile_photon_SCPt_vs_Pt -> Fill(fullinfo.Pt_photon, fullinfo.Pt_photonSC, eventWeight);
 	  h_photon_SCPt_vs_Pt         -> Fill(fullinfo.Pt_photon, fullinfo.Pt_photonSC, eventWeight);
	}
	
	//fill N vertices as a function of eta/pT
	Nvertices[etaBin][ptBin]->Fill(fullinfo.nVtx, eventWeight);		
        if (vertexBin >= 0) {
          vertex_responseBalancing[etaBin][vertexBin] -> Fill(fullinfo.Rbalancing, eventWeight);
          vertex_responseMPF[etaBin][vertexBin]         -> Fill(fullinfo.RMPF, eventWeight);
        }

	//double TotEne = firstRawJet.jet_CHEnF + firstRawJet.jet_NHEnF + firstRawJet.jet_CEmEnF + firstRawJet.jet_NEmEnF + firstRawJet.jet_MuEnF;
	//	if(TotEne > 1 || TotEne < 0) continue;
	//std::cout<< "Tot En = "<< TotEne << std::endl;

	//fill jet energy composition histo vectors
	ChHadronFraction[etaBin][ptBin]->Fill(fullinfo.chargedHadEnFrac_j1, eventWeight);
	NHadronFraction[etaBin][ptBin]->Fill(fullinfo.neutrHadEnFrac_j1, eventWeight);
	CEmFraction[etaBin][ptBin]->Fill(fullinfo.chargedElectromFrac_j1, eventWeight);
	NEmFraction[etaBin][ptBin]->Fill(fullinfo.neutrElectromFrac_j1, eventWeight);
	MuFraction[etaBin][ptBin]->Fill(fullinfo.muEnFract_j1, eventWeight);
	LeptFraction[etaBin][ptBin]->Fill((fullinfo.muEnFract_j1+fullinfo.chargedElectromFrac_j1), eventWeight);

        //fill jet multiplicities histo vectors
        ChHadronMult[etaBin][ptBin]->Fill(fullinfo.chargedMult_j1, eventWeight);
        NHadronMult[etaBin][ptBin]->Fill(fullinfo.neutrMult_j1, eventWeight);
       // ElMult[etaBin][ptBin]->Fill(fullinfo.jet_ElMult, eventWeight);
        PhMult[etaBin][ptBin]->Fill(fullinfo.photonMult_j1, eventWeight);
       // MuonMult[etaBin][ptBin]->Fill(fullinfo.jet_MuonMult, eventWeight);

	//Special case
	if (fabs(firstJet.eta) <1.305) {
	  ChHadronFractionEta013[ptBin]->Fill(fullinfo.chargedHadEnFrac_j1, eventWeight);
	  NHadronFractionEta013[ptBin]->Fill(fullinfo.neutrHadEnFrac_j1, eventWeight);
	  CEmFractionEta013[ptBin]->Fill(fullinfo.chargedElectromFrac_j1, eventWeight);
	  NEmFractionEta013[ptBin]->Fill(fullinfo.neutrElectromFrac_j1, eventWeight);
	  MuFractionEta013[ptBin]->Fill(fullinfo.muEnFract_j1, eventWeight);	    
	  LeptFractionEta013[ptBin]->Fill((fullinfo.muEnFract_j1+fullinfo.chargedElectromFrac_j1), eventWeight);
	  
	  responseBalancingEta013[ptBin]->Fill(fullinfo.Rbalancing, eventWeight);
	  responseBalancingRawEta013[ptBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
	  responseMPFEta013[ptBin]->Fill(fullinfo.RMPF, eventWeight);
	  responseMPFRawEta013[ptBin]->Fill(respMPFRaw, eventWeight);
	}
        responseBalancing[etaBin][ptBin]->Fill(fullinfo.Rbalancing, eventWeight);
        responseBalancingRaw[etaBin][ptBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
        responseMPF[etaBin][ptBin]->Fill(fullinfo.RMPF, eventWeight);
        responseMPFRaw[etaBin][ptBin]->Fill(respMPFRaw, eventWeight);	
	
        // Gen MC values
	if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0) {
	  if (fabs(firstGenJet.eta) <1.305) {
	    responseBalancingGenEta013[ptBinGen] ->Fill(respBalancingGen, eventWeight);
	    responseMPFGenEta013[ptBinGen]        ->Fill(respMPFGen, eventWeight);
	    responseBalancingPhotGammaEta013[ptBinGen]  ->Fill(respPhotGamma, eventWeight);
	    responseBalancingGenPhotEta013[ptBinGen]        ->Fill(respGenPhot, eventWeight);
	  }
	  responseBalancingGen[etaBinGen][ptBinGen]->Fill(respBalancingGen, eventWeight);
	  responseMPFGen[etaBinGen][ptBinGen]->Fill(respMPFGen, eventWeight);
	  responseBalancingPhotGamma[etaBinGen][ptBinGen]  ->Fill(respPhotGamma, eventWeight);
	  responseBalancingGenPhot[etaBinGen][ptBinGen]  ->Fill(respGenPhot, eventWeight);
	}
      } while (false);
      
#if ADD_TREES
      if (! mUncutTrees) {
      
        fullinfoTree->Fill();
      }
#endif
      
      passedEvents++;
      if(mVerbose) std::cout<<"passedEvents"<<std::endl;
      
    }// if secondJetOK

    if (fullinfo.pTAK4_j2 !=0) { //extrapolation if second jet is present
      if(mVerbose) std::cout << "Extrapolating... " << std::endl;
      do {
	
        int extrapBin = mExtrapBinning.getBin(fullinfo.Pt_photon, fullinfo.pTAK4_j2, ptBin);
	if(mIsMC) extrapGenBin = mExtrapBinning.getBin(fullinfo.Pt_photonGEN, fullinfo.pTAK4_j2GEN, ptBin);
	
	do {
          if (extrapBin < 0) {
	    if(mVerbose) std::cout << "No bin found for extrapolation: " << fullinfo.pTAK4_j2 / fullinfo.Pt_photonGEN << std::endl;
	    break;
          }	  
	  
          // Special case 
	  if (fabs(firstJet.eta) < 1.305) {
            extrap_responseBalancingEta013[ptBin][extrapBin]->Fill(fullinfo.Rbalancing, eventWeight);
            extrap_responseMPFEta013[ptBin][extrapBin]->Fill(fullinfo.RMPF, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][extrapBin]->Fill(fullinfo.Rbalancing, eventWeight);
          extrap_responseMPF[etaBin][ptBin][extrapBin]->Fill(fullinfo.RMPF, eventWeight);

	  if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0 && extrapGenBin >= 0) {
	    if (fabs(firstGenJet.eta) < 1.305) {
              extrap_responseBalancingGenEta013[ptBinGen][extrapGenBin] -> Fill(respBalancingGen, eventWeight);
              extrap_responseMPFGenEta013[ptBinGen][extrapGenBin]         -> Fill(respMPFGen, eventWeight);
	      extrap_responseBalancingGenPhotEta013[ptBinGen][extrapGenBin]  -> Fill(respGenPhot, eventWeight);
	    }
	    extrap_responseBalancingGen[etaBinGen][ptBinGen][extrapGenBin] -> Fill(respBalancingGen, eventWeight);
            extrap_responseMPFGen[etaBinGen][ptBinGen][extrapGenBin]         -> Fill(respMPFGen, eventWeight);
	    extrap_responseBalancingGenPhot[etaBinGen][ptBinGen][extrapGenBin]  ->Fill(respGenPhot, eventWeight);
	  }
        } while (false);
	
        int rawExtrapBin = mExtrapBinning.getBin(fullinfo.Pt_photon, fullinfo.pTAK4_j2*(1./fullinfo.jetJecAK4_j2), ptBin); 
	
	
        do {
          if (rawExtrapBin < 0) {
	    if(mVerbose) std::cout << "No bin found for RAW extrapolation: " << fullinfo.pTAK4_j2*(1./fullinfo.jetJecAK4_j2) / fullinfo.Pt_photon << std::endl;
            break;
          }
	  
          // Special case 
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingRawEta013[ptBin][rawExtrapBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
            extrap_responseMPFRawEta013[ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);
          }	  
          extrap_responseBalancingRaw[etaBin][ptBin][rawExtrapBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
          extrap_responseMPFRaw[etaBin][ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);
        } while (false);
	
      } while (false);
      
    } 
     
  }

  std::cout << std::endl;
  std::cout << "Absolute efficiency : related to initial number of event =  " << to-from << std::endl;
  std::cout << "Efficiency for photon/jet cut: " << MAKE_RED << (double) passedPhotonJetCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency for trigger selection: " << MAKE_RED << (double) passedEventsFromTriggers / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  // federico
  //  std::cout << "Selection efficiency for photon requests: " << MAKE_RED << (double) passedEventsFromPhotonRequests / (to - from) * 100 << "%" << RESET_COLOR << std::endl;

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
  // federico
  // std::cout << "Selection efficiency for photon requests: " << MAKE_RED << (double) passedEventsFromPhotonRequests / passedEventsFromTriggers * 100 << "%" << RESET_COLOR << std::endl;
  // std::cout << "Efficiency for Δφ cut: " << MAKE_RED << (double) passedDeltaPhiCut / passedEventsFromPhotonRequests * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for Δφ cut: " << MAKE_RED << (double) passedDeltaPhiCut / passedEventsFromTriggers * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pixel seed veto cut: " << MAKE_RED << (double) passedPixelSeedVetoCut / passedDeltaPhiCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for muons cut: " << MAKE_RED << (double) passedMuonsCut / passedPixelSeedVetoCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for electrons cut: " << MAKE_RED << (double) passedElectronsCut / passedMuonsCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for pT(j1) cut: " << MAKE_RED << (double) passedJetPtCut / passedElectronsCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Efficiency for α cut: " << MAKE_RED << (double) passedAlphaCut / passedJetPtCut * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << std::endl;
  std::cout<< "Histo entries -->    " << passedJetPtCut << std::endl;
  std::cout<< "Histo entries (passedID) -->    " << passedEvents << std::endl;

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



void GammaJetFinalizer::cleanTriggerName(std::string& trigger) {
  boost::replace_first(trigger, "_.*", "");
  boost::replace_first(trigger, ".*", "");
}

//// PU Reweighting
void GammaJetFinalizer::computePUWeight() {

  //  std::cout<< "Using PU Reweighting CODE"<<std::endl;
  
  if (mNoPUReweighting)
    return;
  
  mPUWeight = reweighter->weight(fullinfo.trueInteraction);
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

    TTree* analysis = static_cast<TTree*>(f->Get("rootTupleTree/tree"));
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


// necessary adaptation 
int GammaJetFinalizer::checkTriggerfulltree(std::string& passedTrigger, bool& HLT1, bool& HLT2, bool& HLT3, bool& HLT4, bool& HLT5, bool& HLT6, float& weight) {
if (! mIsMC) {
    const PathVector& mandatoryTriggers = mTriggers->getTriggers(fullinfo.run);
    
    // - With the photon p_t, find the trigger it should pass
    // - Then, look on trigger list if it pass it or not (only for data)
    //    std::cout << "photon.pt  " << photon.pt << std::endl;
    
    const PathData* mandatoryTrigger = nullptr;
    for (auto& path: mandatoryTriggers) {
      if (path.second.range.in(fullinfo.Pt_photon)) {
        mandatoryTrigger = &path;
      }
    }
   
    if (!mandatoryTrigger)
      return TRIGGER_FOUND_BUT_PT_OUT;

    //    weight = mandatoryTrigger->second.weight; // from file
    //    std::cout << "weight  " << weight << std::endl;
    
    // Photon had to pass mandatoryTrigger.first
   // size_t size = fullinfo.trigger_names->size();
 
    
      bool passed1 = HLT1;
      bool passed2 = HLT2;
      bool passed3 = HLT3;
      bool passed4 = HLT4;
      bool passed5 = HLT5;
      bool passed6 = HLT6;
      
      int passedtriggerresult ;
      
      cout<<"test trigger parsing " <<endl;
      
      if (   passed6 == 1 )
       passedtriggerresult = 6 ;
       
       if (   passed5 == 1 )
       passedtriggerresult = 5 ;
       
       if (   passed4 == 1 )
       
        passedtriggerresult = 4 ;
         if (   passed3 == 1 )
       
       passedtriggerresult = 3 ;
        if (   passed2 == 1  )
       
        passedtriggerresult = 2 ;
       
     
        if ( passed1 == 1 ) 
        passedtriggerresult = 1 ;
        

       switch( passedtriggerresult ) { 
         
      case 6:
      if (boost::regex_match("HLT_Photon165_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon165_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
      break;
      
      case 5:
       if (boost::regex_match("HLT_Photon120_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon120_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
      break; 
      
      case 4:
       if (boost::regex_match("HLT_Photon75_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon75_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
       break;
       case 3:
      if (boost::regex_match("HLT_Photon75_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon75_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
      case 2:
       if (boost::regex_match("HLT_Photon50_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon50_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
      break;
      case 1:
      if (boost::regex_match("HLT_Photon30_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = 1; // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << "HLT_Photon30_R9Id90_HE10_IsoM_v.*" << std::endl;
	//if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
	//	std::cout << "Trigger OK"<<std::endl;
	return TRIGGER_OK;
      }
      break;
    }
  } else { // IsMC

    const std::map<Range<float>, std::vector<MCTrigger>>& triggers = mMCTriggers->getTriggers();

    //    std::cout<<photon.pt <<std::endl;
    
    const std::vector<MCTrigger>* mandatoryTrigger = nullptr;
    for (auto& path: triggers) {
      if (path.first.in(fullinfo.Pt_photon)) {
	mandatoryTrigger = &path.second;
      }
    }
    
    if (!mandatoryTrigger)
      return TRIGGER_FOUND_BUT_PT_OUT;;
    
    return TRIGGER_OK;

    
    
  
  

}
 return TRIGGER_NOT_FOUND;
}

int GammaJetFinalizer::checkTrigger(std::string& passedTrigger, float& weight) {

  if (! mIsMC) {
    const PathVector& mandatoryTriggers = mTriggers->getTriggers(fullinfo.run);
    
    // - With the photon p_t, find the trigger it should pass
    // - Then, look on trigger list if it pass it or not (only for data)
    //    std::cout << "photon.pt  " << photon.pt << std::endl;
    
    const PathData* mandatoryTrigger = nullptr;
    for (auto& path: mandatoryTriggers) {
      if (path.second.range.in(fullinfo.Pt_photon)) {
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
      if (! passed)
        continue;
      
      if (boost::regex_match(analysis.trigger_names->at(i), mandatoryTrigger->first)) {
        passedTrigger = mandatoryTrigger->first.str();
	weight = analysis.trigger_prescale->at(i); // prescale from ntupla
	if(mVerbose) std::cout << "Trigger name   " << analysis.trigger_names->at(i) << std::endl;
	if(mVerbose) std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
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
    
    return TRIGGER_OK;

    /*
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
    }*/
    
  }
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
    jetTypes.push_back("puppi");
    TCLAP::ValuesConstraint<std::string> allowedJetTypes(jetTypes);
    TCLAP::ValueArg<std::string> typeArg("", "type", "jet type", true, "pf", &allowedJetTypes, cmd);
    std::vector<std::string> algoTypes;
    algoTypes.push_back("ak4");
    algoTypes.push_back("ak8");
    TCLAP::ValuesConstraint<std::string> allowedAlgoTypes(algoTypes);
    TCLAP::ValueArg<std::string> algoArg("", "algo", "jet algo", true, "ak4", &allowedAlgoTypes, cmd);
    TCLAP::ValueArg<float> chsArg("", "chs", "Use CHS branches", false, true, "bool", cmd);
    // alpha cut
    TCLAP::ValueArg<float> alphaCutArg("", "alpha", "P_t^{second jet} / p_t^{photon} cut (default: 0.3)", false, 0.3, "float", cmd);

    TCLAP::SwitchArg mcArg("", "mc", "MC?", cmd);
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
