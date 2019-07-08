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
#include "parserootPileUpfromJson.h"
#include "parseCSVtoExtractLumiPerHLT.h"
#include <boost/regex.hpp>

#define RESET_COLOR "\033[m"
#define MAKE_RED "\033[31m"
#define MAKE_BLUE "\033[34m"
#define ADD_TREES false
#define DELTAPHI_CUT (2.8)
#define TRIGGER_OK                    0
#define TRIGGER_NOT_FOUND            -1
#define TRIGGER_FOUND_BUT_PT_OUT     -2

boost::shared_ptr<PUReweighter> reweighter30;
boost::shared_ptr<PUReweighter> reweighter50;
boost::shared_ptr<PUReweighter> reweighter75;
boost::shared_ptr<PUReweighter> reweighter90;
boost::shared_ptr<PUReweighter> reweighter120;
boost::shared_ptr<PUReweighter> reweighter165;
boost::shared_ptr<PUReweighter> reweighter200;
TFile* PUFile;
//TFile* EoverP_dataMCRatio_File;
//TH1D *h_test=0;
TFile* EtaPhiCleaning_File;
//TH2D *h_hotjets=0;
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

                //HLTphoton33   

                static std::string puMC30 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();    
                static std::string puData30 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton33%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter30 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData30, puMC30));

                //HLTphoton50


                static std::string puMC50 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();
                static std::string puData50 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton50%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter50 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData50, puMC50));

                //HLTphoton75 


                static std::string puMC75 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();   
                static std::string puData75 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton75%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter75 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData75, puMC75));


                //HLTphoton90 


                static std::string puMC90 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();  
                static std::string puData90 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton90%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter90 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData90, puMC90));

                //HLTphoton120 


                static std::string puMC120 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();
                static std::string puData120 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton120%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter120 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData120, puMC120));

                //HLTphoton165


                static std::string puMC165 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();   
                static std::string puData165 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton165%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter165 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData165, puMC165));

		static std::string puMC200 = TString::Format("%s/computed_mc_files_2018_MC_for_%s_pu_truth_100bins.root", puPrefix.c_str(), mRunera.c_str()).Data();//computed_mc_files_2018_MC_pu_truth_100bins.root", puPrefix.c_str()).Data();   
                static std::string puData200 = TString::Format("%s/pu_truth_data2018_100bins_HLTphoton200%s.root", puPrefix.c_str(), mRunera.c_str()).Data();                                                           
                reweighter200 = boost::shared_ptr<PUReweighter>(new PUReweighter(puData200, puMC200));


                // Trigger
                std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_mc.xml", cmsswBase.c_str()).Data();
                std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
                mMCTriggers      = new MCTriggers( TriggerFile.c_str() ) ;
                
        } else {
                // Get PU
                parsePileUpJSON2();
               // parserootPileUpfromJson();
                parsePrescalefromJson();
                //Trigger
                static std::string cmsswBase = getenv("CMSSW_BASE");
                std::string TriggerFile = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/bin/triggers_data.xml", cmsswBase.c_str()).Data();
                std::cout<< "Trigger File "<< TriggerFile.c_str() << std::endl;
                mTriggers      = new Triggers( TriggerFile.c_str() ) ;

                static std::string Prefix = TString::Format("%s/src/JetMETCorrections/GammaJetFilter/data", cmsswBase.c_str()).Data();
            //    TString EtaPhiCleaning_FileName = TString::Format("%s/hotjets-runBCDEFGH.root", Prefix.c_str()).Data();  //version from Mikko https://github.com/miquork/jecsys/tree/master/rootfiles   
                //EtaPhiCleaning_File = TFile::Open(EtaPhiCleaning_FileName);
               // assert(EtaPhiCleaning_File && !EtaPhiCleaning_File->IsZombie());
               // h_hotjets = (TH2D*)EtaPhiCleaning_File->Get("h2jet"); 
               // assert(h_hotjets);



                

        }

        // Initialization
        mExtrapBinning.initialize(/*mPtBinning, (mJetType == PF) ? "PFlow" : "PUPPI"*/);




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
        std::string outputFile = TString::Format("PhotonJet_%s_%s_%s.root", mDatasetName.c_str(),mRunera.c_str(), postFix.c_str()).Data();
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
  TH1F* h_inst_Lumi = analysisDir.make<TH1F>("Lumi_inst", "luminosity", 100, 0., 2600.);
  TH1F* h_nvertex = analysisDir.make<TH1F>("nvertex", "nvertex", 51, 0., 50.);
  TH1F* h_nvertex_reweighted = analysisDir.make<TH1F>("nvertex_reweighted", "nvertex_reweighted", 51, 0., 50.);
  TH1F* h_ntrue_interactions = analysisDir.make<TH1F>("ntrue_interactions", "ntrue_interactions", 76, 0., 75.);
  TH1F* h_ntrue_interactions_reweighted = analysisDir.make<TH1F>("ntrue_interactions_reweighted", "ntrue_interactions_reweighted", 76, 0., 75.);

  TH1F* h_mPUWeight = analysisDir.make<TH1F>("mPUWeight", "mPUWeight", 50, 0., 5.);
  TH1F* h_generatorWeight = analysisDir.make<TH1F>("generatorWeight", "generatorWeight", 50, 0., 1.);
  TH1F* h_analysis_evtWeightTot = analysisDir.make<TH1F>("analysis_evtWeightTot", "analysis_evtWeightTot", 50, 0., 10.);
  TH1F* h_event_weight_used = analysisDir.make<TH1F>("event_weight_used", "event_weight_used", 150, 0., 150.);

  TH1F* h_ptPhoton_NoCut = analysisDir.make<TH1F>("ptPhoton_NoCut", "ptPhoton NoCut", 50, 40., 1500.);

  TH1F* h_ptPhoton = analysisDir.make<TH1F>("ptPhoton", "ptPhoton", 50, 40., 2000.);
  
  // define spécial binning for pt photon plots.
  
  Double_t BINwidth[76];
  
  BINwidth[0] = 40 ;
  BINwidth[1] = 45 ;
  BINwidth[2] = 50 ;
  BINwidth[3] = 55 ;
  BINwidth[4] = 60 ;
  BINwidth[5] = 66.25 ;
  BINwidth[6] = 72.5 ;
  BINwidth[7] = 78.75 ;
  BINwidth[8] = 85 ;
  BINwidth[9] = 90 ;
  BINwidth[10] = 95 ;
  BINwidth[11] = 100 ;
  BINwidth[12] = 105 ;
  BINwidth[13] = 111.25 ;
  BINwidth[14] = 117.5 ;
  BINwidth[15] = 123.75 ;
  BINwidth[16] = 130 ;
  BINwidth[17] = 141.25 ;
  BINwidth[18] = 152.5 ;  
  BINwidth[19] = 163.75 ;
  BINwidth[20] = 175. ;

  double bin_up = 175.;
  for(int ibin = 0 ; ibin < 55 ; ibin ++){
  	bin_up += 15. ;
  	std::cout<<" bin edge "<< bin_up << std::endl;
  	BINwidth[21 + ibin ] = bin_up;
  
  }
  
  h_ptPhoton->SetBins(75, BINwidth);
  
  
  TH1F* h_ptPhoton_1 = analysisDir.make<TH1F>("ptPhoton_130_175", "ptPhoton_130_175", 50, 130., 175.);
  TH1F* h_ptPhoton_2 = analysisDir.make<TH1F>("ptPhoton_105_130", "ptPhoton_105_130", 50, 105., 130.);
  TH1F* h_ptPhoton_3 = analysisDir.make<TH1F>("ptPhoton_85_105", "ptPhoton_85_105", 50, 85., 105.);
  TH1F* h_ptPhoton_4 = analysisDir.make<TH1F>("ptPhoton_60_85", "ptPhoton_60_85", 50, 60., 85.);
  TH1F* h_ptPhoton_5 = analysisDir.make<TH1F>("ptPhoton_40_60", "ptPhoton_40_60", 50, 40., 60.);
  TH1F* h_EtaPhoton = analysisDir.make<TH1F>("EtaPhoton", "EtaPhoton", 50, -5, 5.);
  TH1F* h_PhiPhoton = analysisDir.make<TH1F>("PhiPhoton", "PhiPhoton", 50, -3.5, 3.5);
  TH1F* h_ptFirstJet = analysisDir.make<TH1F>("ptFirstJet", "ptFirstJet", 50, 15., 1000.);
  TH1F* h_EtaFirstJet = analysisDir.make<TH1F>("EtaFirstJet", "EtaFirstJet", 50, -5, 5.);
  TH1F* h_PhiFirstJet = analysisDir.make<TH1F>("PhiFirstJet", "PhiFirstJet", 30, -3.5, 3.5);// put back 60 bins
  TH1F* h_ptSecondJet = analysisDir.make<TH1F>("ptSecondJet", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet = analysisDir.make<TH1F>("EtaSecondJet", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet = analysisDir.make<TH1F>("PhiSecondJet", "PhiSecondJet", 60, -3.5, 3.5);
  TH1F* h_MET = analysisDir.make<TH1F>("MET", "MET", 50, 0., 250.);
  TH1F* h_alpha = analysisDir.make<TH1F>("alpha", "alpha", 50, 0., 2.);

  TH1F* h_deltaPhi_NoCut = analysisDir.make<TH1F>("deltaPhi_NoCut", "deltaPhi before cut", 60, M_PI / 2, M_PI);
  TH1F* h_deltaPhi = analysisDir.make<TH1F>("deltaPhi", "deltaPhi", 60, M_PI / 2, M_PI);
  TH1F* h_deltaPhi_2ndJet = analysisDir.make<TH1F>("deltaPhi_2ndjet", "deltaPhi of 2nd jet", 60, M_PI / 2., M_PI);
  TH1F* h_deltaPhi_Photon_MET = analysisDir.make<TH1F>("deltaPhi_Photon_MEt", "deltaPhi MET", 60, M_PI / 2., M_PI);

  TH1F* h_rho = analysisDir.make<TH1F>("rho", "rho", 50, 0, 50);
  TH1F* h_hadTowOverEm = analysisDir.make<TH1F>("hadTowOverEm", "hadTowOverEm", 50, 0, 0.035);
  TH1F* h_sigmaIetaIeta = analysisDir.make<TH1F>("sigmaIetaIeta", "sigmaIetaIeta", 50, 0, 0.011);
  TH1F* h_chargedHadronsIsolation = analysisDir.make<TH1F>("chargedHadronsIsolation", "chargedHadronsIsolation", 50, 0, 2);
  TH1F* h_neutralHadronsIsolation = analysisDir.make<TH1F>("neutralHadronsIsolation", "neutralHadronsIsolation", 50, 0, 2);//100, 0, 100);
  TH1F* h_photonIsolation = analysisDir.make<TH1F>("photonIsolation", "photonIsolation", 50, 0, 3);
  
  TH2F* h_deltaPhi_vs_alpha = analysisDir.make<TH2F>("deltaPhi_vs_alpha", "deltaPhi_vs_alpha", 50, 0, 4, 50, 0, 0.5);
  TH2F* h_eta_vs_phi = analysisDir.make<TH2F>("eta_phi", "eta_phi", 100, -5., 5., 100, -3.5, 3.5);

  TH1F* h_deltaPhi_passedID = analysisDir.make<TH1F>("deltaPhi_passedID", "deltaPhi", 40, M_PI / 2, M_PI);

  double ptBins[] = {40, 50, 60, 85, 105, 130, 175, 230, 300, 400, 500, 700, 1000, 3000};
  int  binnum = sizeof(ptBins)/sizeof(double) -1;
  double etaBins[] = {0, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
  int  binnumEta = sizeof(etaBins)/sizeof(double) -1;

  TH1F* h_ptPhoton_Binned = new TH1F("ptPhoton_Binned","ptPhoton", binnum, ptBins); 
  TH1F* h_ptPhoton_passedID_Binned = new TH1F("ptPhoton_passedID_Binned","ptPhoton_passedID", binnum, ptBins);  

  TH1F* h_ptPhoton_passedID = analysisDir.make<TH1F>("ptPhoton_passedID", "ptPhoton", 50, 40., 2000.);
  h_ptPhoton_passedID->Sumw2();
  TH1F* h_EtaPhoton_passedID = analysisDir.make<TH1F>("EtaPhoton_passedID", "EtaPhoton", 50, -5, 5.);
  TH1F* h_PhiPhoton_passedID = analysisDir.make<TH1F>("PhiPhoton_passedID", "PhiPhoton", 50, -3.5, 3.5);
  TH1F* h_ptFirstJet_passedID = analysisDir.make<TH1F>("ptFirstJet_passedID", "ptFirstJet", 50, 0., 1500.);
  TH1F* h_EtaFirstJet_passedID = analysisDir.make<TH1F>("EtaFirstJet_passedID", "EtaFirstJet", 50, -5, 5.);
  TH1F* h_PhiFirstJet_passedID = analysisDir.make<TH1F>("PhiFirstJet_passedID", "PhiFirstJet", 50, -3.5, 3.5);
  TH1F* h_ptSecondJet_passedID = analysisDir.make<TH1F>("ptSecondJet_passedID", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet_passedID = analysisDir.make<TH1F>("EtaSecondJet_passedID", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet_passedID = analysisDir.make<TH1F>("PhiSecondJet_passedID", "PhiSecondJet", 60, -3.5, 3.5);
  TH2F* h_PtEtaSecondJet_passedID = analysisDir.make<TH2F>("PtEtaSecondJet_passedID", "Pt vs Eta SecondJet", 50, -5, 5, 50, 0, 1500);

  TH1F* h_ptSecondJet_2ndJetOK = analysisDir.make<TH1F>("ptSecondJet_2ndJetOK", "ptSecondJet", 20, 0., 200.);
  TH1F* h_EtaSecondJet_2ndJetOK = analysisDir.make<TH1F>("EtaSecondJet_2ndJetOK", "EtaSecondJet", 60, -5, 5.);
  TH1F* h_PhiSecondJet_2ndJetOK = analysisDir.make<TH1F>("PhiSecondJet_2ndJetOK", "PhiSecondJet", 60, -3.5, 3.5);

  TH1F* h_MET_passedID = analysisDir.make<TH1F>("MET_passedID", "MET", 50, 0., 250.);
  TH1F* h_rawMET_passedID = analysisDir.make<TH1F>("rawMET_passedID", "raw MET", 75, 0., 300.);
  TH1F* h_alpha_passedID = analysisDir.make<TH1F>("alpha_passedID", "alpha", 100, 0., 2.);
  
    TH1F* h_MET_parr = analysisDir.make<TH1F>("MET_parr", "MET_{//}", 75, 0., 600.);
    TH1F* h_MET_ortho = analysisDir.make<TH1F>("MET_ortho", "MET_{ #perp } ", 75, 0., 600.);

  TH1F* h_mu = analysisDir.make<TH1F>("mu", "mu", 50, 0, 50);
  TH1F* h_npvGood = analysisDir.make<TH1F>("npvGood", "npvGood", 50, 0, 50);
  TH2F* h_rho_vs_mu = analysisDir.make<TH2F>("rho_vs_mu", "Rho vs mu", 50, 0, 50, 100, 0, 50);
  TH2F* h_npvGood_vs_mu = analysisDir.make<TH2F>("npvGood_vs_mu", "npv_good vs mu", 50, 0, 50, 50, 0, 50);
  
  TH2F *EtaPhiJet_afterhot = analysisDir.make<TH2F>("EtaPhioccupency_after", "EtaPhioccupency_after", 100, -5.2, 5.2 ,100, -3.5, 3.5  ) ;
  
 TH2F *EtaPhiphoton_afterhot = analysisDir.make<TH2F>("EtaPhioccupency_photon", "EtaPhioccupency_photon", 100, -5.2, 5.2 ,100, -3.5, 3.5  ) ;
  
  
  
 //flavor composition
        TFileDirectory FcompositionDir = analysisDir.mkdir("flavorcomposition");
        std::vector<std::vector<TH1F*> > fSumEntries = buildEtaPtVector<TH1F>(FcompositionDir, "fSumEntries", 23, -0.5, 22.5);
        std::vector<std::vector<TH1F*> > fSumWeights = buildEtaPtVector<TH1F>(FcompositionDir, "fSumWeights", 23, -0.5, 22.5);
        std::vector<TH1F*> fSumEntries_0013 = buildPtVector<TH1F>(FcompositionDir, "fSumEntries", "eta0013", 23, -0.5, 22.5);
        std::vector<TH1F*> fSumWeights_0013 = buildPtVector<TH1F>(FcompositionDir, "fSumWeights", "eta0013", 23, -0.5, 22.5); 
  
  
  //plots per HLT for control 
  
  //isolation variables : 
  TFileDirectory HLTDirCHiso = analysisDir.mkdir("HLT_CH_iso");
  std::vector<std::vector<TH1F*> > HLTChHadronIso = buildEtaHLTPtVector<TH1F>(HLTDirCHiso, "ChHadronisoHLT", 50, 0., 2.);
  std::vector<TH1F*>      HLTChHadronIsoEta013    = buildHLTPtVector<TH1F>(HLTDirCHiso, "ChHadronisoHLT", "eta0013", 50, 0., 2.);
  
  TFileDirectory HLTDirNHiso = analysisDir.mkdir("HLT_NH_iso");
  std::vector<std::vector<TH1F*> > HLTNhHadronIso = buildEtaHLTPtVector<TH1F>(HLTDirNHiso, "NhHadronisoHLT", 50, 0., 15.);
  std::vector<TH1F*>      HLTNhHadronIsoEta013    = buildHLTPtVector<TH1F>(HLTDirNHiso, "NhHadronisoHLT", "eta0013", 50, 0., 15.);
  
  TFileDirectory HLTDirphoiso = analysisDir.mkdir("HLT_Photon_iso");
  std::vector<std::vector<TH1F*> > HLTPhotonIso = buildEtaHLTPtVector<TH1F>(HLTDirphoiso, "PhotonisoHLT", 30, 0., 4.);
  std::vector<TH1F*>      HLTPhotonIsoEta013    = buildHLTPtVector<TH1F>(HLTDirphoiso, "PhotonisoHLT", "eta0013", 30, 0., 4.);
  
  // selection variables : 
  TFileDirectory HLTDirsigieta = analysisDir.mkdir("HLT_sigieta");
  std::vector<std::vector<TH1F*> > HLTsigieta = buildEtaHLTPtVector<TH1F>(HLTDirsigieta, "sigmaIetaIetaHLT", 100, 0, 0.03);
  std::vector<TH1F*>      HLTsigietaEta013    = buildHLTPtVector<TH1F>(HLTDirsigieta, "sigmaIetaIetaHLT", "eta0013", 100, 0, 0.03);
  
  
  TFileDirectory HLTDirHoverE = analysisDir.mkdir("HLT_HoverE");
  std::vector<std::vector<TH1F*> > HLTHoverE = buildEtaHLTPtVector<TH1F>(HLTDirHoverE, "HoverEHLT", 50, 0, 0.035);
  std::vector<TH1F*>      HLTHoverEEta013    = buildHLTPtVector<TH1F>(HLTDirHoverE, "HoverEHLT", "eta0013", 50, 0, 0.035);
  
  //other : 
  
  TFileDirectory HLTDirrho = analysisDir.mkdir("HLT_rho");
  std::vector<std::vector<TH1F*> > HLTrho = buildEtaHLTPtVector<TH1F>(HLTDirrho, "rhoHLT",100, 0, 50);
  std::vector<TH1F*>      HLTrhoEta013    = buildHLTPtVector<TH1F>(HLTDirrho, "rhoHLT", "eta0013", 100, 0, 50);
  
  
  TFileDirectory HLTDirmetparr = analysisDir.mkdir("HLT_metparr");
  std::vector<std::vector<TH1F*> > HLTmetparr = buildEtaHLTPtVector<TH1F>(HLTDirmetparr, "metparrHLT",50, 0., 300.);
  std::vector<TH1F*>      HLTmetparrEta013    = buildHLTPtVector<TH1F>(HLTDirmetparr, "metparrHLT", "eta0013", 50, 0., 200.);
  
  
  TFileDirectory HLTDirmetperp = analysisDir.mkdir("HLT_metperp");
  std::vector<std::vector<TH1F*> > HLTmetperp = buildEtaHLTPtVector<TH1F>(HLTDirmetperp, "metperpHLT",50, 0., 200.);
  std::vector<TH1F*>      HLTmetperpEta013    = buildHLTPtVector<TH1F>(HLTDirmetperp, "metperpHLT", "eta0013", 50, 0., 200.);
  
  
  
  TFileDirectory HLTDirmet = analysisDir.mkdir("HLT_met");
  std::vector<std::vector<TH1F*> > HLTmet = buildEtaHLTPtVector<TH1F>(HLTDirmet, "metpHLT",50, 0., 300.);
  std::vector<TH1F*>      HLTmetEta013    = buildHLTPtVector<TH1F>(HLTDirmet, "metpHLT", "eta0013", 50, 0., 300.);
  
  TFileDirectory HLTDirNvertex = analysisDir.mkdir("HLT_Nvertex");
  std::vector<std::vector<TH1F*> > HLTNvertex = buildEtaHLTPtVector<TH1F>(HLTDirNvertex, "NvertexHLT",50, 0., 50.);
  std::vector<TH1F*>      HLTNvertexEta013    = buildHLTPtVector<TH1F>(HLTDirNvertex, "NvertexHLT", "eta0013", 50, 0., 50.);
  
  TFileDirectory HLTDiralpha = analysisDir.mkdir("HLT_alpha");
  std::vector<std::vector<TH1F*> > HLTalpha = buildEtaHLTPtVector<TH1F>(HLTDiralpha, "alphaHLT",50, 0., 0.5);
  std::vector<TH1F*>      HLTalphaEta013    = buildHLTPtVector<TH1F>(HLTDiralpha, "alphaHLT", "eta0013", 50, 0., 0.5);
  
  
  TFileDirectory HLTDirHCCalPrecleaning = analysisDir.mkdir("HLT_Hcalpreclean");
  std::vector<std::vector<TH1F*> > HLTphiprecleaning = buildEtaHLTPtVector<TH1F>(HLTDirHCCalPrecleaning, "etaHLTpreclean",50, -4., 0.);
  std::vector<std::vector<TH1F*> > HLTetaprecleaning = buildEtaHLTPtVector<TH1F>(HLTDirHCCalPrecleaning, "phiHLTpreclean",50, 0., 3.);
  
  TFileDirectory HLTDirHCCalPostcleaning = analysisDir.mkdir("HLT_Hcalpostclean");
  std::vector<std::vector<TH1F*> > HLTphipostcleaning = buildEtaHLTPtVector<TH1F>(HLTDirHCCalPostcleaning, "etaHLTpostclean",50, -4., 0.);
  std::vector<std::vector<TH1F*> > HLTetapostcleaning = buildEtaHLTPtVector<TH1F>(HLTDirHCCalPostcleaning, "phiHLTpostclean",50, 0., 3.);

  
  
  //jet variables
  
  
  TFileDirectory HLTDirjetpt = analysisDir.mkdir("HLT_jetpt");
  std::vector<std::vector<TH1F*> > HLTjetpt = buildEtaHLTPtVector<TH1F>(HLTDirjetpt, "jetptHLT",75, 0., 600.);
  std::vector<TH1F*>      HLTjetptEta013    = buildHLTPtVector<TH1F>(HLTDirjetpt, "jetptHLT", "eta0013", 1000, 0., 2000.);
  
  
  std::vector<std::vector<TH1F*> > HLTjet_2pt = buildEtaHLTPtVector<TH1F>(HLTDirjetpt, "jet_2ptHLT",75, 0., 600.);
  std::vector<TH1F*>      HLTjet_2ptEta013    = buildHLTPtVector<TH1F>(HLTDirjetpt, "jet_2ptHLT", "eta0013", 1000, 0., 2000.);
  
  
  
  //end per HLT plots

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
  TH1F* h_photonIsolation_passedID = analysisDir.make<TH1F>("photonIsolation_passedID", "photonIsolation", 100, 0, 3);

  TH2F* h_METvsfirstJet = analysisDir.make<TH2F>("METvsfirstJet", "MET vs firstJet", 150, 0., 300., 150, 0., 500.);
  TH2F* h_firstJetvsSecondJet = analysisDir.make<TH2F>("firstJetvsSecondJet", "firstJet vs secondJet", 60, 5., 100., 60, 5., 100.);
  
  //Ptgamma
  TFileDirectory PtgammaDir = analysisDir.mkdir("Ptgamma");
  std::vector<std::vector<TH1F*> > Ptgamma = buildEtaPtVector<TH1F>(PtgammaDir, "Pt_gamma", 600, 0., 3000.);
  std::vector<TH1F*> PtgammaEta013       = buildPtVector<TH1F>(PtgammaDir, "Pt_gamma", "eta0013", 600, 0., 3000.);
  
  //Pt1st jet
  TFileDirectory Pt1stjetDir = analysisDir.mkdir("Ptfirstjets");
  std::vector<std::vector<TH1F*> > Pt1st = buildEtaPtVector<TH1F>(Pt1stjetDir, "Pt_1stjet", 600, 0., 3000.);
  std::vector<TH1F*> Pt1stEta013       = buildPtVector<TH1F>(Pt1stjetDir, "Pt_1stjet", "eta0013", 600, 0., 3000.);
  //Pt2nd jet
   TFileDirectory Pt2ndjetDir = analysisDir.mkdir("Ptsecondjets");
   std::vector<std::vector<TH1F*> > Pt2nd = buildEtaPtVector<TH1F>(Pt2ndjetDir, "Pt_2ndjet", 600, 0., 3000.);
   std::vector<TH1F*> Pt2ndEta013       = buildPtVector<TH1F>(Pt2ndjetDir, "Pt_2ndjet", "eta0013", 600, 0., 3000.);
  //MET
   TFileDirectory METDir = analysisDir.mkdir("Met");
   std::vector<std::vector<TH1F*> > Met = buildEtaPtVector<TH1F>(METDir, "met", 250, 0., 250.);
   std::vector<TH1F*> MetEta013       = buildPtVector<TH1F>(METDir, "met", "eta0013", 250, 0., 250.); 
   
    TFileDirectory trueinterDir = analysisDir.mkdir("MUDir");
   std::vector<std::vector<TH1F*> > Mu = buildEtaPtVector<TH1F>(trueinterDir, "mu", 100, 0., 100.);
   std::vector<TH1F*> MuEta013       = buildPtVector<TH1F>(trueinterDir, "mu", "eta0013", 101, 0., 100.);
   
   TFileDirectory NverticeDir = analysisDir.mkdir("nvertices");
   std::vector<std::vector<TH1F*> > Nverticesh = buildEtaPtVector<TH1F>(NverticeDir, "nvertices", 100, 0., 100.);
   std::vector<TH1F*> NverticeshEta013       = buildPtVector<TH1F>(NverticeDir, "nvertices", "eta0013", 101, 0., 100.);
  
  
  //plot fine eta binning.
  mfineEtaBinning.Is_Jer_computation(misJER);
  TFileDirectory fine_eta   = analysisDir.mkdir("fine_eta_binning");
  TFileDirectory balancing_fine  = fine_eta.mkdir("balancing");
  TFileDirectory mpf_fine        = fine_eta.mkdir("mpf");
  
  std::vector<std::vector<TH1F*> > responseBalancing_fine = buildfineEtaPtVector<TH1F>(balancing_fine, "resp_balancing_fine_bining", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseMPF_fine = buildfineEtaPtVector<TH1F>(mpf_fine, "resp_mpf_fine_bining", 150, 0., 2.);
  
  
  
  
  // Balancing
  TFileDirectory balancingDir = analysisDir.mkdir("balancing");
  std::vector<std::vector<TH1F*> > responseBalancing = buildEtaPtVector<TH1F>(balancingDir, "resp_balancing", 100, 0., 2.);//to change to 150 bin
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
   
  
  
  //plot versus ETA
  
  TFileDirectory bal_etaDir = analysisDir.mkdir("balancing_vs_eta");  
  std::vector<TH1F*> responseBalancingPt175  = buildEtaVector<TH1F>(bal_etaDir, "resp_balancing_pt175", 150, 0., 2.);
  
  // MPF
  TFileDirectory mpfDir = analysisDir.mkdir("mpf");
  std::vector<std::vector<TH1F*> > responseMPF = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf", 100, 0., 2.);
   std::vector<std::vector<TH1F*> > responseMPFRaw = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_raw", 150, 0., 2.);
  std::vector<std::vector<TH1F*> > responseMPFGen;
  std::vector<TH1F*> responseMPFEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFRawEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_raw", "eta0013", 150, 0., 2.);
  std::vector<TH1F*> responseMPFGenEta013;
  if (mIsMC) {
    responseMPFGen = buildEtaPtVector<TH1F>(mpfDir, "resp_mpf_gen", 150, 0., 2.);
    responseMPFGenEta013 = buildPtVector<TH1F>(mpfDir, "resp_mpf_gen", "eta0013", 150, 0., 2.);
  }
  
  //plot versus ETA
  
  TFileDirectory mpf_etaDir = analysisDir.mkdir("mpf_vs_eta");
  std::vector<TH1F*> responseMPFPt175  = buildEtaVector<TH1F>(mpf_etaDir, "resp_mpf_pt175", 150, 0., 2.);
  
  
  TFileDirectory pli_etaDir = analysisDir.mkdir("PLI_vs_eta");
  std::vector<TH1F*> PLIPt175 ;
  if (mIsMC) {
  PLIPt175  = buildEtaVector<TH1F>(pli_etaDir, "PLI_pt175", 150, 0., 2.);
  }
  // vs number of vertices
  TFileDirectory vertexDir = analysisDir.mkdir("vertex");
  //check Nvtx vs ptphoton
  std::vector<std::vector<TH1F*> > Nvertices = buildEtaPtVector<TH1F>(vertexDir, "Nvertices", 50, 0., 50.);
  std::vector<std::vector<TH1F*>> vertex_responseBalancing = buildEtaVertexVector<TH1F>(vertexDir, "resp_balancing", 150, 0., 2.);
  std::vector<std::vector<TH1F*>> vertex_responseMPF = buildEtaVertexVector<TH1F>(vertexDir, "resp_mpf", 150, 0., 2.);
  
  //vs run number
  
  TFileDirectory Rundir = analysisDir.mkdir("run");
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT30 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT30",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT50 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT50",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT75 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT75",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT90 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT90",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT120 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT120",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseBalancingHLT165 = buildEtaRunVector<TH1F>(Rundir,"resp_balancing_HLT165",100,0.,2.);
  
  std::vector<std::vector<TH1F*>> run_responseMpfHLT30 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT30",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseMpfHLT50 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT50",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseMpfHLT75 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT75",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseMpfHLT90 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT90",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseMpfHLT120 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT120",100,0.,2.);
  std::vector<std::vector<TH1F*>> run_responseMpfHLT165 = buildEtaRunVector<TH1F>(Rundir,"resp_mpf_HLT165",100,0.,2.);


  std::vector<TH1F*> run_responseBalancing_0013HLT30 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT30", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseBalancing_0013HLT50 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT50", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseBalancing_0013HLT75 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT75", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseBalancing_0013HLT90 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT90", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseBalancing_0013HLT120 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT120","eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseBalancing_0013HLT165 = buildRunVector<TH1F>(Rundir,"resp_balancing_HLT165", "eta0013",100,0.,2.);
  
  std::vector<TH1F*> run_responseMpf_0013HLT30 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT30", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseMpf_0013HLT50 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT50", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseMpf_0013HLT75 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT75", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseMpf_0013HLT90 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT90", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseMpf_0013HLT120 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT120", "eta0013",100,0.,2.);
  std::vector<TH1F*> run_responseMpf_0013HLT165 = buildRunVector<TH1F>(Rundir,"resp_mpf_HLT165", "eta0013",100,0.,2.);
  // TProfile
  // response vs Pt (all eta)
  TProfile* Profile_Bal_vs_Pt = balancingDir.make<TProfile>("Bal_vs_Pt", "Balancing vs Pt", binnum, ptBins, 0, 2);
  TProfile* Profile_MPF_vs_Pt = mpfDir.make<TProfile>("MPF_vs_Pt", "MPF vs Pt", binnum, ptBins, 0, 2);
  
  TProfile* Profile_Pt_gamma_vs_Pt = PtgammaDir.make<TProfile>("Ptgamma_vs_Pt", "Ptgamma vs Pt", binnum, ptBins, 0., 3000.);
  TProfile* Profile_Pt_1stjet_vs_Pt = PtgammaDir.make<TProfile>("Pt_1stjet_vs_Pt", "Pt_1stjet vs Pt", binnum, ptBins, 0., 3000.);
  TProfile* Profile_Pt2nd_jet_vs_Pt = PtgammaDir.make<TProfile>("Pt2nd_jet_vs_Pt", "Pt2nd vs Pt", binnum, ptBins, 0., 3000.);
  TProfile* Profile_Met_vs_Pt = PtgammaDir.make<TProfile>("Met_vs_Pt", "Met vs Pt", binnum, ptBins, 0., 3000.);  
  // response vs Eta (all pT)
  TProfile* Profile_Bal_vs_Eta = balancingDir.make<TProfile>("Bal_vs_Eta", "Balancing vs Eta", binnumEta, etaBins, 0, 2);
  TProfile* Profile_MPF_vs_Eta = mpfDir.make<TProfile>("MPF_vs_Eta", "MPF vs Eta", binnumEta, etaBins, 0, 2);
  
  TProfile* Profile_Pt_gamma_vs_Eta = PtgammaDir.make<TProfile>("Ptgamma_vs_Eta", "Ptgamma vs Eta", binnumEta, etaBins, 0., 3000.);
  // response vs nVtx (all pT and eta)
  TProfile* Profile_Bal_vs_Nvtx = balancingDir.make<TProfile>("Bal_vs_Nvtx", "Balancing vs Nvtx", 10, 0, 40, 0, 2);
  TProfile* Profile_MPF_vs_Nvtx = mpfDir.make<TProfile>("MPF_vs_Nvtx", "MPF vs Nvtx", 10, 0, 40, 0, 2);
  
  TProfile* Profile_Pt_gamma_vs_Nvtx = PtgammaDir.make<TProfile>("Ptgamma_vs_Nvtx", "Ptgamma vs Nvtx", 10, 0, 40, 0., 3000.);
  
  // Phot SC vs Pho
  // |eta| < 1.3 --- for the others? --- vector of TProfile!
  TProfile* Profile_photon_SCPt_vs_Pt = analysisDir.make<TProfile>("PhotSCPt_vs_Pt", "SCPt vs Pt", binnum, ptBins, 0, 2000); 
  TH2F* h_photon_SCPt_vs_Pt = analysisDir.make<TH2F>("PhotonSCPt_vs_Pt", "SCPt vs Pt", binnum, ptBins, 200, 0, 2000);  
 
  //////////////////////////////////////////////////// Extrapolation
  int extrapolationBins = 50;
  double extrapolationMin = 0.;
  double extrapolationMax = 2.;
  TFileDirectory extrapDir = analysisDir.mkdir("extrapolation");
  ExtrapolationVectors<TH1F>::type extrap_responseBalancing_finebin = buildExtrapolationfineEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_fine_bining", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseMPF_finebin = buildExtrapolationfineEtaVector<TH1F>(extrapDir, "extrap_resp_mpf_fine_bining", extrapolationBins, extrapolationMin, extrapolationMax);

  ExtrapolationVectors<TH1F>::type extrap_PLI_fine;  
  if (mIsMC) {
    extrap_PLI_fine  = buildExtrapolationfineEtaVector<TH1F>(extrapDir, "extrap_PLI_fine", extrapolationBins, extrapolationMin, extrapolationMax);
  }
  
  ExtrapolationVectors<TH1F>::type extrap_responseBalancing = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingRaw = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", extrapolationBins, extrapolationMin, extrapolationMax);
  ExtrapolationVectors<TH1F>::type extrap_responseBalancingGen;

  std::vector<std::vector<TH1F*> > extrap_responseBalancingEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingRawEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_raw", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
  std::vector<std::vector<TH1F*> > extrap_responseBalancingGenEta013; 
  
  ExtrapolationVectors<TH1F>::type extrap_PLI;
  std::vector<std::vector<TH1F*> > extrap_PLI_Eta013;
  if (mIsMC) {
    extrap_responseBalancingGen = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_responseBalancingGenEta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_resp_balancing_gen", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_PLI  = buildExtrapolationEtaVector<TH1F>(extrapDir, "extrap_PLI", extrapolationBins, extrapolationMin, extrapolationMax);
    extrap_PLI_Eta013 = buildExtrapolationVector<TH1F>(extrapDir, "extrap_PLI", "eta0013", extrapolationBins, extrapolationMin, extrapolationMax);
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
  

  std::vector<TProfile*> Profile_NpvGood_vs_mu = buildPtVector<TProfile>( trueinterDir,"NpvGood_vs_mu",  "eta0013", 80, 0, 80, 0, 80);
  std::vector<TProfile*> Profile_rho_vs_mu     = buildPtVector<TProfile>( trueinterDir,"rho_vs_mu"    ,  "eta0013"   , 80, 0, 80, 0, 80);
  
  // Luminosity
   Double_t totalluminosity = 0;
  if (! mIsMC) {
    // For data, there's only one file, so open it in order to read the luminosity
    TFile* f = TFile::Open(mInputFiles[0].c_str());        
    analysisDir.make<TParameter<float>>("luminosity", static_cast<TParameter<float>*>(f->Get("totallumi"))->GetVal());
    f->Close();
    delete f;
  }
  
  TFile* f = TFile::Open(mInputFiles[0].c_str());
  
  TH2F * Hdeltaphi_all = static_cast<TH2F*>(f->Get("DeltaPhi_vs_alpha"));
  f->Close();
  delete f;
  
 // TH2F* Hdeltaphi_ = analysisDir.make<TH2F>(*Hdeltaphi_all);
  
 
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
  int triggernotzero = 0;
  int nEvent_rejected = 0 ;
    
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
   
    
    if (EXIT) {
      break;
    }
      
     fullinfo.GetEntry(i);
      

      
    // if you want analyze a run range
    //    if( !mIsMC && analysis.run>274315 ) continue;
      
    //skip all events -- usefull to check the crab output
    

    if(mVerbose) std::cout<< std::endl;          
    if(mVerbose) std::cout<<"Event: "<< i << std::endl;  
    
    //if(fullinfo.Pt_photon < 175.) continue ;        
     
    passedPhotonJetCut++; 
    if(mVerbose)        std::cout<<" passedPhotonJetPresence  " << std::endl;    

    if(mEndcaps && fabs(fullinfo.Eta_photon) <= 1.566 ) continue;
    if(fullinfo.passHLT_Photon120) triggernotzero++;
    
     int checkTriggerResult = 0;
    std::string passedTrigger;
    float triggerWeight = 1.;
    if(mVerbose) std::cout<<"Finding trigger... " << std::endl;
    if ((checkTriggerResult = checkTriggerfulltree(passedTrigger, fullinfo.passHLT_Photon30, fullinfo.passHLT_Photon50, fullinfo.passHLT_Photon75, fullinfo.passHLT_Photon90, fullinfo.passHLT_Photon120, fullinfo.passHLT_Photon165, fullinfo.passHLT_Photon200, fullinfo.phomatchHLT_Photon30, fullinfo.phomatchHLT_Photon50, fullinfo.phomatchHLT_Photon75, fullinfo.phomatchHLT_Photon90, fullinfo.phomatchHLT_Photon120, fullinfo.phomatchHLT_Photon165, triggerWeight)) != TRIGGER_OK) {
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
     // triggerWeight = triggerWeight;
      if(mVerbose)  std::cout<< triggerWeight << std::endl;
    }
      
    // Weights
    
     
    double generatorWeight = (mIsMC) ? fullinfo.weight : 1.;
    if (generatorWeight == 0.)
      generatorWeight = 1.;
    double evtWeightSum = (mIsMC) ? fullinfo.evtWeightTotA : 1.;

    if (evtWeightSum == 0.)
      evtWeightSum = 1.;
    double eventWeight = (mIsMC) ? mPUWeight * generatorWeight * evtWeightSum * ComputePreScaleForMC(fullinfo.Pt_photon) : 1.;//triggerWeight;
      
    if(mVerbose){
      if( mIsMC){
        std::cout << "pt photon  "<< fullinfo.Pt_photon << std::endl; 
	std::cout << "generatorWeight   "<< generatorWeight << std::endl; 
	std::cout << "evtWeightTot   "<<evtWeightSum << std::endl; 
	std::cout << "mPUWeight     "<< mPUWeight << std::endl;
	std::cout << "prescale weight "<<ComputePreScaleForMC(fullinfo.Pt_photon)<<std::endl;
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
      mu = -99;//getAvgPUfromlatest( fullinfo.run, fullinfo.lumi );
    }
   // if(mu!=0.)  std::cout<<"test n true interaction : "<<mu<<std::endl;
    
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
  
  //etaphi cleaning
  
  //  if (!mIsMC && h_hotjets->GetBinContent(h_hotjets->FindBin(fullinfo.etaAK4_j1, fullinfo.phiAK4_j1)) > 0) keepEvent=false; // -10 good, +10 bad
   // if( mIsMC && fullinfo.PassGenmatching == 0) continue;
    if (! keepEvent)
      continue;

    passedElectronsCut++;    
    if(mVerbose) std::cout<<"passedElectronsCut"<<std::endl;
    
    if (fullinfo.pTAK4_j1 < 15)
      continue;

    passedJetPtCut++;
    if(mVerbose) std::cout<<"passedJetPtCut"<<std::endl;
    
     
  // if(fullinfo.Isfakephoton != 0 && mIsMC) continue;
    
    bool secondJetOK =  fullinfo.pTAK4_j2==0 || (fullinfo.pTAK4_j2 < 10 ||  fullinfo.alpha < mAlphaCut );
    
  //  std::cout<<" alpha : "<<fullinfo.alpha<<" secondJetOK "<< secondJetOK << std::endl;
    
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
    
    double dataMCRatio = 1;
   

    TLorentzVector PhotonCorr;
    PhotonCorr.SetPtEtaPhiE( (fullinfo.Pt_photon/dataMCRatio), fullinfo.Eta_photon, fullinfo.Phi_photon, (fullinfo.Energy_photon/dataMCRatio)  );
   
    PhotonCorr.SetPtEtaPhiE( (fullinfo.Pt_photon), fullinfo.Eta_photon, fullinfo.Phi_photon, (fullinfo.Energy_photon)  );
    
   // PhotonCorr.SetPtEtaPhiE( (fullinfo.Pt_photonSC), fullinfo.Eta_photonSC, fullinfo.Phi_photonSC, (fullinfo.Energy_photonSC)  );
    
    
    TLorentzVector Photon;
    Photon.SetPtEtaPhiE( fullinfo.Pt_photon, fullinfo.Eta_photon, fullinfo.Phi_photon, fullinfo.Energy_photon );
  //   Photon.SetPtEtaPhiE( (fullinfo.Pt_photonSC), fullinfo.Eta_photonSC, fullinfo.Phi_photonSC, (fullinfo.Energy_photonSC)  );
    
    TLorentzVector met;
    met.SetPtEtaPhiE( fullinfo.MET_Pt, fullinfo.MET_Eta, fullinfo.MET_Phi, fullinfo.MET);

    TLorentzVector METCorr; 
    METCorr = met + Photon - PhotonCorr;
    
    
    TLorentzVector PhotonGen;
    
    PhotonGen.SetPtEtaPhiE(fullinfo.Pt_photonGEN,fullinfo.Eta_photonGEN,fullinfo.Phi_photonGEN,fullinfo.Energy_photonGEN);
    
    TLorentzVector MetGen;
    MetGen.SetPtEtaPhiE( fullinfo.METGEN_Pt, fullinfo.METGEN_Eta, fullinfo.METGEN_Phi, fullinfo.METGEN);
    
    bool skip_event_Hcalveto = false ;
    
    if(!mIsMC){
    if(fullinfo.run >=272007 && fullinfo.run <=  275376 ){
    
    for(size_t ijet = 0 ; ijet < fullinfo.pT_jets->size() ; ++ijet){
    
    

      if(fullinfo.pT_jets->at(ijet) >= 15. && fullinfo.Eta_jets->at(ijet) >= -2.25 && fullinfo.Eta_jets->at(ijet) <= -1.93 && fullinfo.Phi_jets->at(ijet) >= 2.2 && fullinfo.Phi_jets->at(ijet) <= 2.5){

      
      skip_event_Hcalveto = true ;
      break;
      
      }
    }
    
    }
    
    
    if (fullinfo.run >=275657 && fullinfo.run <=  276283 ){
    
    for(size_t ijet = 0 ; ijet < fullinfo.pT_jets->size() ; ++ijet){
    
    
      if(fullinfo.pT_jets->at(ijet) >= 15. && fullinfo.Eta_jets->at(ijet) >= -3.489 && fullinfo.Eta_jets->at(ijet) <= -3.139 && fullinfo.Phi_jets->at(ijet) >= 2.237 && fullinfo.Phi_jets->at(ijet) <= 2.475){
      
      skip_event_Hcalveto = true ;
      break;
      
      }
    }
    
    }
    
    if(fullinfo.run >=276315 && fullinfo.run <=  276811){
    
    for(size_t ijet = 0 ; ijet < fullinfo.pT_jets->size() ; ++ijet){
    
    
      if(fullinfo.pT_jets->at(ijet) >= 15. && fullinfo.Eta_jets->at(ijet) >= -3.60 && fullinfo.Eta_jets->at(ijet) <= -3.139 && fullinfo.Phi_jets->at(ijet) >= 2.237 && fullinfo.Phi_jets->at(ijet) <= 2.475){
      
      skip_event_Hcalveto = true ;
      break;
      
      }
    }
    
    }
    }
 int HLTptBin = mHLTPtBinning.getPtBin(PhotonCorr.Pt()/*fullinfo.Pt_photon*/);   
 int etaBin = mEtaBinning.getBin(fullinfo.etaAK4_j1);   
 
 //   if(PhotonCorr.Pt()>60.) continue;
    
   //if( /* PhotonCorr.Pt() < 130. ||*/ PhotonCorr.Pt() < 175. ) continue;
   // if( fabs(fullinfo.etaAK4_j1) < 2.85 || fabs(fullinfo.etaAK4_j1) > 2.964 ) continue;   
    
    h_mPUWeight                   ->Fill(mPUWeight);
    h_generatorWeight           ->Fill(generatorWeight);
    h_analysis_evtWeightTot  ->Fill(evtWeightSum);
    h_event_weight_used       ->Fill(eventWeight);
    
    h_nvertex->Fill(fullinfo.nVtx);
    
    h_ntrue_interactions->Fill(mu);
    h_nvertex_reweighted->Fill(fullinfo.nVtx, eventWeight);
    h_ntrue_interactions_reweighted->Fill(mu, eventWeight);
    
    h_ptPhoton               ->Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
    h_ptPhoton_Binned        ->Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
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
   

    if(mVerbose && PhotonCorr.Pt() < 60.){
    std::cout<<"photon phi = "<< fullinfo.Phi_photon << std::endl;
    std::cout<<"MET phi = "<< MET.phi << std::endl;
    std::cout<<"deltaPhi PhotMET = "<< deltaPhi_Photon_MET << std::endl;
    std::cout<<"photon pT = "<< PhotonCorr.Pt()/*fullinfo.Pt_photon*/ << std::endl;
    std::cout<<"MET = "<< fullinfo.MET << std::endl;
    std::cout<<"resp MPF = "<< fullinfo.RMPF << std::endl;
    }
   //Met parallel && Met perp: 
   double Met_para = (METCorr.Pt() * PhotonCorr.Pt() * cos(METCorr.DeltaPhi(PhotonCorr))/PhotonCorr.Pt());
   double Met_perp =  std::sqrt(std::pow(METCorr.Px()-(METCorr.Px()*PhotonCorr.Px()/PhotonCorr.Pt()),2)+std::pow(METCorr.Py()-(METCorr.Py()*PhotonCorr.Py()/PhotonCorr.Pt()),2));


   // deltaPhi_Photon_MET_raw = reco::deltaPhi(photon.phi, rawMET.phi);
    respMPFRaw = fullinfo.RMPFRAW;
    respMPF = 1. + METCorr.Pt() * PhotonCorr.Pt() * cos(METCorr.DeltaPhi(PhotonCorr)) / (PhotonCorr.Pt() * PhotonCorr.Pt());
    if ( mIsMC){
    //  deltaPhi_Photon_MET_gen = reco::deltaPhi(genPhoton.phi, genMET.phi);
      respMPFGen = 1. + MetGen.Pt() * PhotonGen.Pt() * cos(MetGen.DeltaPhi(PhotonGen)) / (PhotonGen.Pt() * PhotonGen.Pt());
    } // true MPF response

    // Balancing
    respBalancing = fullinfo.pTAK4_j1 / PhotonCorr.Pt();/*fullinfo.Rbalancing*/;
  //  respBalancingRaw = firstRawJet.pt / photon.pt;
    if( mIsMC && fullinfo.pTAK4_j1GEN!=0)    respBalancingGen = fullinfo.pTAK4_j1 / fullinfo.pTAK4_j1GEN; // true balancing response 
    if( mIsMC )    respGenPhot = fullinfo.pTAK4_j1GEN / PhotonCorr.Pt()/*fullinfo.Pt_photon*/; // used to constrain extrapolation fits // no more
    if( mIsMC )    respPhotGamma = PhotonCorr.Pt()/*fullinfo.Pt_photon*/ / fullinfo.Pt_photonGEN; // to check photon response

    int ptBin = mPtBinning.getPtBin(PhotonCorr.Pt()/*fullinfo.Pt_photon*/);
    
    if (ptBin < 0) {
      if(mVerbose) std::cout << "Photon pt " << PhotonCorr.Pt()/*fullinfo.Pt_photon*/ << " is not covered by our pt binning. Dumping event." << std::endl;
      continue;
    }
    if ( mIsMC)   ptBinGen = mPtBinning.getPtBin(fullinfo.Pt_photonGEN);
    
   // int etaBin = mEtaBinning.getBin(fullinfo.etaAK4_j1);
    int etafineBin = mfineEtaBinning.getBin(fullinfo.etaAK4_j1);
    int runBinning = 0;
    
    if(!mIsMC) runBinning = mRunBinning.getRunBin(fullinfo.run);
    
    
    
    // time dependence studies
    bool dotimedep  = false ;
    if(!mIsMC && dotimedep){
       if(PhotonCorr.Pt() >= 175.) {
    run_responseBalancingHLT165[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT165[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    run_responseBalancing_0013HLT165[runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpf_0013HLT165[runBinning]->Fill(respMPF,eventWeight);
    }
    }
    
   if(PhotonCorr.Pt() < 175. && PhotonCorr.Pt() >130.) {
    run_responseBalancingHLT120[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT120[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    	run_responseBalancing_0013HLT120[runBinning]->Fill(respBalancing,eventWeight);
        run_responseMpf_0013HLT120[runBinning]->Fill(respMPF,eventWeight);
    }
    }
   
   if(PhotonCorr.Pt() < 130. && PhotonCorr.Pt() >105.) {
    run_responseBalancingHLT90[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT90[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    run_responseBalancing_0013HLT90[runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpf_0013HLT90[runBinning]->Fill(respMPF,eventWeight);
    
    }
    }
   
   if(PhotonCorr.Pt() < 105. && PhotonCorr.Pt() >85.) {
    run_responseBalancingHLT75[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT75[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    run_responseBalancing_0013HLT75[runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpf_0013HLT75[runBinning]->Fill(respMPF,eventWeight);
    
    }
    }
    
   if(PhotonCorr.Pt() < 85. && PhotonCorr.Pt() >60.){
    run_responseBalancingHLT50[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT50[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    
    run_responseBalancing_0013HLT50[runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpf_0013HLT50[runBinning]->Fill(respMPF,eventWeight);
    
    }
    } 
   
   if(PhotonCorr.Pt() < 60. && PhotonCorr.Pt() >40.){
    run_responseBalancingHLT30[etaBin][runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpfHLT30[etaBin][runBinning]->Fill(respMPF,eventWeight);
    
    if(fabs(fullinfo.etaAK4_j1) <1.305) {
    
    run_responseBalancing_0013HLT30[runBinning]->Fill(respBalancing,eventWeight);
    run_responseMpf_0013HLT30[runBinning]->Fill(respMPF,eventWeight);
    
    }
    }
    }
    
    if (fabs(fullinfo.etaAK4_j1) <1.305) {
        Profile_NpvGood_vs_mu[ptBin] -> Fill(fullinfo.nVtx,mu,  eventWeight);
	Profile_rho_vs_mu[ptBin] -> Fill(fullinfo.rho, mu , eventWeight);
//if(fullinfo.nVtx > mu && fullinfo.nVtx > 20) std::cout<<"before alpha cut NPV : "<< fullinfo.nVtx << " rho : "<< fullinfo.rho <<" mu "<<mu<<std::endl;
    }
    
    h_eta_vs_phi                        -> Fill(fullinfo.etaAK4_j1,fullinfo.phiAK4_j1,eventWeight);
    
    if (etaBin < 0) {
      if(mVerbose) std::cout << "Jet Bin " << fullinfo.etaAK4_j1 << " is not covered by our eta binning. Dumping event." << std::endl;
      continue;
    } 
    if ( mIsMC)   etaBinGen = mEtaBinning.getBin(fullinfo.etaAK4_j1GEN);
    
    int vertexBin = mVertexBinning.getVertexBin(fullinfo.nVtx);

    if (secondJetOK) { // ! is_present || pT < 10 || pT < 0.3*pT(pho)
    
       HLTphiprecleaning[etaBin][HLTptBin]->Fill(fullinfo.phiAK4_j1, eventWeight);
       HLTetaprecleaning[etaBin][HLTptBin]->Fill(fullinfo.etaAK4_j1, eventWeight);
      if(mVerbose) std::cout << "Filling histograms passedID"<< std::endl; 

        if(skip_event_Hcalveto && !mIsMC){ 
        nEvent_rejected ++;
      //  continue ;
        }

       
      do {
        
      HLTphipostcleaning[etaBin][HLTptBin]->Fill(fullinfo.phiAK4_j1, eventWeight);
      HLTetapostcleaning[etaBin][HLTptBin]->Fill(fullinfo.etaAK4_j1, eventWeight);
      //  std::cout<<" value of alpha  : "<<fullinfo.alpha<<std::endl;
        
        
      /*  std::cout<<" prescale associated to this event "<<std::endl;
        for(int ihlt = 0; ihlt <6 ; ihlt++){
        std::cout<<" prescale for HLT["<<ihlt+1<<"] "<<getPrescaleperHLT(fullinfo.run, fullinfo.lumi,ihlt+1)<<std::endl;
        }*/
        h_inst_Lumi                         -> Fill(fullinfo.lumi, eventWeight);
        h_MET_ortho                         -> Fill(Met_perp, eventWeight);
        h_MET_parr                          -> Fill(Met_para, eventWeight);
        h_ptPhoton_passedID                 -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
        h_ptPhoton_passedID_Binned    -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
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
	
	EtaPhiphoton_afterhot->Fill(fullinfo.Eta_photon,fullinfo.Phi_photon,eventWeight);
	
	if(fullinfo.pTAK4_j2 > 0){
	EtaPhiJet_afterhot->Fill(fullinfo.etaAK4_j2,fullinfo.phiAK4_j2,eventWeight);}
	EtaPhiJet_afterhot->Fill(fullinfo.etaAK4_j1,fullinfo.phiAK4_j1,eventWeight);
	
	Profile_Bal_vs_Pt     -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, respBalancing/*fullinfo.Rbalancing*/, eventWeight);
	Profile_MPF_vs_Pt   -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, respMPF/* fullinfo.RMPF*/, eventWeight);
	
	Profile_Pt_gamma_vs_Pt   -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.Pt_photon, eventWeight);
	Profile_Pt_1stjet_vs_Pt  -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.pTAK4_j1, eventWeight);
	Profile_Pt2nd_jet_vs_Pt   -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.pTAK4_j2, eventWeight); 
	Profile_Met_vs_Pt  -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.MET, eventWeight);
	Profile_Bal_vs_Eta   -> Fill(fabs(fullinfo.etaAK4_j1), respBalancing/*fullinfo.Rbalancing*/, eventWeight);
	Profile_MPF_vs_Eta -> Fill(fabs(fullinfo.etaAK4_j1),respMPF/* fullinfo.RMPF*/, eventWeight);
	
        Profile_Pt_gamma_vs_Eta   -> Fill(fabs(fullinfo.etaAK4_j1), fullinfo.Pt_photon, eventWeight);
        
	Profile_Bal_vs_Nvtx -> Fill(fullinfo.nVtx, respBalancing/*fullinfo.Rbalancing*/, eventWeight);	
	Profile_MPF_vs_Nvtx -> Fill(fullinfo.nVtx,respMPF/* fullinfo.RMPF*/, eventWeight);	
	
        Profile_Pt_gamma_vs_Nvtx   -> Fill(fullinfo.nVtx, PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
        

 	if (fabs(fullinfo.etaAK4_j1) <1.305) { //only the special case now
 	  Profile_photon_SCPt_vs_Pt -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.Pt_photonSC, eventWeight);
 	  h_photon_SCPt_vs_Pt         -> Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, fullinfo.Pt_photonSC, eventWeight);

	}
	
	//fill N vertices as a function of eta/pT
	Nvertices[etaBin][ptBin]->Fill(fullinfo.nVtx, eventWeight);		
        if (vertexBin >= 0) {
          vertex_responseBalancing[etaBin][vertexBin] -> Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          vertex_responseMPF[etaBin][vertexBin]         -> Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
        }

	

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
       
       if (mIsMC) {
                                        // fill histos for flavor fractions information 
                                        if(fullinfo.PDGIDAK4_j1<-900.)
                                        {       
                                                fSumWeights[etaBin][ptBin]->Fill(0., eventWeight);
                                                fSumEntries[etaBin][ptBin]->Fill(0., 1.);
                                        } else {
                                                fSumWeights[etaBin][ptBin]->Fill(fabs(fullinfo.PDGIDAK4_j1), eventWeight);
                                                fSumEntries[etaBin][ptBin]->Fill(fabs(fullinfo.PDGIDAK4_j1), 1.);
                                        }
                                        
                                        if(fabs(fullinfo.etaAK4_j1) <1.305){
                                        
                                        if(fullinfo.PDGIDAK4_j1<-900.)
                                        {       
                                                fSumWeights_0013[ptBin]->Fill(0., eventWeight);
                                                fSumEntries_0013[ptBin]->Fill(0., 1.);
                                        } else {
                                                fSumWeights_0013[ptBin]->Fill(fabs(fullinfo.PDGIDAK4_j1), eventWeight);
                                                fSumEntries_0013[ptBin]->Fill(fabs(fullinfo.PDGIDAK4_j1), 1.);
                                        }
                                        
                                        }
}
       
       //HLT plots;

       HLTChHadronIso[etaBin][HLTptBin]->Fill(fullinfo.CHiso_photon, eventWeight);
       HLTNhHadronIso[etaBin][HLTptBin]->Fill(fullinfo.NHiso_photon, eventWeight);
       HLTPhotonIso[etaBin][HLTptBin]->Fill(fullinfo.Photoniso_photon, eventWeight);
       HLTsigieta[etaBin][HLTptBin]->Fill(fullinfo.sigmaietaieta_photon, eventWeight);
       HLTHoverE[etaBin][HLTptBin]->Fill(fullinfo.hadTowOverEm, eventWeight);
       HLTrho[etaBin][HLTptBin]->Fill(fullinfo.rho, eventWeight);
       HLTmetparr[etaBin][HLTptBin]->Fill(Met_para, eventWeight);
       HLTmetperp[etaBin][HLTptBin]->Fill(Met_perp, eventWeight);
       HLTjetpt[etaBin][HLTptBin]->Fill(fullinfo.pTAK4_j1, eventWeight);
       HLTmet[etaBin][HLTptBin] ->Fill(fullinfo.MET, eventWeight);
       HLTNvertex[etaBin][HLTptBin] ->Fill(fullinfo.nVtx, eventWeight);
       HLTalpha[etaBin][HLTptBin]->Fill(fullinfo.alpha, eventWeight);
       HLTjet_2pt[etaBin][HLTptBin]->Fill(fullinfo.pTAK4_j2, eventWeight);
       // special case 
       if(fabs(fullinfo.etaAK4_j1) <1.305){
        
         HLTChHadronIsoEta013[HLTptBin] ->Fill(fullinfo.CHiso_photon, eventWeight);
         HLTNhHadronIsoEta013[HLTptBin] ->Fill(fullinfo.NHiso_photon, eventWeight);
         HLTPhotonIsoEta013[HLTptBin] ->Fill(fullinfo.Photoniso_photon, eventWeight);
         HLTsigietaEta013[HLTptBin] ->Fill(fullinfo.sigmaietaieta_photon, eventWeight);
         HLTHoverEEta013[HLTptBin] ->Fill(fullinfo.hadTowOverEm, eventWeight);
         HLTrhoEta013[HLTptBin] ->Fill(fullinfo.rho, eventWeight);
         HLTmetparrEta013[HLTptBin] ->Fill(Met_para, eventWeight);
         HLTmetperpEta013[HLTptBin] ->Fill(Met_perp, eventWeight);
         HLTjetptEta013[HLTptBin] ->Fill(fullinfo.pTAK4_j1, eventWeight);
         HLTmetEta013[HLTptBin] ->Fill(fullinfo.MET, eventWeight);
         HLTNvertexEta013[HLTptBin] ->Fill(fullinfo.nVtx, eventWeight);
         HLTalphaEta013[HLTptBin]->Fill(fullinfo.alpha, eventWeight);
         HLTjet_2ptEta013[HLTptBin]->Fill(fullinfo.pTAK4_j2, eventWeight);
       }
  // if(PhotonCorr.Pt() >= 175.) h_ptPhoton ->Fill(PhotonCorr.Pt(), eventWeight);
   if(PhotonCorr.Pt() < 175. && PhotonCorr.Pt() >130.) h_ptPhoton_1->Fill(PhotonCorr.Pt(), eventWeight);
   if(PhotonCorr.Pt() < 130. && PhotonCorr.Pt() >105.) h_ptPhoton_2 ->Fill(PhotonCorr.Pt(), eventWeight);
   if(PhotonCorr.Pt() < 105. && PhotonCorr.Pt() >85.) h_ptPhoton_3 ->Fill(PhotonCorr.Pt(), eventWeight);
   if(PhotonCorr.Pt() < 85. && PhotonCorr.Pt() >60.) h_ptPhoton_4 ->Fill(PhotonCorr.Pt(), eventWeight);
   if(PhotonCorr.Pt() < 60. && PhotonCorr.Pt() >40.) h_ptPhoton_5 ->Fill(PhotonCorr.Pt(), eventWeight);
     
     
     
       
    
    
    
       
       
       //versus eta plots pt > 175  : 
        if( PhotonCorr.Pt() > 175.){
        
        responseBalancingPt175[etaBin]->Fill(respBalancing, eventWeight);
        responseMPFPt175[etaBin]->Fill(respMPF, eventWeight);
         if(mIsMC){
         
             PLIPt175[etaBin]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
         }

        }
  //Special case
                        if (fabs(fullinfo.etaAK4_j1) <1.305) {
                                ChHadronFractionEta013[ptBin]->Fill(fullinfo.chargedHadEnFrac_j1, eventWeight);
                                NHadronFractionEta013[ptBin]->Fill(fullinfo.neutrHadEnFrac_j1, eventWeight);
                                CEmFractionEta013[ptBin]->Fill(fullinfo.chargedElectromFrac_j1, eventWeight);
                                NEmFractionEta013[ptBin]->Fill(fullinfo.neutrElectromFrac_j1, eventWeight);
                                MuFractionEta013[ptBin]->Fill(fullinfo.muEnFract_j1, eventWeight);        
                                LeptFractionEta013[ptBin]->Fill((fullinfo.muEnFract_j1+fullinfo.chargedElectromFrac_j1), eventWeight);

                                responseBalancingEta013[ptBin]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
                                responseBalancingRawEta013[ptBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
                                responseMPFEta013[ptBin]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
                                responseMPFRawEta013[ptBin]->Fill(respMPFRaw, eventWeight);
                                PtgammaEta013[ptBin]->Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
                                Pt1stEta013[ptBin]->Fill(fullinfo.pTAK4_j1, eventWeight);
                                Pt2ndEta013[ptBin]->Fill(fullinfo.pTAK4_j2, eventWeight);
                                MetEta013[ptBin]->Fill(fullinfo.MET, eventWeight);
                                MuEta013[ptBin]->Fill(mu, eventWeight);
                                NverticeshEta013[ptBin]->Fill(fullinfo.nVtx, eventWeight);


                        }

                        responseBalancing[etaBin][ptBin]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
                        responseBalancingRaw[etaBin][ptBin]->Fill((fullinfo.Rbalancing)*(1./fullinfo.jetJecAK4_j1), eventWeight);
                        responseMPF[etaBin][ptBin]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
                        responseMPFRaw[etaBin][ptBin]->Fill(respMPFRaw, eventWeight);   
                        Ptgamma[etaBin][ptBin]->Fill(PhotonCorr.Pt()/*fullinfo.Pt_photon*/, eventWeight);
                        Pt1st[etaBin][ptBin]->Fill(fullinfo.pTAK4_j1, eventWeight);
                        Pt2nd[etaBin][ptBin]->Fill(fullinfo.pTAK4_j2, eventWeight);
                        Met[etaBin][ptBin]->Fill(fullinfo.MET, eventWeight);
                        Mu[etaBin][ptBin]->Fill(mu, eventWeight);
                        Nverticesh[etaBin][ptBin]->Fill(fullinfo.nVtx, eventWeight);

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

      //plot fine eta bin baland mpf
    responseBalancing_fine [etafineBin][ptBin]->Fill(respBalancing,eventWeight);
    responseMPF_fine[etafineBin][ptBin]->Fill(respMPF,eventWeight);
      passedEvents++;
      if(mVerbose) std::cout<<"passedEvents"<<std::endl;
      
    }// if secondJetOK
  //  if(skip_event_Hcalveto) continue ;
    if (fullinfo.pTAK4_j2 > 0.) { //extrapolation if second jet is present -> nor for exclusive alpha binning
      if(mVerbose) std::cout << "Extrapolating... " << std::endl;
      do {
	
        int extrapBin = mExtrapBinning.getBin(fullinfo.Pt_photon, fullinfo.pTAK4_j2, ptBin);
        
        
	if(mIsMC) extrapGenBin = mExtrapBinning.getBin(fullinfo.Pt_photonGEN, fullinfo.pTAK4_j2GEN, ptBin);
	
	
	do {
          if (extrapBin < 0) {
	    if(mVerbose) std::cout << "No bin found for extrapolation: " << fullinfo.pTAK4_j2 / fullinfo.Pt_photon << std::endl;
	    break;
          }	  
	  
	  // Ugly code to make inclusive extrapolation bin for response and PLI
	  
	  //bin 9 0 to 0.3 :
       //bin 9 0 to 0.3 :
       // /*if(fullinfo.pTAK4_j2/fullinfo.Pt_photon > 0.2) */ std::cout<<" value of alpha extrap : "<<fullinfo.pTAK4_j2/fullinfo.Pt_photon<<" extrapBin "<<extrapBin<<std::endl;
	  if((fullinfo.pTAK4_j2/fullinfo.Pt_photon ) < 0.3 /*&& (fullinfo.pTAK4_j2/fullinfo.Pt_photon) >= 0.25*/){
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingEta013[ptBin][4]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
            extrap_responseMPFEta013[ptBin][4]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][4]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF[etaBin][ptBin][4]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          extrap_responseBalancing_finebin[etafineBin][ptBin][4]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF_finebin[etafineBin][ptBin][4]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);

          if (mIsMC){
	    extrap_PLI[etaBin][ptBin][4]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    extrap_PLI_fine[etafineBin][ptBin][4]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    if(fabs(fullinfo.etaAK4_j1) < 1.305){
	    extrap_PLI_Eta013[ptBin][4]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    
	    }}}

	  
	   //bin 9 0 to 0.25 :

	  if((fullinfo.pTAK4_j2/fullinfo.Pt_photon ) < 0.25/* && (fullinfo.pTAK4_j2/fullinfo.Pt_photon) >= 0.2*/){
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingEta013[ptBin][3]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
            extrap_responseMPFEta013[ptBin][3]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][3]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF[etaBin][ptBin][3]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          extrap_responseBalancing_finebin[etafineBin][ptBin][3]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF_finebin[etafineBin][ptBin][3]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          if (mIsMC){
	    extrap_PLI[etaBin][ptBin][3]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    extrap_PLI_fine[etafineBin][ptBin][3]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    if(fabs(fullinfo.etaAK4_j1) < 1.305){
	    extrap_PLI_Eta013[ptBin][3]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    }
	  }}
	  
	  
	  //bin 9 0 to 0.20 :

	  if((fullinfo.pTAK4_j2/ fullinfo.Pt_photon) < 0.2 /*&& (fullinfo.pTAK4_j2/fullinfo.Pt_photon) >= 0.15*/){
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingEta013[ptBin][2]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
            extrap_responseMPFEta013[ptBin][2]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][2]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF[etaBin][ptBin][2]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          extrap_responseBalancing_finebin[etafineBin][ptBin][2]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF_finebin[etafineBin][ptBin][2]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          if (mIsMC){
	    extrap_PLI[etaBin][ptBin][2]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    
	    extrap_PLI_fine[etafineBin][ptBin][2]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    
	    if(fabs(fullinfo.etaAK4_j1) < 1.305){
	    extrap_PLI_Eta013[ptBin][2]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    }
	  }}
	  
	  
	  
	  //bin 9 0 to 0.15 :
	  

	  if((fullinfo.pTAK4_j2/fullinfo.Pt_photon) < 0.15/* && (fullinfo.pTAK4_j2/fullinfo.Pt_photon) >= 0.1*/){
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingEta013[ptBin][1]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
            extrap_responseMPFEta013[ptBin][1]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][1]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF[etaBin][ptBin][1]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          extrap_responseBalancing_finebin[etafineBin][ptBin][1]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF_finebin[etafineBin][ptBin][1]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          if (mIsMC){
	    extrap_PLI[etaBin][ptBin][1]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    
	    extrap_PLI_fine[etafineBin][ptBin][1]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    
	    if(fabs(fullinfo.etaAK4_j1) < 1.305){
	    extrap_PLI_Eta013[ptBin][1]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    }
	  }}
	  
	  //bin 9 0 to 0.1 :

	  if((fullinfo.pTAK4_j2/fullinfo.Pt_photon) < 0.1 /*&& (fullinfo.pTAK4_j2/fullinfo.Pt_photon) >= 0.05 */){
	  if (fabs(fullinfo.etaAK4_j1) < 1.305) {
            extrap_responseBalancingEta013[ptBin][0]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
            extrap_responseMPFEta013[ptBin][0]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
	  }
          extrap_responseBalancing[etaBin][ptBin][0]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF[etaBin][ptBin][0]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          extrap_responseBalancing_finebin[etafineBin][ptBin][0]->Fill(respBalancing/*fullinfo.Rbalancing*/, eventWeight);
          extrap_responseMPF_finebin[etafineBin][ptBin][0]->Fill(respMPF/* fullinfo.RMPF*/, eventWeight);
          
          if (mIsMC){
	    extrap_PLI[etaBin][ptBin][0]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    extrap_PLI_fine[etafineBin][ptBin][0]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    if(fabs(fullinfo.etaAK4_j1) < 1.305){
	    extrap_PLI_Eta013[ptBin][0]->Fill((fullinfo.pTAK4_j1GEN)/(PhotonGen.Pt()), eventWeight);
	    }
}}

	  
	
         
	    

	  if (mIsMC && ptBinGen >= 0 && etaBinGen >= 0 && extrapGenBin >= 0) {
	    if (fabs(fullinfo.etaAK4_j1GEN) < 1.305) {
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
            extrap_responseBalancingRawEta013[ptBin][rawExtrapBin]->Fill((respBalancing/*fullinfo.Rbalancing*/)*(1./fullinfo.jetJecAK4_j1), eventWeight);
            extrap_responseMPFRawEta013[ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);
          }	  
          extrap_responseBalancingRaw[etaBin][ptBin][rawExtrapBin]->Fill((respBalancing/*fullinfo.Rbalancing*/)*(1./fullinfo.jetJecAK4_j1), eventWeight);
          extrap_responseMPFRaw[etaBin][ptBin][rawExtrapBin]->Fill(respMPFRaw, eventWeight);
        } while (false);
	
      } while (false);
      
    }//else{
    
    // put back bin zero of extrap here
    
    
    
   // } 
     
  }

  std::cout << std::endl;
  std::cout << "Absolute efficiency : related to initial number of event =  " << to-from << std::endl;
  //std::cout << "Efficiency for photon/jet cut: " << MAKE_RED << (double) passedPhotonJetCut / (to - from) * 100 << "%" << RESET_COLOR << std::endl;
  std::cout << "Selection efficiency for trigger selection: " << MAKE_RED << (double) passedEventsFromTriggers  << RESET_COLOR << std::endl;

  std::cout << std::endl;
  std::cout<< "Histo entries -->    " << passedJetPtCut << std::endl;
  std::cout<< "Histo entries (passedID) -->    " << passedEvents << std::endl;

  std::cout << std::endl;
  std::cout << "Rejected events because trigger was not found: " << MAKE_RED << (double) rejectedEventsTriggerNotFound  << RESET_COLOR << std::endl;
  std::cout << "Rejected events because pT was out of range: " << MAKE_RED << (double) rejectedEventsPtOut / (rejectedEventsFromTriggers) * 100 << "%" << RESET_COLOR << std::endl;

  std::cout << "Rejected events because HCAL cleaning: " << MAKE_RED << (double) nEvent_rejected << RESET_COLOR << std::endl;

  

}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

        bool appendText = (xMin >= -0.5 && xMax >= 0);
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
std::vector<T*> GammaJetFinalizer::buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax, double yMin, double yMax) {

        bool appendText = (xMin >= -0.5 && xMax >= 0);
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

                T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax, yMin, yMax);
                vector.push_back(object);
        }

        return vector;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax, double yMin, double yMax) {
        return buildPtVector<T>(dir, branchName + "_" + etaName, nBins, xMin, xMax, yMin, yMax);
}




template<typename T>
std::vector<T*> GammaJetFinalizer::buildHLTPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

      //  bool appendText = (xMin >= 0 && xMax >= 0);
        std::vector<T*> vector;
        size_t ptBinningSize = mHLTPtBinning.size();
        for (size_t j = 0; j < ptBinningSize; j++) {

                const std::pair<float, float> bin = mHLTPtBinning.getBinValue(j);
                std::stringstream ss;
             //   if (appendText)
                        ss << branchName << "_HLTptPhot_" << (int) bin.first << "_" << (int) bin.second;
             //   else
             //           ss << branchName << "_" << (int) bin.first << "_" << (int) bin.second;

               /* if (!appendText) {
                        xMin = bin.first;
                }

                if (!appendText) {
                        xMax = bin.second;
                }
*/
                T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
                vector.push_back(object);
        }

        return vector;
}

template<typename T>
std::vector<T*> GammaJetFinalizer::buildHLTPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax) {
        return buildHLTPtVector<T>(dir, branchName + "_" + etaName, nBins, xMin, xMax);
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
std::vector<std::vector<T*> > GammaJetFinalizer::buildfineEtaPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
        size_t etaBinningSize = mfineEtaBinning.size();
        std::vector<std::vector<T*> > etaBinning;

        for (size_t i = 0; i < etaBinningSize; i++) {
                const std::string etaName = mfineEtaBinning.getBinName(i);
                etaBinning.push_back(buildPtVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
        }

        return etaBinning;
}


template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildEtaHLTPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {
        size_t etaBinningSize = mEtaBinning.size();
        std::vector<std::vector<T*> > etaBinning;

        for (size_t i = 0; i < etaBinningSize; i++) {
                const std::string etaName = mEtaBinning.getBinName(i);
                etaBinning.push_back(buildHLTPtVector<T>(dir, branchName, etaName, nBins, xMin, xMax));
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
std::vector<T*> GammaJetFinalizer::buildRunVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax){


        std::vector<T*> vector;
        size_t runBinningSize = mRunBinning.size();
        for (size_t j = 0; j < runBinningSize; j++) {

                const std::pair<int, int> bin = mRunBinning.getBinValue(j);
                std::stringstream ss;
                ss << branchName << "_" << etaName << "_Run_" << bin.first << "_" << bin.second;

                T* object = dir.make<T>(ss.str().c_str(), ss.str().c_str(), nBins, xMin, xMax);
                vector.push_back(object);
        }

        return vector;

}


template<typename T>
std::vector<std::vector<T*> > GammaJetFinalizer::buildEtaRunVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax){

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










template<typename T>
std::vector<std::vector<std::vector<T*> > > GammaJetFinalizer::buildExtrapolationfineEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax) {

        size_t etaBinningSize = mfineEtaBinning.size();
        std::vector<std::vector<std::vector<T*> > > etaBinning;

        for (size_t i = 0; i < etaBinningSize; i++) {
                const std::string etaName = mfineEtaBinning.getBinName(i);

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

        if(fullinfo.Pt_photon >= 40 && fullinfo.Pt_photon < 60)             mPUWeight = reweighter30->weight(fullinfo.trueInteraction);  
        if(fullinfo.Pt_photon >= 60 && fullinfo.Pt_photon < 85)             mPUWeight = reweighter50->weight(fullinfo.trueInteraction);  
        if(fullinfo.Pt_photon >= 85 && fullinfo.Pt_photon < 100)            mPUWeight = reweighter75->weight(fullinfo.trueInteraction);  
        if(fullinfo.Pt_photon >= 100 && fullinfo.Pt_photon < 130)           mPUWeight = reweighter90->weight(fullinfo.trueInteraction);  
        if(fullinfo.Pt_photon >= 130 && fullinfo.Pt_photon < 175)           mPUWeight = reweighter120->weight(fullinfo.trueInteraction);  
        if(fullinfo.Pt_photon >= 175 && fullinfo.Pt_photon < 210 )          mPUWeight = reweighter165->weight(fullinfo.trueInteraction);
        if(fullinfo.Pt_photon >= 210 && fullinfo.Pt_photon < 5000 )         mPUWeight = reweighter200->weight(fullinfo.trueInteraction);

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
int GammaJetFinalizer::checkTriggerfulltree(std::string& passedTrigger, double& HLT1, double& HLT2, double& HLT3, double& HLT4, double& HLT5, double& HLT6, double& HLT7, double& triggHLT1,double& triggHLT2, double& triggHLT3, double& triggHLT4, double& triggHLT5, double& triggHLT6, float& weight) {
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


                double passed1 = HLT1;
                double passed2 = HLT2;
                double passed3 = HLT3;
                double passed4 = HLT4;
                double passed5 = HLT5;
                double passed6 = HLT6;
                double passed7 = HLT7;


                //to change only for 2017 datas
                /* 
                   double passed1 = HLT6;
                   double passed2 = HLT1;
                   double passed3 = HLT2;
                   double passed4 = HLT3;
                   double passed5 = HLT4;
                   double passed6 = HLT5;
                   */ 
                double triggpassed1 = triggHLT1;//1.;//
                double triggpassed2 = triggHLT2;//1.;//
                double triggpassed3 = triggHLT3;//1.;//
                double triggpassed4 = triggHLT4;//1.;//
                double triggpassed5 = triggHLT5;//1.;//
                double triggpassed6 = triggHLT6;//1.;//
               // double triggpassed7 = triggHLT7;

                //   int passedtriggerresult ;

                if (  passed7 == 1. && boost::regex_match("HLT_Photon200_v.*"                  , mandatoryTrigger->first)) return TRIGGER_OK;
                if (triggpassed6 == 1. &&  passed6 == 1. && boost::regex_match("HLT_Photon165_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)){
              //   weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,6);
                 return TRIGGER_OK;}
                if (triggpassed5 == 1. &&  passed5 == 1. && boost::regex_match("HLT_Photon120_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)){
                // weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,5);
                 return TRIGGER_OK;}
                if (triggpassed4 == 1. &&  passed4 == 1. && boost::regex_match("HLT_Photon90_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)){
                 //weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,4);
                 return TRIGGER_OK;}
                if (triggpassed3 == 1. &&  passed3 == 1. && boost::regex_match("HLT_Photon75_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)){
                 //weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,3);
                 return TRIGGER_OK;}
                if (triggpassed2 == 1. &&  passed2 == 1. && boost::regex_match("HLT_Photon50_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->first)){
                 //weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,2);
                 return TRIGGER_OK;}
                if (triggpassed1 == 1. &&  passed1 == 1. && boost::regex_match("HLT_Photon33_v.*", mandatoryTrigger->first)){
                 //weight = getPrescaleperHLT(fullinfo.run, fullinfo.lumi,1);
                 return TRIGGER_OK;}



                return TRIGGER_NOT_FOUND;

        } else { // IsMC

                const std::map<Range<float>, std::vector<MCTrigger>>& triggers = mMCTriggers->getTriggers();

                //    std::cout<<photon.pt <<std::endl;

                const std::vector<MCTrigger>* mandatoryTrigger = nullptr;
                for (auto& path: triggers) {
                        if (path.first.in(fullinfo.Pt_photon)) {
                                mandatoryTrigger = &path.second;
                        }
                }

                if (!mandatoryTrigger) return TRIGGER_FOUND_BUT_PT_OUT;;

                double passed1 = HLT1;
                double passed2 = HLT2;
                double passed3 = HLT3;
                double passed4 = HLT4;
                double passed5 = HLT5;
                double passed6 = HLT6;
                double passed7 = HLT7;

                double triggpassed1 = triggHLT1;//1.;//
                double triggpassed2 = triggHLT2;//1.;//
                double triggpassed3 = triggHLT3;//1.;//
                double triggpassed4 = triggHLT4;//1.;//
                double triggpassed5 = triggHLT5;//1.;//
                double triggpassed6 = triggHLT6;//1.;//
		//double triggpassed7 = triggHLT7;

                //   int passedtriggerresult ;
                std::string p = mandatoryTrigger->at(0).name.str() ;
                std::string H1 = "HLT_Photon165_R9Id90_HE10_IsoM_v.*";
                std::string H2 = "HLT_Photon120_R9Id90_HE10_IsoM_v.*";
                std::string H3 = "HLT_Photon90_R9Id90_HE10_IsoM_v.*";
                std::string H4 = "HLT_Photon75_R9Id90_HE10_IsoM_v.*";
                std::string H5 = "HLT_Photon50_R9Id90_HE10_IsoM_v.*";
                std::string H6 = "HLT_Photon30_R9Id90_HE10_IsoM_v.*";
		// return TRIGGER_OK;
                if (passed7 == 1.  && boost::regex_match("HLT_Photon200_v.*"                , mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed6 == 1. && */  passed6 == 1. && boost::regex_match("HLT_Photon165_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed5 == 1. && */  passed5 == 1. && boost::regex_match("HLT_Photon120_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed4 == 1. && */  passed4 == 1. && boost::regex_match("HLT_Photon90_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed3 == 1. && */  passed3 == 1. && boost::regex_match("HLT_Photon75_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed2 == 1. && */  passed2 == 1. && boost::regex_match("HLT_Photon50_R9Id90_HE10_IsoM_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;
                if (/*triggpassed1 == 1. && */  passed1 == 1. && boost::regex_match("HLT_Photon33_v.*", mandatoryTrigger->at(0).name)) return TRIGGER_OK;

                /* 

                   std::string p = mandatoryTrigger->at(0).name.str() ; 
                   std::cout<<" trigger vector size "<<mandatoryTrigger->size()<<" mandatory trigger first "<< p.c_str() << " pt photon " << fullinfo.Pt_photon <<std::endl;
                   */
                return TRIGGER_NOT_FOUND;

                //return TRIGGER_OK; 




        }
        return TRIGGER_NOT_FOUND;
}
double GammaJetFinalizer::ComputePreScaleForMC(double Pt_photon) {

	double preScale =1.;
	// TODO redefine for 2018
	// Without prescale=1 for Pt_photon<189
	// return preScale;
	if 	(Pt_photon<60) 	{ preScale= 0.0000401153758676; }
	else if (Pt_photon<85) 	{ preScale= 0.0039473720141; }
	else if (Pt_photon<105) { preScale= 0.0156656382257; }
	else if (Pt_photon<130) { preScale= 0.0312899931745; }
	else if (Pt_photon<175) { preScale= 0.125036122867; }
	else if (Pt_photon<230) { preScale= 0.250030962458; }
	else 			{ preScale= 1.; }
	
	/*
	if 	(Pt_photon<60) 	{ preScale= 0.00209; }
	else if (Pt_photon<85) 	{ preScale= 0.00955769; }
	else if (Pt_photon<105) { preScale= 0.04606247; }
	else if (Pt_photon<130) { preScale= 0.090614; }
	else if (Pt_photon<175) { preScale= 0.255642; }
	else 			{ preScale= 0.674888; }
	*/
    
	return preScale;
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
                                //  std::cout << "Trigger OK"<<std::endl;
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
                //  std::cout << "Trigger prescale   " << analysis.trigger_prescale->at(i) << std::endl;
                //  weight = analysis.trigger_prescale->at(i); // prescale from ntupla
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
                        TCLAP::ValueArg<std::string> runarg("r", "runera", "Run era", true, "", "string", cmd);
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
                        TCLAP::SwitchArg endcapArg("", "endcap", "endcap?", cmd);
                        TCLAP::SwitchArg jerArg("", "JER", "JER?", cmd);
                        
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
                        finalizer.setRunera(runarg.getValue());
                        finalizer.setEndcap(endcapArg.getValue());
                        finalizer.setJER(jerArg.getValue());
                        //    if (totalJobsArg.isSet() && currentJobArg.isSet()) {
                        //   finalizer.setBatchJob(currentJobArg.getValue(), totalJobsArg.getValue());
                        //  }

                        finalizer.runAnalysis();

                } catch (TCLAP::ArgException &e) {
                        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
                        return 1;
                }
        }
