// system include files
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <unordered_map>
#include <memory>
#include <string>
#include "Math/GenVector/LorentzVector.h"
// Boost
#include "boost/shared_ptr.hpp"
// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/PtComparator.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/GammaJetFilter/interface/json/json.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include <TParameter.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <boost/regex.hpp>

// class declaration
//
enum JetAlgorithm {
  AK4,
  AK8,
  PUPPI
};

struct JetInfos {
  JetAlgorithm algo;
  edm::EDGetTokenT<pat::JetCollection> inputTag;

};

#define FOREACH(x) for (std::vector<std::string>::const_iterator it = x.begin(); it != x.end(); ++it)

class GammaJetFilter : public edm::EDFilter {
public:
  explicit GammaJetFilter(const edm::ParameterSet&);
  ~GammaJetFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  
  void correctJets(pat::JetCollection& jets, edm::Event& iEvent, const JetAlgorithm algo, const edm::EventSetup& iSetup);
  void extractRawJets(pat::JetCollection& jets);
  void processJets(pat::Photon* photon, pat::JetCollection& jets, const JetAlgorithm algo, edm::Handle<edm::ValueMap<float>>& qgTag, const edm::Handle<pat::JetCollection>& handleForRef, std::vector<TTree*>& trees);
  
  void correctMETWithTypeI(pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, const JetAlgorithm algo, edm::Event& event, pat::Photon& photon);
  void correctMETWithFootprintAndTypeI(pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, const JetAlgorithm algo,  edm::Event& event, pat::Photon& photon);
   
  bool isValidPhotonEB(const pat::PhotonRef& photonRef, edm::Event& event, double generatorWeight);
  bool isValidPhotonEB_DiPhotonSelection(const pat::PhotonRef& photonRef, edm::Event& event, double generatorWeight);
  bool isValidJet(const pat::Jet& jet);
   
  // ----------member data ---------------------------
  bool mVerbose = false;

  bool mIsMC;   
  bool mRedoTypeI; 
  bool mDoFootprint;
  bool mDoJEC;  
  bool mApplyL2Res;  
  bool mApplyL2L3Res;  
  bool mJECFromRaw;
  GreaterByPt<pat::Jet> mSorter;
  
  bool mFirstJetPtCut;
  double mFirstJetThreshold;

  std::vector<std::string> mJetCollections;
  std::map<std::string, JetInfos> mJetCollectionsData;
  
  edm::EDGetTokenT<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > phoChargedIsolationToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoNeutralHadronIsolationToken_; 
  edm::EDGetTokenT<edm::ValueMap<float> > phoPhotonIsolationToken_; 
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonsToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsToken_;
  edm::EDGetTokenT<pat::JetCollection> jetsAK8Token_;
  edm::EDGetTokenT<pat::JetCollection> jetsPUPPIToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::METCollection> metPUPPIToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonsToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUInfoToken_;
  edm::EDGetTokenT<edm::ValueMap<float> > qgToken_; 
  
  // Events Counter
  int Event_Initial =0 ;
  int Event_VtxCut =0 ;
  int Event_NPhotons=0 ;
  int Event_AfterJets =0;
  int Event_Final =0;
  
  int N_Photons_Total = 0;
  int N_Photons_Eta13 = 0;
  int N_Photons_PtRange = 0;
  
  int NoGoodPhotons = 0;
  int OnlyOneGoodPhotons = 0;
  int MoreGoodPhotons =0;
  // Counter inside the photon method
  int N_pho_Initial =0 ;
  int N_pho_genPhoton =0 ;
  int N_pho_HE =0 ;
  int N_pho_R9 =0 ;
  int N_pho_Sigma =0 ;
  int N_pho_ElecVeto =0 ;
  int N_pho_ChargedIsolation =0 ;
  int N_pho_NeutralIsolation =0 ;
  int N_pho_PhotonIsolation =0 ;
  int N_pho_Isolation =0 ;
     
  // Trees
  void createTrees(const std::string& rootName, TFileService& fs);
  TTree* mPhotonTree;
  TTree* mPhotonGenTree;
  TTree* mAnalysisTree;
  TTree* mElectronsTree;
  TTree* mMuonsTree;
  std::map<std::string, std::vector<TTree*> > mJetTrees;
  std::map<std::string, std::vector<TTree*> > mMETTrees;
  std::map<std::string, TTree*>               mMiscTrees;
  TParameter<double>*    mTotalLuminosity;  
  TParameter<bool>*             mJECRedone;
  TParameter<bool>*             mJECFromRawParameter;
  TParameter<bool>*             mFirstJetPtCutParameter;
  TParameter<double>*          mFirstJetThresholdParameter;
  TParameter<long long>* mProcessedEvents;
  TParameter<long long>* mSelectedEvents;   
  double crossSection;
  
  // to store sum of event weights
  TH1F* h_sumW;
  // DEBUG
  TH1F* EventCounter;
  TH1F* EventCounter_Raw;
  TH1F* EventCounterPhoton;
  TH1F* EventCounterPhoton_Raw;
  TH1F* mSelectedFirstJetIndex;
  TH1F* mSelectedSecondJetIndex;
  TH1F* mSelectedFirstJetPhotonDeltaPhi;
  TH1F* mSelectedFirstJetPhotonDeltaR;  
  TH1F* mSelectedSecondJetPhotonDeltaPhi;
  TH1F* mSelectedSecondJetPhotonDeltaR;  

  // For B / C jets neutrinos
  TClonesArray* mNeutrinos;
  TClonesArray* mNeutrinosPDG;
    
  void updateBranch(TTree* tree, void* address, const std::string& name, const std::string& type = "F");
  template<typename U>
  void updateBranch(TTree* tree, std::vector<U>*& address, const std::string& name);
  void updateBranchArray(TTree* tree, void* address, const std::string& name, const std::string& size, const std::string& type = "F");
  void particleToTree(const reco::Candidate* particle, TTree* t, std::vector<boost::shared_ptr<void> >& addresses);  
  void photonToTree(const pat::PhotonRef& photonRef, pat::Photon& photon, const edm::Event& event);
  void metsToTree(const pat::MET& met, const pat::MET& rawMet, const std::vector<TTree*>& trees);
  void metToTree(const pat::MET* met, TTree* tree, TTree* genTree);
  void jetsToTree(const pat::Jet* firstJet, const pat::Jet* secondJet, const std::vector<TTree*>& trees);
  void jetToTree(const pat::Jet* jet, bool findNeutrinos, TTree* tree, TTree* genTree);
  void electronsToTree(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& pv);
  void muonsToTree(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& pv);  
  int getMotherIndex(const edm::Handle<reco::GenParticleCollection>& genParticles, const reco::Candidate* mother);
  void genParticlesToTree(const edm::Handle<reco::GenParticleCollection>& genParticles);

  //FactorizedJetCorrector
  FactorizedJetCorrector *jetCorrector;
  FactorizedJetCorrector *jetCorrectorForTypeI;
  std::vector<JetCorrectorParameters> vPar;
  std::vector<JetCorrectorParameters> vParTypeI;

  FactorizedJetCorrector *PUPPIjetCorrector;
  FactorizedJetCorrector *PUPPIjetCorrectorForTypeI;
  std::vector<JetCorrectorParameters> vParPUPPI;
  std::vector<JetCorrectorParameters> vParTypeIPUPPI;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GammaJetFilter::GammaJetFilter(const edm::ParameterSet& iConfig):
  mIsMC(false),
  full5x5SigmaIEtaIEtaMapToken_(consumes <edm::ValueMap<float> >
				(iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))),
  phoChargedIsolationToken_(consumes <edm::ValueMap<float> >
			    (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
  phoNeutralHadronIsolationToken_(consumes <edm::ValueMap<float> >
				  (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
  phoPhotonIsolationToken_(consumes <edm::ValueMap<float> >
			   (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation"))),
  triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>
			 (iConfig.getParameter<edm::InputTag>("prescalesTag"))),
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResultsTag"))),
  pfToken_(consumes<pat::PackedCandidateCollection>
	   (iConfig.getParameter<edm::InputTag>("pfCands"))),
  generatorToken_(consumes<GenEventInfoProduct>
		  (iConfig.getParameter<edm::InputTag>("generatorTag"))),
  vertexToken_(consumes<reco::VertexCollection>
	       (iConfig.getParameter<edm::InputTag>("vertexTag"))),
  photonsToken_(consumes<pat::PhotonCollection>
		(iConfig.getParameter<edm::InputTag>("photonsTag"))),
  jetsToken_(consumes<pat::JetCollection>
	     (iConfig.getParameter<edm::InputTag>("jetsTag"))),
  jetsAK8Token_(consumes<pat::JetCollection>
		(iConfig.getParameter<edm::InputTag>("jetsAK8Tag"))),
  jetsPUPPIToken_(consumes<pat::JetCollection>
		  (iConfig.getParameter<edm::InputTag>("jetsPUPPITag"))),
  metToken_(consumes<pat::METCollection>
	    (iConfig.getParameter<edm::InputTag>("metTag"))),
  metPUPPIToken_(consumes<pat::METCollection>
		 (iConfig.getParameter<edm::InputTag>("metPUPPITag"))),
  electronsToken_(consumes<pat::ElectronCollection>
		  (iConfig.getParameter<edm::InputTag>("electronsTag"))),
  muonsToken_(consumes<pat::MuonCollection>
	      (iConfig.getParameter<edm::InputTag>("muonsTag"))),
  rhoToken_(consumes<double>
	    (iConfig.getParameter<edm::InputTag>("rhoTag"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo>>
	       (iConfig.getParameter<edm::InputTag>("PUInfoTag")))
{
  
  mIsMC = iConfig.getUntrackedParameter<bool>("isMC", "false");
  qgToken_ = consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"));
  
  mDoJEC             = iConfig.getUntrackedParameter<bool>("doJetCorrection", false);
  mApplyL2Res     = iConfig.getUntrackedParameter<bool>("applyL2Res", false);
  mApplyL2L3Res = iConfig.getUntrackedParameter<bool>("applyL2L3Res", false);
  mRedoTypeI       = iConfig.getUntrackedParameter<bool>("redoTypeIMETCorrection", false);
  mDoFootprint    = iConfig.getUntrackedParameter<bool>("doFootprintMETCorrection", false);
  
  if (mDoJEC){ 
    mJECFromRaw = iConfig.getUntrackedParameter<bool>("correctJecFromRaw", false);
  }
  
  mFirstJetPtCut = iConfig.getUntrackedParameter<bool>("firstJetPtCut", false);
  mFirstJetThreshold = iConfig.getUntrackedParameter<double>("firstJetThreshold", 0.3);
    
  bool runOnPFAK4    = iConfig.getUntrackedParameter<bool>("runOnPFAK4", false);
  bool runOnPFAK8    = iConfig.getUntrackedParameter<bool>("runOnPFAK8", false);
  bool runOnPUPPIAK4    = iConfig.getUntrackedParameter<bool>("runOnPUPPIAK4", false);

  // Load the corrections files
  if (! mIsMC) { // DATA
    // PF AK4 chs
    edm::FileInPath    L1corr_DATA_          = iConfig.getParameter<edm::FileInPath>("L1corr_DATA");
    edm::FileInPath    L2corr_DATA_          = iConfig.getParameter<edm::FileInPath>("L2corr_DATA");
    edm::FileInPath    L3corr_DATA_          = iConfig.getParameter<edm::FileInPath>("L3corr_DATA");
    edm::FileInPath    L1RCcorr_DATA_      = iConfig.getParameter<edm::FileInPath>("L1RCcorr_DATA");
    edm::FileInPath    L2Rescorr_DATA_     = iConfig.getParameter<edm::FileInPath>("L2Rescorr_DATA");
    edm::FileInPath    L2L3Rescorr_DATA_ = iConfig.getParameter<edm::FileInPath>("L2L3Rescorr_DATA");
    JetCorrectorParameters *L1JetPar               = new JetCorrectorParameters(L1corr_DATA_.fullPath());
    JetCorrectorParameters *L2JetPar               = new JetCorrectorParameters(L2corr_DATA_.fullPath());
    JetCorrectorParameters *L3JetPar               = new JetCorrectorParameters(L3corr_DATA_.fullPath());
    JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(L1RCcorr_DATA_.fullPath());
    JetCorrectorParameters *L2ResJetPar         = new JetCorrectorParameters(L2Rescorr_DATA_.fullPath());
    JetCorrectorParameters *ResJetPar             = new JetCorrectorParameters(L2L3Rescorr_DATA_.fullPath());
    // PUPPI jets
    edm::FileInPath    L1PUPPIcorr_DATA_          = iConfig.getParameter<edm::FileInPath>("L1PUPPIcorr_DATA");
    edm::FileInPath    L2PUPPIcorr_DATA_          = iConfig.getParameter<edm::FileInPath>("L2PUPPIcorr_DATA");
    edm::FileInPath    L3PUPPIcorr_DATA_          = iConfig.getParameter<edm::FileInPath>("L3PUPPIcorr_DATA");
    edm::FileInPath    L1RCPUPPIcorr_DATA_      = iConfig.getParameter<edm::FileInPath>("L1RCPUPPIcorr_DATA");
    edm::FileInPath    L2ResPUPPIcorr_DATA_     = iConfig.getParameter<edm::FileInPath>("L2ResPUPPIcorr_DATA");
    edm::FileInPath    L2L3ResPUPPIcorr_DATA_ = iConfig.getParameter<edm::FileInPath>("L2L3ResPUPPIcorr_DATA");
    JetCorrectorParameters *L1PUPPIJetPar               = new JetCorrectorParameters(L1PUPPIcorr_DATA_.fullPath());
    JetCorrectorParameters *L2PUPPIJetPar               = new JetCorrectorParameters(L2PUPPIcorr_DATA_.fullPath());
    JetCorrectorParameters *L3PUPPIJetPar               = new JetCorrectorParameters(L3PUPPIcorr_DATA_.fullPath());
    JetCorrectorParameters *L1PUPPIJetParForTypeI = new JetCorrectorParameters(L1RCPUPPIcorr_DATA_.fullPath());
    JetCorrectorParameters *L2ResPUPPIJetPar         = new JetCorrectorParameters(L2ResPUPPIcorr_DATA_.fullPath());
    JetCorrectorParameters *ResPUPPIJetPar             = new JetCorrectorParameters(L2L3ResPUPPIcorr_DATA_.fullPath());

    // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!    
    // PF AK4 chs
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    if(mApplyL2Res){
      std::cout<<"Applying L2Res" << std::endl ;
      vPar.push_back(*L2ResJetPar);
    }else if(mApplyL2L3Res){
      std::cout<<"Applying L2L3Res" << std::endl ;
      vPar.push_back(*ResJetPar);
    }
    jetCorrector = new FactorizedJetCorrector(vPar);
    //vPar for MET typeI -- only RC
    vParTypeI.push_back(*L1JetParForTypeI);
    jetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeI);
    // PUPPI jets
    vParPUPPI.push_back(*L1PUPPIJetPar);
    vParPUPPI.push_back(*L2PUPPIJetPar);
    vParPUPPI.push_back(*L3PUPPIJetPar);
    if(mApplyL2Res){
      std::cout<<"Applying L2Res" << std::endl ;
      vParPUPPI.push_back(*L2ResPUPPIJetPar);
    }else if(mApplyL2L3Res){
      std::cout<<"Applying L2L3Res" << std::endl ;
      vParPUPPI.push_back(*ResPUPPIJetPar);
    }
    PUPPIjetCorrector = new FactorizedJetCorrector(vParPUPPI);
    //vPar for MET typeI -- only RC
    vParTypeIPUPPI.push_back(*L1PUPPIJetParForTypeI);
    PUPPIjetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeIPUPPI);

    delete L1JetPar;    
    delete L2JetPar;
    delete L3JetPar;
    delete L1JetParForTypeI;
    delete L2ResJetPar;
    delete ResJetPar;
    delete L1PUPPIJetPar;    
    delete L2PUPPIJetPar;
    delete L3PUPPIJetPar;
    delete L1PUPPIJetParForTypeI;
    delete L2ResPUPPIJetPar;
    delete ResPUPPIJetPar;
  } else {  // MC
    // PF AK4 CHS 
    edm::FileInPath L1corr_MC_ = iConfig.getParameter<edm::FileInPath>("L1corr_MC");
    edm::FileInPath L2corr_MC_ = iConfig.getParameter<edm::FileInPath>("L2corr_MC");
    edm::FileInPath L3corr_MC_ = iConfig.getParameter<edm::FileInPath>("L3corr_MC");
    edm::FileInPath L1RCcorr_MC_ = iConfig.getParameter<edm::FileInPath>("L1RCcorr_MC");
    JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(L1corr_MC_.fullPath());
    JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(L2corr_MC_.fullPath());
    JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(L3corr_MC_.fullPath());
    JetCorrectorParameters *L1JetParForTypeI = new JetCorrectorParameters(L1RCcorr_MC_.fullPath());
    // PUPPI jets
    edm::FileInPath L1PUPPIcorr_MC_ = iConfig.getParameter<edm::FileInPath>("L1PUPPIcorr_MC");
    edm::FileInPath L2PUPPIcorr_MC_ = iConfig.getParameter<edm::FileInPath>("L2PUPPIcorr_MC");
    edm::FileInPath L3PUPPIcorr_MC_ = iConfig.getParameter<edm::FileInPath>("L3PUPPIcorr_MC");
    edm::FileInPath L1RCPUPPIcorr_MC_ = iConfig.getParameter<edm::FileInPath>("L1RCPUPPIcorr_MC");
    JetCorrectorParameters *L1PUPPIJetPar = new JetCorrectorParameters(L1PUPPIcorr_MC_.fullPath());
    JetCorrectorParameters *L2PUPPIJetPar = new JetCorrectorParameters(L2PUPPIcorr_MC_.fullPath());
    JetCorrectorParameters *L3PUPPIJetPar = new JetCorrectorParameters(L3PUPPIcorr_MC_.fullPath());
    JetCorrectorParameters *L1PUPPIJetParForTypeI = new JetCorrectorParameters(L1RCPUPPIcorr_MC_.fullPath());

    // Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!!
    // PF AK4 chs
    vPar.push_back(*L1JetPar);
    vPar.push_back(*L2JetPar);
    vPar.push_back(*L3JetPar);
    jetCorrector = new FactorizedJetCorrector(vPar);
    // for MET typeI -- only RC
    vParTypeI.push_back(*L1JetParForTypeI);
    jetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeI);
    // PUPPI jet
    vParPUPPI.push_back(*L1PUPPIJetPar);
    vParPUPPI.push_back(*L2PUPPIJetPar);
    vParPUPPI.push_back(*L3PUPPIJetPar);
    PUPPIjetCorrector = new FactorizedJetCorrector(vParPUPPI);
    //vPar for MET typeI -- only RC
    vParTypeIPUPPI.push_back(*L1PUPPIJetParForTypeI);
    PUPPIjetCorrectorForTypeI = new FactorizedJetCorrector(vParTypeIPUPPI);

    delete L1JetPar;
    delete L2JetPar;
    delete L3JetPar;
    delete L1JetParForTypeI;
    delete L1PUPPIJetPar;
    delete L2PUPPIJetPar;
    delete L3PUPPIJetPar;
    delete L1PUPPIJetParForTypeI;
  }
  
  edm::Service<TFileService> fs;
  mPhotonTree = fs->make<TTree>("photon", "photon tree");
  
  if (mIsMC){
    mPhotonGenTree = fs->make<TTree>("photon_gen", "photon gen tree");
  }  else{
    mPhotonGenTree = nullptr;
  }
  mAnalysisTree = fs->make<TTree>("analysis", "analysis tree");
  mMuonsTree = fs->make<TTree>("muons", "muons tree");
  mElectronsTree = fs->make<TTree>("electrons", "electrons tree");

  mTotalLuminosity = fs->make<TParameter<double> >("total_luminosity", 1.);  

  crossSection = -1.;
  
  if (mIsMC) {
    crossSection = iConfig.getParameter<double>("crossSection");    
  }
  
  if (runOnPFAK4) {
    mJetCollections.push_back("PFlowAK4chs");
    mJetCollectionsData["PFlowAK4chs"] = {AK4, jetsToken_};
  }
  if (runOnPFAK8) {
    mJetCollections.push_back("PFlowAK8chs");
    mJetCollectionsData["PFlowAK8chs"] = {AK8, jetsAK8Token_};
  }
  if (runOnPUPPIAK4) {
    mJetCollections.push_back("PUPPIAK4");
    mJetCollectionsData["PUPPIAK4"] = {PUPPI, jetsPUPPIToken_};
  }
  
  FOREACH(mJetCollections) {
    createTrees(*it, *fs);
  }
  
  mProcessedEvents = fs->make<TParameter<long long> >("total_events", 0);
  mSelectedEvents = fs->make<TParameter<long long> >("passed_events", 0);
  
  mJECRedone = fs->make<TParameter<bool> >("jec_redone", mDoJEC, '*');
  mFirstJetPtCutParameter = fs->make<TParameter<bool> >("cut_on_first_jet_pt", mFirstJetPtCut, '*');
  if (mDoJEC) {
    mJECFromRawParameter = fs->make<TParameter<bool> >("jec_from_raw_jet", mJECFromRaw, '*');
  }
  
  if (mFirstJetPtCut) {
    mFirstJetThresholdParameter = fs->make<TParameter<double> >("cut_on_first_jet_treshold", mFirstJetThreshold);
  }
  
  if (mIsMC){
    // to store the sum of weights
    h_sumW = fs->make<TH1F>("h_sumW", "h_sumW", 1, -0.5, 5.5);
    h_sumW->Sumw2();
  }  else{
    h_sumW = nullptr;
  }
  
  EventCounter = fs->make<TH1F>("EventCounter", "EventCounter", 7, -0.5, 6.5);
  EventCounter_Raw = fs->make<TH1F>("EventCounter_Raw", "EventCounter_Raw", 7, -0.5, 6.5);
  EventCounterPhoton = fs->make<TH1F>("EventCounterPhoton", "EventCounterPhoton", 11, -0.5, 10.5);
  EventCounterPhoton_Raw = fs->make<TH1F>("EventCounterPhoton_Raw", "EventCounterPhoton_Raw", 11, -0.5, 10.5);

  mSelectedFirstJetIndex = fs->make<TH1F>("selectedFirstJetIndex", "selectedFirstJetIndex", 20, 0, 20);
  mSelectedFirstJetPhotonDeltaPhi = fs->make<TH1F>("selectedFirstJetPhotonDeltaPhi", "selectedFirstJetPhotonDeltaPhi", 50, 0., M_PI);
  mSelectedFirstJetPhotonDeltaR = fs->make<TH1F>("selectedFirstJetPhotonDeltaR", "selectedFirstJetPhotonDeltaR", 80, 0, 10);
  mSelectedSecondJetIndex = fs->make<TH1F>("selectedSecondJetIndex", "selectedSecondJetIndex", 20, 0, 20);
  mSelectedSecondJetPhotonDeltaPhi = fs->make<TH1F>("selectedSecondJetPhotonDeltaPhi", "selectedSecondJetPhotonDeltaPhi", 50, 0., M_PI);
  mSelectedSecondJetPhotonDeltaR = fs->make<TH1F>("selectedSecondJetPhotonDeltaR", "selectedSecondJetPhotonDeltaR", 80, 0, 10);
  
  mNeutrinos = NULL;
  mNeutrinosPDG = NULL;
  if (mIsMC) {
    mNeutrinos = new TClonesArray("TLorentzVector", 3);
    mNeutrinosPDG = new TClonesArray("TParameter<int>", 3);
  }
}

GammaJetFilter::~GammaJetFilter()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.) 
  delete mNeutrinos;
  delete mNeutrinosPDG;
}

void GammaJetFilter::createTrees(const std::string& rootName, TFileService& fs) {
  
  TFileDirectory dir = fs.mkdir(rootName);
  std::vector<TTree*>& trees = mJetTrees[rootName];
  
  trees.push_back(dir.make<TTree>("first_jet", "first jet tree"));
  trees.push_back(dir.make<TTree>("second_jet", "second jet tree"));
  
  trees.push_back(dir.make<TTree>("first_jet_raw", "first raw jet tree"));
  trees.push_back(dir.make<TTree>("second_jet_raw", "second raw jet tree"));
  
  if (mIsMC) {
    trees.push_back(dir.make<TTree>("first_jet_gen", "first gen jet tree"));
    trees.push_back(dir.make<TTree>("second_jet_gen", "second gen jet tree"));
  } else {
    trees.push_back(nullptr);
    trees.push_back(nullptr);
  }
  
  // MET
  std::vector<TTree*>& met = mMETTrees[rootName];
  met.push_back(dir.make<TTree>("met", "met tree"));
  met.push_back(dir.make<TTree>("met_raw", "met raw tree"));
  
  if (mIsMC)
    met.push_back(dir.make<TTree>("met_gen", "met gen tree"));
  else
    met.push_back(nullptr);
  
  // Misc
  mMiscTrees[rootName] = dir.make<TTree>("misc", "misc tree");
}

void GammaJetFilter::updateBranch(TTree* tree, void* address, const std::string& name, const std::string& type/* = "F"*/) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), address, std::string(name + "/" + type).c_str()); 
  } else {
    branch->SetAddress(address);
  }
}

template<typename U> void GammaJetFilter::updateBranch(TTree* tree, std::vector<U>*& address, const std::string& name) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), &address); 
  } else {
    branch->SetAddress(&address);
  }
}

void GammaJetFilter::updateBranchArray(TTree* tree, void* address, const std::string& name, const std::string& size, const std::string& type/* = "F"*/) {
  TBranch* branch = tree->GetBranch(name.c_str());
  if (branch == NULL) {
    branch = tree->Branch(name.c_str(), address, std::string(name + "[" + size + "]/" +type).c_str()); 
  } else {
    branch->SetAddress(address);
  }
}

//
// member functions
//
// ------------ method called on each new Event  ------------
bool GammaJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  mProcessedEvents->SetVal(mProcessedEvents->GetVal() + 1);  
  
  double generatorWeight = 1.;
  
  if (mIsMC)
    {
      edm::Handle<GenEventInfoProduct> eventInfos;
      iEvent.getByToken( generatorToken_, eventInfos);
      
      generatorWeight = eventInfos->weight();
      if (generatorWeight == 0.) {
	generatorWeight = 1.;
      }
      h_sumW->Fill(0.,generatorWeight); 
    }
    
  Event_Initial++;
  EventCounter -> AddBinContent(1, generatorWeight );
  
  // Vertex
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken( vertexToken_, vertices);
  
  // Keep events with at least one vertex
  if (!vertices.isValid() || vertices->size() == 0 || vertices->front().isFake())
    return false;
  
  const reco::Vertex& primaryVertex = vertices->at(0);
  
  Event_VtxCut++;
  EventCounter -> AddBinContent(2, generatorWeight );
  
  int nPVGood =0;  
  reco::VertexCollection::const_iterator i_pv, endpv = vertices->end();
  for (i_pv = vertices->begin();  i_pv != endpv;  ++i_pv) {
    if ( !i_pv->isFake() && i_pv->ndof() > 4){
      nPVGood++;
    }
  }
  
  edm::Handle<pat::PhotonCollection> photonsHandle;
  iEvent.getByToken( photonsToken_, photonsHandle);
  pat::PhotonCollection photons = *photonsHandle;
  pat::PhotonCollection::iterator it = photons.begin();
  std::vector<pat::Photon> photonsVec;
  pat::Photon pho;
  
  uint32_t index = 0;
  uint32_t goodPhoIndex = -1;
  for (; it != photons.end(); ++it, index++) { 
    pho =*it;
    
    N_Photons_Total++ ;
    EventCounterPhoton -> AddBinContent(1, generatorWeight );
    
    if (fabs(it->eta()) <= 1.3) {       
      N_Photons_Eta13++;          
      EventCounterPhoton -> AddBinContent(2, generatorWeight );
      
      N_Photons_PtRange++;      
      EventCounterPhoton -> AddBinContent(3, generatorWeight );
      
      pat::PhotonRef PhotonReftmp(photonsHandle, index);
      
      //      if (isValidPhotonEB_DiPhotonSelection(PhotonReftmp, iEvent, generatorWeight)) {
      if (isValidPhotonEB(PhotonReftmp, iEvent, generatorWeight)) {
	photonsVec.push_back(*it);
	goodPhoIndex=index;
      }
    }
  }
  
  if(photonsVec.size() == 0)    NoGoodPhotons++;
  if(photonsVec.size() == 1)    OnlyOneGoodPhotons++;
  if(photonsVec.size()    > 1)    MoreGoodPhotons++;
  
  // Only one good photon per event
  if (photonsVec.size() != 1)    return false;
  
  Event_NPhotons++;   
  EventCounterPhoton -> AddBinContent(10, generatorWeight );
  EventCounter -> AddBinContent(3, generatorWeight );
  
  pat::Photon photon = photonsVec[0];
  pat::PhotonRef GoodphotonRef(photonsHandle, goodPhoIndex);

  // Process jets
  edm::Handle<pat::JetCollection> jetsHandle;
  
  FOREACH(mJetCollections) {

    if(mVerbose) std::cout << std::endl;
    
    JetInfos infos = mJetCollectionsData[*it];    
    iEvent.getByToken(infos.inputTag, jetsHandle);
    pat::JetCollection jets = *jetsHandle;

    if(mVerbose){    
      if(infos.algo == PUPPI){
	std::cout<< "PUPPI jets"<< std::endl;
      }else{
	std::cout<< "PF AK4 jets"<< std::endl;
      }
    }
    for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
      pat::Jet& jet = *it;
      
      pat::Jet rawJet = jet.correctedJet("Uncorrected");
      rawJet.addUserFloat("qgLikelihood", -10); //quark/gluon discriminator not performed on rawJet
      jet.addUserData("rawJet", rawJet, true); // Store raw jet inside our jet. This allow us to correctly sort the resulting collection
      
    }
    
    if (mDoJEC) {
      correctJets(jets, iEvent, infos.algo, iSetup);
    } else {
      extractRawJets(jets);
    }
    
    edm::Handle<edm::ValueMap<float>> qgTagHandle;
    iEvent.getByToken(qgToken_, qgTagHandle);
    
    processJets(&photon, jets, infos.algo, qgTagHandle, jetsHandle, mJetTrees[*it]);
    
    Event_AfterJets++;     
    EventCounter -> AddBinContent(4, generatorWeight );
    
    // MET
    edm::Handle<pat::METCollection> metsHandle;
    edm::Handle<pat::METCollection> rawMetsHandle;
    if(infos.algo == PUPPI){
      iEvent.getByToken( metPUPPIToken_, metsHandle); 
      iEvent.getByToken( metPUPPIToken_,  rawMetsHandle);    
    }else{
      iEvent.getByToken( metToken_, metsHandle); 
      iEvent.getByToken( metToken_,  rawMetsHandle);    
    }
    pat::METCollection mets = *metsHandle;
    pat::MET& met = mets[0];

    if(mVerbose) std::cout<<"Default Met: "<< met.pt() << " "<< met.eta() << " " << met.phi() << " " << met.et() <<std::endl;     

    pat::METCollection rawMets = *rawMetsHandle;
    pat::MET& rawMet = rawMets[0];
    rawMet.setP4( met.uncorP4() );

    if (mDoJEC || mRedoTypeI) { // authomatic done if mDoJEC is done
      if (mDoFootprint && infos.algo != PUPPI) {
	correctMETWithFootprintAndTypeI(rawMet, met, jets, infos.algo, iEvent, photon);
      }else {
	correctMETWithTypeI(rawMet, met, jets, infos.algo, iEvent, photon);
      }
    }
    
    if (rawMetsHandle.isValid())
      metsToTree(met, rawMet, mMETTrees[*it]);
    else {
      pat::MET emptyRawMet = pat::MET();
      metsToTree(met, emptyRawMet, mMETTrees[*it]);
    }
    
    // Rho
    edm::Handle<double> rhos;
    iEvent.getByToken( rhoToken_, rhos); 
    double rho = *rhos;
    updateBranch(mMiscTrees[*it], &rho, "rho", "D");
    
    mMiscTrees[*it]->Fill();
  }//FOREACH(mJetsCollection)
  
  // Number of vertices for pu reweighting
  edm::Handle<std::vector<PileupSummaryInfo> > puInfos;
  iEvent.getByToken( PUInfoToken_  , puInfos);
  
  float nTrueInteractions = -1;
  int nPUVertex = -1;
  unsigned int nVertex = vertices->size();
  
  edm::EventID eventId = iEvent.id();
  EventNumber_t event = eventId.event();
  RunNumber_t run = eventId.run();
  LuminosityBlockNumber_t lumiBlock = eventId.luminosityBlock();
  
  if (mIsMC) {
    nTrueInteractions = puInfos->at(1).getTrueNumInteractions();
  } else {
    // for Data: recipe to get the PU info
    nTrueInteractions = -1;
  }
  
  // Triggers
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken( triggerResultsToken_, triggerResults);
  
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);
  
  std::vector<std::string>* trigNames = new std::vector<std::string>();
  std::vector<bool>* trigResults = new std::vector<bool>();
  std::vector<double>* trigPrescale = new std::vector<double>();
  
  if (triggerResults.isValid()) {
    static std::vector<boost::regex> validTriggers = { boost::regex("HLT_Photon.*_v.*", boost::regex_constants::icase) };
    const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResults); 
    size_t size = triggerResults->size();
    
    for (size_t i = 0; i < size; i++) {
      std::string triggerName = triggerNames.triggerName(i);
      bool isValid = false;
      
      for (boost::regex& validTrigger: validTriggers) {
	if (boost::regex_match(triggerName, validTrigger)) {
	  isValid = true;
	  break;
	}
      }
      
      if (!isValid)
	continue;
      
      unsigned int index = triggerNames.triggerIndex(triggerName);
      bool passed = triggerResults->accept(index);
      double prescale = triggerPrescales->getPrescaleForIndex(index);
      
      trigPrescale -> push_back(prescale);	
      trigResults->push_back(passed);
      trigNames->push_back(triggerName);
    }
  }
   
  updateBranch(mAnalysisTree, &run, "run", "i");
  updateBranch(mAnalysisTree, &lumiBlock, "lumi_block", "i");
  updateBranch(mAnalysisTree, &event, "event", "i");
  updateBranch(mAnalysisTree, &crossSection, "crossSection", "D");
  updateBranch(mAnalysisTree, &nTrueInteractions, "ntrue_interactions");
  updateBranch(mAnalysisTree, &nVertex, "nvertex", "i");
  updateBranch(mAnalysisTree, &nPVGood, "nvertexGood", "i");
  updateBranch(mAnalysisTree, &nPUVertex, "pu_nvertex", "I");
  updateBranch(mAnalysisTree, &generatorWeight, "generator_weight", "D");
  updateBranch(mAnalysisTree, trigNames, "trigger_names");
  updateBranch(mAnalysisTree, trigResults, "trigger_results");
  updateBranch(mAnalysisTree, trigPrescale, "trigger_prescale");
  
  mAnalysisTree->Fill();
    
  delete trigNames;
  delete trigResults;
  delete trigPrescale;
  
  photonToTree(GoodphotonRef, photon, iEvent);
  
  // Electrons
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken( electronsToken_, electrons);
  electronsToTree(electrons, primaryVertex);
  
  // Muons
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken( muonsToken_, muons);
  muonsToTree(muons, primaryVertex);
  
  mSelectedEvents->SetVal(mSelectedEvents->GetVal() + 1);
  
  Event_Final ++ ;  
  EventCounter -> AddBinContent(5, generatorWeight );
  
  return true;
}// end filter

void GammaJetFilter::correctJets(pat::JetCollection& jets, edm::Event& iEvent, const JetAlgorithm algo, const edm::EventSetup& iSetup) {
  
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it)  {
    pat::Jet& jet = *it;   
    
    double toRaw = jet.jecFactor("Uncorrected"); // factor =1/ correction
    jet.setP4(jet.p4() * toRaw); // It's now a raw jet
    
    if(mVerbose) std::cout<<"Raw Jet   "<< jet.pt() <<std::endl;
    
    double corrections =1.;
    
    // Corrections from txt files for both DATA and MC
    edm::Handle<double> rho_;
    iEvent.getByToken( rhoToken_, rho_);
    
    if(algo == PUPPI){
      PUPPIjetCorrector->setJetEta(jet.eta());
      PUPPIjetCorrector->setJetPt(jet.pt());
      PUPPIjetCorrector->setJetA(jet.jetArea());
      PUPPIjetCorrector->setRho(*rho_);
      corrections = PUPPIjetCorrector->getCorrection();     
    }else{
      jetCorrector->setJetEta(jet.eta());
      jetCorrector->setJetPt(jet.pt());
      jetCorrector->setJetA(jet.jetArea());
      jetCorrector->setRho(*rho_);
      corrections = jetCorrector->getCorrection();     
    }
    
    jet.scaleEnergy(corrections) ; // L1L2L3 + Res (only for data)
    
    if(mVerbose) std::cout<<"Corrected Jet   "<< jet.pt() <<std::endl;
  }
  
  // Sort collection by pt
  std::sort(jets.begin(), jets.end(), mSorter);
}


void GammaJetFilter::correctMETWithTypeI(pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, const JetAlgorithm algo, edm::Event& event, pat::Photon& photon) {
  
  if(mVerbose){
    if(algo == PUPPI){
      std::cout<< "PUPPI"<< std::endl; 
    }else{
      std::cout<< "PF"<< std::endl; 
    }
  }
  
  if(mVerbose) std::cout<< " Correct MET only Type I "<< std::endl; 
  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4    
  
  if(mVerbose) std::cout<<"rawMet: "<< rawMet.pt() << " "<< rawMet.eta() << " " << rawMet.phi() << " " << rawMet.et() <<std::endl; 
  
  if(algo==PUPPI){ 
    float rawMetPx = rawMet.px();
    float rawMetPy = rawMet.py();
    // Re-adding  photon
    rawMetPx += -1.* photon.px();
    rawMetPy += -1.* photon.py();   
    double rawMetPt = sqrt( rawMetPx * rawMetPx + rawMetPy * rawMetPy );
    rawMet.setP4(reco::Candidate::LorentzVector(rawMetPx, rawMetPy, 0., rawMetPt));
  }
  
  if(mVerbose) std::cout<<"NEW rawMet: "<< rawMet.pt() << " "<< rawMet.eta() << " " << rawMet.phi() << " " << rawMet.et() <<std::endl; 
  
  double deltaPx = 0., deltaPy = 0.;
  
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) {
    
    const pat::Jet* rawJet = it->userData<pat::Jet>("rawJet");
    
    double corrs = 1.;
    double corrsForTypeI = 1.;
    edm::Handle<double> rho_;
    event.getByToken( rhoToken_, rho_);
    
    if( algo == PUPPI){      
      PUPPIjetCorrectorForTypeI->setJetEta(rawJet->eta());
      PUPPIjetCorrectorForTypeI->setJetPt(rawJet->pt());
      PUPPIjetCorrectorForTypeI->setJetA(rawJet->jetArea());
      PUPPIjetCorrectorForTypeI->setRho(*rho_);
      corrsForTypeI = PUPPIjetCorrectorForTypeI->getCorrection(); //only RC
      
      PUPPIjetCorrector ->setJetEta(rawJet->eta());
      PUPPIjetCorrector ->setJetPt(rawJet->pt());
      PUPPIjetCorrector ->setJetA(rawJet->jetArea());
      PUPPIjetCorrector ->setRho(*rho_);
      corrs = PUPPIjetCorrector->getCorrection(); // L1L2L3
    }else{
      jetCorrectorForTypeI->setJetEta(rawJet->eta());
      jetCorrectorForTypeI->setJetPt(rawJet->pt());
      jetCorrectorForTypeI->setJetA(rawJet->jetArea());
      jetCorrectorForTypeI->setRho(*rho_);
      corrsForTypeI = jetCorrectorForTypeI->getCorrection(); //only RC

      jetCorrector ->setJetEta(rawJet->eta());
      jetCorrector ->setJetPt(rawJet->pt());
      jetCorrector ->setJetA(rawJet->jetArea());
      jetCorrector ->setRho(*rho_);
      corrs = jetCorrector->getCorrection(); // L1L2L3
    }

    pat::Jet jetRC = *rawJet;
    jetRC.scaleEnergy(corrsForTypeI); // only RC
      
    pat::Jet jet = *rawJet;
    jet.scaleEnergy(corrs); // L1L2L3
      
    double dR = reco::deltaR(photon, jet);
      
    if(jet.pt() > 15 && dR > 0.25) {
	
      double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
      if (emEnergyFraction > 0.90)
	continue;
	
      deltaPx += (jet.px() - jetRC.px());
      deltaPy += (jet.py() - jetRC.py());
    } // jet.pt() && dR
  }//loop over jets
    
  double correctedMetPx = rawMet.px() - deltaPx;
  double correctedMetPy = rawMet.py() - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
    
  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));

  if(mVerbose) std::cout<<"MET: "<< met.pt() << " "<< met.eta() << " " << met.phi() << " " << met.et() <<std::endl; 

}

void GammaJetFilter::correctMETWithFootprintAndTypeI(pat::MET& rawMet, pat::MET& met, const pat::JetCollection& jets, const JetAlgorithm algo, edm::Event& event, pat::Photon& photon) {
  
  if(mVerbose) std::cout<< " Correct MET With FootPrint "<< std::endl; 
  
  edm::Handle<pat::PackedCandidateCollection> pfs;
  event.getByToken(pfToken_, pfs);
  
  float FootprintMEx = 0;
  float FootprintMEy = 0;
  
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0, n = photon.numberOfSourceCandidatePtrs(); i < n; ++i) {
    footprint.push_back(photon.sourceCandidatePtr(i) );
  }
  // now loop on pf candidates
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    // pfcandidate-based footprint removal
    if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
      continue;
    }
    // federico
    //    if(pf.fromPV() == 0) std::cout<<"Dropped PF candidates"<<std::endl;
    if(pf.fromPV() == 0) continue;

    FootprintMEx += -1.* pf.px();
    FootprintMEy += -1.* pf.py();    
  }// loop over pfCand
  
  // Re-adding  photon but reco 
  FootprintMEx += -1.* photon.px();
  FootprintMEy += -1.* photon.py();
  
  double FootprintMEPt = sqrt(FootprintMEx * FootprintMEx + FootprintMEy * FootprintMEy);   
  
  rawMet.setP4(reco::Candidate::LorentzVector(FootprintMEx, FootprintMEy, 0., FootprintMEPt));

  if(mVerbose) std::cout<<"FromPF --> rawMet: "<< rawMet.pt() << " "<< rawMet.eta() << " " << rawMet.phi() << " " << rawMet.et() <<std::endl; 
  
  /////////////// Propagate the JEC to MET 
  // See https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=174324 slide 4
  
  double deltaPx = 0., deltaPy = 0.;
  
  for (pat::JetCollection::const_iterator it = jets.begin(); it != jets.end(); ++it) { 
    
    const pat::Jet* rawJet = it->userData<pat::Jet>("rawJet");
    
    double corrs = 1.;
    double corrsForTypeI = 1.;
    edm::Handle<double> rho_;
    event.getByToken( rhoToken_, rho_);

    if( algo == PUPPI){      
      PUPPIjetCorrectorForTypeI->setJetEta(rawJet->eta());
      PUPPIjetCorrectorForTypeI->setJetPt(rawJet->pt());
      PUPPIjetCorrectorForTypeI->setJetA(rawJet->jetArea());
      PUPPIjetCorrectorForTypeI->setRho(*rho_);
      corrsForTypeI = PUPPIjetCorrectorForTypeI->getCorrection(); //only RC

      PUPPIjetCorrector ->setJetEta(rawJet->eta());
      PUPPIjetCorrector ->setJetPt(rawJet->pt());
      PUPPIjetCorrector ->setJetA(rawJet->jetArea());
      PUPPIjetCorrector ->setRho(*rho_);
      corrs = PUPPIjetCorrector->getCorrection(); // L1L2L3
    }else{
      jetCorrectorForTypeI->setJetEta(rawJet->eta());
      jetCorrectorForTypeI->setJetPt(rawJet->pt());
      jetCorrectorForTypeI->setJetA(rawJet->jetArea());
      jetCorrectorForTypeI->setRho(*rho_);
      corrsForTypeI = jetCorrectorForTypeI->getCorrection(); //only RC

      jetCorrector ->setJetEta(rawJet->eta());
      jetCorrector ->setJetPt(rawJet->pt());
      jetCorrector ->setJetA(rawJet->jetArea());
      jetCorrector ->setRho(*rho_);
      corrs = jetCorrector->getCorrection(); // L1L2L3
    }    
    
    pat::Jet jetRC = *rawJet;
    jetRC.scaleEnergy(corrsForTypeI); //only RC

    pat::Jet jet = *rawJet;
    jet.scaleEnergy(corrs); //L1L2L3
    
    double dR = reco::deltaR(photon, jet);
    
    if (jet.pt()>15 && dR>0.25) {
      
      double emEnergyFraction = rawJet->chargedEmEnergyFraction() + rawJet->neutralEmEnergyFraction();
      if (emEnergyFraction > 0.90)
	continue;
      
      deltaPx += (jet.px() - jetRC.px());
      deltaPy += (jet.py() - jetRC.py());
    }
  }
  
  double correctedMetPx = FootprintMEx  - deltaPx;
  double correctedMetPy = FootprintMEy  - deltaPy;
  double correctedMetPt = sqrt(correctedMetPx * correctedMetPx + correctedMetPy * correctedMetPy);
  
  met.setP4(reco::Candidate::LorentzVector(correctedMetPx, correctedMetPy, 0., correctedMetPt));
  
  if(mVerbose) std::cout<<"MET: "<< met.pt() << " "<< met.eta() << " " << met.phi() << " " << met.et() <<std::endl; 
  
} 

void GammaJetFilter::extractRawJets(pat::JetCollection& jets) {
  
  for (pat::JetCollection::iterator it = jets.begin(); it != jets.end(); ++it) {
    pat::Jet& jet = *it;
    
    const pat::Jet rawJet = jet.correctedJet("Uncorrected");
    jet.addUserData("rawJet", rawJet, true);
  } 
}

void GammaJetFilter::processJets(pat::Photon* photon, pat::JetCollection& jets, const JetAlgorithm algo, edm::Handle<edm::ValueMap<float>>& qgTag, const edm::Handle<pat::JetCollection>& handleForRef, std::vector<TTree*>& trees) {
  
  bool verbose = false;
  
  if(verbose)  std::cout<< "NUOVO EVENTO" << std::endl;

  pat::JetCollection selectedJets; 
  pat::JetCollection::iterator it = jets.begin();
  
  if(verbose)  std::cout<< "numero di jets: "<< jets.size() << std::endl;
  
  uint32_t index = 0;
  uint32_t goodJetIndex = -1;
  for (; it != jets.end(); ++it, index++) {
    
    if(verbose)    std::cout<<"jet numero " <<index<< std::endl; 
    if (! isValidJet(*it) && verbose) std::cout<<"Questo non e' valido "<< std::endl; 
    
    if (! isValidJet(*it)) 
      continue;
    
    goodJetIndex++;
    
    if(verbose)    std::cout<<"goodJetIndex "<< goodJetIndex << std::endl; 

    if( algo == AK4){    
      pat::JetRef jetRef(handleForRef, index);
      it->addUserFloat("qgLikelihood", (*qgTag)[jetRef]);
    }else{ // qg discriminator not performed on AK8 jets
      it->addUserFloat("qgLikelihood", -10);
    }
    
    const double deltaR_threshold = (algo == AK4) ? 0.4 : 0.8;
    
    if(verbose) std::cout<<"Algo = "<< deltaR_threshold << std::endl; 
    
    if(verbose)    std::cout<<"selectedJets.size = "<< selectedJets.size() << std::endl; 
    if(verbose)    std::cout<<"index = "<< index << std::endl; 
    
    if (selectedJets.size() == 0) {
      // First jet selection
      
      if (index > 1) {
      	// It's the third jet of the event. We only want to consider the first two jets for our leading jet,
	// so, throw this event
	break;
      }
      
      const double deltaPhi = reco::deltaPhi(*photon, *it);
      if(verbose)    std::cout<<"deltaPhi phot Jet = "<< deltaPhi << std::endl; 
      if (fabs(deltaPhi) < M_PI / 2.)
	continue; // Only back 2 back event are interesting
      
      const double deltaR = reco::deltaR(*photon, *it);
      if(verbose)    std::cout<<"deltaR phot Jet = "<< deltaR << std::endl; 
      if (deltaR < deltaR_threshold) // This jet is inside the photon. This is probably the photon mis-reconstructed as a jet
      	continue;
      
      // Jet are ordered by pt value.
      // Events are supposed to be balanced between Jet and Gamma
      // If the leading jet has less than 30% of the Photon pt,
      // dump the event as it's not interesting
      if (mFirstJetPtCut && (it->pt() < photon->pt() * mFirstJetThreshold))
      	break;
      
      mSelectedFirstJetIndex->Fill(goodJetIndex);
      selectedJets.push_back(*it);
      
    } else {
      // Second jet selection
      
      const double deltaR = reco::deltaR(*photon, *it);
      if(verbose)    std::cout<<"deltaR phot secondJet = "<< deltaR << std::endl;       
      if (deltaR > deltaR_threshold) {
	mSelectedSecondJetIndex->Fill(goodJetIndex);
	selectedJets.push_back(*it);
      } else {
      	continue;
      }     
      break;
    }
    
  }// end loop over jets
  
  const pat::Jet* firstJet = NULL;
  const pat::Jet* secondJet = NULL;
  
  if (selectedJets.size() > 0) {    
    firstJet = &selectedJets[0];
    mSelectedFirstJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *firstJet)));
    mSelectedFirstJetPhotonDeltaR->Fill(reco::deltaR(*photon, *firstJet));
    if(verbose)    std::cout << "process Jet#1 deltaR= "<< reco::deltaR(*photon,*firstJet) <<std::endl;
    
    if (selectedJets.size() > 1) {
      secondJet = &selectedJets[1];
      mSelectedSecondJetPhotonDeltaPhi->Fill(fabs(reco::deltaPhi(*photon, *secondJet)));
      mSelectedSecondJetPhotonDeltaR->Fill(reco::deltaR(*photon, *secondJet));
      if(verbose)      std::cout << "process Jet#2 deltaR= "<< reco::deltaR(*photon,*secondJet) <<std::endl;
    }
  }
  
  jetsToTree(firstJet, secondJet, trees);
  
  return;
}
////////////////////////////////////////

// ------------ method called once each job just before starting event loop  ------------
void GammaJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void GammaJetFilter::endJob() {
  
  EventCounter -> GetXaxis()->SetBinLabel(1,"All events");
  EventCounter -> GetXaxis()->SetBinLabel(2,"Vertex Cut");
  EventCounter -> GetXaxis()->SetBinLabel(3,"N photons Cut");
  EventCounter -> GetXaxis()->SetBinLabel(4,"After process Jets");
  EventCounter -> GetXaxis()->SetBinLabel(5,"Final events");
  
  EventCounter_Raw -> GetXaxis()->SetBinLabel(1,"All events");
  EventCounter_Raw -> GetXaxis()->SetBinLabel(2,"Vertex Cut");
  EventCounter_Raw -> GetXaxis()->SetBinLabel(3,"N photons Cut");
  EventCounter_Raw -> GetXaxis()->SetBinLabel(4,"After process Jets");
  EventCounter_Raw -> GetXaxis()->SetBinLabel(5,"Final events");
  EventCounter_Raw -> SetBinContent(1,Event_Initial);
  EventCounter_Raw -> SetBinContent(2,Event_VtxCut);
  EventCounter_Raw -> SetBinContent(3,Event_NPhotons);
  EventCounter_Raw -> SetBinContent(4,Event_AfterJets);
  EventCounter_Raw -> SetBinContent(5,Event_Final);
  
  EventCounterPhoton -> GetXaxis()->SetBinLabel(1,"All photons");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(2,"Photons Eta Cut");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(3,"Photons Pt Range");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(4,"! Photon Gen");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(5,"HE  cut");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(6,"Sigma cut");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(7,"Isolation Cuts");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(8,"Electron Veto");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(9,"R9  cut");
  EventCounterPhoton -> GetXaxis()->SetBinLabel(10,"Only One Good Photon");
  
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(1,"All photons");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(2,"Photons Eta Cut");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(3,"Photons Pt Range");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(4,"! Photon Gen");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(5,"HE  cut");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(6,"Sigma cut");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(7,"Isolation Cuts");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(8,"Electron Veto");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(9,"R9  cut");
  EventCounterPhoton_Raw -> GetXaxis()->SetBinLabel(10,"Only One Good Photon");  
  EventCounterPhoton_Raw -> SetBinContent(1,N_Photons_Total);
  EventCounterPhoton_Raw -> SetBinContent(2,N_Photons_Eta13);
  EventCounterPhoton_Raw -> SetBinContent(3,N_Photons_PtRange);
  EventCounterPhoton_Raw -> SetBinContent(4,N_pho_genPhoton);
  EventCounterPhoton_Raw -> SetBinContent(5,N_pho_HE);
  EventCounterPhoton_Raw -> SetBinContent(6,N_pho_Sigma);
  EventCounterPhoton_Raw -> SetBinContent(7,N_pho_Isolation);
  EventCounterPhoton_Raw -> SetBinContent(8,N_pho_ElecVeto);
  EventCounterPhoton_Raw -> SetBinContent(9,N_pho_R9);
  EventCounterPhoton_Raw -> SetBinContent(10,OnlyOneGoodPhotons);
    
  std::cout << " " << std::endl;
  std::cout<< "Event_Initial =  " << Event_Initial << std::endl ;
  std::cout<< "Event_VtxCut = " <<  Event_VtxCut << std::endl ;
  std::cout<< "Event_NPhotons = " << Event_NPhotons << std::endl ;
  std::cout<< "Event_AfterJets = " << Event_AfterJets << std::endl ;
  std::cout<< "Event_Final = " << Event_Final << std::endl ;
  
  std::cout<< "N_Photons_Total = " << N_Photons_Total << std::endl ;
  std::cout<< "N_Photons_Eta13 = " << N_Photons_Eta13 << std::endl ;
  std::cout<< "N_Photons_PtRange = " << N_Photons_PtRange << std::endl ;
  std::cout<< "N_pho_genPhoton =  " << N_pho_genPhoton << std::endl ;
  std::cout<< "N_pho_HE =  " << N_pho_HE << std::endl ;
  std::cout<< "N_pho_Sigma =  " << N_pho_Sigma << std::endl ;
  std::cout<< "N_pho_Isolation =  " << N_pho_Isolation << std::endl ;  
  std::cout<< "N_pho_ElecVeto =  " << N_pho_ElecVeto << std::endl ;
  std::cout<< "N_pho_R9 =  " << N_pho_R9 << std::endl ;
  std::cout<< "OnlyOneGoodPhotons = " << OnlyOneGoodPhotons << std::endl ;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GammaJetFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

bool GammaJetFilter::isValidJet(const pat::Jet& jet) {
  // First, check if this pat::Jet has a gen jet
  if (mIsMC && !jet.genJet()) {
    return false;
  }
  
  //  if (jet.isPFJet()) {    
  // Jet ID
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data    
  // Jet ID works on uncorrected jets. *EnergyFraction take that into account when calculating the fraction,
  // so there's *NO* need to use an uncorrected jet
  bool isValid = true;

  double chf    = jet.chargedHadronEnergyFraction();
  double nhf   = jet.neutralHadronEnergyFraction();
  double cemf = jet.chargedEmEnergyFraction();
  double nemf = jet.neutralEmEnergyFraction();
  //    double muf = jet.muonEnergyFraction();  // not used in Loose Jet ID    
  int NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity(); 
  int neMult = jet.neutralMultiplicity(); 
  int chMult = jet.chargedMultiplicity(); 
    
  if (fabs(jet.eta()) < 3.0) {
    isValid &= nhf < 0.99;
    isValid &= nemf < 0.99;
    isValid &= NumConst > 1;
    if (fabs(jet.eta()) < 2.4) {
      isValid &= chf > 0.;
      isValid &= chMult > 0;
      isValid &= cemf < 0.99;
    }
  }else{
    isValid &= nemf < 0.90;
    isValid &= neMult > 10; 
  }  
  return isValid;
  //  }else{
  // throw cms::Exception("UnsupportedJetType")
  // << "Only PF are supported at this time" << std::endl;
  //}
  // return false;
}

enum class IsolationType {
  CHARGED_HADRONS,
    NEUTRAL_HADRONS,
    PHOTONS
    };

//updated effective areas for SPRING25
float getEffectiveArea(float eta, IsolationType type) {
  eta = fabs(eta);
  switch (type) {
  case IsolationType::CHARGED_HADRONS:
    if (eta < 1.0)
      return 0.0158;
    else if (eta < 1.479)
      return 0.0143;
    else if (eta < 2.0)
      return 0.0115;
    else if (eta < 2.2)
      return 0.0094;
    else if (eta < 2.3)
      return 0.0095;
    else if (eta < 2.4)
      return 0.0068;
    else
      return 0.0053;
    break;
    
  case IsolationType::NEUTRAL_HADRONS:
    if (eta < 1.0)
      return 0.0599;
    else if (eta < 1.479)
      return 0.0819;
    else if (eta < 2.0)
      return 0.0696;
    else if (eta < 2.2)
      return 0.0360;
    else if (eta < 2.3)
      return 0.0360;
    else if (eta < 2.4)
      return 0.0462;
    else
      return 0.0656;
    break;

    //Official  
  case IsolationType::PHOTONS:
    if (eta < 1.0)
      return 0.1271;
    else if (eta < 1.479)
      return 0.1101;
    else if (eta < 2.0)
      return 0.0756;
    else if (eta < 2.2)
      return 0.1175;
    else if (eta < 2.3)
      return 0.1498;
    else if (eta < 2.4)
      return 0.1857;
    else
      return 0.2183;
    break;
    /*
  case IsolationType::PHOTONS:
    if (eta <= 0.9)
      return 0.17;
    else if (eta <= 1.5)
      return 0.14;
    else if (eta <= 2.0)
      return 0.11;
    else if (eta <= 2.2)
      return 0.14;
    else if (eta <= 2.5)
      return 0.22;
    else
      return 0.22;
    break;
*/
  }
  
  return -1;
}

double getCorrectedPFIsolation(double isolation, double rho, float eta, IsolationType type) {
  float effectiveArea = getEffectiveArea(eta, type); 
  //  return std::max(isolation - rho*effectiveArea, 0.);
  return isolation - rho*effectiveArea;
}

// See https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedPhotonIdentificationRun2 -- tight WP
bool GammaJetFilter::isValidPhotonEB(const pat::PhotonRef& photonRef, edm::Event& event, double generatorWeight) {

  bool isValid = true;
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);  
  
  //  real photon matching a gen level photon 
  if (mIsMC && !photonRef->genPhoton())   
    return false;
  EventCounterPhoton -> AddBinContent(4, generatorWeight );
  N_pho_genPhoton++;
  
  // #1: H/E
  isValid &= photonRef->hadTowOverEm() < 0.05;
  if (! isValid)
    return false;
  EventCounterPhoton -> AddBinContent(5, generatorWeight );
  N_pho_HE++;

  //#2: sigma ietaieta
  isValid &= (*full5x5SigmaIEtaIEtaMap)[photonRef] < 0.0100; //Official    
  if (! isValid)
    return false;
  EventCounterPhoton -> AddBinContent(6, generatorWeight );
  N_pho_Sigma++;

  edm::Handle<double> rhos;
  event.getByToken( rhoToken_, rhos);
  double rho = *rhos;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
  // #5
  isValid &= (*phoChargedIsolationMap)[photonRef] < 0.76;
  isValid &= getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::NEUTRAL_HADRONS) < (0.97 + 0.014*photonRef->pt()+0.000019*(photonRef->pt()*photonRef->pt() ) );
  isValid &= getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) < (0.08 + 0.0053*photonRef->pt());
 
  EventCounterPhoton -> AddBinContent(7, generatorWeight );
  N_pho_Isolation++;
  
  isValid &= photonRef->passElectronVeto();
  if (! isValid)
    return false;
  
  EventCounterPhoton -> AddBinContent(8, generatorWeight );
  N_pho_ElecVeto++;

  // added to emule trigger
  isValid &= photonRef->r9() >0.90;
  if (! isValid)
    return false;
  
  EventCounterPhoton -> AddBinContent(9, generatorWeight );
  N_pho_R9++;
  
  return isValid;
  
}
/////////////
bool GammaJetFilter::isValidPhotonEB_DiPhotonSelection(const pat::PhotonRef& photonRef, edm::Event& event, double generatorWeight) {

  bool isValid = true;
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);  
  
  //  real photon matching a gen level photon 
  if (mIsMC && !photonRef->genPhoton())   
    return false;
  EventCounterPhoton -> AddBinContent(4, generatorWeight );
  N_pho_genPhoton++;
  
  // #1: H/E   OK
  isValid &= photonRef->hadTowOverEm() < 0.05;
  if (! isValid)
    return false;
  EventCounterPhoton -> AddBinContent(5, generatorWeight );
  N_pho_HE++;
  
  //#2: sigma ietaieta    OK
  isValid &= (*full5x5SigmaIEtaIEtaMap)[photonRef] < 1.05e-02;    
  if (! isValid)
    return false;
  //#3: New
  isValid &= (*full5x5SigmaIEtaIEtaMap)[photonRef] > 0.001;
  if (! isValid)
    return false;
  EventCounterPhoton -> AddBinContent(6, generatorWeight );
  N_pho_Sigma++;
  
  edm::Handle<double> rhos;
  event.getByToken( rhoToken_, rhos);
  double rho = *rhos;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);
  // #5
  isValid &= (*phoChargedIsolationMap)[photonRef] < 5.;
  if (! isValid)
    return false;
  // #6
  isValid &= 2.5 + getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photonRef->eta(), IsolationType::PHOTONS) - 0.0045*photonRef->pt() < 2.75;
  if (! isValid)
    return false;  
  
  EventCounterPhoton -> AddBinContent(7, generatorWeight );
  N_pho_Isolation++;
  
  EventCounterPhoton -> AddBinContent(8, generatorWeight );
  N_pho_ElecVeto++;
  
  EventCounterPhoton -> AddBinContent(9, generatorWeight );
  N_pho_R9++;
  
  return isValid;
  
}
//---------------------

void GammaJetFilter::particleToTree(const reco::Candidate* particle, TTree* t, std::vector<boost::shared_ptr<void> >& addresses) {
  addresses.clear();
  
  addresses.push_back(boost::shared_ptr<void>(new int((particle) ? 1 : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->et() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->pt() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->eta() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->phi() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->px() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->py() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->pz() : 0)));
  addresses.push_back(boost::shared_ptr<void>(new float((particle) ? particle->energy() : 0)));

  updateBranch(t, addresses[0].get(), "is_present", "I");
  updateBranch(t, addresses[2].get(), "pt");
  updateBranch(t, addresses[3].get(), "eta");
  updateBranch(t, addresses[4].get(), "phi");
  updateBranch(t, addresses[1].get(), "et");
  updateBranch(t, addresses[5].get(), "px");
  updateBranch(t, addresses[6].get(), "py");
  updateBranch(t, addresses[7].get(), "pz");
  updateBranch(t, addresses[8].get(), "e");
}


void GammaJetFilter::photonToTree(const pat::PhotonRef& photonRef, pat::Photon& photon, const edm::Event& event) {
  std::vector<boost::shared_ptr<void> > addresses;

  // std::cout<<"Photon Pt = "<< photon.pt() << std::endl;
  // std::cout<<"Photon Eta = "<< photon.eta() << std::endl;
  // std::cout<<"Photon Phi = "<< photon.phi() << std::endl;
  // std::cout<<"Photon En = "<< photon.energy() << std::endl;

  int pho_is_present =1;
  float pho_et = photon.et();
  float pho_pt = photon.pt();
  float pho_eta = photon.eta();
  float pho_phi = photon.phi();
  float pho_px = photon.px();
  float pho_py = photon.py();
  float pho_pz = photon.pz();
  float pho_e = photon.energy();
  bool hasPixelSeed = photon.hasPixelSeed();
  float hadTowOverEm = photon.hadTowOverEm();
  edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
  event.getByToken(full5x5SigmaIEtaIEtaMapToken_, full5x5SigmaIEtaIEtaMap);  
  float sigmaIetaIeta = (*full5x5SigmaIEtaIEtaMap)[photonRef] ;  
  float r9 = photon.r9();  
  edm::Handle<double> rhos;
  event.getByToken( rhoToken_, rhos);
  float rho = *rhos;
  bool hasMatchedPromptElectron ; 
  if (photon.passElectronVeto())
    hasMatchedPromptElectron = false;
  else 
    hasMatchedPromptElectron = true;

  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  event.getByToken(phoChargedIsolationToken_, phoChargedIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  event.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  event.getByToken(phoPhotonIsolationToken_, phoPhotonIsolationMap);  
  //  float chargedHadronsIsolation = getCorrectedPFIsolation((*phoChargedIsolationMap)[photonRef], rho, photon.eta(), IsolationType::CHARGED_HADRONS);
  float chargedHadronsIsolation = (*phoChargedIsolationMap)[photonRef];
  float neutralHadronsIsolation  = getCorrectedPFIsolation((*phoNeutralHadronIsolationMap)[photonRef], rho, photon.eta(), IsolationType::NEUTRAL_HADRONS);
  float photonIsolation               = getCorrectedPFIsolation((*phoPhotonIsolationMap)[photonRef], rho, photon.eta(), IsolationType::PHOTONS);     
  // variable for the trigger
  float trkSumPtHollowConeDR03  = photon.trkSumPtHollowConeDR03();
  float ecalPFClusterIso = photon.ecalPFClusterIso();
  float hcalPFClusterIso = photon.hcalPFClusterIso();
  // std::cout<<"trkSumPtHollowConeDR03 = "<< trkSumPtHollowConeDR03  << std::endl;
  // std::cout<<"ecalPFClusterIso = "<< ecalPFClusterIso  << std::endl;
  // std::cout<<"hcalPFClusterIso = "<< hcalPFClusterIso  << std::endl;

  // superCluster Info
  float phoSC_pt = photon.superCluster()->rawEnergy()/ cosh(photon.superCluster()->eta());
  float phoSC_eta = photon.superCluster()->eta();
  float phoSC_phi = photon.superCluster()->phi();
  float phoSC_e = photon.superCluster()->rawEnergy();

  // std::cout<<"SC Photon Pt = "<< photon.superCluster()->rawEnergy()/ cosh(photon.superCluster()->eta()) << std::endl;
  // std::cout<<"SC Photon Eta = "<< photon.superCluster()->eta() << std::endl;
  // std::cout<<"SC Photon Phi = "<< photon.superCluster()->phi() << std::endl;
  // std::cout<<"SC Photon Raw En = "<< photon.superCluster()->rawEnergy() << std::endl;
  // std::cout<<"SC Photon En = "<< photon.superCluster()->energy() << std::endl;
  
  updateBranch(mPhotonTree,&pho_is_present,"is_present","I");
  updateBranch(mPhotonTree,&pho_pt,"pt");
  updateBranch(mPhotonTree,&pho_eta,"eta");
  updateBranch(mPhotonTree,&pho_phi,"phi"); 
  updateBranch(mPhotonTree,&pho_et,"et");
  updateBranch(mPhotonTree,&pho_px,"px");
  updateBranch(mPhotonTree,&pho_py,"py");
  updateBranch(mPhotonTree,&pho_pz,"pz");
  updateBranch(mPhotonTree,&pho_e,"e");
  updateBranch(mPhotonTree, &hasPixelSeed, "has_pixel_seed", "O");
  updateBranch(mPhotonTree, &hadTowOverEm, "hadTowOverEm");  
  updateBranch(mPhotonTree, &sigmaIetaIeta, "sigmaIetaIeta");
  updateBranch(mPhotonTree, &r9, "r9");
  updateBranch(mPhotonTree, &rho, "rho");  
  updateBranch(mPhotonTree, &hasMatchedPromptElectron, "hasMatchedPromptElectron", "O");  
  updateBranch(mPhotonTree, &chargedHadronsIsolation, "chargedHadronsIsolation");
  updateBranch(mPhotonTree, &neutralHadronsIsolation, "neutralHadronsIsolation");
  updateBranch(mPhotonTree, &photonIsolation, "photonIsolation");
  // info for the trigger
  updateBranch(mPhotonTree, &trkSumPtHollowConeDR03, "trkSumPtHollowConeDR03");
  updateBranch(mPhotonTree, &ecalPFClusterIso, "ecalPFClusterIso");
  updateBranch(mPhotonTree, &hcalPFClusterIso, "hcalPFClusterIso");
  // added sc info
  updateBranch(mPhotonTree, &phoSC_pt, "SC_pt");
  updateBranch(mPhotonTree, &phoSC_eta, "SC_eta");
  updateBranch(mPhotonTree, &phoSC_phi, "SC_phi");
  updateBranch(mPhotonTree, &phoSC_e, "SC_e");

  mPhotonTree->Fill();
  
  if (mIsMC) {
    particleToTree(photon.genPhoton(), mPhotonGenTree, addresses);
    mPhotonGenTree->Fill();
  }
}

void GammaJetFilter::jetsToTree(const pat::Jet* firstJet, const pat::Jet* secondJet, const std::vector<TTree*>& trees) {
  jetToTree(firstJet, mIsMC, trees[0], trees[4]);
  jetToTree(secondJet, false, trees[1], trees[5]);
  
  // Raw jets
  const pat::Jet* rawJet = (firstJet) ? firstJet->userData<pat::Jet>("rawJet") : NULL; 
  jetToTree(rawJet, false, trees[2], NULL);
  rawJet = (secondJet) ? secondJet->userData<pat::Jet>("rawJet") : NULL;
  jetToTree(rawJet, false, trees[3], NULL);
}

void findNeutrinos(const reco::Candidate* parent, std::vector<const reco::Candidate*>& neutrinos) {
  
  int pdg_id = abs(parent->pdgId());
  if (pdg_id == 12 || pdg_id == 14 || pdg_id == 16) {
    if (std::find_if(neutrinos.begin(), neutrinos.end(), [parent] (const reco::Candidate* candidate) -> bool {
          static double EPSILON = 0.0001;
          bool same = (candidate->pdgId() == parent->pdgId());
          same &= (fabs(candidate->px() - candidate->px()) < EPSILON);
          same &= (fabs(candidate->py() - candidate->py()) < EPSILON);
          same &= (fabs(candidate->pz() - candidate->pz()) < EPSILON);
	  
	  return same;
	  
	}) == neutrinos.end()) {
      neutrinos.push_back(parent);
    }
    return;
  }
  
  for (unsigned int i = 0; i < parent->numberOfDaughters(); i++) {
    const reco::Candidate* d = parent->daughter(i);
    findNeutrinos(d, neutrinos);
  }
}

void GammaJetFilter::jetToTree(const pat::Jet* jet, bool _findNeutrinos, TTree* tree, TTree* genTree) {
  std::vector<boost::shared_ptr<void> > addresses;
  particleToTree(jet, tree, addresses);
  
  if (mIsMC) {
    mNeutrinos->Clear("C");
    mNeutrinosPDG->Clear("C");
  }
  
  if (jet) {
    float area = jet->jetArea();
    float qgTagLikelihood = jet->userFloat("qgLikelihood");
    ///     
    float EnergyTot = jet->chargedHadronEnergy() + jet->neutralHadronEnergy() + jet->chargedEmEnergy() + jet->neutralEmEnergy() + jet->muonEnergy();
    float jetCHEnF   = jet->chargedHadronEnergy() / EnergyTot; 
    float jetNHEnF   = jet->neutralHadronEnergy() / EnergyTot; 
    float jetCEmEnF = jet->chargedEmEnergy() / EnergyTot; 
    float jetNEmEnF = jet->neutralEmEnergy() / EnergyTot; 
    float jetMuEnF   = jet->muonEnergy() / EnergyTot; 
    /*
    float SumEnFraction_New   = jetCHEnF_New + jetNHEnF_New + jetCEmEnF_New + jetNEmEnF_New + jetMuEnF_New ;
    float jecFactor = jet->jecFactor(0);
    std::cout<<"Jet Eta: " << jet->eta() << std::endl;
    std::cout<<"Jet Energy: " << jet->energy() << std::endl;
    const pat::Jet* rawJet = jet->userData<pat::Jet>("rawJet");
    std::cout<<"RawJet Energy: " << rawJet->energy() << std::endl;
    std::cout<<"JEC Factor: " << jecFactor << std::endl;
    std::cout<<"Energy Tot : " << EnergyTot << std::endl;
    std::cout<<" CH Energy : " << jet->chargedHadronEnergy() << std::endl;
    std::cout<<" NH Energy : " << jet->neutralHadronEnergy() << std::endl;
    std::cout<<" C Em Energy : " << jet->chargedEmEnergy() << std::endl;
    std::cout<<" N Em Energy : " << jet->neutralEmEnergy() << std::endl;
    std::cout<<" Mu Energy : " << jet->muonEnergy() << std::endl;
    std::cout<<"New CHEnF : " << jetCHEnF_New << std::endl;
    std::cout<<"New NHEnF : " << jetNHEnF_New << std::endl;
    std::cout<<"New CEmEnF : " << jetCEmEnF_New << std::endl;
    std::cout<<"New NEmEnF : " << jetNEmEnF_New << std::endl;
    std::cout<<"New MuEnF : " << jetMuEnF_New << std::endl;
    std::cout<<"New SumEnFraction : " << SumEnFraction_New << std::endl;
    
    //jet energy composition
    float jetCHEnF_Old   = jet->chargedHadronEnergyFraction(); 
    float jetNHEnF_Old   = jet->neutralHadronEnergyFraction(); 
    float jetCEmEnF_Old = jet->chargedEmEnergyFraction();
    float jetNEmEnF_Old = jet->neutralEmEnergyFraction();
    float jetMuEnF_Old   = jet->muonEnergyFraction();
    float SumEnFraction_Old   = jetCHEnF + jetNHEnF + jetCEmEnF + jetNEmEnF + jetMuEnF ; 
    std::cout<<"Old CHEnF : " << jetCHEnF_Old << std::endl;
    std::cout<<"Old NHEnF : " << jetNHEnF_Old << std::endl;
    std::cout<<"Old CEmEnF : " << jetCEmEnF_Old << std::endl;
    std::cout<<"Old NEmEnF : " << jetNEmEnF_Old << std::endl;
    std::cout<<"Old MuEnF : " << jetMuEnF_Old << std::endl;
    std::cout<<"Old SumEnFraction : " << SumEnFraction_Old << std::endl;
    */
    //jet constituents multiplicities
    int jetCHMult = jet->chargedHadronMultiplicity();
    int jetNHMult = jet->neutralHadronMultiplicity();
    int jetElMult   = jet->electronMultiplicity();
    int jetPhMult  = jet->photonMultiplicity();
    int jetMuonMult = jet->muonMultiplicity();
    int jetNeutralMult  = jet->neutralMultiplicity(); 
    int jetChargedMult = jet->chargedMultiplicity(); 

    updateBranch(tree, &area, "jet_area");
    updateBranch(tree, &qgTagLikelihood, "qg_tag_likelihood");
    updateBranch(tree, &jetCHEnF, "jet_CHEnF");
    updateBranch(tree, &jetNHEnF, "jet_NHEnF");
    updateBranch(tree, &jetCEmEnF, "jet_CEmEnF");
    updateBranch(tree, &jetNEmEnF, "jet_NEmEnF");
    updateBranch(tree, &jetMuEnF, "jet_MuEnF");
    updateBranch(tree, &jetCHMult, "jet_CHMult", "I");
    updateBranch(tree, &jetNHMult, "jet_NHMult", "I");
    updateBranch(tree, &jetPhMult, "jet_PhMult", "I");
    updateBranch(tree, &jetElMult, "jet_ElMult", "I");
    updateBranch(tree, &jetMuonMult, "jet_MuonMult", "I");
    updateBranch(tree, &jetChargedMult, "jet_ChargedMult", "I");
    updateBranch(tree, &jetNeutralMult, "jet_NeutralMult", "I");
    
    tree->Fill(); // This Fill() must be called inside the {} block, otherwise it'll crash. Don't move it!
  } else {
    tree->Fill();
  }
  
  if (genTree) {
    particleToTree((jet) ? jet->genJet() : NULL, genTree, addresses);
    
    if (_findNeutrinos) {
      static bool init = false;
      if (! init) {
        genTree->Branch("neutrinos", &mNeutrinos, 32000, 0);
        genTree->Branch("neutrinos_pdg_id", &mNeutrinosPDG, 32000, 0);
        init = true;
      }
    }
    
    if (jet && _findNeutrinos) {
      const reco::Candidate* parton = (jet) ? jet->genParton() : NULL;
      
      if (parton) {
        if (abs(parton->pdgId()) == 5 || abs(parton->pdgId()) == 4) {
	  
          std::vector<const reco::Candidate*> neutrinos;
          findNeutrinos(parton, neutrinos);
	  
          if (neutrinos.size() > 0) {
            // Build TCloneArray of TLorentzVector
            unsigned int index = 0;
            for (const reco::Candidate* neutrino: neutrinos) {
              TLorentzVector* p4 = (TLorentzVector*) mNeutrinos->ConstructedAt(index);
              p4->SetPxPyPzE(neutrino->px(), neutrino->py(), neutrino->pz(), neutrino->energy());
	      
              TParameter<int>* pdg_id = (TParameter<int>*) mNeutrinosPDG->ConstructedAt(index++);
              pdg_id->SetVal(neutrino->pdgId());
            }
          }
        }
      }
    }
    
    // Add parton id and pt
    const reco::Candidate* parton = (jet) ? jet->genParton() : NULL;
    int pdgId = (parton) ? parton->pdgId() : 0;
    updateBranch(genTree, &pdgId, "parton_pdg_id", "I");
    
    TLorentzVector parton_p4;
    if (parton) {
      parton_p4.SetPxPyPzE(parton->px(), parton->py(), parton->pz(), parton->energy());
    }
    TLorentzVector* p_parton_p4 = &parton_p4;
    
    TBranch* branch = genTree->GetBranch("parton_p4");
    if (branch == NULL) {
      branch = genTree->Branch("parton_p4", &p_parton_p4);
    } else {
      branch->SetAddress(&p_parton_p4);
    }
    
    int flavour = (jet) ? jet->partonFlavour() : 0;
    updateBranch(genTree, &flavour, "parton_flavour", "I");
    
    genTree->Fill();
  }
}

void GammaJetFilter::metsToTree(const pat::MET& met, const pat::MET& rawMet, const std::vector<TTree*>& trees) {
  metToTree(&met, trees[0], trees[2]);
  metToTree(&rawMet, trees[1], NULL);
}

void GammaJetFilter::metToTree(const pat::MET* met, TTree* tree, TTree* genTree) {
  std::vector<boost::shared_ptr<void> > addresses;
  particleToTree(met, tree, addresses);
  
  tree->Fill();
  
  if (genTree) {
    particleToTree((met) ? met->genMET() : NULL, genTree, addresses);
    genTree->Fill();
  }
}


void GammaJetFilter::electronsToTree(const edm::Handle<pat::ElectronCollection>& electrons, const reco::Vertex& pv) {
  
  int n = electrons->size();
  static int   id[30];
  static float isolation[30];
  static float pt[30];
  static float px[30];
  static float py[30];
  static float pz[30];
  static float eta[30];
  static float phi[30];
  static int   charge[30];
  
  int i = 0;
  for (pat::ElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); ++it, i++) {
    const pat::Electron& electron = *it;
    
    if (i >= 30)
      break;
    
    // See https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_ele
    bool elecID = fabs(pv.z() - it->vertex().z()) < 1.;
    elecID     &= it->et() > 30.;
    elecID     &= fabs(it->eta()) < 2.5 && (it->superCluster()->eta() > 1.4442 && it->superCluster()->eta() < 1.5660);
    elecID     &= it->dB() < 0.02;
    elecID     &= ((int) it->electronID("eidLoose") & 0x1);
    
    float iso     = (it->dr03TkSumPt() + it->dr03EcalRecHitSumEt() + it->dr03HcalTowerSumEt()) / it->et();
    
    id[i]         = elecID;
    isolation[i]  = iso;
    pt[i]         = electron.pt();
    px[i]         = electron.px();
    py[i]         = electron.py();
    pz[i]         = electron.pz();
    eta[i]        = electron.eta();
    phi[i]        = electron.phi();
    charge[i]     = electron.charge();
  }

  updateBranch(mElectronsTree, &n, "n", "I");
  updateBranchArray(mElectronsTree, id, "id", "n", "I");
  updateBranchArray(mElectronsTree, isolation, "isolation", "n");
  updateBranchArray(mElectronsTree, pt, "pt", "n");
  updateBranchArray(mElectronsTree, px, "px", "n");
  updateBranchArray(mElectronsTree, py, "py", "n");
  updateBranchArray(mElectronsTree, pz, "pz", "n");
  updateBranchArray(mElectronsTree, eta, "eta", "n");
  updateBranchArray(mElectronsTree, phi, "phi", "n");
  updateBranchArray(mElectronsTree, charge, "charge", "n", "I");

  mElectronsTree->Fill();  

}

void GammaJetFilter::muonsToTree(const edm::Handle<pat::MuonCollection>& muons, const reco::Vertex& pv) {
  
  int n = muons->size();
  static int   id[30];
  static float isolation[30];
  static float delta_beta_isolation[30];
  static int isLooseMuon[30];
  int nLooseMuon = 0;
  static float pt[30];
  static float px[30];
  static float py[30];
  static float pz[30];
  static float eta[30];
  static float phi[30];
  static int   charge[30];
  
  int i = 0;
  for (pat::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); ++it, i++) {
    const pat::Muon& muon = *it;
    
    if (i >= 30)
      break;
    
    // See https://twiki.cern.ch/twiki/bin/view/CMS/TopLeptonPlusJetsRefSel_mu
    bool muonID = it->isGlobalMuon();
    //FIXME: reco::Tracks need to be keept in PF2PAT.
    //It's not the case right now, so muon ID will be incorrect
    if (it->globalTrack().isNull() || it->innerTrack().isNull() || it->muonBestTrack().isNull() || it->track().isNull()) {
      muonID = false;
    } else {
      muonID     &= it->globalTrack()->normalizedChi2() < 10.;
      muonID     &= it->globalTrack()->hitPattern().numberOfValidMuonHits() > 0;
      muonID     &= it->numberOfMatchedStations() > 1;
      muonID     &= it->dB() < 0.2;
      muonID     &= fabs(it->muonBestTrack()->dz(pv.position())) < 0.5;
      muonID     &= it->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
      muonID     &= it->track()->hitPattern().trackerLayersWithMeasurement() > 5;
    }
    
    float relIso = (it->chargedHadronIso() + it->neutralHadronIso() + it->photonIso()) / it->pt();
    float deltaBetaRelIso = (it->chargedHadronIso() + std::max((it->neutralHadronIso() + it->photonIso()) - 0.5 * it->puChargedHadronIso(), 0.0)) / it->pt();
    
    id[i]          = muonID;
    isolation[i]   = relIso;
    delta_beta_isolation[i] = deltaBetaRelIso;
    isLooseMuon[i]          = muon.isLooseMuon();
    pt[i]          = muon.pt();
    px[i]          = muon.px();
    py[i]          = muon.py();
    pz[i]          = muon.pz();
    eta[i]         = muon.eta();
    phi[i]         = muon.phi();
    charge[i]      = muon.charge();
    
    if(pt[i] > 5 && isLooseMuon[i]) nLooseMuon++ ;
  }

  updateBranch(mMuonsTree, &n, "n", "I");
  updateBranchArray(mMuonsTree, id, "id", "n", "I");
  updateBranchArray(mMuonsTree, isolation, "relative_isolation", "n");
  updateBranchArray(mMuonsTree, delta_beta_isolation, "delta_beta_relative_isolation", "n");
  updateBranchArray(mMuonsTree, isLooseMuon, "isLooseMuon", "n", "I");
  updateBranch(mMuonsTree, &nLooseMuon, "nLooseMuon", "I");
  updateBranchArray(mMuonsTree, pt, "pt", "n");
  updateBranchArray(mMuonsTree, px, "px", "n");
  updateBranchArray(mMuonsTree, py, "py", "n");
  updateBranchArray(mMuonsTree, pz, "pz", "n");
  updateBranchArray(mMuonsTree, eta, "eta", "n");
  updateBranchArray(mMuonsTree, phi, "phi", "n");
  updateBranchArray(mMuonsTree, charge, "charge", "n", "I");
  
  mMuonsTree->Fill();  
}

//define this as a plug-in
DEFINE_FWK_MODULE(GammaJetFilter);
