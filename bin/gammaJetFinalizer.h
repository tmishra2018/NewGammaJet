#pragma once

#include <TRandom3.h>

#include "Tree/AnalysisTree.h"
#include "Tree/PhotonTree.h"
#include "Tree/JetTree.h"
#include "Tree/GenJetTree.h"
#include "Tree/METTree.h"
#include "Tree/MiscTree.h"
#include "Tree/ElectronTree.h"
#include "Tree/MuonTree.h"
#include "Tree/rootNTupleClass.h"


#include "etaBinning.h"
#include "ptBinning.h"
#include "HLTptBinning.h"
#include "vertexBinning.h"
#include "runBinning.h"
#include "extrapBinning.h"
#include "triggers.h"
#include "GaussianProfile.h"
#include "fineetaBinning.h"

#include <vector>
#include <memory>
#include <unordered_map>

namespace fwlite {
  class TFileService;
}

namespace pat {
  class Photon;
  class Jet;
}

namespace edm {
  class EventBase;
}

class TTree;
class TChain;
class TFileDirectory;

enum JetAlgo {
  AK4,
  AK8
};

enum JetType {
  PF,
  PUPPI
};

template <typename T>
struct ExtrapolationVectors {
  typedef std::vector<std::vector<std::vector<T*> > > type;
};


class PUReweighter;

class GammaJetFinalizer
{
 public:
  GammaJetFinalizer();
  ~GammaJetFinalizer();
  
    void setInputFiles(const std::vector<std::string>& files) {
      mInputFiles = files;
      checkInputFiles();
    }

    void setDatasetName(const std::string& name) {
      mDatasetName = name;
    }
    
    void setRunera(const std::string& runera) {
      mRunera = runera;
    }

    void setJetAlgo(const std::string& jetType, const std::string& jetAlgo) {
      if (jetType == "pf") {
        mJetType = PF;
      } else {
        mJetType = PUPPI;
      }
      
      if (jetAlgo == "ak4") {
        mJetAlgo = AK4;
      } else {
        mJetAlgo = AK8;
      }
    }

    void setMC(bool isMC) {
      mIsMC = isMC;
    }

    void setUseExternalJEC(bool useExternalJEC) {
      mUseExternalJECCorrecion = useExternalJEC;
    }
    /*
    void setBatchJob(int currentJob, int totalJobs) {
      mIsBatchJob = true;
      mCurrentJob = currentJob;
      mTotalJobs = totalJobs;
    }
    */
    void setAlphaCut(float alphaCut) {
      mAlphaCut = alphaCut;
    }

    void setCHS(bool chs) {
      mUseCHS = chs;
    }

    void setVerbose(bool verbose) {
      mVerbose = verbose;
    }

    void setUncutTrees(bool uncutTrees) {
      mUncutTrees = uncutTrees;
    }

    void runAnalysis();

  private:
    void checkInputFiles();
    void loadFiles(TChain& chain);

    //bool passTrigger(const TRegexp& regexp) const;
    int checkTrigger(std::string& passedTrigger, float& weight);
    int checkTriggerfulltree(std::string& passedTrigger, double& HLT1, double& HLT2, double& HLT3, double& HLT4, double& HLT5, double& HLT6, double& HLT7, double& triggHLT1, double& triggHLT2, double& triggHLT3, double& triggHLT4, double& triggHLT5, double& triggHLT6, float& weight);

    void cleanTriggerName(std::string& trigger);
    // new RD PU reweighting
    // void computePUWeight(const std::string& passedTrigger, int run_period);
    void computePUWeight();

    void computePUWeight_NVtxBased(double ptPhot, int nvertex);

    template<typename T>
      std::vector<std::vector<T*> > buildEtaPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
    template<typename T>
      std::vector<std::vector<T*> > buildfineEtaPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
    template<typename T>
      std::vector<std::vector<T*> > buildEtaHLTPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
      
    template<typename T>
      std::vector<T*> buildPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildHLTPtVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
      template<typename T>
      std::vector<T*> buildHLTPtVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
    template<typename T>
      std::vector<T*> buildEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);

    template<typename T>
      std::vector<std::vector<T*> > buildEtaVertexVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<T*> buildVertexVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<std::vector<T*> > > buildExtrapolationEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
      
   template<typename T>
      std::vector<std::vector<std::vector<T*> > > buildExtrapolationfineEtaVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<T*> > buildExtrapolationVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    // federico
    template<typename T>
      std::vector<T*> buildRunVector(TFileDirectory dir, const std::string& branchName, const std::string& etaName, int nBins, double xMin, double xMax);
    template<typename T>
      std::vector<std::vector<T*> > buildEtaRunVector(TFileDirectory dir, const std::string& branchName, int nBins, double xMin, double xMax);
    

    void cloneTree(TTree* from, TTree*& to);

    std::string buildPostfix();

    // Datas from step 2
    rootNtupleClass fullinfo;
    
    AnalysisTree analysis;
    PhotonTree photon;
    GenTree genPhoton;
    MuonTree muons;
    ElectronTree electrons;

    JetTree firstJet;
    JetTree firstRawJet;
    GenJetTree firstGenJet;

    JetTree secondJet;
    JetTree secondRawJet;
    GenJetTree secondGenJet;

    METTree MET;
    GenTree genMET;
    METTree rawMET;

    MiscTree misc;

    EtaBinning mEtaBinning;
    fineEtaBinning mfineEtaBinning;
    PtBinning mPtBinning;
    HLTPtBinning mHLTPtBinning;
    VertexBinning mVertexBinning;
    RunBinning mRunBinning;
    ExtrapBinning mExtrapBinning;

    std::vector<std::string> mInputFiles;
    std::string mDatasetName;
    std::string mRunera;
    JetType mJetType;
    JetAlgo mJetAlgo;
    bool mNoPUReweighting;
    bool mIsMC;
    bool mIsBatchJob;
    int mTotalJobs;
    int mCurrentJob;
    bool mUseExternalJECCorrecion;

    float  mAlphaCut;
    bool   mDoMCComparison;
    bool   mUseCHS;
    bool   mVerbose;
    bool   mUncutTrees;

//new RD PU rweighting
    std::map<std::pair<std::string,int>, boost::shared_ptr<PUReweighter>> mLumiReweighting;
//old wrong S10 PU reweigting
//    std::unordered_map<std::string, boost::shared_ptr<PUReweighter>> mLumiReweighting;

    float mPUWeight;

    int ptBinGen;
    int etaBinGen;    
    int extrapGenBin;    

    //MPF
    float deltaPhi_Photon_MET;
    float respMPF;
    float deltaPhi_Photon_MET_raw;
    float respMPFRaw;
    float deltaPhi_Photon_MET_gen;
    float respMPFGen;
            
    // Balancing
    float respBalancing;
    float respBalancingRaw;
    float respBalancingGen;

    float respPhotGamma;
    float respGenPhot;
     
     // Triggers on data
     Triggers* mTriggers;
     MCTriggers* mMCTriggers;
     TRandom3 mRandomGenerator;
};
