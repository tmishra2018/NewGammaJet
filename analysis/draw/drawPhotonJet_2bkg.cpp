#include <stdlib.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"

#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>


bool useMCassoc_ = false;
bool ONEVTX = false;
bool OUTPUT_GRAPHS = true;

int main(int argc, char* argv[]) {

  if (argc != 7 && argc != 8) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')] [flags=\"\"]" << std::endl;
    exit(23);
  }

  std::string data_dataset(argv[1]);
  std::string mc_photonjet(argv[2]);
  std::string mc_QCD(argv[3]);
  std::string recoType(argv[4]);
  std::string jetAlgo(argv[5]);
  std::string norm(argv[6]);
  if (norm != "LUMI" && norm != "SHAPE") {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }
  std::string flags = "";
  if (argc == 8) {
    std::string flags_str(argv[7]);
    flags = flags_str;
  }


  //  std::string algoType;
  //  if (recoType == "calo") {
  //    algoType = jetAlgo;
  //  } else { //pf
  //    algoType = recoType + jetAlgo;
  //  }
  //  if (recoType == "jpt" && jetAlgo == "akt4") {
  //    algoType = "jptak4";
  //  }

  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  std::string postFix = recoType + jetAlgo;

  postFix += "chs";
  

  //Float_t etamax = 3.;
  //bool sameEvents = false; //until njets histos have no overflows... or maybe use GetEntries instead of integral?

  drawBase* db = new drawBase("PhotonJet", recoType, jetAlgo, OUTPUT_GRAPHS);
  db->set_pdf_aussi((bool)false);
  db->set_flags(flags);
  db->set_isCMSArticle(false);

  std::cout << "flags set." << std::endl;

  TString dataFileName;
  if (flags.length() > 0) {
    dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  }

  TFile* dataFile = TFile::Open(dataFileName);

  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
    db->add_dataFile(dataFile, data_dataset);
  }

  TString mc1FileName;
  if (flags.length() > 0) {
    mc1FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_photonjet.c_str(), postFix.c_str(), flags.c_str());
  } else {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened GJet file '" << mc1FileName << "'." << std::endl;

  if (mcPhotonJetFile) {
    db->add_mcFile(mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", TColor::GetColor(192, 41, 66));
  }

  if (mc_QCD == " ") std::cout<<"Not opened QCD file "<<std::endl;

  if (mc_QCD != " ") {
    TString mc2FileName;
    if (flags.length() > 0) {
      mc2FileName = TString::Format("PhotonJet_%s_%s_%s.root", mc_QCD.c_str(), postFix.c_str(), flags.c_str());
    } else {
      mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc_QCD.c_str(), postFix.c_str());
    }
    TFile* mcQCDFile = TFile::Open(mc2FileName);
    std::cout << "Opened mc QCD file '" << mc2FileName << "'." << std::endl;

    if (mcQCDFile && mc_QCD != mc_photonjet) {
      db->add_mcFile(mcQCDFile, mc_QCD, "QCD MC", TColor::GetColor("#ECD078"));
    }
  }

  // MC should already be normalized to a lumi of 1 pb-1
  // Read luminosity
  double dLumi = 1e6;

  //Federico --> restored Lumi
  if (dataFile) {
    TParameter<double>* lumi = static_cast<TParameter<double>*>(dataFile->Get("analysis/luminosity"));
    dLumi = lumi->GetVal();
  
  }

  std::cout<< "Lumi "<< dLumi << std::endl;

  //  db->set_lumi(dLumi * 1e-6);
  db->set_lumi(dLumi);
  //db->set_lumi(3000);
  if (norm == "LUMI") {
    db->set_lumiNormalization();
  } else {
    db->set_shapeNormalization();
  }

  db->setFolder("analysis");
  db->set_outputdir();

  bool log = true;
  gErrorIgnoreLevel = kWarning;

  db->setOutputGraphs(true);

  db->set_rebin(1);

  // Data / MC comparison
  db->drawHisto("ptPhoton", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("EtaPhoton", "Photon Eta", " ", "Events" , log);
  db->drawHisto("PhiPhoton", "Photon Phi", " ", "Events" , false);
  db->drawHisto("ptFirstJet", "Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("EtaFirstJet", "FirstJet Eta", " ", "Events" , log);
  db->drawHisto("PhiFirstJet", "FirstJet Phi", " ", "Events" , false);
  db->drawHisto("ptSecondJet", "2nd Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("EtaSecondJet", "2nd Jet Eta", "", "Events" , log);
  db->drawHisto("PhiSecondJet", "2nd Jet Phi", " ", "Events" , false);
  db->drawHisto("MET", "Missing E_{T}", "GeV", "Events", log);
  db->drawHisto("alpha", "#alpha", "", "Events", log);
  db->drawHisto("deltaPhi_NoCut", "#Delta#phi", "", "Events", log);


  //  db->drawHisto("responsePhotGamma", "Response Photon (Reco-Gen)", " ", "Events" , false);
  //  db->drawHisto("responseBalancingGen", "Response Jet (Reco-Gen)", " ", "Events" , false);

  db->drawHisto("ptPhoton_passedID", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("EtaPhoton_passedID", "Photon Eta", " ", "Events" , log);
  db->drawHisto("PhiPhoton_passedID", "Photon Phi", " ", "Events" , false);
  db->drawHisto("ptFirstJet_passedID", "Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("EtaFirstJet_passedID", "FirstJet Eta", " ", "Events" , log);
  db->drawHisto("PhiFirstJet_passedID", "FirstJet Phi", " ", "Events" , false);
  db->drawHisto("ptSecondJet_passedID", "2nd Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("EtaSecondJet_passedID", "2nd Jet Eta", " ", "Events" , log);
  db->drawHisto("PhiSecondJet_passedID", "2nd Jet Phi", " ", "Events" , false);

  db->drawHisto("ptSecondJet_2ndJetOK", "2nd Jet Transverse Momentum", "GeV", "Events", log);
  db->drawHisto("EtaSecondJet_2ndJetOK", "2nd Jet Eta", " ", "Events" , log);
  db->drawHisto("PhiSecondJet_2ndJetOK", "2nd Jet Phi", " ", "Events" , false);

  db->drawHisto("MET_passedID", "Missing E_{T}", "GeV", "Events", log);
  db->drawHisto("alpha_passedID", "#alpha", "", "Events", log);

  //db->drawHisto("etaPhot", "#eta", "", "", log);

  db->set_rebin(1);
  db->drawHisto("nvertex", "Number of Reconstructed Vertices", "", "Events", log);
  db->drawHisto("nvertex_reweighted", "Number of Reconstructed Vertices (after reweighting)", "", "Events", log);

  db->set_rebin(2);
  db->drawHisto("deltaPhi_passedID", "#Delta #varphi", "", "Events", log);

  db->drawHisto("hadTowOverEm", "", "", "Events", log);
  db->drawHisto("sigmaIetaIeta", "", "", "Events", log);
  db->drawHisto("rho", "#rho", "", "Events", log);
  db->drawHisto("chargedHadronsIsolation", "Charged hadrons isolation", "", "Events", log);
  db->drawHisto("neutralHadronsIsolation", "Neutral hadrons isolation", "", "Events", log);
  db->drawHisto("photonIsolation", "Photon isolation", "", "Events", log);

  db->drawHisto("hadTowOverEm_passedID", "", "", "Events", log);
  db->drawHisto("sigmaIetaIeta_passedID", "", "", "Events", log);
  db->drawHisto("rho_passedID", "#rho", "", "Events", log);
  db->drawHisto("chargedHadronsIsolation_passedID", "Charged hadrons isolation", "", "Events", log);
  db->drawHisto("neutralHadronsIsolation_passedID", "Neutral hadrons isolation", "", "Events", log);
  db->drawHisto("photonIsolation_passedID", "Photon isolation", "", "Events", log);


  db->set_rebin(5);

  db->setOutputGraphs(OUTPUT_GRAPHS);

  PtBinning ptBinning;
  std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();

  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();

     //giulia --- not useful
//  db->set_rebin(2);
////Jet energy composition
//  db->setFolder("analysis/ecomposition");
//  for (size_t i = 0; i < etaBinningSize; i++) {
//    db->set_legendTitle(etaBinning.getBinTitle(i));
//
//    TString histoName = TString::Format("ChHadronEnergy_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Charged Hadron Energy", "", "Events", false);
///*    histoName = TString::Format("NHadronEnergy_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Neutral Hadron Energy", "", "Events", log);
//    histoName = TString::Format("ElEnergy_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Electron Energy", "", "Events", log);
//    histoName = TString::Format("PhEnergy_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Photon Energy", "", "Events", log);
//    histoName = TString::Format("MuEnergy_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Muon Energy", "", "Events", log);
////multiplicities
//    histoName = TString::Format("ChHadronMult_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Charged Hadron Multiplicity", "", "Events", log);
//    histoName = TString::Format("NHadronMult_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Neutral Hadron Multiplicity", "", "Events", log);
//    histoName = TString::Format("ElMult_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Electron Multiplicity", "", "Events", log);
//    histoName = TString::Format("PhMult_%s", etaBinning.getBinName(i).c_str());
//    db->drawHisto_vs_pt(ptBins, histoName.Data(), "Electron Multiplicity", "", "Events", log);
//*/
//  }
//
//


  // Balancing
  db->setFolder("analysis/balancing");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
  
    //cout << "after set_legendTitle" << endl << endl;
  
    TString responseName = TString::Format("resp_balancing_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, responseName.Data(), "Balancing Response", "", "Events", log);

    // Raw jets
    responseName = TString::Format("resp_balancing_raw_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, responseName.Data(), "Balancing Response (raw jets)", "", "Events", log);
    
    //  -- federico    
    //Balancing Reco - Gen JET
    //    responseName = TString::Format("resp_balancing_gen_%s", etaBinning.getBinName(i).c_str());
    //    db->drawHisto_vs_pt(ptBins, responseName.Data(), "Balancing Response (Reco-Gen)", "", "Events", log);
    //  -- federico    
    //Balancing Reco - Gen PHOTON
    //    responseName = TString::Format("resp_photGamma_%s", etaBinning.getBinName(i).c_str());
    //    db->drawHisto_vs_pt(ptBins, responseName.Data(), "Balancing Response (Reco-Gen)", "", "Events", log);
    
  }
  
  // Special case eta < 1.3

  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(ptBins, "resp_balancing_eta0013", "Balancing Response", "", "Events", log);
  db->drawHisto_vs_pt(ptBins, "resp_balancing_raw_eta0013", "Balancing Response (raw jets)", "", "Events", log);

  // MPF
  db->setFolder("analysis/mpf");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_mpf_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, responseName.Data(), "MPF Response", "", "Events", log);

    // Raw jets
    responseName = TString::Format("resp_mpf_raw_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, responseName.Data(), "MPF Response (raw ME_{T})", "", "Events", log);

  }
  // Special case eta < 1.3

  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(ptBins, "resp_mpf_eta0013", "MPF Response", "", "Events", log);
  db->drawHisto_vs_pt(ptBins, "resp_mpf_raw_eta0013", "MPF Response (raw ME_{T})", "", "Events", log);

  delete db;
  db = NULL;

  return 0;

}
