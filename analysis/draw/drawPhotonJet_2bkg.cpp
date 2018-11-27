#include <stdlib.h>

#include "drawBase.h"
#include "TParameter.h"
#include "TError.h"
#include "fitTools.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include <TColor.h>
#include "etaBinning.h"
#include "ptBinning.h"
#include "HLTptBinning.h"
#include "fineetaBinning.h"


bool OUTPUT_GRAPHS = true;
bool RAW = false;

int main(int argc, char* argv[]) {

  if (argc != 7) {
    std::cout << "USAGE: ./drawPhotonJet [data_dataset] [mc_SIGNAL_dataset] [mc_BG_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')]" << std::endl;
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
  
  recoType = (recoType == "pf") ? "PFlow" : "PUPPI";
  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  std::string postFix = recoType + jetAlgo;  
  if(recoType == "PFlow") postFix += "chs";
  
  drawBase* db = new drawBase("PhotonJet", recoType, jetAlgo, OUTPUT_GRAPHS);
  db->set_pdf_aussi((bool)false);
  db->set_isCMSArticle(false);

  TString dataFileName;
  dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  TFile* dataFile = TFile::Open(dataFileName);

  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
    db->add_dataFile(dataFile, data_dataset);
  }else{
    std::cout << "Impossible open data file '" << dataFileName << "'." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  TString mc1FileName;
  mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);

  if (mcPhotonJetFile) {
    std::cout << "Opened GJet file '" << mc1FileName << "'." << std::endl;
    db->add_mcFile(mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", TColor::GetColor(192, 41, 66));
  }else{
    std::cout << "Impossible open GJet file '" << mc1FileName << "'." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  if (mc_QCD != mc_photonjet){ 
    TString mc2FileName;
    mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc_QCD.c_str(), postFix.c_str());
    TFile* mcQCDFile = TFile::Open(mc2FileName);
 
    if (mcQCDFile) {
      std::cout << "Opened mc QCD file '" << mc2FileName << "'." << std::endl;
      db->add_mcFile(mcQCDFile, mc_QCD, "QCD MC", TColor::GetColor("#ECD078"));
    }else{
      std::cout << "Impossible open QCD file '" << mc2FileName << "'." << std::endl;
      std::cout << "Exiting." << std::endl;
      exit(9811);
    }   
  }else{
    std::cout<<"Not passed any QCD file "<<std::endl;
  }
  
  // MC should already be normalized to a lumi of 1 pb-1
  // Read luminosity
  double dLumi=1;
  TParameter<float>* lumi = static_cast<TParameter<float>*>(dataFile->Get("analysis/luminosity"));
  dLumi = lumi->GetVal();   
  std::cout<< "Lumi "<< dLumi << std::endl;  
  db->set_lumi(dLumi);
  
  
  
  
  
  
  
  if (norm == "LUMI") {
    db->set_lumiNormalization();
  } else {
    db->set_shapeNormalization();
  }
  
  bool log = true;
  gErrorIgnoreLevel = kWarning;  
  db->setFolder("analysis");
  db->set_outputdir();
  db->setOutputGraphs(true);
  db->set_rebin(1);

  // Data/MC comparison
  db->drawHisto("ptPhoton", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_NoCut", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
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
  db->drawHisto("MET_parr", "Missing E_{T //} ", "GeV", "Events", log);
  db->drawHisto("MET_ortho", "Missing E_{T #perp}", "GeV", "Events", log);
  db->drawHisto("alpha_passedID", "#alpha", "", "Events", log);
  //  db->set_rebin(1);
  db->drawHisto("nvertex", "Number of Reconstructed Vertices", "", "Events", log);
  db->drawHisto("nvertex_reweighted", "Number of Reconstructed Vertices (after reweighting)", "", "Events", log);
 // db->drawHisto("ntrue_interactions", "Number of True interaction", "", "Events", log);
  //db->drawHisto("ntrue_interactions_reweighted", "Number of True interaction (after reweighting)", "", "Events", log);
  
  //  db->set_rebin(2);
  db->drawHisto("deltaPhi_passedID", "#Delta #varphi", "", "Events", log);  
  db->drawHisto("hadTowOverEm", "H/E", "", "Events", log);
  db->drawHisto("sigmaIetaIeta", "#sigma_{i#eta i#eta}", "", "Events", log);
  db->drawHisto("rho", "#rho", "", "Events", log);
  db->drawHisto("chargedHadronsIsolation", "Charged hadrons isolation", "", "Events", log);
  db->drawHisto("neutralHadronsIsolation", "Neutral hadrons isolation", "", "Events", log);
  db->drawHisto("photonIsolation", "Photon isolation", "", "Events", log);
  db->drawHisto("hadTowOverEm_passedID", "H/E", "", "Events", log);
  db->drawHisto("sigmaIetaIeta_passedID", "#sigma_{i#eta i#eta}", "", "Events", log);
  db->drawHisto("rho_passedID", "#rho", "", "Events", log);
  db->drawHisto("chargedHadronsIsolation_passedID", "Charged hadrons isolation", "", "Events", log);
  db->drawHisto("neutralHadronsIsolation_passedID", "Neutral hadrons isolation", "", "Events", log);
  db->drawHisto("photonIsolation_passedID", "Photon isolation", "", "Events", log);

  db->drawProfile("PhotSCPt", "Pt", "p_{T} [GeV]", "p_{T}(SC) [GeV]", log, 0);
  db->drawHisto2D("PhotonSCPt_vs_Pt", "p_{T} [GeV]", "p_{T}(SC) [GeV]","", log);
  db->drawHisto2D("deltaPhi_vs_alpha", "#delta #phi ", "#alpha","", log);
 
 // db->set_rebin(5);
  db->setOutputGraphs(OUTPUT_GRAPHS);

  TH1D *ptPhot = (TH1D*)dataFile->Get("analysis/ptPhoton_passedID");
  PtBinning ptBinning;
  size_t ptBinningSize = ptBinning.size();
  std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();
  std::vector<float> ptMean;
  for( size_t i = 0 ; i< ptBinningSize ; i++){ 
    std::pair<float, float> currentBin = ptBinning.getBinValue(i);
    ptPhot ->GetXaxis()->SetRangeUser(currentBin.first, currentBin.second);
    double Mean = ptPhot->GetMean();
    std::cout<< "Bin " << currentBin.first<< "-"<<currentBin.second<<" -> Mean  "<< Mean << std::endl; 
    ptMean.push_back(Mean);
  }

  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();
  std::vector<EtaBin > etaBins = etaBinning.getBinning();
  
  fineEtaBinning fineetaBinning;
  size_t fineetaBinningSize = fineetaBinning.size();
  std::vector<fineEtaBin > fineetaBins = fineetaBinning.getBinning();
  
  
  HLTPtBinning hltptBinning;
  size_t hltptBinningSize = hltptBinning.size();
  std::vector<std::pair<float, float> > hltptBins = hltptBinning.getBinning();
  std::vector<float> hltptMean;
  for( size_t i = 0 ; i< hltptBinningSize ; i++){ 
    std::pair<float, float> currentBin = hltptBinning.getBinValue(i);
    ptPhot ->GetXaxis()->SetRangeUser(currentBin.first, currentBin.second);
    double Mean = ptPhot->GetMean();
    std::cout<< "Bin " << currentBin.first<< "-"<<currentBin.second<<" -> Mean  "<< Mean << std::endl; 
    hltptMean.push_back(Mean);
  }
  
  
  // HLT plots
  /*
  db->setFolder("analysis/HLT_CH_iso");

  
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "ChHadronisoHLT_eta0013", "Charge hadron isolation", "", "Events", true);
  
  db->setFolder("analysis/HLT_NH_iso");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("NhHadronisoHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "Neutral hadron isolation", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "NhHadronisoHLT_eta0013", "Neutral hadron isolation", "", "Events", false);
  
  db->setFolder("analysis/HLT_Photon_iso");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("PhotonisoHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "Photon isolation", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "PhotonisoHLT_eta0013", "Photon isolation", "", "Events", false);
  
  db->setFolder("analysis/HLT_sigieta");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("sigmaIetaIetaHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "#sigma _{i#etai#eta}", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "sigmaIetaIetaHLT_eta0013", "#sigma _{i#etai#eta}", "", "Events", false);
  
  db->setFolder("analysis/HLT_HoverE");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("HoverEHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "#frac{H}{E}", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "HoverEHLT_eta0013", "#frac{H}{E}", "", "Events", false);
  
  db->setFolder("analysis/HLT_rho");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("rhoHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "rho", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "rhoHLT_eta0013", "rho", "", "Events", false);
  
  
  db->setFolder("analysis/HLT_metparr");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("metparrHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "metparr", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "metparrHLT_eta0013", "metparr", "", "Events", false);
  
  db->setFolder("analysis/HLT_metperp");

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("metperpHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "metperp", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "metperpHLT_eta0013", "metperp", "", "Events", false);
  
  
  db->setFolder("analysis/HLT_jetpt");
   for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("jetptHLT_%s", etaBinning.getBinName(i).c_str());
    TString responseName2 = TString::Format("jet_2ptHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "jet pt", "", "Events", false);
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName2.Data(), "jet pt", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "jetptHLT_eta0013", "jet pt", "", "Events", false);
  db->drawHisto_vs_pt(hltptBins, hltptMean, "jet_2ptHLT_eta0013", "jet pt", "", "Events", false);
  
  
  db->setFolder("analysis/HLT_met");
   for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("metpHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "MET ", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "metpHLT_eta0013", "MET", "", "Events", false);
  
  db->setFolder("analysis/HLT_Nvertex");
   for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("NvertexHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "npv ", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "NvertexHLT_eta0013", "npv", "", "Events", false);
  
  db->setFolder("analysis/HLT_alpha");
   for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("alphaHLT_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(hltptBins, hltptMean, responseName.Data(), "#alpha ", "", "Events", false);
    
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(hltptBins, hltptMean, "alphaHLT_eta0013", "#alpha", "", "Events", false);
  */
  
  // fine binning eta plots
  
  db->setFolder("analysis/fine_eta_binning/balancing");
  
  for (size_t i = 0; i < fineetaBinningSize; i++) {
    db->set_legendTitle(fineetaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_balancing_fine_bining_%s", fineetaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Balancing Response", "", "Events", false);
   
  }
  
  db->setFolder("analysis/fine_eta_binning/mpf");
  
  for (size_t i = 0; i < fineetaBinningSize; i++) {
    db->set_legendTitle(fineetaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_mpf_fine_bining_%s", fineetaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "MPF Response", "", "Events", false);
   
  }
  
  
  
  
  // Balancing
  db->setFolder("analysis/balancing");
  db->drawProfile("Bal", "Pt", "p_{T} [GeV]", "Balancing Response", log, 0);
  db->drawProfile("Bal", "Eta", "#eta", "Balancing Response", false, 0); 
  db->drawProfile("Bal", "Nvtx", "N_{vertex}", "Balancing Response", false, 0);
  
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("resp_balancing_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Balancing Response", "", "Events", false);
    if(RAW){ // Raw jets
      responseName = TString::Format("resp_balancing_raw_%s", etaBinning.getBinName(i).c_str());
      db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Balancing Response (raw jets)", "", "Events", false);
    }
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(ptBins, ptMean, "resp_balancing_eta0013", "Balancing Response", "", "Events", false);
  if(RAW) db->drawHisto_vs_pt(ptBins, ptMean, "resp_balancing_raw_eta0013", "Balancing Response (raw jets)", "", "Events", false);
  
  
  
  
  
  
  // MPF
  db->setFolder("analysis/mpf");
  db->drawProfile("MPF", "Pt", "p_{T} [GeV]", "MPF Response", log, 0);
  db->drawProfile("MPF", "Eta", "#eta", "MPF Response", false, 0);
  db->drawProfile("MPF", "Nvtx", "N_{vertex}", "MPF Response", false, 0); 

  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
   
    TString responseName = TString::Format("resp_mpf_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "MPF Response", "", "Events", false);
    if(RAW){ // Raw jets
      responseName = TString::Format("resp_mpf_raw_%s", etaBinning.getBinName(i).c_str());
      db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "MPF Response (raw MET)", "", "Events", false);
    }
  }
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(ptBins, ptMean, "resp_mpf_eta0013", "MPF Response", "", "Events", false);
  if(RAW) db->drawHisto_vs_pt(ptBins, ptMean, "resp_mpf_raw_eta0013", "MPF Response (raw MET)", "", "Events", false);
  
  db->setFolder("analysis");
  db->set_shapeNormalization();
 // db->drawHisto("ptPhoton_175_3000", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_130_175", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_105_130", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_85_105", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_60_85", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  db->drawHisto("ptPhoton_40_60", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
  
  /*
  db->set_lumiNormalization();
  
  db->setFolder("analysis/balancing_vs_eta");
  db->set_legendTitle("p_{T}^{#gamma} > 175 GeV ");
  db->drawHisto_vs_eta(etaBins, "resp_balancing_pt175", "Balancing Response", "", "Events", false);
  
  db->setFolder("analysis/mpf_vs_eta");
  db->set_legendTitle("p_{T}^{#gamma} > 175 GeV ");
  db->drawHisto_vs_eta(etaBins, "resp_mpf_pt175", "MPF Response", "", "Events", false);
  */
  
  

  delete db;
  db = NULL;

  return 0;

}


/*
  
  // Pt 1st
  db->setFolder("analysis/Ptfirstjets");
  db->drawProfile("Pt_1stjet", "Pt", "p_{T} [GeV]", "Pt_1stjet", log, 0);
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("Pt_1stjet_%s", etaBinning.getBinName(i).c_str());
   // bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Pt_1stjet", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "Pt_1stjet_eta0013", "Pt_1stjet", "", "Events", false);
  
  // Pt 2nd
  db->setFolder("analysis/Ptsecondjets");
  db->drawProfile("Pt2nd_jet", "Pt", "p_{T} [GeV]", "Pt2nd_jet", log, 0);
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("Pt_2ndjet_%s", etaBinning.getBinName(i).c_str());
   // bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Pt_2ndjet", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "Pt_2ndjet_eta0013", "met", "", "Events", false);
  ;
  // Pt met
  db->setFolder("analysis/Met");
  db->drawProfile("Met", "Pt", "p_{T} [GeV]", "Met", log, 0);
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("met_%s", etaBinning.getBinName(i).c_str());
   // bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "met", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "met_eta0013", "met", "", "Events", false);
  
  
  
  // Pt gamma
  db->setFolder("analysis/Ptgamma");
  db->drawProfile("Ptgamma", "Pt", "p_{T} [GeV]", "Pt gamma", log, 0);
  db->drawProfile("Ptgamma", "Eta", "#eta", "Pt gamma", false, 0); 
  db->drawProfile("Ptgamma", "Nvtx", "N_{vertex}", "Pt gamma", false, 0);
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("Pt_gamma_%s", etaBinning.getBinName(i).c_str());
   // bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Pt gamma", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "Pt_gamma_eta0013", "Pt gamma", "", "Events", false);

  
  
  // Mu
  db->setFolder("analysis/MUDir");
  
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("mu_%s", etaBinning.getBinName(i).c_str());
  //  bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "mu", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "mu_eta0013", "mu", "", "Events", false);
  
  
  // N vertice 
  db->setFolder("analysis/nvertices");
  
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    
    TString responseName = TString::Format("nvertices_%s", etaBinning.getBinName(i).c_str());
  //  bool isgen = TString(etaBinning.getBinName(i).c_str()).Contains("gen", TString::kIgnoreCase);

    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "nvertices", "", "Events", false);

  }

  db->set_legendTitle("|#eta| < 1.3");

  db->drawHisto_vs_pt(ptBins, ptMean, "nvertices_eta0013", "mu", "", "Events", false);
  */
