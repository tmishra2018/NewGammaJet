#include <stdlib.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"


#include "etaBinning.h"
#include "ptBinning.h"

#include <TColor.h>


bool useMCassoc_ = false;
bool ONEVTX = false;
bool OUTPUT_GRAPHS = true;

int main(int argc, char* argv[]) {

  if (argc != 5) {
    std::cout << "USAGE: ./drawPhotonJet [mc_SIGNAL_dataset] [recoType] [jetAlgo] [norm ('LUMI' or 'SHAPE')]" << std::endl;
    exit(23);
  }

  std::string mc_photonjet(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);
  std::string norm(argv[4]);
  if (norm != "LUMI" && norm != "SHAPE") {
    std::cout << "'" << norm << "' normalization not implemented yet." << std::endl;
    std::cout << "Only 'LUMI' and 'SHAPE' currently supported." << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(9811);
  }

  recoType = (recoType == "pf") ? "PFlow" : "PUPPI";
  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  std::string postFix = recoType + jetAlgo;
  if(recoType == "PFlow")  postFix += "chs";
  
  drawBase* db = new drawBase("PhotonJet", recoType, jetAlgo, OUTPUT_GRAPHS);
  db->set_pdf_aussi((bool)false);
  db->set_isCMSArticle(false);

  TString mc1FileName;
  mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened GJet file '" << mc1FileName << "'." << std::endl;
  if (mcPhotonJetFile) {
    db->add_mcFile(mcPhotonJetFile, mc_photonjet,  "MC", TColor::GetColor(192, 41, 66));
  }


  double dLumi = 1.;
  std::cout<< "Lumi "<< dLumi << std::endl;
  db->set_lumi(dLumi);
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
  /*
  db->drawHisto("ptPhoton_Binned", "Photon Transverse Momentum", "GeV", "Events", log, 1, "", false, 50);
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
  // db->set_rebin(1);
  db->drawHisto("nvertex", "Number of Reconstructed Vertices", "", "Events", log);
  db->drawHisto("nvertex_reweighted", "Number of Reconstructed Vertices (after reweighting)", "", "Events", log);
  db->drawHisto("ntrue_interactions", "N true interactions", "", "Events", log);
  db->drawHisto("ntrue_interactions_reweighted", "N true interactions (after reweighting)", "", "Events", log);
  //db->set_rebin(2);
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
  */
  db->setOutputGraphs(OUTPUT_GRAPHS);

  TH1D *ptPhot = (TH1D*)mcPhotonJetFile->Get("analysis/ptPhoton_Binned");
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

  // Balancing
  db->setFolder("analysis/balancing");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    //TString responseName = TString::Format("resp_balancing_%s", etaBinning.getBinName(i).c_str());
    // db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "Balancing Response", "", "Events", log);
    // Balancing Reco - Gen JET
    //    responseName = TString::Format("resp_balancing_gen_%s", etaBinning.getBinName(i).c_str());
    //    db->drawHisto_vs_pt(ptBins, responseName.Data(), "Balancing Response (Reco-Gen)", "", "Events", log);
    // Balancing Reco - Gen PHOTON
    TString responseName = TString::Format("resp_balancing_photGamma_%s", etaBinning.getBinName(i).c_str());
    db->drawHistoGen_vs_pt(ptBins, ptMean, responseName.Data(), "Balancing #gamma Response (Reco/Gen)", "", "Events", false);    
  }
  
  // Special case eta < 1.3
  db->set_legendTitle("|#eta| < 1.3");
  //db->drawHisto_vs_pt(ptBins, ptMean, "resp_balancing_eta0013", "Balancing Response", "", "Events", log);
  // db->drawHisto_vs_pt(ptBins, ptMean, "resp_balancing_raw_eta0013", "Balancing Response (raw jets)", "", "Events", log);
  db->drawHistoGen_vs_pt(ptBins, ptMean, "resp_balancing_photGamma_eta0013", "Balancing #gamma Response (Reco/Gen)", "", "Events", log);

  /*
  // MPF
  db->setFolder("analysis/mpf");
  for (size_t i = 0; i < etaBinningSize; i++) {
    db->set_legendTitle(etaBinning.getBinTitle(i));
    TString responseName = TString::Format("resp_mpf_%s", etaBinning.getBinName(i).c_str());
    db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "MPF Response", "", "Events", log);
    // Raw jets
    // responseName = TString::Format("resp_mpf_raw_%s", etaBinning.getBinName(i).c_str());
    // db->drawHisto_vs_pt(ptBins, ptMean, responseName.Data(), "MPF Response (raw ME_{T})", "", "Events", log);
  }
  // Special case eta < 1.3
  db->set_legendTitle("|#eta| < 1.3");
  db->drawHisto_vs_pt(ptBins, ptMean, "resp_mpf_eta0013", "MPF Response", "", "Events", log);
  // db->drawHisto_vs_pt(ptBins, ptMean, "resp_mpf_raw_eta0013", "MPF Response (raw ME_{T})", "", "Events", log);
*/

  delete db;
  db = NULL;
  return 0;

}
