#include <stdlib.h>
#include <sys/stat.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"

#include "etaBinning.h"
#include "ptBinning.h"
#include "runBinning.h"

#include <TColor.h>

bool useMCassoc_ = false;
bool ONEVTX = false;
bool OUTPUT_GRAPHS = true;

#define BALANCING TColor::GetColor(217, 91, 67)
#define MPF TColor::GetColor(192, 41, 66)


int main(int argc, char* argv[]) {

  // std::cout<< "Comincio programma "<< std::endl;

  if (argc != 4 && argc != 5) {
    std::cout << "USAGE: ./draw_vs_npv [data_dataset] [recoType] [jetAlgo]" << std::endl;
    exit(23);
  }

  gROOT->SetBatch();

  std::string data_dataset(argv[1]);
  std::string recoType(argv[2]);
  std::string jetAlgo(argv[3]);
  std::string flags = "";
  if (argc == 5) {
    std::string flags_str(argv[4]);
    flags = flags_str;
  }


  std::string algoType;
  if (recoType == "calo") {
    algoType = jetAlgo;
  } else {
    algoType = recoType + jetAlgo;
  }
  if (recoType == "jpt" && jetAlgo == "akt4") {
    algoType = "jptak4";
  }

  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  recoType = (recoType == "pf") ? "PFlow" : "Calo";
  std::string postFix = recoType + jetAlgo;

  postFix += "chs";


  TString dataFileName;
  if (flags.length() > 0) {
    dataFileName = TString::Format("PhotonJet_%s_%s_%s.root", data_dataset.c_str(), postFix.c_str(), flags.c_str());
  } else {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  }

  TFile* dataFile = TFile::Open(dataFileName);

  if (dataFile) {
    std::cout << "Opened data file '" << dataFileName << "'." << std::endl;
  }

  //  TDirectory* dir  = dataFile->cd("analysis/run")


  //  db->setFolder("analysis");
  std::string outputDir = "vs_run";
  mkdir(outputDir.c_str(), 0777); 

  gErrorIgnoreLevel = kWarning;

  RunBinning runBinning;
  size_t runBinningSize = runBinning.size();
  EtaBinning etaBinning;
  size_t etaBinningSize = etaBinning.size();

  std::vector<double> etaBinsError;                                                                                                                                
  etaBinsError.push_back( 0.3915 );                                                                                                                                
  etaBinsError.push_back( 0.261 );                                                                                                                                
  etaBinsError.push_back( 0.3125 );                                                                                                                               
  etaBinsError.push_back( 0.285 );                                                                                                                               
  etaBinsError.push_back( 0.232 );                                                                                                                               
  etaBinsError.push_back( 0.118 );                                                                                                                               
  etaBinsError.push_back( 0.9955 );    
  
  TLegend *leg;
  TGraphErrors* gr_response_vs_run[runBinningSize] ;
  TCanvas *c = new TCanvas("c","",800,800);
  
  // methods: Bal or MPF
  for (size_t kk = 0; kk < 2; kk++) {

    if(kk == 0)   std::cout<<"Balancing method "<< std::endl;
    if(kk == 1)   std::cout<<"MPF method "<< std::endl;
    
    c->cd();
    
    leg=new TLegend(0.15,0.15,0.45,0.45);
    leg->SetFillColor(0);
    //  leg->SetTextFont(42);                                                                                                                                                    
    leg->SetBorderSize(0);                                                                                                                                                   
    leg->SetFillColor(kWhite);                                                                                                                                               
    leg->SetFillStyle(0);                                                                                                                                                    
    //  leg->SetTextSize(0.036);                                                                                                                                                 
    
    
    // Balancing && MPF
    for (size_t jj = 0; jj < runBinningSize; jj++) {
      const std::pair<int, int>  runBin = runBinning.getBinValue(jj);
      
      if(jj == runBinningSize -1 ) continue;      

      gr_response_vs_run[jj] = new TGraphErrors(0);
      

      std::cout<<"Bin run number "<< jj+1 << std::endl;
      
      int run_low  = runBin.first ;
      int run_high = runBin.second ;

      std::cout<<"Run range "<< run_low <<"-"<< run_high << std::endl;
      
      for (size_t i = 0; i < etaBinningSize; i++) {
	
	std::pair<float, float> currentBin =etaBinning.getBinValue(i) ;
	float etaMean = (currentBin.first + currentBin.second) / 2.;
	
	TString responseName;
	
	if(kk==0)  responseName = TString::Format("resp_balancing_%s_run_%i_%i", etaBinning.getBinName(i).c_str(), run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_%s_run_%i_%i", etaBinning.getBinName(i).c_str(), run_low, run_high );
	
	TH1F *h = (TH1F*)dataFile->Get("analysis/run/"+responseName);
	
	double Mean = h->GetMean();
	double MeanError = h->GetMeanError();
	
	std::cout<< "Set point "<< i << " :   "<< etaMean << "  " << Mean<< std::endl;      
	
	gr_response_vs_run[jj]-> SetPoint(i, etaMean, Mean);
	gr_response_vs_run[jj]-> SetPointError(i, etaBinsError.at(i), MeanError);
	
      }//for eta
      
      std::cout<< "Drawing graph"<< std::endl; 
      if(jj < 4){
      gr_response_vs_run[jj] -> SetLineColor(jj+1);
      gr_response_vs_run[jj] -> SetMarkerColor(jj+1);
      }else{
      gr_response_vs_run[jj] -> SetLineColor(jj+2);
      gr_response_vs_run[jj] -> SetMarkerColor(jj+2);
      }
      
      leg->AddEntry(gr_response_vs_run[jj] ,Form("Run %d - %d", run_low, run_high) ,"l");
      
      
      if(jj == 0){
	//      std::cout << jj << std::endl;
	gr_response_vs_run[jj] -> SetTitle();
	gr_response_vs_run[jj] -> GetHistogram()-> SetXTitle("#eta (jet)");  
	if(kk ==0)      gr_response_vs_run[jj] -> GetHistogram()-> SetYTitle("Balancing Response");  
	if(kk ==1)      gr_response_vs_run[jj] -> GetHistogram()-> SetYTitle("MPF Response");  
	gr_response_vs_run[jj] -> GetYaxis()-> SetTitleOffset(1.4);  
	gr_response_vs_run[jj] -> GetYaxis()->SetRangeUser(0.3, 1.2);  
	gr_response_vs_run[jj] -> GetXaxis()->SetRangeUser(0., 5.2);  
	gr_response_vs_run[jj] -> Draw("ZAP");
      }else{
	gr_response_vs_run[jj] ->Draw("PSAME");
    }  
      
    }// for run
    
    leg -> Draw() ;
     
    //  c->SaveAs("vs_run/Balancing_vs_run.png");
    if(kk == 0)  c->SaveAs("vs_run/Balancing_vs_run.png");
    if(kk == 1)  c->SaveAs("vs_run/MPF_vs_run.png");
  }


}


