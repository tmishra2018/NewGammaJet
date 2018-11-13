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
  TGraphErrors* gr_response_vs_run[6] ;
  
  
  TCanvas *c = new TCanvas("c","",800,800);
  
  // methods: Bal or MPF
  for (size_t kk = 0; kk < 2; kk++) {

    if(kk == 0)   std::cout<<"Balancing method "<< std::endl;
    if(kk == 1)   std::cout<<"MPF method "<< std::endl;
    
    c->cd();
    
    leg=new TLegend(0.2,0.2,0.8,0.3);
    leg->SetFillColor(0);
    //  leg->SetTextFont(42);                                                                                                                                                    
    leg->SetBorderSize(0);                                                                                                                                                   
    leg->SetFillColor(kWhite);                                                                                                                                               
    leg->SetFillStyle(0);
    leg->SetNColumns(6);                                                                                                                                                    
    //  leg->SetTextSize(0.036);                                                                                                                                                 
    
    TH2D* h2_axes = new TH2D("axes_again", "", 31, 0., 30., 10, 0.6, 1.1);
   // h2_axes->SetXTitle("p_{T}(#gamma) [GeV/c]");
    if(kk == 0){
      h2_axes->SetYTitle("p_{T} Balance");
     }else{
    
    h2_axes->SetYTitle("MPF response");    
    }
  h2_axes->GetYaxis()->SetRangeUser(0.8,1.05);  
  h2_axes->SetXTitle("Run Ranges");
  h2_axes->SetStats(0);
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleSize(0.045);
  //h2_axes->GetXaxis()->SetMoreLogLabels();
  //h2_axes->GetXaxis()->SetNoExponent();
    h2_axes->GetXaxis()->SetLabelSize(0.015);
    
    TLine* BCD_iov = new TLine(4,0.8,4,1 ) ;
    BCD_iov->SetVertical();
    BCD_iov->SetLineStyle (9);
    BCD_iov->SetLineColor(kBlue);
    
    TLine* EF_iov = new TLine(13,0.8,13,1 ) ;
    EF_iov->SetVertical();
    EF_iov->SetLineStyle (9);
    EF_iov->SetLineColor(kBlue);
    
    TLine* FG_iov = new TLine(18,0.8,18,1 ) ;
    FG_iov->SetVertical();
    FG_iov->SetLineStyle (9);
    FG_iov->SetLineColor(kBlue);
    
    TLine* E_iov = new TLine(26,0.8,26,1 ) ;
    E_iov->SetVertical();
    E_iov->SetLineStyle (9);
    E_iov->SetLineColor(kBlue);
    
    
    
    
    // Balancing && MPF
    for(int ii = 0 ; ii < 6 ; ii++){
    int HLTnumber = 0 ;
    gr_response_vs_run[ii] = new TGraphErrors(0);
    for (size_t jj = 0; jj < runBinningSize; jj++) {
      const std::pair<int, int>  runBin = runBinning.getBinValue(jj);
      
     // if(jj == runBinningSize -1 ) continue;      

      int run_low = runBin.first;
      int run_high = runBin.second;
      
      TString runBinName;
      runBinName = TString::Format("%i-%i", run_low,run_high);
      h2_axes->GetXaxis()->SetBinLabel(jj + 1 , runBinName);
      
      std::cout<<"Run range "<< run_low <<"-"<< run_high << std::endl;
	
	TString responseName;
	if(ii == 0){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT30_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT30_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber = 30;
	}
	if(ii == 1){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT50_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT50_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber =50;
	}
	if(ii == 2){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT75_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT75_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber = 75 ;
	}
	if(ii == 3){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT90_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT90_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber = 90 ; 
	}
	if(ii == 4){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT120_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT120_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber = 120 ;
	}
	if(ii == 5){
	if(kk==0)  responseName = TString::Format("resp_balancing_HLT165_eta0013_Run_%i_%i", run_low, run_high );
	if(kk==1)  responseName = TString::Format("resp_mpf_HLT165_eta0013_Run_%i_%i", run_low, run_high );
	HLTnumber = 165 ;
	}
	
	
	
	TH1F *h = (TH1F*)dataFile->Get("analysis/run/"+responseName);
	std::cout<<"histo name "<< responseName << std::endl;
	Double_t Mean = 0;
	Double_t MeanError = 0 ;
	
	if ( h->GetEntries() != 0 ) Mean = h->GetMean();
	if ( h->GetEntries() != 0 ) MeanError = h->GetMeanError();
	Double_t highrun =  run_high;
	Double_t lowrun  =  run_low;
	Double_t meanrunBin = (highrun - lowrun)/2. + lowrun ;
	
	std::cout<<" graph number "<<ii<< " Set point "<< jj   << " :   "<< meanrunBin << "  " << Mean<< " +- "<<MeanError<< std::endl;      
	
	gr_response_vs_run[ii]-> SetPoint(jj , jj+0.5 , Mean);
	gr_response_vs_run[ii]-> SetPointError(jj , 0., MeanError);
	

        
      
    }// for run
    
    if(ii > 0){
    std::cout<< "Drawing graph  "<<gr_response_vs_run[ii]->GetN()<< std::endl; 
    leg->AddEntry(gr_response_vs_run[ii] ,Form("HLT%d ", HLTnumber ) ,"p");
    }
    }
  
 
       h2_axes->Draw();
    
    for(size_t ii = 0; ii<6  ;ii++){
    
      if(ii < 4){
      gr_response_vs_run[ii] -> SetLineColor(ii+1);
      gr_response_vs_run[ii] -> SetMarkerColor(ii+1);
      gr_response_vs_run[ii] -> SetMarkerStyle(ii+20);
      gr_response_vs_run[ii] -> SetMarkerSize(2);
      }else{
      gr_response_vs_run[ii] -> SetLineColor(ii+2);
      gr_response_vs_run[ii] -> SetMarkerColor(ii+2);
      gr_response_vs_run[ii] -> SetMarkerStyle(ii+20);
      gr_response_vs_run[ii] -> SetMarkerSize(1);
      }
      
      gr_response_vs_run[5] -> SetLineColor(1);
      gr_response_vs_run[5] -> SetMarkerColor(1);
      gr_response_vs_run[5] -> SetMarkerStyle(8);
      gr_response_vs_run[5] -> SetMarkerSize(1);
      
      
      
	//      std::cout << jj << std::endl;
	//gr_response_vs_run[ii] -> SetTitle();
	//gr_response_vs_run[ii] -> GetHistogram()-> SetXTitle("Run Number");  
	//if(kk ==0)      gr_response_vs_run[ii] -> GetHistogram()-> SetYTitle("Balancing Response");  
	//if(kk ==1)      gr_response_vs_run[ii] -> GetHistogram()-> SetYTitle("MPF Response");  
	//gr_response_vs_run[ii] -> GetYaxis()-> SetTitleOffset(1.4);  
	//gr_response_vs_run[ii] -> GetYaxis()->SetRangeUser(0.8, 0.95);  
	//gr_response_vs_run[ii] -> GetXaxis()->SetRangeUser(297045, 300575);  
	if(ii > 0) gr_response_vs_run[ii] -> Draw("ZPSAME");
	BCD_iov->Draw("DSAME");
	EF_iov->Draw("DSAME");
	FG_iov->Draw("DSAME");
	E_iov->Draw("DSAME");
      
    }
    leg -> Draw("P") ;
     


    if(kk == 0)  c->SaveAs("vs_run/Balancing_vs_run.png");
    if(kk == 1)  c->SaveAs("vs_run/MPF_vs_run.png");
  }


}


