#include <stdint.h>

#include "drawExtrap.h"
#include "fitTools.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include <TColor.h>

#define BALANCING_COLOR TColor::GetColor(217, 91, 67)
#define MPF_COLOR TColor::GetColor(192, 41, 66)
#define LINE_COLOR TColor::GetColor("#0B486B")
#define OUTPUT_GRAPHS true

drawExtrap::drawExtrap(const std::string& analysisType, const std::string& recoType, const std::string& jetAlgo, bool outputGraphs, const std::string& flags) : drawBase(analysisType, recoType, jetAlgo, outputGraphs, flags) {
  
  FIT_RMS_ = "FIT";
  NOQ_ = false;
  INTPERC_ = 95.;
  FIXM_ = false;
  EXCLUDE_FIRST_POINT_ = false;
  
  if(recoType == "PFlow") std::cout<<"recoType PFlow" << std::endl;
  mExtrapBinning.initialize();
}


void drawExtrap::drawResponseExtrap(std::vector<float> ptMeanVec, std::vector<float> alphaMeanVec, const std::string& etaRegion, const std::string& etaRegion_str, bool rawJets) {
  
  int MC_color = MPF_COLOR;
  
  std::vector<std::pair<float, float> > ptPhot_binning = mPtBinning.getBinning();
  //Response
  //DATA
  TGraphErrors* gr_DATAResp_vs_pt = new TGraphErrors(0);
  gr_DATAResp_vs_pt->SetName("gr_DATAResp_vs_pt");
  gr_DATAResp_vs_pt->SetMarkerStyle(25);
  gr_DATAResp_vs_pt->SetMarkerColor(kBlack);
  gr_DATAResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors* gr_DATARespMPF_vs_pt = new TGraphErrors(0);
  gr_DATARespMPF_vs_pt->SetName("gr_DATARespMPF_vs_pt");
  gr_DATARespMPF_vs_pt->SetMarkerStyle(25);
  gr_DATARespMPF_vs_pt->SetMarkerColor(kBlack);
  gr_DATARespMPF_vs_pt->SetMarkerSize(1.5);
  //MC
  TGraphErrors* gr_extrapResp_vs_pt = new TGraphErrors(0);
  gr_extrapResp_vs_pt->SetName("gr_extrapResp_vs_pt");
  gr_extrapResp_vs_pt->SetMarkerStyle(25);
  gr_extrapResp_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors*  gr_extrapRespMPF_vs_pt = new TGraphErrors(0);
  gr_extrapRespMPF_vs_pt->SetName("gr_extrapRespMPF_vs_pt");
  gr_extrapRespMPF_vs_pt->SetMarkerStyle(25);
  gr_extrapRespMPF_vs_pt->SetMarkerColor(kBlack);
  gr_extrapRespMPF_vs_pt->SetMarkerSize(1.5);
  // DATA / MC ratio 
  TGraphErrors* gr_RatioResp_vs_pt = new TGraphErrors(0);
  gr_RatioResp_vs_pt->SetName("gr_RatioResp_vs_pt");
  gr_RatioResp_vs_pt->SetMarkerStyle(25);
  gr_RatioResp_vs_pt->SetMarkerColor(kBlack);
  gr_RatioResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors* gr_RatioMPFResp_vs_pt = new TGraphErrors(0);
  gr_RatioMPFResp_vs_pt->SetName("gr_RatioMPFResp_vs_pt");
  gr_RatioMPFResp_vs_pt->SetMarkerStyle(25);
  gr_RatioMPFResp_vs_pt->SetMarkerColor(kBlack);
  gr_RatioMPFResp_vs_pt->SetMarkerSize(1.5);


  // Particle Level Imbalance
  TGraphErrors* gr_extrapPLI_vs_pt = new TGraphErrors(0);
  gr_extrapPLI_vs_pt->SetName("gr_extrapPLI_vs_pt");
  gr_extrapPLI_vs_pt->SetMarkerStyle(25);
  gr_extrapPLI_vs_pt->SetMarkerColor(kBlack);
  gr_extrapPLI_vs_pt->SetMarkerSize(1.5);
  
  
  
  TGraphErrors* gr_extrapPLIreso_vs_pt = new TGraphErrors(0);
  gr_extrapPLIreso_vs_pt->SetName("gr_extrapPLIreso_vs_pt");
  gr_extrapPLIreso_vs_pt->SetMarkerStyle(25);
  gr_extrapPLIreso_vs_pt->SetMarkerColor(kBlack);
  gr_extrapPLIreso_vs_pt->SetMarkerSize(1.5);


  // Resolution
  // DATA
  TGraphErrors* gr_DATAReso_vs_pt = new TGraphErrors(0);
  gr_DATAReso_vs_pt->SetName("gr_DATAReso_vs_pt");
  gr_DATAReso_vs_pt->SetMarkerStyle(25);
  gr_DATAReso_vs_pt->SetMarkerColor(kBlack);
  gr_DATAReso_vs_pt->SetMarkerSize(1.5);
  TGraphErrors* gr_DATAResoMPF_vs_pt = new TGraphErrors(0);
  gr_DATAResoMPF_vs_pt->SetName("gr_DATAResoMPF_vs_pt");
  gr_DATAResoMPF_vs_pt->SetMarkerStyle(25);
  gr_DATAResoMPF_vs_pt->SetMarkerColor(kBlack);
  gr_DATAResoMPF_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_DATAReso_vs_pt_PLI_not_sub = new TGraphErrors(0);
  gr_DATAReso_vs_pt_PLI_not_sub->SetName("gr_DATAReso_vs_pt_PLI_not_sub");
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerStyle(25);
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerColor(kBlack);
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerSize(1.5);
  
  
  //MC
  TGraphErrors* gr_extrapReso_vs_pt = new TGraphErrors(0);
  gr_extrapReso_vs_pt->SetName("gr_extrapReso_vs_pt");
  gr_extrapReso_vs_pt->SetMarkerStyle(25);
  gr_extrapReso_vs_pt->SetMarkerColor(kBlack);
  gr_extrapReso_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_extrapReso_vs_pt_PLI_not_sub = new TGraphErrors(0);
  gr_extrapReso_vs_pt_PLI_not_sub->SetName("gr_extrapReso_vs_pt_PLI_not_sub");
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerStyle(25);
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerColor(kBlack);
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerSize(1.5);
  
  
  TGraphErrors* gr_extrapResoMPF_vs_pt = new TGraphErrors(0);
  gr_extrapResoMPF_vs_pt->SetName("gr_extrapResoMPF_vs_pt");
  gr_extrapResoMPF_vs_pt->SetMarkerStyle(25);
  gr_extrapResoMPF_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResoMPF_vs_pt->SetMarkerSize(1.5);
  //DATA/MC ratio
  TGraphErrors* gr_RatioReso_vs_pt = new TGraphErrors(0);
  gr_RatioReso_vs_pt->SetName("gr_MCRatioReso_vs_pt");
  gr_RatioReso_vs_pt->SetMarkerStyle(25);
  gr_RatioReso_vs_pt->SetMarkerColor(kBlack);
  gr_RatioReso_vs_pt->SetMarkerSize(1.5);

  std::string suffix = get_fullSuffix();
  std::string graphFileName = "PhotonJetExtrapGraphs_" + suffix + "_" + etaRegion + ((rawJets) ? "RAW" : "") + "_" + FIT_RMS_;

  if (NOQ_) graphFileName += "_NOQ";
  if (FIXM_) graphFileName += "_FIXM";
  if (EXCLUDE_FIRST_POINT_) graphFileName += "_NOFIRSTP";

  graphFileName += ".root";

  TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
  graphFile->cd();

  // To draw also the last PTbins -> ptPhot_binning.size  ---- with -XX not draws last XX bins
  for (uint32_t iPtBin = 0; iPtBin < (ptPhot_binning.size()); //-3 instead of -1 (extrap reaches up to ~2 less bins in pt wrt balancing)
       ++iPtBin) {

    std::pair<float, float> currentBin = mPtBinning.getBinValue(iPtBin);
    float ptMin = currentBin.first;
    float ptMax = currentBin.second;

    std::string rawPostfix = (rawJets) ? "_raw" : "";

    char projName[100];
    sprintf(projName, "projection_%d", iPtBin);

    float ptPhotReco_thisBin =  ptMeanVec.at(iPtBin);
    float ptPhotReco_err_thisBin = 0;

    // npoints in alpha : scan in alpha
    int nPoints = mExtrapBinning.size();
    std::cout<<"binning size : "<<nPoints<<std::endl;
    float x[nPoints];
    float x_err[nPoints];
    getXPoints(iPtBin, x, x_err);
    
    for (int il = 0 ; il < nPoints ; il++){
    x[il]=alphaMeanVec.at(il);
    x_err[il]=0.;
    
    }
    
    
   /* x[0]=0.1;
    x_err[0]=0.0;
    
    x[1]=0.15;
    x_err[1]=0.0;
    
    x[2]=0.2;
    x_err[2]=0.0;
    
    x[3]=0.25;
    x_err[3]=0.0;
    
    x[4]=0.3;
    x_err[4]=0.0;*/
    Float_t y_resp_DATA[nPoints];
    Float_t y_resp_err_DATA[nPoints];

    Float_t y_resp_MC[nPoints];
    Float_t y_resp_MC_err[nPoints];

    Float_t y_resp_MPFDATA[nPoints];
    Float_t y_resp_err_MPFDATA[nPoints];

    Float_t y_resp_MPFMC[nPoints];
    Float_t y_resp_MPFMC_err[nPoints];
    
    //PLI
    Float_t y_PLI[nPoints];
    Float_t y_PLI_err[nPoints];
    
    Float_t y_reso_PLI[nPoints];
    Float_t y_reso_PLI_err[nPoints];

    //Resolution
    Float_t y_reso_DATA[nPoints];
    Float_t y_reso_err_DATA[nPoints];

    Float_t y_reso_MC[nPoints];
    Float_t y_reso_MC_err[nPoints];

    Float_t y_reso_MPFDATA[nPoints];
    Float_t y_reso_err_MPFDATA[nPoints];

    Float_t y_reso_MPFMC[nPoints];
    Float_t y_reso_MPFMC_err[nPoints];
    
    Float_t y_Nevents[nPoints];
    Float_t y_Nevents_MC[nPoints];
    Float_t y_NeventsMPF[nPoints];
    Float_t y_Nevents_MCMPF[nPoints];
    
     Float_t y_NeventsPLI[nPoints];

    TString yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_resp_balancing%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    
    // infinite loop for amc at NLO start
    getYPoints(get_dataFile(0), yHistoName, nPoints, y_resp_DATA, y_resp_err_DATA,  y_reso_DATA, y_reso_err_DATA,y_Nevents);
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_resp_MC, y_resp_MC_err,  y_reso_MC, y_reso_MC_err,y_Nevents_MC);
    // infinite loop for amc at NLO stop
    yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_resp_mpf%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    getYPoints(get_dataFile(0), yHistoName, nPoints, y_resp_MPFDATA, y_resp_err_MPFDATA,  y_reso_MPFDATA, y_reso_err_MPFDATA,y_NeventsMPF);
    
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_resp_MPFMC, y_resp_MPFMC_err,  y_reso_MPFMC, y_reso_MPFMC_err,y_Nevents_MCMPF);

     //PLI
     yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_PLI%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_PLI, y_PLI_err,  y_reso_PLI, y_reso_PLI_err,y_NeventsPLI);
    
    int is_empty = 0;
    for(int i = 1; i < nPoints ; i++){
          if(y_resp_DATA[i]== -1){
             is_empty++;
          }      
    }
    
    if(iPtBin < 3){
    	 x[0]= (alphaMeanVec.at(0)*y_Nevents[0]+alphaMeanVec.at(1)*y_Nevents[1]+alphaMeanVec.at(3)*y_Nevents[3])/(y_Nevents[0]+y_Nevents[1]); // (mExtrapBinning.getBinValue(is_empty).second) / 2. ;
    	 y_resp_DATA[0] = (y_resp_DATA[0]*y_Nevents[0] + y_resp_DATA[1]*y_Nevents[1] + y_resp_DATA[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_resp_MC[0] = (y_resp_MC[0]*y_Nevents_MC[0] + y_resp_MC[1]*y_Nevents_MC[1] + y_resp_MC[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_reso_DATA[0] = (y_reso_DATA[0]*y_Nevents[0] + y_reso_DATA[1]*y_Nevents[1] + y_reso_DATA[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_reso_MC[0] = (y_reso_MC[0]*y_Nevents_MC[0] + y_reso_MC[1]*y_Nevents_MC[1] + y_reso_MC[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_resp_err_DATA[0] = (y_resp_err_DATA[0]*y_Nevents[0] + y_resp_err_DATA[1]*y_Nevents[1] + y_resp_err_DATA[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_resp_MC_err[0] = (y_resp_MC_err[0]*y_Nevents_MC[0] + y_resp_MC_err[1]*y_Nevents_MC[1] + y_resp_MC_err[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_reso_err_DATA[0] = (y_reso_err_DATA[0]*y_Nevents[0] + y_reso_err_DATA[1]*y_Nevents[1] + y_reso_err_DATA[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_reso_MC_err[0] = (y_reso_MC_err[0]*y_Nevents_MC[0] + y_reso_MC_err[1]*y_Nevents_MC[1] + y_reso_MC_err[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_PLI[0] = (y_PLI[0]*y_Nevents_MC[0] + y_PLI[1]*y_Nevents_MC[1] + y_PLI[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;   	 
    	 y_PLI_err[0] = (y_PLI_err[0]*y_Nevents[0] + y_PLI_err[1]*y_Nevents[1] + y_PLI_err[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 
    	 y_reso_PLI[0] = (y_reso_PLI[0]*y_Nevents_MC[0] + y_reso_PLI[1]*y_Nevents_MC[1] + y_reso_PLI[2]*y_Nevents_MC[2])/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;   	 
    	 y_reso_PLI_err[0] = (y_reso_PLI_err[0]*y_Nevents[0] + y_reso_PLI_err[1]*y_Nevents[1] + y_reso_PLI_err[2]*y_Nevents[2])/(y_Nevents[0]+y_Nevents[1]) ;
    	 
    	 
    	 y_resp_DATA[1] = -1 ;
    	 y_resp_MC[1] = -1 ;
    	 
    	 y_reso_DATA[1] = -1 ;
    	 y_reso_MC[1] = -1 ;
    	 
    	 y_resp_err_DATA[1] = 100000000 ;
    	 y_resp_MC_err[1] = 1000000000 ;
    	 
    	 y_reso_err_DATA[0] = 100000000 ;
    	 y_reso_MC_err[0] = 10000000 ;
    	 
    	 y_PLI[1] = -1;   	 
    	 y_PLI_err[1] = -1;
    	 
    	 y_reso_PLI[1] = 10000000 ;   	 
    	 y_reso_PLI_err[1] = 1000000;
    	 
    	  y_resp_DATA[2] = -1 ;
    	 y_resp_MC[2] = -1 ;
    	 
    	 y_reso_DATA[2] = -1 ;
    	 y_reso_MC[2] = -1 ;
    	 
    	 y_resp_err_DATA[2] = 100000000 ;
    	 y_resp_MC_err[2] = 1000000000 ;
    	 
    	 y_reso_err_DATA[2] = 100000000 ;
    	 y_reso_MC_err[2] = 10000000 ;
    	 
    	 y_PLI[2] = -1;   	 
    	 y_PLI_err[2] = -1;
    	 
    	 y_reso_PLI[2] = 10000000 ;   	 
    	 y_reso_PLI_err[2] = 1000000;
    	 
    }
    //draw response histograms:
    TGraphErrors* gr_resp_DATA = new TGraphErrors(nPoints, x, y_resp_DATA, x_err, y_resp_err_DATA);
    gr_resp_DATA->SetMarkerStyle(20);
    gr_resp_DATA->SetMarkerColor(MC_color);
    gr_resp_DATA->SetLineColor(MC_color);

    TGraphErrors* gr_resp_MC = new TGraphErrors(nPoints, x, y_resp_MC, x_err, y_resp_MC_err);
    gr_resp_MC->SetMarkerStyle(24);
    gr_resp_MC->SetMarkerColor(MC_color);
    gr_resp_MC->SetLineColor(MC_color);

    TGraphErrors* gr_resp_MPFDATA = new TGraphErrors(nPoints, x, y_resp_MPFDATA, x_err, y_resp_err_MPFDATA);
    gr_resp_MPFDATA->SetMarkerStyle(20);
    gr_resp_MPFDATA->SetMarkerColor(MC_color);
    gr_resp_MPFDATA->SetLineColor(MC_color);
  
    TGraphErrors* gr_resp_MPF = new TGraphErrors(nPoints, x, y_resp_MPFMC, x_err, y_resp_MPFMC_err);
    gr_resp_MPF->SetMarkerStyle(24);
    gr_resp_MPF->SetMarkerColor(MC_color);
    gr_resp_MPF->SetLineColor(MC_color);
    
    
    //PLI 
    
    TGraphErrors* gr_PLI = new TGraphErrors(nPoints, x, y_PLI, x_err, y_PLI_err);
    gr_PLI->SetMarkerStyle(20);
    gr_PLI->SetMarkerColor(MC_color);
    gr_PLI->SetLineColor(MC_color);

    Float_t lastX = x[nPoints - 1];
    Float_t xMax_fit = lastX + 2. / 100.;
    Float_t xMax_axis;
    xMax_axis = 30. / 100.;
  
    std::string xTitle = "p_{T}^{2^{nd} Jet} / p_{T}^{#gamma}";

    std::string fitFunct_name;
   // fitFunct_name = "[0] + x*[1]";
    
    fitFunct_name = "[0] + x*[1] ";


    // MC Balancing
    TF1* fit_resp = new TF1("fit_resp", fitFunct_name.c_str());
    fit_resp->SetRange(0., xMax_fit);
    fit_resp->SetLineColor(2);
    fit_resp->SetLineColor(MC_color);
    fit_resp->SetLineStyle(2);
    fit_resp->SetLineWidth(1.);
    gr_resp_MC->Fit(fit_resp, "RQ");

    // DATA Balancing
    TF1* fit_resp_DATA = new TF1("fit_resp_DATA", fitFunct_name.c_str());
    fit_resp_DATA->SetRange(0., xMax_fit);
    fit_resp_DATA->SetLineColor(MC_color);
    fit_resp_DATA->SetLineWidth(1.);
    gr_resp_DATA->Fit(fit_resp_DATA, "RQ");

    const std::string lineFunction = "[0] + x*[1] ";
    //lineFunction = "[0] + x*[1] + x*x*[2]";
    // MC MPF
    TF1* fit_respMPF = new TF1("fit_respMPF", lineFunction.c_str());
    fit_respMPF->SetRange(0., xMax_fit);
    fit_respMPF->SetParameter(1, 0);
    fit_respMPF->SetLineColor(2);
    fit_respMPF->SetLineColor(MC_color);
    fit_respMPF->SetLineStyle(2);
    fit_respMPF->SetLineWidth(1.);
    gr_resp_MPF->Fit(fit_respMPF, "RQ");
     
    // DATA MPF
    TF1* fit_resp_MPFDATA = new TF1("fit_resp_MPFDATA", lineFunction.c_str());
    fit_resp_MPFDATA->SetRange(0., xMax_fit);
    //    fit_resp_MPFDATA->SetParameter(0, fit_resp_genMPF->GetParameter(0));
    fit_resp_MPFDATA->SetParameter(1, 0);
    fit_resp_MPFDATA->SetLineColor(MC_color);
    fit_resp_MPFDATA->SetLineWidth(1.);
    gr_resp_MPFDATA->Fit(fit_resp_MPFDATA, "RQ");
    
    //PLI 
    
    TF1* fit_PLI = new TF1("fit_PLI", lineFunction.c_str());
    fit_PLI->SetRange(0., xMax_fit);
    //    fit_resp_MPFDATA->SetParameter(0, fit_resp_genMPF->GetParameter(0));
    fit_PLI->SetParameter(1, 0);
    fit_PLI->SetLineColor(MC_color);
    fit_PLI->SetLineWidth(1.);
    gr_PLI->Fit(fit_PLI, "RQ");
    
   
       std::cout<< " +++++++ Extrapolated Response (alpha = 0) " <<  std::endl;
    // set response graph points:
    if (fit_resp_DATA->GetParameter(0) > 0) {
      std::cout<< "DATA Balancing: " <<  fit_resp_DATA->GetParameter(0)<< std::endl;
      gr_DATAResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp_DATA->GetParameter(0));
      gr_DATAResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp_DATA->GetParError(0));
    }

    if (fit_resp->GetParameter(0) > 0) {
      std::cout<< "MC Balancing: " <<  fit_resp->GetParameter(0)<< std::endl;
      gr_extrapResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp->GetParameter(0));
      gr_extrapResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp->GetParError(0));
    }

    if (fit_resp_MPFDATA->GetParameter(0) > 0) {
      std::cout<< "DATA MPF: " <<  fit_resp_MPFDATA->GetParameter(0)<< std::endl;
      gr_DATARespMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp_MPFDATA->GetParameter(0));
      gr_DATARespMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp_MPFDATA->GetParError(0));
    }

    if (fit_respMPF->GetParameter(0) > 0) {
      std::cout<< "MC MPF: " <<  fit_respMPF->GetParameter(0)<< std::endl;
      gr_extrapRespMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_respMPF->GetParameter(0));
      gr_extrapRespMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_respMPF->GetParError(0));
    }
   //pli graph
    if (fit_PLI->GetParameter(0) > 0) {
      std::cout<< "PLI: " <<  fit_PLI->GetParameter(0)<< std::endl;
      gr_extrapPLI_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_PLI->GetParameter(0));
      gr_extrapPLI_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_PLI->GetParError(0));
    }
    

    Float_t yMin_axis;
    if (this->get_recoType() == "calo") {
      yMin_axis = (currentBin.first < 80.) ? 0.3 : 0.5;
      if (currentBin.first < 20.) yMin_axis = 0.2;
    } else {
      yMin_axis = (currentBin.first < 80.) ? 0.7 : 0.7;
      if (currentBin.first < 20.) yMin_axis = 0.6;
    }

    float yMax_resp = (! rawJets) ? 1.3 : 1.2;

    TH2D* h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "p_{T} Balancing (with raw jets)" : "p_{T} Balancing");
    //h2_axes_resp->GetXaxis()->SetTitleOffset(1.1);
    //h2_axes_resp->GetYaxis()->SetTitleOffset(1.5);
    
    LegendBox legbox;
    legbox.xMin = 0.6;
    legbox.yMin = 0.67;
    legbox.xMax = 0.92;
    legbox.yMax = 0.92;
    
    TLegend* legend_resp;
    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_resp_DATA, "Data", "PL");
    legend_resp->AddEntry(gr_resp_MC, "MC", "PL");

    TString labelText = TString::Format("%d < p_{T}^{#gamma} < %d GeV", (int) ptMin, (int) ptMax);
    TPaveText* label_resp = new TPaveText(0.24, 0.18, 0.4, 0.21, "brNDC");
    label_resp->SetTextFont(42);
    label_resp->SetFillColor(kWhite);
    label_resp->SetTextSize(0.035);
    label_resp->AddText(labelText);

    TPaveText* label_algo = this->get_labelAlgo(4);
    TPaveText* label_cms = this->get_labelCMS(0);
    TPaveText* label_sqrt = this->get_labelSqrt(0);

    //TLatex* label_CMStop = this->get_labelCMStop();

    float markerSize = 1.4;
    gr_resp_MC->SetMarkerSize(markerSize);
    gr_resp_DATA->SetMarkerSize(markerSize);
    gr_resp_MPFDATA->SetMarkerSize(markerSize);
    gr_resp_MPF->SetMarkerSize(markerSize);

    TCanvas* c1_resp = new TCanvas("c1_resp", "c1_resp", 600, 600);
    //c1_resp->SetLeftMargin(0.12);
    //c1_resp->SetBottomMargin(0.12);
    c1_resp->cd();
    h2_axes_resp->Draw();
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");
    gr_resp_MC->Draw("Psame");
    gr_resp_DATA->Draw("Psame");
    
    rawPostfix = (rawJets) ? "RAW" : "";

    // Fit parameters
    TPaveText* label_fit_ = new TPaveText(0.45, 0.60, 0.4, 0.54, "blNDC");
    label_fit_->SetFillColor(kWhite);
    label_fit_->SetTextSize(0.030);

    TString line1_ = TString::Format("MC: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp->GetParameter(0), fit_resp->GetParError(0), fit_resp->GetParameter(1), fit_resp->GetParError(1));
    TString line2_ = TString::Format("Data: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp_DATA->GetParameter(0), fit_resp_DATA->GetParError(0), fit_resp_DATA->GetParameter(1), fit_resp_DATA->GetParError(1));

    label_fit_->SetTextAlign(11);
    label_fit_->AddText(line2_);
    label_fit_->AddText(line1_);
    label_fit_->Draw("same");

    char canvasName_resp_pdf[500];
    char canvasName_resp_png[500];

    sprintf(canvasName_resp_pdf, "%s/response%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/response%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
    if (OUTPUT_GRAPHS) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }

    delete legend_resp;
    /////////////////////////////////////////////

    TGraphErrors* gr_Ratio = fitTools::get_graphRatio(gr_resp_DATA, gr_resp_MC);
    gr_Ratio->SetName("gr_RatioResp");
    gr_Ratio->SetMarkerStyle(20);
    gr_Ratio->SetMarkerColor(MC_color);
    gr_Ratio->SetLineColor(MC_color);
    gr_Ratio->SetMarkerSize(markerSize);
    
    for (int iPointDATA = 0; iPointDATA < gr_Ratio->GetN(); ++iPointDATA) {
      Double_t x, y;
      gr_Ratio->GetPoint(iPointDATA, x, y);
    }

    TF1* fit_Ratio_ = new TF1("fit_Ratio_", fitFunct_name.c_str());
    fit_Ratio_->SetRange(0., xMax_fit);
    fit_Ratio_->SetParameter(0, 1.);
    fit_Ratio_->SetParameter(1, 0.);
    fit_Ratio_->SetLineColor(MC_color);
    fit_Ratio_->SetLineWidth(1.);
    gr_Ratio->Fit(fit_Ratio_, "RQ");
    
    gr_RatioResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_Ratio_->GetParameter(0));
    gr_RatioResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_Ratio_->GetParError(0));
    
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95,1.05 );//yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle("Data/MC");
    h2_axes_resp->Draw();

    TPaveText* label_fit_Ratio = new TPaveText(0.45, 0.65, 0.4, 0.59, "blNDC");
    label_fit_Ratio->SetFillColor(kWhite);
    label_fit_Ratio->SetTextSize(0.030);

    TString line1_Ratio = TString::Format("Ratio: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_Ratio_->GetParameter(0), fit_Ratio_->GetParError(0), fit_Ratio_->GetParameter(1), fit_Ratio_->GetParError(1));

    label_fit_Ratio ->SetTextAlign(11);
    label_fit_Ratio ->AddText(line1_Ratio);
    label_fit_Ratio ->Draw("same");

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_Ratio, "Data/MC", "P");
    
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");
    
    gr_Ratio->Draw("Psame");
    
    if (OUTPUT_GRAPHS) {
      sprintf(canvasName_resp_png, "%s/dataMC_ratio_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
      c1_resp->SaveAs(canvasName_resp_png);
    }
    
    /////////////////////////////////////
    // MPF Extrap
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "MPF Response (with raw MET)" : "MPF Response");
    h2_axes_resp->Draw();

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_resp_MPFDATA, "Data", "PL");
    legend_resp->AddEntry(gr_resp_MPF, "MC", "PL");
    gr_resp_MPF->Draw("Psame");
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");
    gr_resp_MPFDATA->Draw("Psame");

    // Fit parameters
    TPaveText* label_fit = new TPaveText(0.45, 0.40, 0.4, 0.34, "blNDC");
    label_fit->SetFillColor(kWhite);
    label_fit->SetTextSize(0.030);

    TString line1 = TString::Format("MC: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_respMPF->GetParameter(0), fit_respMPF->GetParError(0), fit_respMPF->GetParameter(1), fit_respMPF->GetParError(1));
    TString line2 = TString::Format("Data: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp_MPFDATA->GetParameter(0), fit_resp_MPFDATA->GetParError(0), fit_resp_MPFDATA->GetParameter(1), fit_resp_MPFDATA->GetParError(1));

    label_fit->SetTextAlign(11);
    label_fit->AddText(line2);
    label_fit->AddText(line1);
    label_fit->Draw("same");
  
    std::string MPF = "MPF";

    sprintf(canvasName_resp_pdf, "%s/response%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), MPF.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/response%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), MPF.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);

    if (OUTPUT_GRAPHS || ptMin == 155) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }
  
    TGraphErrors* gr_RatioMPF = fitTools::get_graphRatio(gr_resp_MPFDATA, gr_resp_MPF);
    gr_RatioMPF -> SetName("gr_RatioMPFResp");
    gr_RatioMPF -> SetMarkerStyle(20);
    gr_RatioMPF -> SetMarkerColor(MC_color);
    gr_RatioMPF -> SetLineColor(MC_color);
    gr_RatioMPF -> SetMarkerSize(markerSize);
    
    for (int iPointDATA = 0; iPointDATA < gr_RatioMPF->GetN(); ++iPointDATA) {
      Double_t x, y;
      gr_RatioMPF->GetPoint(iPointDATA, x, y);
    }

    TF1* fit_RatioMPF_ = new TF1("fit_RatioMPF_", fitFunct_name.c_str());
    fit_RatioMPF_->SetRange(0., xMax_fit);
    fit_RatioMPF_->SetParameter(0, 1.);
    fit_RatioMPF_->SetParameter(1, 0.);
    fit_RatioMPF_->SetLineColor(MC_color);
    fit_RatioMPF_->SetLineWidth(1.);
    gr_RatioMPF->Fit(fit_RatioMPF_, "RQ");
    
    gr_RatioMPFResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_RatioMPF_->GetParameter(0));
    gr_RatioMPFResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_RatioMPF_->GetParError(0));
    
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95, 1.05 );//yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle("Data/MC");
    h2_axes_resp->Draw();

    TPaveText* label_fit_RatioMPF = new TPaveText(0.45, 0.65, 0.4, 0.59, "blNDC");
    label_fit_RatioMPF->SetFillColor(kWhite);
    label_fit_RatioMPF->SetTextSize(0.030);

    TString line1_RatioMPF = TString::Format("Fit: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_RatioMPF_->GetParameter(0), fit_RatioMPF_->GetParError(0), fit_RatioMPF_->GetParameter(1), fit_RatioMPF_->GetParError(1));

    label_fit_RatioMPF ->SetTextAlign(11);
    label_fit_RatioMPF ->AddText(line1_RatioMPF);
    label_fit_RatioMPF ->Draw("same");

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_Ratio, "Data/MC", "P");
    
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");
    
    gr_RatioMPF->Draw("Psame");
    
    if (OUTPUT_GRAPHS) {
      sprintf(canvasName_resp_png, "%s/dataMC_ratioMPF_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
      c1_resp->SaveAs(canvasName_resp_png);
    }



  // PLI  

    c1_resp->Clear();
    c1_resp->cd();
    
    
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "PLI (with raw MET)" : "PLI");
    h2_axes_resp->Draw();

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_PLI, "PLI", "PL");
    gr_PLI->Draw("Psame");
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");

    // Fit parameters
    TPaveText* label_fit_pli = new TPaveText(0.45, 0.40, 0.4, 0.34, "blNDC");
    label_fit_pli->SetFillColor(kWhite);
    label_fit_pli->SetTextSize(0.030);

     line1 = TString::Format("PLI: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_PLI->GetParameter(0), fit_PLI->GetParError(0), fit_PLI->GetParameter(1), fit_PLI->GetParError(1));

    label_fit_pli->SetTextAlign(11);
    label_fit_pli->AddText(line1);
    label_fit_pli->Draw("same");
  
    std::string pli = "PLI";

    sprintf(canvasName_resp_pdf, "%s/PLI%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), pli.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/PLI%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), pli.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);

    if (OUTPUT_GRAPHS || ptMin == 155) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }
    
    c1_resp->Clear();
  

  // end PLI
  delete gr_Ratio;
  delete fit_Ratio_;
  delete legend_resp;
  delete c1_resp;
  delete h2_axes_resp;

  ////////////////////////////////////////////////
  /////// RESOLUTION
  //draw resolution histograms:
  // PLI resolution to subscract .
  
  TGraphErrors* gr_reso_PLI = new TGraphErrors(nPoints, x, y_reso_PLI, x_err, y_reso_PLI_err);
  gr_reso_PLI->SetMarkerStyle(20);
  gr_reso_PLI->SetMarkerColor(8);
  gr_reso_PLI->SetLineColor(MC_color);
  // take out points with reso=0:
  for (int iPointDATA = 0; iPointDATA < gr_reso_PLI->GetN(); ++iPointDATA) {
    Double_t x, y;
    gr_reso_PLI->GetPoint(iPointDATA, x, y);
    Double_t yerr = gr_reso_PLI->GetErrorY(iPointDATA);
    if (y < 0.00000001 || yerr == 0.00000000001) gr_reso_PLI->RemovePoint(iPointDATA);
  }
  gr_reso_PLI->SetLineColor(8);
  gr_reso_PLI->SetLineWidth(1.);
  
  TF1* fit_extrapToZero_sqrt_PLI = new TF1("fit_extrapToZero_sqrt_PLI", "[0] + [1]*x");
  fit_extrapToZero_sqrt_PLI->SetRange(0., xMax_fit);
  fit_extrapToZero_sqrt_PLI->SetParameter(0, 1.);
  fit_extrapToZero_sqrt_PLI->SetParameter(1, 0.);
  fit_extrapToZero_sqrt_PLI->SetLineColor(8);

  gr_reso_PLI->Fit(fit_extrapToZero_sqrt_PLI, "RQ");

  if (EXCLUDE_FIRST_POINT_) {
    gr_reso_PLI->RemovePoint(0);
  }
  
  gr_extrapPLIreso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin,  fit_extrapToZero_sqrt_PLI->GetParameter(0));
  gr_extrapPLIreso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt_PLI->GetParError(0));
  
  

  TGraphErrors* gr_reso_DATA = new TGraphErrors(nPoints, x, y_reso_DATA, x_err, y_reso_err_DATA);
  gr_reso_DATA->SetMarkerStyle(20);
  gr_reso_DATA->SetMarkerColor(MC_color);
  gr_reso_DATA->SetLineColor(MC_color);
  // take out points with reso=0:
  for (int iPointDATA = 0; iPointDATA < gr_reso_DATA->GetN(); ++iPointDATA) {
    Double_t x, y;
    gr_reso_DATA->GetPoint(iPointDATA, x, y);
    Double_t yerr = gr_reso_DATA->GetErrorY(iPointDATA);
    if (y < 0.00000001 || yerr == 0.00000000001) gr_reso_DATA->RemovePoint(iPointDATA);
  }
  gr_reso_DATA->SetLineColor(MC_color);
  gr_reso_DATA->SetLineWidth(1.);

  TGraphErrors* gr_reso_MPFDATA = new TGraphErrors(nPoints, x, y_reso_MPFDATA, x_err, y_reso_err_MPFDATA);
  gr_reso_MPFDATA->SetMarkerStyle(20);
  gr_reso_MPFDATA->SetMarkerColor(MC_color);
  gr_reso_MPFDATA->SetLineColor(MC_color);

  TGraphErrors* gr_reso_MC = new TGraphErrors(nPoints, x, y_reso_MC, x_err, y_reso_MC_err);
  gr_reso_MC->SetMarkerStyle(24);
  gr_reso_MC->SetMarkerColor(MC_color);
  gr_reso_MC->SetLineColor(MC_color);
  gr_reso_MC->SetLineStyle(2);
  gr_reso_MC->SetLineWidth(1.);

  TGraphErrors* gr_reso_MPF = new TGraphErrors(nPoints, x, y_reso_MPFMC, x_err, y_reso_MPFMC_err);
  gr_reso_MPF->SetMarkerStyle(24);
  gr_reso_MPF->SetMarkerColor(MC_color);
  gr_reso_MPF->SetLineColor(MC_color);

  //  Double_t x1, x2, y1, y2;
  //  Double_t x1_reco, y1_reco;
  //  gr_reso_MC->GetPoint(1, x1_reco, y1_reco);
  //  Double_t x2_reco, y2_reco;
  //  gr_reso_MC->GetPoint(2, x2_reco, y2_reco);

  //TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
  //TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "[0] + [1]*x + [2]*x*x");
  //fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
 // fit_extrapToZero_sqrt->SetParLimits(0, 0.001, 0.3);
  //  fit_extrapToZero_sqrt->SetParameter(0, c);
  //  fit_extrapToZero_sqrt->FixParameter(1, q);
  //  fit_extrapToZero_sqrt->SetParameter(2, m);
  
  TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "[0] + [1]*x");
  fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
  fit_extrapToZero_sqrt->SetParameter(0, 1.);
  fit_extrapToZero_sqrt->SetParameter(1, 0.);
 // fit_extrapToZero_sqrt->SetLineColor(MC_color);
 // fit_extrapToZero_sqrt->SetLineWidth(1.);
  

  fit_extrapToZero_sqrt->SetLineStyle(2);
  fit_extrapToZero_sqrt->SetLineColor(MC_color);
  fit_extrapToZero_sqrt->SetLineWidth(1.);

  TF1* fit_extrapToZero_sqrt_DATA = new TF1(*fit_extrapToZero_sqrt);
  fit_extrapToZero_sqrt_DATA->SetLineStyle(1);
  fit_extrapToZero_sqrt_DATA->SetLineWidth(1.);

  gr_reso_MC->Fit(fit_extrapToZero_sqrt, "RQ");

  if (EXCLUDE_FIRST_POINT_) {
    gr_reso_DATA->RemovePoint(0);
  }
  gr_reso_DATA->Fit(fit_extrapToZero_sqrt_DATA, "RQ");

  gr_DATAReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)));
  
  gr_DATAReso_vs_pt_PLI_not_sub->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt_DATA->GetParameter(0));
  gr_DATAReso_vs_pt_PLI_not_sub->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt_DATA->GetParError(0));

  
  double errorextrapol_data =(0.5/(std::sqrt(std::pow(fit_extrapToZero_sqrt_DATA->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*std::sqrt(std::pow(2.*fit_extrapToZero_sqrt_DATA->GetParameter(0)*fit_extrapToZero_sqrt_DATA->GetParError(0),2) + std::pow(2.*fit_extrapToZero_sqrt_PLI->GetParameter(0)*fit_extrapToZero_sqrt_PLI->GetParError(0),2));
  
  
  gr_DATAReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, errorextrapol_data /*fit_extrapToZero_sqrt_DATA->GetParError(0)*/);
  
  gr_extrapReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin,   std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)) );
  gr_extrapReso_vs_pt_PLI_not_sub->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt->GetParameter(0));
  gr_extrapReso_vs_pt_PLI_not_sub->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt->GetParError(0));

 double errorextrapol_mc =(0.5/(std::sqrt(std::pow(fit_extrapToZero_sqrt->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*std::sqrt(std::pow(2.*fit_extrapToZero_sqrt->GetParameter(0)*fit_extrapToZero_sqrt->GetParError(0),2) + std::pow(2.*fit_extrapToZero_sqrt_PLI->GetParameter(0)*fit_extrapToZero_sqrt_PLI->GetParError(0),2)); 
  
  gr_extrapReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, /*fit_extrapToZero_sqrt->GetParError(0)*/errorextrapol_mc);
  //extrapReso = sqrt(fit_extrapToZero_sqrt->GetParameter(0));
  //extrapReso_err = ( 1. / (2.*sqrt(extrapReso)) * fit_extrapToZero_sqrt->GetParError(0)); //error propag
  
  
  std::cout<<"reso MC : "<<fit_extrapToZero_sqrt->GetParameter(0)<<" PLI : "<<fit_extrapToZero_sqrt_PLI->GetParameter(0)<<" resolution - PLI  MC : "<<std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)) << " error : "<<errorextrapol_mc<<std::endl;
  std::cout<<"reso data : "<<fit_extrapToZero_sqrt_DATA->GetParameter(0)<<" PLI : "<<fit_extrapToZero_sqrt_PLI->GetParameter(0)<<" resolution - PLI  DATA : "<<std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))<< " error : "<<errorextrapol_data<<std::endl;
  // MPF
  TF1* fit_reso_MPFDATA = new TF1("fit_reso_MPFDATA", "[0] + [1]*x ", 0, xMax_fit);
  fit_reso_MPFDATA->SetLineColor(MC_color);
  fit_reso_MPFDATA->SetLineWidth(1.);
  //  fit_reso_MPFDATA->SetParameter(0, fit_reso_genMPF->GetParameter(0));
  fit_reso_MPFDATA->SetParameter(1, 0);
  gr_reso_MPFDATA->Fit(fit_reso_MPFDATA, "QR");

  TF1* fit_reso_MPF = new TF1("fit_reso_MPF", "[0] + [1]*x ", 0, xMax_fit);
  fit_reso_MPF->SetLineColor(MC_color);
  fit_reso_MPF->SetLineWidth(1.);
  //  fit_reso_MPF->SetParameter(0, fit_reso_genMPF->GetParameter(0));
  fit_reso_MPF->SetParameter(1, 0);
  gr_reso_MPF->Fit(fit_reso_MPF, "QR");

  gr_DATAResoMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_reso_MPFDATA->GetParameter(0));
  gr_DATAResoMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_reso_MPFDATA->GetParError(0));

  gr_extrapResoMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_reso_MPF->GetParameter(0));
  gr_extrapResoMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_reso_MPF->GetParError(0));

  Float_t ymax;
  if (currentBin.first < 40.) ymax = (get_recoType() == "pf") ? 0.7 : 0.9;
  else if (currentBin.first < 80.) ymax = (get_recoType() == "pf") ? 0.6 : 0.8;
  else ymax = (get_recoType() == "pf") ? 0.4 : 0.6;

  TH2D* h2_axes_reso = new TH2D("axes_reso", "", 10, 0., xMax_axis, 10, 0., ymax);
  h2_axes_reso->SetXTitle(xTitle.c_str());
  if (! rawJets)
    h2_axes_reso->SetYTitle("Jet p_{T} Resolution");
  else
    h2_axes_reso->SetYTitle("Raw Jet p_{T} Resolution");

  Float_t minLegend = 0.2;
  TLegend* legend_reso;
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetTextFont(42);
  legend_reso->SetBorderSize(0);
  legend_reso->AddEntry(gr_reso_DATA, "Data", "PL");
  legend_reso->AddEntry(gr_reso_MC, "MC", "PL");
  legend_reso->AddEntry(gr_reso_PLI, "PLI", "PL");

  TPaveText* label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_MC->SetMarkerSize(markerSize);
  gr_reso_DATA->SetMarkerSize(markerSize);

  TCanvas* c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  //c1_reso->SetLeftMargin(0.12);
  //c1_reso->SetBottomMargin(0.12);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  //label_CMStop->Draw("same");
  label_algo->Draw("same");
  gr_reso_MC->Draw("Psame");
  gr_reso_DATA->Draw("Psame");
  gr_reso_PLI->Draw("Psame");

  char canvasName_reso[500];
  sprintf(canvasName_reso, "%s/resolution%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  std::string canvasName_reso_str(canvasName_reso);
  std::string canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  std::string canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS || ptMin == 155) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;

  // MPF Resolution
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetBorderSize(0);
  legend_reso->SetTextFont(42);
  legend_reso->AddEntry(gr_reso_MPFDATA, "Data", "L");
  legend_reso->AddEntry(gr_reso_MPF, "MC", "PL");

  delete label_reso;
  label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_MPF->SetMarkerSize(markerSize);
  gr_reso_MPFDATA->SetMarkerSize(markerSize);

  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_reso_MPF->Draw("Psame");
  gr_reso_MPFDATA->Draw("Psame");

  sprintf(canvasName_reso, "%s/resolutionMPF%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
  canvasName_reso_str = canvasName_reso;
  canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;
  
  
  //--- PLI resolution----
  
  
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetBorderSize(0);
  legend_reso->SetTextFont(42);
  legend_reso->AddEntry(gr_reso_PLI, "PLI", "L");


  delete label_reso;
  label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_PLI->SetMarkerSize(markerSize);


  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_reso_PLI->Draw("Psame");


  sprintf(canvasName_reso, "%s/resolutionPLI%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
  canvasName_reso_str = canvasName_reso;
  canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;
  

  // for Balancing
  TGraphErrors* gr_RatioReso = this->get_graphRatio(gr_reso_DATA, gr_reso_MC);
  gr_RatioReso->SetName("gr_RatioReso");
  gr_RatioReso->SetMarkerStyle(20);
  gr_RatioReso->SetMarkerColor(MC_color);
  gr_RatioReso->SetLineColor(MC_color);
  gr_RatioReso->SetMarkerSize(markerSize);

  // genPhot: y = mx + q
  // recoGen: y = c
  // [0] = c; [1] = q; [2] = m
  /*TF1* fit_RatioReso = new TF1("fit_RatioReso", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
    fit_RatioReso->SetRange(0., xMax_fit);
    fit_RatioReso->SetParameter(0, 1);
  //if( NOQ_ )
  //  fit_RatioReso->FixParameter(1, 0.5*q); //to evaluate syst
  //else
  fit_RatioReso->FixParameter(1, 1); //fixed
  if( FIXM_ ) {
  fit_RatioReso->FixParameter(2, 1);
  } else {
  fit_RatioReso->SetParameter(2, 1);
  //fit_extrapToZero_sqrt->SetParLimits(2, 0., 0.05);
  fit_RatioReso->SetParLimits(2, 0., 0.05*100.);
  }*/
  TF1* fit_RatioReso = new TF1("fit_RatioReso", "[0] + [1]*x");
  fit_RatioReso->SetParameter(0, 1);
  fit_RatioReso->SetParameter(1, 0);
  fit_RatioReso->SetLineStyle(2);
  fit_RatioReso->SetLineColor(MC_color);
  fit_RatioReso->SetLineWidth(1.);
  gr_RatioReso->Fit(fit_RatioReso, "RQ");

  gr_RatioReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_RatioReso->GetParameter(0));
  gr_RatioReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_RatioReso->GetParError(0));

  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->Clear();
  c1_reso->cd();
  h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95, 1.05 );//yMin_axis, yMax_resp);
  h2_axes_resp->SetYTitle("Data/MC");
  h2_axes_resp->Draw();

  if (etaRegion_str != "")
    legend_reso = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
  legend_reso->SetTextSize(0.035);
  legend_reso->SetTextFont(42);
  //legend_resp->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->AddEntry(gr_RatioReso, "Data/MC", "P");
  legend_reso->Draw("same");
  label_resp->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_RatioReso->Draw("Psame");

  if (OUTPUT_GRAPHS || ptMin == 155) {
    char output_name[500];
    sprintf(output_name, "%s/reso_dataMC_ratio_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    c1_reso->SaveAs(output_name);
  }

  delete gr_RatioReso;
  delete fit_RatioReso;
  delete legend_reso;
  delete c1_reso;
  delete h2_axes_reso;
  delete h2_axes_resp;
  delete label_reso;
  delete label_resp;
  delete gr_reso_MC;
  delete gr_reso_PLI;
 


  } //for iPtBin

  gr_DATAResp_vs_pt->Write();
  gr_extrapResp_vs_pt->Write();
  gr_DATARespMPF_vs_pt->Write();
  gr_extrapRespMPF_vs_pt->Write();
  gr_RatioResp_vs_pt->Write();
  gr_RatioMPFResp_vs_pt->Write();

  gr_DATAReso_vs_pt->Write();
  gr_extrapReso_vs_pt->Write();
  gr_DATAResoMPF_vs_pt->Write();
  gr_extrapResoMPF_vs_pt->Write();
  gr_RatioReso_vs_pt->Write();
  gr_extrapPLIreso_vs_pt->Write();
  gr_extrapPLI_vs_pt->Write();
  
  
  graphFile->Close();



} //drawExtrap







void drawExtrap::drawResponseExtrapfine(std::vector<float> ptMeanVec, std::vector<float> alphaMeanVec, const std::string& etaRegion, const std::string& etaRegion_str, bool rawJets) {
  
  int MC_color = MPF_COLOR;
  
  std::vector<std::pair<float, float> > ptPhot_binning = mPtBinning.getBinning();
  //Response
  //DATA
  TGraphErrors* gr_DATAResp_vs_pt = new TGraphErrors(0);
  gr_DATAResp_vs_pt->SetName("gr_DATAResp_vs_pt");
  gr_DATAResp_vs_pt->SetMarkerStyle(25);
  gr_DATAResp_vs_pt->SetMarkerColor(kBlack);
  gr_DATAResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors* gr_DATARespMPF_vs_pt = new TGraphErrors(0);
  gr_DATARespMPF_vs_pt->SetName("gr_DATARespMPF_vs_pt");
  gr_DATARespMPF_vs_pt->SetMarkerStyle(25);
  gr_DATARespMPF_vs_pt->SetMarkerColor(kBlack);
  gr_DATARespMPF_vs_pt->SetMarkerSize(1.5);
  //MC
  TGraphErrors* gr_extrapResp_vs_pt = new TGraphErrors(0);
  gr_extrapResp_vs_pt->SetName("gr_extrapResp_vs_pt");
  gr_extrapResp_vs_pt->SetMarkerStyle(25);
  gr_extrapResp_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors*  gr_extrapRespMPF_vs_pt = new TGraphErrors(0);
  gr_extrapRespMPF_vs_pt->SetName("gr_extrapRespMPF_vs_pt");
  gr_extrapRespMPF_vs_pt->SetMarkerStyle(25);
  gr_extrapRespMPF_vs_pt->SetMarkerColor(kBlack);
  gr_extrapRespMPF_vs_pt->SetMarkerSize(1.5);
  // DATA / MC ratio 
  TGraphErrors* gr_RatioResp_vs_pt = new TGraphErrors(0);
  gr_RatioResp_vs_pt->SetName("gr_RatioResp_vs_pt");
  gr_RatioResp_vs_pt->SetMarkerStyle(25);
  gr_RatioResp_vs_pt->SetMarkerColor(kBlack);
  gr_RatioResp_vs_pt->SetMarkerSize(1.5);
  TGraphErrors* gr_RatioMPFResp_vs_pt = new TGraphErrors(0);
  gr_RatioMPFResp_vs_pt->SetName("gr_RatioMPFResp_vs_pt");
  gr_RatioMPFResp_vs_pt->SetMarkerStyle(25);
  gr_RatioMPFResp_vs_pt->SetMarkerColor(kBlack);
  gr_RatioMPFResp_vs_pt->SetMarkerSize(1.5);


  // Particle Level Imbalance
  TGraphErrors* gr_extrapPLI_vs_pt = new TGraphErrors(0);
  gr_extrapPLI_vs_pt->SetName("gr_extrapPLI_vs_pt");
  gr_extrapPLI_vs_pt->SetMarkerStyle(25);
  gr_extrapPLI_vs_pt->SetMarkerColor(kBlack);
  gr_extrapPLI_vs_pt->SetMarkerSize(1.5);
  
  
  
  TGraphErrors* gr_extrapPLIreso_vs_pt = new TGraphErrors(0);
  gr_extrapPLIreso_vs_pt->SetName("gr_extrapPLIreso_vs_pt");
  gr_extrapPLIreso_vs_pt->SetMarkerStyle(25);
  gr_extrapPLIreso_vs_pt->SetMarkerColor(kBlack);
  gr_extrapPLIreso_vs_pt->SetMarkerSize(1.5);


  // Resolution
  // DATA
  TGraphErrors* gr_DATAReso_vs_pt = new TGraphErrors(0);
  gr_DATAReso_vs_pt->SetName("gr_DATAReso_vs_pt");
  gr_DATAReso_vs_pt->SetMarkerStyle(25);
  gr_DATAReso_vs_pt->SetMarkerColor(kBlack);
  gr_DATAReso_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_DATAkfr_vs_pt = new TGraphErrors(0);
  gr_DATAkfr_vs_pt->SetName("gr_DATAkfr_vs_pt");
  gr_DATAkfr_vs_pt->SetMarkerStyle(25);
  gr_DATAkfr_vs_pt->SetMarkerColor(kBlack);
  gr_DATAkfr_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_DATAReso_vs_pt_PLI_not_sub = new TGraphErrors(0);
  gr_DATAReso_vs_pt_PLI_not_sub->SetName("gr_DATAReso_vs_pt_PLI_not_sub");
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerStyle(25);
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerColor(kBlack);
  gr_DATAReso_vs_pt_PLI_not_sub->SetMarkerSize(1.5);
  
  TGraphErrors* gr_DATAResoMPF_vs_pt = new TGraphErrors(0);
  gr_DATAResoMPF_vs_pt->SetName("gr_DATAResoMPF_vs_pt");
  gr_DATAResoMPF_vs_pt->SetMarkerStyle(25);
  gr_DATAResoMPF_vs_pt->SetMarkerColor(kBlack);
  gr_DATAResoMPF_vs_pt->SetMarkerSize(1.5);
  //MC
  TGraphErrors* gr_extrapkfr_vs_pt = new TGraphErrors(0);
  gr_extrapkfr_vs_pt->SetName("gr_extrapkfr_vs_pt");
  gr_extrapkfr_vs_pt->SetMarkerStyle(25);
  gr_extrapkfr_vs_pt->SetMarkerColor(kBlack);
  gr_extrapkfr_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_extrapReso_vs_pt = new TGraphErrors(0);
  gr_extrapReso_vs_pt->SetName("gr_extrapReso_vs_pt");
  gr_extrapReso_vs_pt->SetMarkerStyle(25);
  gr_extrapReso_vs_pt->SetMarkerColor(kBlack);
  gr_extrapReso_vs_pt->SetMarkerSize(1.5);
  
  TGraphErrors* gr_extrapReso_vs_pt_PLI_not_sub = new TGraphErrors(0);
  gr_extrapReso_vs_pt_PLI_not_sub->SetName("gr_extrapReso_vs_pt_PLI_not_sub");
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerStyle(25);
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerColor(kBlack);
  gr_extrapReso_vs_pt_PLI_not_sub->SetMarkerSize(1.5);
  
  
  TGraphErrors* gr_extrapResoMPF_vs_pt = new TGraphErrors(0);
  gr_extrapResoMPF_vs_pt->SetName("gr_extrapResoMPF_vs_pt");
  gr_extrapResoMPF_vs_pt->SetMarkerStyle(25);
  gr_extrapResoMPF_vs_pt->SetMarkerColor(kBlack);
  gr_extrapResoMPF_vs_pt->SetMarkerSize(1.5);
  //DATA/MC ratio
  TGraphErrors* gr_RatioReso_vs_pt = new TGraphErrors(0);
  gr_RatioReso_vs_pt->SetName("gr_MCRatioReso_vs_pt");
  gr_RatioReso_vs_pt->SetMarkerStyle(25);
  gr_RatioReso_vs_pt->SetMarkerColor(kBlack);
  gr_RatioReso_vs_pt->SetMarkerSize(1.5);

  std::string suffix = get_fullSuffix();
  std::string graphFileName = "PhotonJetExtrapGraphs_fineetabin_" + suffix + "_" + etaRegion + ((rawJets) ? "RAW" : "") + "_" + FIT_RMS_;

  if (NOQ_) graphFileName += "_NOQ";
  if (FIXM_) graphFileName += "_FIXM";
  if (EXCLUDE_FIRST_POINT_) graphFileName += "_NOFIRSTP";

  graphFileName += ".root";

  TFile* graphFile = new TFile(graphFileName.c_str(), "recreate");
  graphFile->cd();

  // To draw also the last PTbins -> ptPhot_binning.size  ---- with -XX not draws last XX bins
  for (uint32_t iPtBin = 0; iPtBin < (ptPhot_binning.size()); //-3 instead of -1 (extrap reaches up to ~2 less bins in pt wrt balancing)
       ++iPtBin) {

    std::pair<float, float> currentBin = mPtBinning.getBinValue(iPtBin);
    float ptMin = currentBin.first;
    float ptMax = currentBin.second;

    std::string rawPostfix = (rawJets) ? "_raw" : "";

    char projName[100];
    sprintf(projName, "projection_%d", iPtBin);
    float ptPhotReco_thisBin =  ptMeanVec.at(iPtBin);
    float ptPhotReco_err_thisBin = 0;

    // npoints in alpha : scan in alpha
    int nPoints = mExtrapBinning.size() ; 
    std::cout<<"binning size : "<<nPoints<<std::endl;
    float x[nPoints];
    float x_err[nPoints];
   // getXPoints(iPtBin, x, x_err);
    // std::cout<<" test vector out of range size vec alpha mean : "<<alphaMeanVec.size() <<std::endl; 
    for (int il = 0 ; il < nPoints ; il++){
    x[il]=alphaMeanVec.at(il);
    x_err[il]=0.;
    
    }
    
   
    /*
    x[0]=0.0;
    x_err[0]=0.;
    
    x[1]=0.1;
    x_err[1]=0.0;
    
    x[2]=0.15;
    x_err[2]=0.0;
    
    x[3]=0.2;
    x_err[3]=0.0;
    
    x[4]=0.25;
    x_err[4]=0.0;
    
    x[5]=0.3;
    x_err[5]=0.0;*/
    
    
   /* x[0]=0.1;
    x_err[0]=0.0;
    
    x[1]=0.15;
    x_err[1]=0.0;
    
    x[2]=0.2;
    x_err[2]=0.0;
    
    x[3]=0.25;
    x_err[3]=0.0;
    
    x[4]=0.3;
    x_err[4]=0.0;*/
    Float_t y_resp_DATA[nPoints];
    Float_t y_resp_err_DATA[nPoints];

    Float_t y_resp_MC[nPoints];
    Float_t y_resp_MC_err[nPoints];

    Float_t y_resp_MPFDATA[nPoints];
    Float_t y_resp_err_MPFDATA[nPoints];

    Float_t y_resp_MPFMC[nPoints];
    Float_t y_resp_MPFMC_err[nPoints];
    
    //PLI
    Float_t y_PLI[nPoints];
    Float_t y_PLI_err[nPoints];
    
    Float_t y_reso_PLI[nPoints];
    Float_t y_reso_PLI_err[nPoints];

    //Resolution
    Float_t y_reso_DATA[nPoints];
    Float_t y_reso_err_DATA[nPoints];

    Float_t y_reso_MC[nPoints];
    Float_t y_reso_MC_err[nPoints];

    Float_t y_reso_MPFDATA[nPoints];
    Float_t y_reso_err_MPFDATA[nPoints];

    Float_t y_reso_MPFMC[nPoints];
    Float_t y_reso_MPFMC_err[nPoints];
    Float_t y_Nevents[nPoints];
    Float_t y_Nevents_MC[nPoints];
    Float_t y_NeventsMPF[nPoints];
    Float_t y_Nevents_MCMPF[nPoints];
    
     Float_t y_NeventsPLI[nPoints];


    TString yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_resp_balancing_fine_bining%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    
    // infinite loop for amc at NLO start
    getYPoints(get_dataFile(0), yHistoName, nPoints, y_resp_DATA, y_resp_err_DATA,  y_reso_DATA, y_reso_err_DATA,y_Nevents);
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_resp_MC, y_resp_MC_err,  y_reso_MC, y_reso_MC_err,y_Nevents_MC);
    // infinite loop for amc at NLO stop
    yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_resp_mpf_fine_bining%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    getYPoints(get_dataFile(0), yHistoName, nPoints, y_resp_MPFDATA, y_resp_err_MPFDATA,  y_reso_MPFDATA, y_reso_err_MPFDATA,y_NeventsMPF);
    
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_resp_MPFMC, y_resp_MPFMC_err,  y_reso_MPFMC, y_reso_MPFMC_err,y_Nevents_MCMPF);

     //PLI
     yHistoName = TString::Format("analysis/extrapolation/extrap_ptPhot_%d_%d/extrap_PLI_fine%s_%s", (int) currentBin.first, (int) currentBin.second, rawPostfix.c_str(), etaRegion.c_str());
    getYPoints(get_mcFile(0), yHistoName, nPoints, y_PLI, y_PLI_err,  y_reso_PLI, y_reso_PLI_err,y_NeventsPLI);
    
    int is_empty = 0;
    for(int i = 1; i < nPoints ; i++){
          if(y_resp_DATA[i]== -1){
             is_empty++;
          }      
    }
    
    
    	 
    	 if(iPtBin < 3){
    	 x[0]= (alphaMeanVec.at(0)*y_Nevents[0]+alphaMeanVec.at(1)*y_Nevents[1])/(y_Nevents[0]+y_Nevents[1]); // (mExtrapBinning.getBinValue(is_empty).second) / 2. ;
    	 y_resp_DATA[0] = (y_resp_DATA[0]*y_Nevents[0] + y_resp_DATA[1]*y_Nevents[1] )/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_resp_MC[0] = (y_resp_MC[0]*y_Nevents_MC[0] + y_resp_MC[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_reso_DATA[0] = (y_reso_DATA[0]*y_Nevents[0] + y_reso_DATA[1]*y_Nevents[1] )/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_reso_MC[0] = (y_reso_MC[0]*y_Nevents_MC[0] + y_reso_MC[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_resp_err_DATA[0] = (y_resp_err_DATA[0]*y_Nevents[0] + y_resp_err_DATA[1]*y_Nevents[1] )/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_resp_MC_err[0] = (y_resp_MC_err[0]*y_Nevents_MC[0] + y_resp_MC_err[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_reso_err_DATA[0] = (y_reso_err_DATA[0]*y_Nevents[0] + y_reso_err_DATA[1]*y_Nevents[1] )/(y_Nevents[0]+y_Nevents[1]) ;
    	 y_reso_MC_err[0] = (y_reso_MC_err[0]*y_Nevents_MC[0] + y_reso_MC_err[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 y_PLI[0] = (y_PLI[0]*y_Nevents_MC[0] + y_PLI[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;   	 
    	 y_PLI_err[0] = (y_PLI_err[0]*y_Nevents[0] + y_PLI_err[1]*y_Nevents[1])/(y_Nevents[0]+y_Nevents[1]) ;
    	 
    	 y_reso_PLI[0] = (y_reso_PLI[0]*y_Nevents_MC[0] + y_reso_PLI[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;   	 
    	 y_reso_PLI_err[0] = (y_reso_PLI_err[0]*y_Nevents_MC[0] + y_reso_PLI_err[1]*y_Nevents_MC[1] )/(y_Nevents_MC[0]+y_Nevents_MC[1]) ;
    	 
    	 
    	 y_resp_DATA[1] = -1 ;
    	 y_resp_MC[1] = -1 ;
    	 
    	 y_reso_DATA[1] = -1 ;
    	 y_reso_MC[1] = -1 ;
    	 
    	 y_resp_err_DATA[1] = 100000000 ;
    	 y_resp_MC_err[1] = 1000000000 ;
    	 
    	 y_reso_err_DATA[1] = 100000000 ;
    	 y_reso_MC_err[1] = 10000000 ;
    	 
    	 y_PLI[1] = -1;   	 
    	 y_PLI_err[1] = -1;
    	 
    	 y_reso_PLI[1] = 10000000 ;   	 
    	 y_reso_PLI_err[1] = 1000000;
    	 
    	 
    	 
    }
    	 
    
    
    
    
    //draw response histograms:
    TGraphErrors* gr_resp_DATA = new TGraphErrors(nPoints, x, y_resp_DATA, x_err, y_resp_err_DATA);
    gr_resp_DATA->SetMarkerStyle(20);
    gr_resp_DATA->SetMarkerColor(MC_color);
    gr_resp_DATA->SetLineColor(MC_color);

    TGraphErrors* gr_resp_MC = new TGraphErrors(nPoints, x, y_resp_MC, x_err, y_resp_MC_err);
    gr_resp_MC->SetMarkerStyle(24);
    gr_resp_MC->SetMarkerColor(MC_color);
    gr_resp_MC->SetLineColor(MC_color);

    TGraphErrors* gr_resp_MPFDATA = new TGraphErrors(nPoints, x, y_resp_MPFDATA, x_err, y_resp_err_MPFDATA);
    gr_resp_MPFDATA->SetMarkerStyle(20);
    gr_resp_MPFDATA->SetMarkerColor(MC_color);
    gr_resp_MPFDATA->SetLineColor(MC_color);
  
    TGraphErrors* gr_resp_MPF = new TGraphErrors(nPoints, x, y_resp_MPFMC, x_err, y_resp_MPFMC_err);
    gr_resp_MPF->SetMarkerStyle(24);
    gr_resp_MPF->SetMarkerColor(MC_color);
    gr_resp_MPF->SetLineColor(MC_color);
    
    
    //PLI 
    
    TGraphErrors* gr_PLI = new TGraphErrors(nPoints, x, y_PLI, x_err, y_PLI_err);
    gr_PLI->SetMarkerStyle(20);
    gr_PLI->SetMarkerColor(MC_color);
    gr_PLI->SetLineColor(MC_color);

    Float_t lastX = x[nPoints - 1];
    Float_t xMax_fit = lastX + 2. / 100.;
    Float_t xMax_axis;
    xMax_axis = 30. / 100.;
  
    std::string xTitle = "p_{T}^{2^{nd} Jet} / p_{T}^{#gamma}";

    std::string fitFunct_name;
    fitFunct_name = "sqrt(pow([0],2)  + x*x*pow([1],2))";


    // MC Balancing
    TF1* fit_resp = new TF1("fit_resp", fitFunct_name.c_str());
    fit_resp->SetRange(0., xMax_fit);
    fit_resp->SetLineColor(2);
    fit_resp->SetLineColor(MC_color);
    fit_resp->SetLineStyle(2);
    fit_resp->SetLineWidth(1.);
    fit_resp->SetParameter(0,1.);
    gr_resp_MC->Fit(fit_resp, "RQ");

    // DATA Balancing
    TF1* fit_resp_DATA = new TF1("fit_resp_DATA", fitFunct_name.c_str());
    fit_resp_DATA->SetRange(0., xMax_fit);
    fit_resp_DATA->SetLineColor(MC_color);
    fit_resp_DATA->SetLineWidth(1.);
    fit_resp_DATA->SetParameter(0,1.);
    gr_resp_DATA->Fit(fit_resp_DATA, "RQ");

    const std::string lineFunction = "[0] + x*[1] ";//"sqrt(pow([0],2)  + x*x*pow([1],2))";
    // MC MPF
    TF1* fit_respMPF = new TF1("fit_respMPF", lineFunction.c_str());
    fit_respMPF->SetRange(0., xMax_fit);
    fit_respMPF->SetParameter(1, 0);
    fit_respMPF->SetLineColor(2);
    fit_respMPF->SetLineColor(MC_color);
    fit_respMPF->SetLineStyle(2);
    fit_respMPF->SetLineWidth(1.);
    fit_respMPF->SetParameter(0,1.);
    gr_resp_MPF->Fit(fit_respMPF, "RQ");
     
    // DATA MPF
    TF1* fit_resp_MPFDATA = new TF1("fit_resp_MPFDATA", lineFunction.c_str());
    fit_resp_MPFDATA->SetRange(0., xMax_fit);
    //    fit_resp_MPFDATA->SetParameter(0, fit_resp_genMPF->GetParameter(0));
    fit_resp_MPFDATA->SetParameter(1, 0);
    fit_resp_MPFDATA->SetLineColor(MC_color);
    fit_resp_MPFDATA->SetLineWidth(1.);
    fit_resp_MPFDATA->SetParameter(0,1.);
    gr_resp_MPFDATA->Fit(fit_resp_MPFDATA, "RQ");
    
    //PLI 
    const std::string lineFunction_PLI = "sqrt(pow([0],2)  + x*x*pow([1],2))";
    
    TF1* fit_PLI = new TF1("fit_PLI", lineFunction_PLI.c_str());
    fit_PLI->SetRange(0., xMax_fit);
    //    fit_resp_MPFDATA->SetParameter(0, fit_resp_genMPF->GetParameter(0));
    fit_PLI->SetParameter(1, 0);
    fit_PLI->SetLineColor(MC_color);
    fit_PLI->SetLineWidth(1.);
    fit_PLI->SetParameter(0,1.);
    gr_PLI->Fit(fit_PLI, "RQ");
    
   
       std::cout<< " +++++++ Extrapolated Response (alpha = 0) " <<  std::endl;
    // set response graph points:
    if (fit_resp_DATA->GetParameter(0) > 0) {
      std::cout<< "DATA Balancing: " <<  fit_resp_DATA->GetParameter(0)<< std::endl;
      gr_DATAResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp_DATA->GetParameter(0));
      gr_DATAResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp_DATA->GetParError(0));
    }

    if (fit_resp->GetParameter(0) > 0) {
      std::cout<< "MC Balancing: " <<  fit_resp->GetParameter(0)<< std::endl;
      gr_extrapResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp->GetParameter(0));
      gr_extrapResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp->GetParError(0));
    }

    if (fit_resp_MPFDATA->GetParameter(0) > 0) {
      std::cout<< "DATA MPF: " <<  fit_resp_MPFDATA->GetParameter(0)<< std::endl;
      gr_DATARespMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_resp_MPFDATA->GetParameter(0));
      gr_DATARespMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_resp_MPFDATA->GetParError(0));
    }

    if (fit_respMPF->GetParameter(0) > 0) {
      std::cout<< "MC MPF: " <<  fit_respMPF->GetParameter(0)<< std::endl;
      gr_extrapRespMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_respMPF->GetParameter(0));
      gr_extrapRespMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_respMPF->GetParError(0));
    }
   //pli graph
    if (fit_PLI->GetParameter(0) > 0) {
      std::cout<< "PLI: " <<  fit_PLI->GetParameter(0)<< std::endl;
      if(fit_PLI->GetParameter(0) > 0){
      gr_extrapPLI_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_PLI->GetParameter(0));
      gr_extrapPLI_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_PLI->GetParError(0));
      }else{
      gr_extrapPLI_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, 0.);
      gr_extrapPLI_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, 0.);
      
      }
      
    }
    

    Float_t yMin_axis;
    if (this->get_recoType() == "calo") {
      yMin_axis = (currentBin.first < 80.) ? 0.3 : 0.5;
      if (currentBin.first < 20.) yMin_axis = 0.2;
    } else {
      yMin_axis = (currentBin.first < 80.) ? 0.7 : 0.7;
      if (currentBin.first < 20.) yMin_axis = 0.6;
    }

    float yMax_resp = (! rawJets) ? 1.3 : 1.2;

    TH2D* h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "p_{T} Balancing (with raw jets)" : "p_{T} Balancing");
    //h2_axes_resp->GetXaxis()->SetTitleOffset(1.1);
    //h2_axes_resp->GetYaxis()->SetTitleOffset(1.5);
    
    LegendBox legbox;
    legbox.xMin = 0.6;
    legbox.yMin = 0.67;
    legbox.xMax = 0.92;
    legbox.yMax = 0.92;
    
    TLegend* legend_resp;
    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_resp_DATA, "Data", "PL");
    legend_resp->AddEntry(gr_resp_MC, "MC", "PL");

    TString labelText = TString::Format("%d < p_{T}^{#gamma} < %d GeV", (int) ptMin, (int) ptMax);
    TPaveText* label_resp = new TPaveText(0.24, 0.18, 0.4, 0.21, "brNDC");
    label_resp->SetTextFont(42);
    label_resp->SetFillColor(kWhite);
    label_resp->SetTextSize(0.035);
    label_resp->AddText(labelText);

    TPaveText* label_algo = this->get_labelAlgo(4);
    TPaveText* label_cms = this->get_labelCMS(0);
    TPaveText* label_sqrt = this->get_labelSqrt(0);

    //TLatex* label_CMStop = this->get_labelCMStop();

    float markerSize = 1.4;
    gr_resp_MC->SetMarkerSize(markerSize);
    gr_resp_DATA->SetMarkerSize(markerSize);
    gr_resp_MPFDATA->SetMarkerSize(markerSize);
    gr_resp_MPF->SetMarkerSize(markerSize);

    TCanvas* c1_resp = new TCanvas("c1_resp", "c1_resp", 600, 600);
    //c1_resp->SetLeftMargin(0.12);
    //c1_resp->SetBottomMargin(0.12);
    c1_resp->cd();
    h2_axes_resp->Draw();
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");
    gr_resp_MC->Draw("Psame");
    gr_resp_DATA->Draw("Psame");
    
    rawPostfix = (rawJets) ? "RAW" : "";

    // Fit parameters
    TPaveText* label_fit_ = new TPaveText(0.45, 0.60, 0.4, 0.54, "blNDC");
    label_fit_->SetFillColor(kWhite);
    label_fit_->SetTextSize(0.030);

    TString line1_ = TString::Format("MC: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp->GetParameter(0), fit_resp->GetParError(0), fit_resp->GetParameter(1), fit_resp->GetParError(1));
    TString line2_ = TString::Format("Data: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp_DATA->GetParameter(0), fit_resp_DATA->GetParError(0), fit_resp_DATA->GetParameter(1), fit_resp_DATA->GetParError(1));

    label_fit_->SetTextAlign(11);
    label_fit_->AddText(line2_);
    label_fit_->AddText(line1_);
    label_fit_->Draw("same");

    char canvasName_resp_pdf[500];
    char canvasName_resp_png[500];

    sprintf(canvasName_resp_pdf, "%s/response%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/response%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
    if (OUTPUT_GRAPHS) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }

    delete legend_resp;
    /////////////////////////////////////////////

    TGraphErrors* gr_Ratio = fitTools::get_graphRatio(gr_resp_DATA, gr_resp_MC);
    gr_Ratio->SetName("gr_RatioResp");
    gr_Ratio->SetMarkerStyle(20);
    gr_Ratio->SetMarkerColor(MC_color);
    gr_Ratio->SetLineColor(MC_color);
    gr_Ratio->SetMarkerSize(markerSize);
    
    for (int iPointDATA = 0; iPointDATA < gr_Ratio->GetN(); ++iPointDATA) {
      Double_t x, y;
      gr_Ratio->GetPoint(iPointDATA, x, y);
    }

    TF1* fit_Ratio_ = new TF1("fit_Ratio_", fitFunct_name.c_str());
    fit_Ratio_->SetRange(0., xMax_fit);
    fit_Ratio_->SetParameter(0, 1.);
    fit_Ratio_->SetParameter(1, 0.);
    fit_Ratio_->SetLineColor(MC_color);
    fit_Ratio_->SetLineWidth(1.);
    gr_Ratio->Fit(fit_Ratio_, "RQ");
    
    gr_RatioResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_Ratio_->GetParameter(0));
    gr_RatioResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_Ratio_->GetParError(0));
    
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95,1.05 );//yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle("Data/MC");
    h2_axes_resp->Draw();

    TPaveText* label_fit_Ratio = new TPaveText(0.45, 0.65, 0.4, 0.59, "blNDC");
    label_fit_Ratio->SetFillColor(kWhite);
    label_fit_Ratio->SetTextSize(0.030);

    TString line1_Ratio = TString::Format("Ratio: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_Ratio_->GetParameter(0), fit_Ratio_->GetParError(0), fit_Ratio_->GetParameter(1), fit_Ratio_->GetParError(1));

    label_fit_Ratio ->SetTextAlign(11);
    label_fit_Ratio ->AddText(line1_Ratio);
    label_fit_Ratio ->Draw("same");

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_Ratio, "Data/MC", "P");
    
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");
    
    gr_Ratio->Draw("Psame");
    
    if (OUTPUT_GRAPHS) {
      sprintf(canvasName_resp_png, "%s/dataMC_ratio_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
      c1_resp->SaveAs(canvasName_resp_png);
    }
    
    /////////////////////////////////////
    // MPF Extrap
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "MPF Response (with raw MET)" : "MPF Response");
    h2_axes_resp->Draw();

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_resp_MPFDATA, "Data", "PL");
    legend_resp->AddEntry(gr_resp_MPF, "MC", "PL");
    gr_resp_MPF->Draw("Psame");
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");
    gr_resp_MPFDATA->Draw("Psame");

    // Fit parameters
    TPaveText* label_fit = new TPaveText(0.45, 0.40, 0.4, 0.34, "blNDC");
    label_fit->SetFillColor(kWhite);
    label_fit->SetTextSize(0.030);

    TString line1 = TString::Format("MC: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_respMPF->GetParameter(0), fit_respMPF->GetParError(0), fit_respMPF->GetParameter(1), fit_respMPF->GetParError(1));
    TString line2 = TString::Format("Data: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_resp_MPFDATA->GetParameter(0), fit_resp_MPFDATA->GetParError(0), fit_resp_MPFDATA->GetParameter(1), fit_resp_MPFDATA->GetParError(1));

    label_fit->SetTextAlign(11);
    label_fit->AddText(line2);
    label_fit->AddText(line1);
    label_fit->Draw("same");
  
    std::string MPF = "MPF";

    sprintf(canvasName_resp_pdf, "%s/response%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), MPF.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/response%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), MPF.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);

    if (OUTPUT_GRAPHS || ptMin == 155) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }
  
    TGraphErrors* gr_RatioMPF = fitTools::get_graphRatio(gr_resp_MPFDATA, gr_resp_MPF);
    gr_RatioMPF -> SetName("gr_RatioMPFResp");
    gr_RatioMPF -> SetMarkerStyle(20);
    gr_RatioMPF -> SetMarkerColor(MC_color);
    gr_RatioMPF -> SetLineColor(MC_color);
    gr_RatioMPF -> SetMarkerSize(markerSize);
    
    for (int iPointDATA = 0; iPointDATA < gr_RatioMPF->GetN(); ++iPointDATA) {
      Double_t x, y;
      gr_RatioMPF->GetPoint(iPointDATA, x, y);
    }

    TF1* fit_RatioMPF_ = new TF1("fit_RatioMPF_", fitFunct_name.c_str());
    fit_RatioMPF_->SetRange(0., xMax_fit);
    fit_RatioMPF_->SetParameter(0, 1.);
    fit_RatioMPF_->SetParameter(1, 0.);
    fit_RatioMPF_->SetLineColor(MC_color);
    fit_RatioMPF_->SetLineWidth(1.);
    gr_RatioMPF->Fit(fit_RatioMPF_, "RQ");
    
    gr_RatioMPFResp_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_RatioMPF_->GetParameter(0));
    gr_RatioMPFResp_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_RatioMPF_->GetParError(0));
    
    c1_resp->Clear();
    c1_resp->cd();
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95, 1.05 );//yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle("Data/MC");
    h2_axes_resp->Draw();

    TPaveText* label_fit_RatioMPF = new TPaveText(0.45, 0.65, 0.4, 0.59, "blNDC");
    label_fit_RatioMPF->SetFillColor(kWhite);
    label_fit_RatioMPF->SetTextSize(0.030);

    TString line1_RatioMPF = TString::Format("Fit: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_RatioMPF_->GetParameter(0), fit_RatioMPF_->GetParError(0), fit_RatioMPF_->GetParameter(1), fit_RatioMPF_->GetParError(1));

    label_fit_RatioMPF ->SetTextAlign(11);
    label_fit_RatioMPF ->AddText(line1_RatioMPF);
    label_fit_RatioMPF ->Draw("same");

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_Ratio, "Data/MC", "P");
    
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_algo->Draw("same");
    
    gr_RatioMPF->Draw("Psame");
    
    if (OUTPUT_GRAPHS) {
      sprintf(canvasName_resp_png, "%s/dataMC_ratioMPF_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
      c1_resp->SaveAs(canvasName_resp_png);
    }



  // PLI  

    c1_resp->Clear();
    c1_resp->cd();
    
    
    h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, yMin_axis, yMax_resp);
    h2_axes_resp->SetXTitle(xTitle.c_str());
    h2_axes_resp->SetYTitle((rawJets) ? "PLI (with raw MET)" : "PLI");
    h2_axes_resp->Draw();

    if (etaRegion_str != "")
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
    else
      legend_resp = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
    legend_resp->SetTextSize(0.035);
    legend_resp->SetTextFont(42);
    legend_resp->SetBorderSize(0);
    //legend_resp->SetFillStyle(0);
    legend_resp->SetFillColor(kWhite);
    legend_resp->AddEntry(gr_PLI, "PLI", "PL");
    gr_PLI->Draw("Psame");
    legend_resp->Draw("same");
    label_resp->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    //label_CMStop->Draw("same");
    label_algo->Draw("same");

    // Fit parameters
    TPaveText* label_fit_pli = new TPaveText(0.45, 0.40, 0.4, 0.34, "blNDC");
    label_fit_pli->SetFillColor(kWhite);
    label_fit_pli->SetTextSize(0.030);

     line1 = TString::Format("PLI: (%.3lf #pm %.3lf)  + (%.3lf #pm %.3lf)x", fit_PLI->GetParameter(0), fit_PLI->GetParError(0), fit_PLI->GetParameter(1), fit_PLI->GetParError(1));

    label_fit_pli->SetTextAlign(11);
    label_fit_pli->AddText(line1);
    label_fit_pli->Draw("same");
  
    std::string pli = "PLI";

    sprintf(canvasName_resp_pdf, "%s/PLI%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), pli.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    sprintf(canvasName_resp_png, "%s/PLI%s%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), pli.c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);

    if (OUTPUT_GRAPHS || ptMin == 155) {
      //      c1_resp->SaveAs(canvasName_resp_pdf);
      c1_resp->SaveAs(canvasName_resp_png);
    }
    
    c1_resp->Clear();
  

  // end PLI
  delete gr_Ratio;
  delete fit_Ratio_;
  delete legend_resp;
  delete c1_resp;
  delete h2_axes_resp;

  ////////////////////////////////////////////////
  /////// RESOLUTION
  //draw resolution histograms:
  // PLI resolution to subscract .
  
  TGraphErrors* gr_reso_PLI = new TGraphErrors(nPoints, x, y_reso_PLI, x_err, y_reso_PLI_err);
  gr_reso_PLI->SetMarkerStyle(20);
  gr_reso_PLI->SetMarkerColor(8);
  gr_reso_PLI->SetLineColor(8);
  gr_reso_PLI->SetMarkerSize(markerSize);
  // take out points with reso=0:
  for (int iPointDATA = 0; iPointDATA < gr_reso_PLI->GetN(); ++iPointDATA) {
    Double_t x, y;
    gr_reso_PLI->GetPoint(iPointDATA, x, y);
    Double_t yerr = gr_reso_PLI->GetErrorY(iPointDATA);
    if (y < 0.00000001 || yerr == 0.00000000001) gr_reso_PLI->RemovePoint(iPointDATA);
  }
  gr_reso_PLI->SetLineColor(8);
  gr_reso_PLI->SetLineWidth(1.);
  
  TF1* fit_extrapToZero_sqrt_PLI = new TF1("fit_extrapToZero_sqrt_PLI", "sqrt(pow([0],2) + x*x*[1])");
  fit_extrapToZero_sqrt_PLI->SetRange(0., xMax_fit);
  fit_extrapToZero_sqrt_PLI->SetParameter(0, 1.);
 // fit_extrapToZero_sqrt_PLI->SetParameter(1, 0.);
  fit_extrapToZero_sqrt_PLI->SetLineColor(8);

  gr_reso_PLI->Fit(fit_extrapToZero_sqrt_PLI, "RQ");

  if (EXCLUDE_FIRST_POINT_) {
    gr_reso_PLI->RemovePoint(0);
  }
  if(fit_extrapToZero_sqrt_PLI->GetParameter(0) > 0){
  gr_extrapPLIreso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin,  fit_extrapToZero_sqrt_PLI->GetParameter(0));
  gr_extrapPLIreso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt_PLI->GetParError(0));
  }else{
  
    gr_extrapPLIreso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin,  0.);
    gr_extrapPLIreso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, 0.);
  
  }
 

  TGraphErrors* gr_reso_DATA = new TGraphErrors(nPoints, x, y_reso_DATA, x_err, y_reso_err_DATA);
  gr_reso_DATA->SetMarkerStyle(20);
  gr_reso_DATA->SetMarkerColor(MC_color);
  gr_reso_DATA->SetLineColor(MC_color);
  // take out points with reso=0:
  for (int iPointDATA = 0; iPointDATA < gr_reso_DATA->GetN(); ++iPointDATA) {
    Double_t x, y;
    gr_reso_DATA->GetPoint(iPointDATA, x, y);
    Double_t yerr = gr_reso_DATA->GetErrorY(iPointDATA);
    if (y < 0.00000001 || yerr == 0.00000000001) gr_reso_DATA->RemovePoint(iPointDATA);
  }
  gr_reso_DATA->SetLineColor(MC_color);
  gr_reso_DATA->SetLineWidth(1.);

  TGraphErrors* gr_reso_MPFDATA = new TGraphErrors(nPoints, x, y_reso_MPFDATA, x_err, y_reso_err_MPFDATA);
  gr_reso_MPFDATA->SetMarkerStyle(20);
  gr_reso_MPFDATA->SetMarkerColor(MC_color);
  gr_reso_MPFDATA->SetLineColor(MC_color);

  TGraphErrors* gr_reso_MC = new TGraphErrors(nPoints, x, y_reso_MC, x_err, y_reso_MC_err);
  gr_reso_MC->SetMarkerStyle(24);
  gr_reso_MC->SetMarkerColor(MC_color);
  gr_reso_MC->SetLineColor(MC_color);
  gr_reso_MC->SetLineStyle(2);
  gr_reso_MC->SetLineWidth(1.);

  TGraphErrors* gr_reso_MPF = new TGraphErrors(nPoints, x, y_reso_MPFMC, x_err, y_reso_MPFMC_err);
  gr_reso_MPF->SetMarkerStyle(24);
  gr_reso_MPF->SetMarkerColor(MC_color);
  gr_reso_MPF->SetLineColor(MC_color);

  //  Double_t x1, x2, y1, y2;
  //  Double_t x1_reco, y1_reco;
  //  gr_reso_MC->GetPoint(1, x1_reco, y1_reco);
  //  Double_t x2_reco, y2_reco;
  //  gr_reso_MC->GetPoint(2, x2_reco, y2_reco);

  //TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
  //TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt", "[0] + [1]*x + [2]*x*x");
  //fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
 // fit_extrapToZero_sqrt->SetParLimits(0, 0.001, 0.3);
  //  fit_extrapToZero_sqrt->SetParameter(0, c);
  //  fit_extrapToZero_sqrt->FixParameter(1, q);
  //  fit_extrapToZero_sqrt->SetParameter(2, m);
  
  TF1* fit_extrapToZero_sqrt = new TF1("fit_extrapToZero_sqrt",  "sqrt(pow([0],2) + x*x*[1])");
  fit_extrapToZero_sqrt->SetRange(0., xMax_fit);
  fit_extrapToZero_sqrt->SetParameter(0, 1.);
 // fit_extrapToZero_sqrt->SetParameter(1, 0.);
 // fit_extrapToZero_sqrt->SetLineColor(MC_color);
 // fit_extrapToZero_sqrt->SetLineWidth(1.);
  

  fit_extrapToZero_sqrt->SetLineStyle(2);
  fit_extrapToZero_sqrt->SetLineColor(MC_color);
  fit_extrapToZero_sqrt->SetLineWidth(1.);

  TF1* fit_extrapToZero_sqrt_DATA = new TF1(*fit_extrapToZero_sqrt);
  fit_extrapToZero_sqrt_DATA->SetLineStyle(1);
  fit_extrapToZero_sqrt_DATA->SetLineWidth(1.);

  gr_reso_MC->Fit(fit_extrapToZero_sqrt, "RQ");

  if (EXCLUDE_FIRST_POINT_) {
    gr_reso_DATA->RemovePoint(0);
  }
  gr_reso_DATA->Fit(fit_extrapToZero_sqrt_DATA, "RQ");

  gr_DATAReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)));
/*  
  double errorextrapol_data =(0.5/(std::sqrt(std::pow(fit_extrapToZero_sqrt_DATA->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))* std::sqrt(std::pow( fabs((2*fit_extrapToZero_sqrt_DATA->GetParameter(0) - std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2 ))*(1/std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*fit_extrapToZero_sqrt_DATA->GetParError(0),2) + std::pow(fabs((std::pow(fit_extrapToZero_sqrt_DATA->GetParameter(0),2 ) - 2*fit_extrapToZero_sqrt_PLI->GetParameter(0))*(1/std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*fit_extrapToZero_sqrt_PLI->GetParError(0),2));
  */
  
  double errorextrapol_data =(0.5/(std::sqrt(std::pow(fit_extrapToZero_sqrt_DATA->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*std::sqrt(std::pow(2.*fit_extrapToZero_sqrt_DATA->GetParameter(0)*fit_extrapToZero_sqrt_DATA->GetParError(0),2) + std::pow(2.*fit_extrapToZero_sqrt_PLI->GetParameter(0)*fit_extrapToZero_sqrt_PLI->GetParError(0),2));
  
  
  gr_DATAReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, errorextrapol_data /*fit_extrapToZero_sqrt_DATA->GetParError(0)*/);
  
  gr_extrapReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin,   std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)) );
  
  
  gr_DATAReso_vs_pt_PLI_not_sub->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt_DATA->GetParameter(1));
  gr_DATAReso_vs_pt_PLI_not_sub->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt_DATA->GetParError(0));
 /* 
  double errorextrapol_mc =(1./(std::sqrt(std::pow(fit_extrapToZero_sqrt->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2)))) * std::sqrt(std::pow( fabs((2*fit_extrapToZero_sqrt->GetParameter(0) - std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2 ))*(1/std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*fit_extrapToZero_sqrt->GetParError(0),2) + std::pow(fabs((std::pow(fit_extrapToZero_sqrt->GetParameter(0),2 ) - 2*fit_extrapToZero_sqrt_PLI->GetParameter(0))*(1/std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*fit_extrapToZero_sqrt_PLI->GetParError(0),2))  ;
 */
 double errorextrapol_mc =(0.5/(std::sqrt(std::pow(fit_extrapToZero_sqrt->GetParameter(0),2)-std::pow(fit_extrapToZero_sqrt_PLI->GetParameter(0),2))))*std::sqrt(std::pow(2.*fit_extrapToZero_sqrt->GetParameter(0)*fit_extrapToZero_sqrt->GetParError(0),2) + std::pow(2.*fit_extrapToZero_sqrt_PLI->GetParameter(0)*fit_extrapToZero_sqrt_PLI->GetParError(0),2)); 
  
  gr_extrapReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, /*fit_extrapToZero_sqrt->GetParError(0)*/errorextrapol_mc);
  //extrapReso = sqrt(fit_extrapToZero_sqrt->GetParameter(0));
  //extrapReso_err = ( 1. / (2.*sqrt(extrapReso)) * fit_extrapToZero_sqrt->GetParError(0)); //error propag
  
   gr_extrapReso_vs_pt_PLI_not_sub->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt->GetParameter(0));
   gr_extrapReso_vs_pt_PLI_not_sub->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt->GetParError(0));
  
  std::cout<<"reso MC : "<<fit_extrapToZero_sqrt->GetParameter(0)<<" PLI : "<<fit_extrapToZero_sqrt_PLI->GetParameter(0)<<" resolution - PLI  MC : "<<std::sqrt(std::pow( fit_extrapToZero_sqrt->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2)) << " error : "<<errorextrapol_mc<<std::endl;
  std::cout<<"reso data : "<<fit_extrapToZero_sqrt_DATA->GetParameter(0)<<" PLI : "<<fit_extrapToZero_sqrt_PLI->GetParameter(0)<<" resolution - PLI  DATA : "<<std::sqrt(std::pow( fit_extrapToZero_sqrt_DATA->GetParameter(0),2) - std::pow(  fit_extrapToZero_sqrt_PLI->GetParameter(0),2))<< " error : "<<errorextrapol_data<<std::endl;
  //kfr vs pt 
  gr_DATAkfr_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt_DATA->GetParameter(1));
  gr_DATAkfr_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt_DATA->GetParError(1));
  
  gr_extrapkfr_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_extrapToZero_sqrt->GetParameter(1));
  gr_extrapkfr_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_extrapToZero_sqrt->GetParError(1));
  // MPF
  TF1* fit_reso_MPFDATA = new TF1("fit_reso_MPFDATA",  "sqrt(pow([0],2) + x*x*[1])", 0, xMax_fit);
  fit_reso_MPFDATA->SetLineColor(MC_color);
  fit_reso_MPFDATA->SetLineWidth(1.);
  //  fit_reso_MPFDATA->SetParameter(0, fit_reso_genMPF->GetParameter(0));
  fit_reso_MPFDATA->SetParameter(0, 1);
  gr_reso_MPFDATA->Fit(fit_reso_MPFDATA, "QR");

  TF1* fit_reso_MPF = new TF1("fit_reso_MPF", "sqrt(pow([0],2) + x*x*[1])", 0, xMax_fit);
  fit_reso_MPF->SetLineColor(MC_color);
  fit_reso_MPF->SetLineWidth(1.);
  //  fit_reso_MPF->SetParameter(0, fit_reso_genMPF->GetParameter(0));
  fit_reso_MPF->SetParameter(0, 1);
  gr_reso_MPF->Fit(fit_reso_MPF, "QR");

  gr_DATAResoMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_reso_MPFDATA->GetParameter(0));
  gr_DATAResoMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_reso_MPFDATA->GetParError(0));

  gr_extrapResoMPF_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_reso_MPF->GetParameter(0));
  gr_extrapResoMPF_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_reso_MPF->GetParError(0));

  Float_t ymax;
  if (currentBin.first < 40.) ymax = (get_recoType() == "pf") ? 0.7 : 0.9;
  else if (currentBin.first < 80.) ymax = (get_recoType() == "pf") ? 0.6 : 0.8;
  else ymax = (get_recoType() == "pf") ? 0.4 : 0.6;

  TH2D* h2_axes_reso = new TH2D("axes_reso", "", 10, 0., xMax_axis, 10, 0., ymax);
  h2_axes_reso->SetXTitle(xTitle.c_str());
  if (! rawJets)
    h2_axes_reso->SetYTitle("Jet p_{T} Resolution");
  else
    h2_axes_reso->SetYTitle("Raw Jet p_{T} Resolution");

  Float_t minLegend = 0.2;
  TLegend* legend_reso;
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetTextFont(42);
  legend_reso->SetBorderSize(0);
  legend_reso->AddEntry(gr_reso_DATA, "Data", "PL");
  legend_reso->AddEntry(gr_reso_MC, "MC", "PL");
  legend_reso->AddEntry(gr_reso_PLI, "PLI", "PL");

  TPaveText* label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_MC->SetMarkerSize(markerSize);
  gr_reso_DATA->SetMarkerSize(markerSize);

  TCanvas* c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  //c1_reso->SetLeftMargin(0.12);
  //c1_reso->SetBottomMargin(0.12);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  //label_CMStop->Draw("same");
  label_algo->Draw("same");
  gr_reso_MC->Draw("Psame");
  gr_reso_DATA->Draw("Psame");
  gr_reso_PLI->Draw("Psame");

  char canvasName_reso[500];
  sprintf(canvasName_reso, "%s/resolution%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  std::string canvasName_reso_str(canvasName_reso);
  std::string canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  std::string canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS || ptMin == 155) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;

  // MPF Resolution
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetBorderSize(0);
  legend_reso->SetTextFont(42);
  legend_reso->AddEntry(gr_reso_MPFDATA, "Data", "L");
  legend_reso->AddEntry(gr_reso_MPF, "MC", "PL");

  delete label_reso;
  label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_MPF->SetMarkerSize(markerSize);
  gr_reso_MPFDATA->SetMarkerSize(markerSize);

  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_reso_MPF->Draw("Psame");
  gr_reso_MPFDATA->Draw("Psame");

  sprintf(canvasName_reso, "%s/resolutionMPF%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
  canvasName_reso_str = canvasName_reso;
  canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;
  
  
  //--- PLI resolution----
  
  
  if (etaRegion_str != "")
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(minLegend, 0.6, 0.55, 0.85);
  legend_reso->SetTextSize(0.04);
  //legend_reso->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->SetBorderSize(0);
  legend_reso->SetTextFont(42);
  legend_reso->AddEntry(gr_reso_PLI, "PLI", "L");


  delete label_reso;
  label_reso = new TPaveText(0.25, 0.85, 0.47, 0.9, "brNDC");
  label_reso->SetFillColor(kWhite);
  label_reso->SetTextSize(0.035);
  label_reso->AddText(labelText);
  label_reso->SetTextFont(42);

  gr_reso_PLI->SetMarkerSize(markerSize);


  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->cd();
  h2_axes_reso->Draw();
  legend_reso->Draw("same");
  label_reso->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_reso_PLI->Draw("Psame");


  sprintf(canvasName_reso, "%s/resolutionPLI%s_%s_ptPhot_%d_%d", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
  
  canvasName_reso_str = canvasName_reso;
  canvasName_reso_pdf = canvasName_reso_str + ".pdf";
  canvasName_reso_png = canvasName_reso_str + ".pdf";
  if (OUTPUT_GRAPHS) {
    //    c1_reso->SaveAs(canvasName_reso_pdf.c_str());
    c1_reso->SaveAs(canvasName_reso_png.c_str());
  }

  delete legend_reso;
  delete c1_reso;
  

  // for Balancing
  TGraphErrors* gr_RatioReso = this->get_graphRatio(gr_reso_DATA, gr_reso_MC);
  gr_RatioReso->SetName("gr_RatioReso");
  gr_RatioReso->SetMarkerStyle(20);
  gr_RatioReso->SetMarkerColor(MC_color);
  gr_RatioReso->SetLineColor(MC_color);
  gr_RatioReso->SetMarkerSize(markerSize);

  // genPhot: y = mx + q
  // recoGen: y = c
  // [0] = c; [1] = q; [2] = m
  /*TF1* fit_RatioReso = new TF1("fit_RatioReso", "sqrt([0]*[0] + [1]*[1] + 2.*[1]*[2]*x + [2]*[2]*x*x)");
    fit_RatioReso->SetRange(0., xMax_fit);
    fit_RatioReso->SetParameter(0, 1);
  //if( NOQ_ )
  //  fit_RatioReso->FixParameter(1, 0.5*q); //to evaluate syst
  //else
  fit_RatioReso->FixParameter(1, 1); //fixed
  if( FIXM_ ) {
  fit_RatioReso->FixParameter(2, 1);
  } else {
  fit_RatioReso->SetParameter(2, 1);
  //fit_extrapToZero_sqrt->SetParLimits(2, 0., 0.05);
  fit_RatioReso->SetParLimits(2, 0., 0.05*100.);
  }*/
  TF1* fit_RatioReso = new TF1("fit_RatioReso", "[0] + [1]*x");
  fit_RatioReso->SetParameter(0, 1);
  fit_RatioReso->SetParameter(1, 0);
  fit_RatioReso->SetLineStyle(2);
  fit_RatioReso->SetLineColor(MC_color);
  fit_RatioReso->SetLineWidth(1.);
  gr_RatioReso->Fit(fit_RatioReso, "RQ");

  gr_RatioReso_vs_pt->SetPoint(iPtBin, ptPhotReco_thisBin, fit_RatioReso->GetParameter(0));
  gr_RatioReso_vs_pt->SetPointError(iPtBin, ptPhotReco_err_thisBin, fit_RatioReso->GetParError(0));

  c1_reso = new TCanvas("c1_reso", "c1_reso", 600, 600);
  c1_reso->Clear();
  c1_reso->cd();
  h2_axes_resp = new TH2D("axes_resp", "", 10, 0., xMax_axis, 10, 0.95, 1.05 );//yMin_axis, yMax_resp);
  h2_axes_resp->SetYTitle("Data/MC");
  h2_axes_resp->Draw();

  if (etaRegion_str != "")
    legend_reso = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax, etaRegion_str.c_str());
  else
    legend_reso = new TLegend(legbox.xMin, legbox.yMin, legbox.xMax, legbox.yMax);
  legend_reso->SetTextSize(0.035);
  legend_reso->SetTextFont(42);
  //legend_resp->SetFillStyle(0);
  legend_reso->SetFillColor(kWhite);
  legend_reso->AddEntry(gr_RatioReso, "Data/MC", "P");
  legend_reso->Draw("same");
  label_resp->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  gr_RatioReso->Draw("Psame");

  if (OUTPUT_GRAPHS || ptMin == 155) {
    char output_name[500];
    sprintf(output_name, "%s/reso_dataMC_ratio_%s_%s_ptPhot_%d_%d.pdf", get_outputdir().c_str(), rawPostfix.c_str(), etaRegion.c_str(), (int)ptMin, (int)ptMax);
    c1_reso->SaveAs(output_name);
  }

  delete gr_RatioReso;
  delete fit_RatioReso;
  delete legend_reso;
  delete c1_reso;
  delete h2_axes_reso;
  delete h2_axes_resp;
  delete label_reso;
  delete label_resp;
  delete gr_reso_MC;
  delete gr_reso_PLI;
 


  } //for iPtBin

  gr_DATAResp_vs_pt->Write();
  gr_extrapResp_vs_pt->Write();
  gr_DATARespMPF_vs_pt->Write();
  gr_extrapRespMPF_vs_pt->Write();
  gr_RatioResp_vs_pt->Write();
  gr_RatioMPFResp_vs_pt->Write();

  gr_DATAReso_vs_pt->Write();
  gr_extrapReso_vs_pt->Write();
  gr_DATAResoMPF_vs_pt->Write();
  gr_extrapResoMPF_vs_pt->Write();
  gr_RatioReso_vs_pt->Write();
  gr_extrapPLIreso_vs_pt->Write();
  gr_extrapPLI_vs_pt->Write();
  gr_DATAkfr_vs_pt->Write();
  gr_extrapkfr_vs_pt->Write();
  
  
  graphFile->Close();



}












void drawExtrap::getXPoints(int ptBin, Float_t* x, Float_t* x_err) const {


  /*for (int i = 0; i < nPoints; ++i) {
    char fullName[100];
    sprintf(fullName, "%s_%d", xHistoName,  i);
    TH1D* h1_pt2ndJetMean = (TH1D*)file->Get(fullName);
    x[i] = h1_pt2ndJetMean->GetMean();
    x_err[i] =  h1_pt2ndJetMean->GetRMS() / sqrt((Float_t)h1_pt2ndJetMean->GetEntries());
    }*/
  std::pair<float, float> bin = mExtrapBinning.getBinValue(0);

  float minPt = 0.075;
  float maxPt = 0;
  float stepPt = 0.05;
  size_t numPoints = mExtrapBinning.size();

  //  std::cout<<"ptBin = "<< ptBin<< std::endl;
  //  std::cout<<"n points = "<< numPoints<< std::endl;
  //  std::cout<<"bin.first = "<< bin.first<< std::endl;
  //  std::cout<<"bin.second = "<< bin.second<< std::endl;

  for (size_t i = 0; i < numPoints; i++) {
    
    //    std::cout<<"Bin starts from: "<< minPt<< std::endl;
    //    std::cout<<"by step of: "<< stepPt<< std::endl;
    
    maxPt = minPt + stepPt;
    //    std::cout<<"Bin ends to: "<< maxPt<< std::endl;

    x[i] = (minPt + maxPt) / 2.;
    //    std::cout<<"points x= "<< x[i]<< std::endl;
    x_err[i] = 0.;

    minPt += stepPt;
  }

} //getxpoints

void drawExtrap::getYPointsVector(TFile *file, const std::string& yHistoName, Int_t nPoints, std::vector<Float_t>& y_resp, std::vector<Float_t>& y_resp_err, std::vector<Float_t>& y_reso, std::vector<Float_t>& y_reso_err) const {

  int min = -1, max = 0;
  for (int i = 0; i < nPoints; i++) {
    TString fullName = TString::Format("%s_%d", yHistoName.c_str(), i);
    TH1* h1_r = static_cast<TH1*>(file->Get(fullName));

    if (! h1_r) {
      std::cout << "Didn't find " << yHistoName << " in file " << file->GetName() << std::endl;
      return;
    }

    if (h1_r->GetMean() > 0.4) {
      if (min < 0)
        min = i;

      max = i;
    } else if (min > 0)
      break;
  }

  for (int i = 0; i < nPoints; ++i) {
    if (i < min)
      continue;
    if (i > max)
      return;

    TString fullName = TString::Format("%s_%d", yHistoName.c_str(), i);
    TH1* h1_r = static_cast<TH1*>(file->Get(fullName));

    if (! h1_r) {
      std::cout << "Didn't find " << yHistoName << " in file " << file->GetName() << std::endl;
      return;
    }

    //Float_t mean, mean_err, rms, rms_err;
    //fitTools::getProjectionMeanAndRMS(h1_r, mean, mean_err, rms, rms_err, 1., 0.9);

    if (FIT_RMS_ == "FIT") {

      TF1* gaussian = new TF1("gaussian", "gaus");
      gaussian->SetLineColor(kRed);
      TH1* newhisto = NULL;
      fitTools::fitProjection_sameArea(h1_r, gaussian, &newhisto, (Float_t)(INTPERC_ / 100.), "RQ");
      //fitTools::fitProjection(h1_r, gaussian, 2., "RQ");

      Float_t mu = gaussian->GetParameter(1);
      Float_t sigma = gaussian->GetParameter(2);
      //TF1* gaussian_chi = new TF1("gaussian_chi", "gaus");
      //fitTools::fitProjection(h1_r, gaussian_chi, 1.5, "RQN");
      Float_t mu_err = gaussian->GetParError(1);
      //Float_t sigma_err = gaussian->GetParError(2);
      //Float_t sigma_err = (LUMI>0) ? sigma/sqrt((Float_t)LUMI*h1_r->Integral()*(Float_t)INT_PERC_/100.) : h1_r->GetRMSError();
      Float_t sigma_err = h1_r->GetRMSError();

      y_resp.push_back(mu);
      y_resp_err.push_back(mu_err);

      y_reso.push_back(sigma / mu);
      y_reso_err.push_back(sqrt(sigma_err * sigma_err / (mu * mu) + sigma * sigma * mu_err * mu_err / (mu * mu * mu * mu)));

    } else if (FIT_RMS_ == "RMS99") {

      Float_t mean99, mean99_err, rms99, rms99_err;
      fitTools::getTruncatedMeanAndRMS(h1_r, mean99, mean99_err, rms99, rms99_err, 0.99, 0.99);

      y_resp.push_back(mean99);
      y_resp_err.push_back(mean99_err);

      y_reso.push_back(rms99 / mean99);
      y_reso_err.push_back(sqrt(rms99_err * rms99_err / (mean99 * mean99) + rms99 * rms99 * mean99_err * mean99_err / (mean99 * mean99 * mean99 * mean99)));

    } else {

      std::cout << "WARNING!! FIT_RMS type '" << FIT_RMS_ << "' currently not supported. Exiting." << std::endl;
      exit(66);

    }

  } //for
}

void drawExtrap::getYPoints(TFile * file, const char* yHistoName, Int_t nPoints, Float_t* y_resp, Float_t* y_resp_err,  Float_t* y_reso, Float_t* y_reso_err, Float_t* y_nevents) const {

  for (int i = 0; i < nPoints ; ++i) {

    TString fullName = TString::Format("%s_%d", yHistoName, i);
    TH1* h1_r = (TH1*) file->Get(fullName);

    if (! h1_r) {
      std::cout << "Didn't find " << fullName << " in file " << file->GetName() << std::endl;
      return;
    }

    //this ugly fix saves from empty relative pt binning plots
    
    
    if (h1_r->GetEntries() == 0) {
      y_resp[i] = -1.;
      y_resp_err[i] = 1000000000.;
      y_reso[i] = -1.;
      y_reso_err[i] = 100000000.;

      continue;
    }

    float mean = h1_r->GetMean();
    float rms = h1_r->GetRMS();
    float mean_err = h1_r->GetMeanError();
    float rms_err = h1_r->GetRMSError();

    float rmsFactor = -1;

    if (FIT_RMS_ == "FIT") {

      TF1* gaussian = new TF1("gaussian", "gaus");
      gaussian->SetLineColor(kRed);
      TH1* newhisto = NULL;
      fitTools::fitProjection_sameArea(h1_r, gaussian, &newhisto, (Float_t)(INTPERC_ / 100.), "RQ");

      float mu = gaussian->GetParameter(1);
      float sigma = gaussian->GetParameter(2);
      float mu_err = gaussian->GetParError(1);
      float sigma_err = h1_r->GetRMSError();

      y_resp[i] = mu;
      y_resp_err[i] = mu_err;

      y_reso[i] = sigma / mu;
      y_reso_err[i] = sqrt(sigma_err * sigma_err / (mu * mu) + sigma * sigma * mu_err * mu_err / (mu * mu * mu * mu));
      y_nevents[i] = h1_r->Integral();
      continue;

    } else if (FIT_RMS_ == "RMS70") {
      rmsFactor = 0.70;
    } else if (FIT_RMS_ == "RMS95") {
      rmsFactor = 0.95;
    } else if (FIT_RMS_ == "RMS99") {
      rmsFactor = 0.99;
    }else if (FIT_RMS_ == "RMS96") {
      rmsFactor = 0.96;
    }else if (FIT_RMS_ == "RMS97") {
      rmsFactor = 0.97;
    } else if (FIT_RMS_ == "RMS98") {
      rmsFactor = 0.98;
    }else if (FIT_RMS_ == "RMS100") {
      rmsFactor = 1.;
    }else if (FIT_RMS_ == "RMS985") {
      rmsFactor = 0.985;
    } else {
      std::cout << "WARNING!! FIT_RMS type '" << FIT_RMS_ << "' currently not supported. Exiting." << std::endl;
      exit(66);
    }

    if (rmsFactor > 0  && h1_r -> Integral()  > 0. ) {
      

      fitTools::getTruncatedMeanAndRMS(h1_r, mean, mean_err, rms, rms_err, rmsFactor, rmsFactor);

    }

    y_resp[i] = mean;
    y_resp_err[i] = mean_err;


    y_reso[i] = rms / mean;
    y_reso_err[i] = sqrt(rms_err * rms_err / (mean * mean) + rms * rms * mean_err * mean_err / (mean * mean * mean * mean));
    y_nevents[i] = h1_r->Integral();
  } //for

} //getYPoints







