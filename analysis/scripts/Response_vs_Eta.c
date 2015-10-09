#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include "TPaveText.h"

#define LIGHT_RED TColor::GetColor(0xcf, 0xa0, 0xa1)

using namespace std;

struct Point {
  double x;
  double y;
  
  double x_err;
  double y_err;
  };

TGraphErrors* get_graphRatio(TGraphErrors* gr_data, TGraphErrors* gr_MC) {

    TGraphErrors* gr_ratio = new TGraphErrors(0);

    // First, build a vector of points
    std::vector<Point> data_points;
    std::vector<Point> mc_points;

    for (int i = 0; i < gr_MC->GetN(); ++i) {
      Point p;
      gr_MC->GetPoint(i, p.x, p.y);

      p.x_err = gr_MC->GetErrorX(i);
      p.y_err = gr_MC->GetErrorY(i);      

      mc_points.push_back(p);
    }

    for (int i = 0; i < gr_data->GetN(); ++i) {
      Point p;
      gr_data->GetPoint(i, p.x, p.y);

      p.x_err = gr_data->GetErrorX(i);
      p.y_err = gr_data->GetErrorY(i);      

      data_points.push_back(p);
    }

    int index = 0;
    for (Point data: data_points) {

      // Find corresponding point inside other array
      Point mc;
      bool found = false;
      for (Point foo: mc_points) {
        if (abs(foo.x - data.x) < 1e-6) {
          found = true;
          mc = foo;

          break;
        }
      }

      if (! found)
        continue;

      double ratio_x = data.x;
      double ratio_x_err = data.x_err;

      double ratio_y = data.y / mc.y;
      double ratio_y_err = sqrt(data.y_err * data.y_err / (mc.y * mc.y) + data.y * data.y * mc.y_err * mc.y_err / (mc.y * mc.y * mc.y * mc.y));

      if (std::isnan(ratio_y))
        continue;

      gr_ratio->SetPoint(index, ratio_x, ratio_y);
      gr_ratio->SetPointError(index, ratio_x_err, ratio_y_err);

      index++;
    }

    return gr_ratio;
  }


void DrawRatioAndSave( const char *nameFile,  TGraphErrors* data, TGraphErrors* mc, int xmin, int xmax, const char *XTitle, const char *YTitle){
  char * nameSaved ;
  nameSaved = new char[strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, nameFile);
  
  gStyle -> SetOptStat(kFALSE);
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);
  c1->cd();
  
  // Data / MC comparison
  TPad* pad_hi = new TPad("pad_hi", "", 0., 0.3, 1., 1.);
  pad_hi->Draw();
  // pad_hi->SetLogx();
  pad_hi->SetLeftMargin(0.15);
  pad_hi->SetBottomMargin(0.015);
  
  // Data / MC ratio
  TPad* pad_lo = new TPad("pad_lo", "", 0., 0.05, 1, 0.305);
  pad_lo->Draw();
  //pad_lo->SetLogx();
  pad_lo->SetLeftMargin(0.15);
  pad_lo->SetTopMargin(1.);
  pad_lo->SetBottomMargin(0.3);
  
  TGraphErrors* gr_resp_ratio = 0;
  
  TH2* h2_axes_lo_resp = NULL;
  
  TLine* line_one = new TLine(xmin, 1., xmax, 1.);
  TLine* line_plus_resp = new TLine(xmin, 1.05, xmax, 1.05);
  TLine* line_minus_resp = new TLine(xmin, 0.95, xmax, 0.95);
  
  pad_lo->cd();
  
  //  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 0, 5, 10, 0.86, 1.14);
  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, xmin, xmax, 10, 0.5, 1.2);
  
  h2_axes_lo_resp->SetXTitle(XTitle);     //"|#eta (jet)|");
  h2_axes_lo_resp->SetYTitle("Data / MC");
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.55);
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();
  //h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.);
  h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);
  h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.09);
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
  h2_axes_lo_resp->GetYaxis()->SetNdivisions(7, kTRUE);
  h2_axes_lo_resp->Draw("");
  
  line_one->Draw("same");
   
  //line_plus_resp->SetLineColor(46);
  line_plus_resp->SetLineWidth(2);
  line_plus_resp->SetLineStyle(2);  
  //line_minus_resp->SetLineColor(46);
  line_minus_resp->SetLineWidth(2);
  line_minus_resp->SetLineStyle(2);
  
  //this->
  gr_resp_ratio = get_graphRatio(data, mc);
  gr_resp_ratio->SetName("response_ratio");
  gr_resp_ratio->SetMarkerStyle(20);
  gr_resp_ratio->SetMarkerSize(1.5);
  gr_resp_ratio->SetMarkerColor(kBlue - 6);  

  // TF1* ratioFit = new TF1("ratioFit", "[0]", 0, 5);
  //ratioFit->SetParameter(0, 0.);
  //ratioFit->SetLineColor(46);
  //ratioFit->SetLineWidth(2);
  // gr_resp_ratio->Fit(ratioFit, "RQ");
  //std::cout << "-> ChiSquare: " << constline->GetChisquare() << "   NDF: " << constline->GetNDF() << std::endl;
  
  //double fitValue = ratioFit->GetParameter(0);
  //double fitError = ratioFit->GetParError(0);
  
  //  TH1D* errors = new TH1D("errors", "errors", 100, 0, 5);
  //(TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors, 0.68);
  //errors->SetStats(false);
  //errors->SetFillColor(LIGHT_RED);
  //errors->SetLineColor(46);

  //  TPaveText* label = new TPaveText(0.55, 0.77, 0.88, 0.83, "brNDC");
  //fitlabel->SetTextSize(0.08);
  //fitlabel->SetFillColor(0);
  // TString fitLabelText = TString::Format("Fit: %.3f #pm %.3f", fitValue, fitError);
  //fitlabel->AddText(fitLabelText);
  //  fitlabel->Draw("same");
  
  line_plus_resp->Draw("same");
  line_minus_resp->Draw("same");
  
  //  errors->Draw("e3 same");
  //ratioFit->Draw("same");
  
  gr_resp_ratio->Draw("Psame");
  
  gPad->RedrawAxis();
  
  pad_hi->cd();
  
  TH2D* h2_axes = new TH2D("axes_again", "", 10, xmin, xmax, 10, 0.4, 1.25);
  h2_axes->SetYTitle("Jet p_{T} response");
  //h2_axes->SetYTitle("< p_{T}^{jet} / p_{T}^{#gamma} >");
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);
  h2_axes->GetYaxis()->SetTitleSize(0.045);
  //h2_axes->GetXaxis()->SetMoreLogLabels();
  //h2_axes->GetXaxis()->SetNoExponent();


  h2_axes->Draw();
  
  //  TPaveText* label = new TPaveText(0.55, 0.77, 0.88, 0.83, "brNDC");
  //  label->SetTextSize(0.08);
  //  label->SetFillColor(0);
  //  TString fitLabelText = TString::Format("|#eta|< 1.3");
  //label->AddText("|#eta|< 1.3");
  //label->Draw("same");


  //  Float_t labelTextSize = 0.035;
  //  TPaveText* label_algo = get_labelAlgo(2);

  TLegend* legend = new TLegend(0.85, 0.85, 0.75, 0.75);
  legend->SetTextFont(42);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.036);
  //  legend->SetHeader("3.0 <|#eta|< 5.0");
  legend->AddEntry(data, "Data", "P");
  legend->AddEntry(mc, "MC (#gamma+Jet)", "P");
 
  legend->Draw("same");

  mc->SetMarkerStyle(24);
  mc->SetMarkerSize(1.5);
  mc->SetMarkerColor(kBlue - 6);
  mc->SetLineColor(kBlue - 6);
  mc->Draw("Psame");

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.5);
  data->SetMarkerColor(kBlue - 6);
  data->Draw("Psame");
   
  gPad->RedrawAxis();

  c1->SaveAs(nameSaved);
  c1->Destructor();
}


int main(int argc, char* argv[]) {
  
  string * inputFile1 = new string(argv[1]);
  string * inputFile2 = new string(argv[2]);

  char * input_file1 = new char [inputFile1->length()+1];
  strcpy (input_file1, inputFile1->c_str());

  char * input_file2 = new char [inputFile2->length()+1];
  strcpy (input_file2, inputFile2->c_str());
  
  TFile file1(input_file1);
  TFile file2(input_file2);

  //  TFile file1( );
  //  TFile file2("PhotonJet_GJet_Pt-15To6000_reduced_vs_Eta_PFlowAK4chs.root");


  TGraphErrors* gr_response_vs_eta = new TGraphErrors(0);
  gr_response_vs_eta->SetName("response_vs_eta");
  TGraphErrors* gr_responseMC_vs_eta = new TGraphErrors(0);
  gr_responseMC_vs_eta->SetName("responseMC_vs_eta");
  //----------------------------------------------------------------------------------------------------------

  TGraphErrors* gr_response_vs_npv_Eta013 = new TGraphErrors(0);
  gr_response_vs_npv_Eta013->SetName("response_vs_npv_Eta013");
  TGraphErrors* gr_responseMC_vs_npv_Eta013 = new TGraphErrors(0);
  gr_responseMC_vs_npv_Eta013->SetName("responseMC_vs_npv_Eta013");

  TGraphErrors* gr_response_vs_npv_Eta030 = new TGraphErrors(0);
  gr_response_vs_npv_Eta030->SetName("response_vs_npv_Eta030");
  TGraphErrors* gr_responseMC_vs_npv_Eta030 = new TGraphErrors(0);
  gr_responseMC_vs_npv_Eta030->SetName("responseMC_vs_npv_Eta030");

  TGraphErrors* gr_response_vs_npv_Eta050 = new TGraphErrors(0);
  gr_response_vs_npv_Eta050->SetName("response_vs_npv_Eta050");
  TGraphErrors* gr_responseMC_vs_npv_Eta050 = new TGraphErrors(0);
  gr_responseMC_vs_npv_Eta050->SetName("responseMC_vs_npv_Eta050");


  ///////////////////////////////////////////////////

  // double etaMean = 0.4;
  
  std::vector<std::pair<float, float> > etaBins;
  
  etaBins.push_back(std::make_pair(0., 0.8) );
  etaBins.push_back(std::make_pair(0.8, 1.3) );
  etaBins.push_back(std::make_pair(1.3, 1.9) );
  etaBins.push_back(std::make_pair(1.9, 2.5) );
  etaBins.push_back(std::make_pair(2.5, 3.0) );
  etaBins.push_back(std::make_pair(3.0, 3.2) );
  etaBins.push_back(std::make_pair(3.2, 5.2) );
  
  for( int iplot =0 ; iplot<7 ; iplot++){

    std::pair<float, float> currentBin = etaBins[iplot];
    float etaMean = (currentBin.first + currentBin.second) / 2.;

    double dataResponse=0;
    double dataResponseErr=0;
    double mcResponse=0;
    double mcResponseErr=0;

    TH1D *h_data;
    TH1D *h_mc;

    if(iplot ==0){  
      h_data     = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_0_8");
      h_mc       = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_0_8");
    }
    if(iplot ==1){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_8_13");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_8_13");
    }
    if(iplot ==2){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_13_19");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_13_19");
    }
    if(iplot ==3){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_19_25");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_19_25");
    }
    if(iplot ==4){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_25_30");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_25_30");   
    }
    if(iplot ==5){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_30_32");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_30_32");
    }
    if(iplot ==6){  
      h_data    = (TH1D*)file1.Get("analysis/responseEtaBinned_passedID_Eta_32_52");
      h_mc      = (TH1D*)file2.Get("analysis/responseEtaBinned_passedID_Eta_32_52");
    }

    dataResponse      = h_data ->GetMean();
    dataResponseErr = h_data->GetMeanError();
    mcResponse        = h_mc->GetMean();
    mcResponseErr   = h_mc->GetMeanError();
    
    cout<< "dataResponse  "<< dataResponse << endl;
    cout<< "mcResponse  "<< mcResponse << endl;
    
    gr_response_vs_eta->SetPoint(iplot, etaMean, dataResponse);
    gr_response_vs_eta->SetPointError(iplot, 0., dataResponseErr);
    
    gr_responseMC_vs_eta->SetPoint(iplot, etaMean, mcResponse);
    gr_responseMC_vs_eta->SetPointError(iplot, 0., mcResponseErr);
    
  }

  //void DrawRatioAndSave( const char *nameFile,  TGraphErrors* data, TGraphErrors* mc, int xmin, int xmax, const char *XTitle, const char *YTitle){
  
  DrawRatioAndSave("Response_vs_eta.png", gr_response_vs_eta, gr_responseMC_vs_eta, 0, 5, "|#eta (jet)|", "Jet p_{T} Response"); //works

   /////////////////////////////////////////////////////////////////////// end vs_eta

  std::vector<std::pair<float, float> > mVertexBins;
 
  mVertexBins.push_back(std::make_pair(0, 5));
   mVertexBins.push_back(std::make_pair(5, 8));
  mVertexBins.push_back(std::make_pair(8, 11));
  mVertexBins.push_back(std::make_pair(11, 13));
  mVertexBins.push_back(std::make_pair(13, 15));
  mVertexBins.push_back(std::make_pair(15, 18));
  mVertexBins.push_back(std::make_pair(18, 21));
  mVertexBins.push_back(std::make_pair(21, 23));
  mVertexBins.push_back(std::make_pair(23, 35));
  

  for( int range = 0 ; range <3 ; range++){
    for( int iplot =0 ; iplot<9 ; iplot++){
      
      std::pair<float, float> currentBin = mVertexBins[iplot];
      float vertexMean = (currentBin.first + currentBin.second) / 2.;
      
      double dataResponse=0;
      double dataResponseErr=0;
      double mcResponse=0;
      double mcResponseErr=0;
      
      TH1D *h_data;
      TH1D *h_mc;
      
      TString HistoName;

      if(range ==0)    HistoName = TString::Format("analysis/vertex/resp_balancing_eta013_nvertex_%d_%d", (int) currentBin.first, (int) currentBin.second);
      if(range ==1)    HistoName = TString::Format("analysis/vertex/resp_balancing_eta030_nvertex_%d_%d", (int) currentBin.first, (int) currentBin.second);
      if(range ==2)    HistoName = TString::Format("analysis/vertex/resp_balancing_eta050_nvertex_%d_%d", (int) currentBin.first, (int) currentBin.second);
      
      cout<< HistoName <<endl;
      
      h_data     = (TH1D*)file1.Get(HistoName);
      h_mc       = (TH1D*)file2.Get(HistoName);
      
      dataResponse      = h_data ->GetMean();
      dataResponseErr = h_data->GetMeanError();
      mcResponse        = h_mc->GetMean();
      mcResponseErr   = h_mc->GetMeanError();
      
      cout<< "dataResponse  "<< dataResponse << endl;
      cout<< "mcResponse  "<< mcResponse << endl;

      if(range ==0){
	gr_response_vs_npv_Eta013      ->SetPoint(iplot, vertexMean, dataResponse);
	gr_responseMC_vs_npv_Eta013 ->SetPoint(iplot, vertexMean, mcResponse);
      }
      if(range ==1){
	gr_response_vs_npv_Eta030      ->SetPoint(iplot, vertexMean, dataResponse);
	gr_responseMC_vs_npv_Eta030 ->SetPoint(iplot, vertexMean, mcResponse);
      }
      if(range ==2){
	gr_response_vs_npv_Eta050      ->SetPoint(iplot, vertexMean, dataResponse);
	gr_responseMC_vs_npv_Eta050 ->SetPoint(iplot, vertexMean, mcResponse);
      }      
    }
  }
  
  DrawRatioAndSave("Response_vs_npv_Eta013.png", gr_response_vs_npv_Eta013, gr_responseMC_vs_npv_Eta013, 0, 35, "n_{pv}", "Jet p_{T} Response"); //works
  DrawRatioAndSave("Response_vs_npv_Eta030.png", gr_response_vs_npv_Eta030, gr_responseMC_vs_npv_Eta030, 0, 35, "n_{pv}", "Jet p_{T} Response"); //works
  DrawRatioAndSave("Response_vs_npv_Eta050.png", gr_response_vs_npv_Eta050, gr_responseMC_vs_npv_Eta050, 0, 35, "n_{pv}", "Jet p_{T} Response"); //works


}






