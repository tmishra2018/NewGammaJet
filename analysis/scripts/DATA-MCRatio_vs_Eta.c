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
  
  TGraphErrors* gr_ratio_vs_eta = new TGraphErrors(0);
  gr_ratio_vs_eta->SetName("ratio_vs_eta");

  //----------------------------------------------------------------------------------------------------------

  std::vector<std::pair<float, float> > etaBins;  
  //  etaBins.push_back(std::make_pair(0., 0.8) );
  etaBins.push_back(std::make_pair(0, 1.3) );
  etaBins.push_back(std::make_pair(1.3, 2.0) );
  etaBins.push_back(std::make_pair(2.0, 2.5) );
  etaBins.push_back(std::make_pair(2.5, 3.0) );
  etaBins.push_back(std::make_pair(3.0, 3.2) );
  etaBins.push_back(std::make_pair(3.2, 5.2) );
  
  std::vector<float> etaRange;  
  //  etaRange.push_back(std::make_pair(0., 0.8) );
  etaRange.push_back(0.65);
  etaRange.push_back(0.35);
  etaRange.push_back(0.25);
  etaRange.push_back(0.25);
  etaRange.push_back(0.1);
  etaRange.push_back(1.0);


  for(int plot=0 ; plot <6; plot++){

    std::vector<float> fitValue;  
    std::vector<float> fitValueError;        

    if(plot == 0){      
      fitValue.push_back(0.972);
      fitValue.push_back(0.936);
      fitValue.push_back(0.947);
      fitValue.push_back(0.885);
      fitValue.push_back(0.751);
      fitValue.push_back(0.846);

      fitValueError.push_back(0.005);
      fitValueError.push_back(0.009);
      fitValueError.push_back(0.012);
      fitValueError.push_back(0.019);
      fitValueError.push_back(0.070);
      fitValueError.push_back(0.090);   
    }

    if(plot == 1){      
      fitValue.push_back(0.958);
      fitValue.push_back(0.919);
      fitValue.push_back(0.910);
      fitValue.push_back(0.848);
      fitValue.push_back(0.819);
      fitValue.push_back(0.713);
      
      fitValueError.push_back(0.004);
      fitValueError.push_back(0.007);
      fitValueError.push_back(0.008);
      fitValueError.push_back(0.012);
      fitValueError.push_back(0.081);
      fitValueError.push_back(0.048);   
    }

    if(plot == 2){      
      fitValue.push_back(0.998);
      fitValue.push_back(0.928);
      fitValue.push_back(0.924);
      fitValue.push_back(0.946);
      fitValue.push_back(0.384);
      fitValue.push_back(0.541);

      fitValueError.push_back(0.007);
      fitValueError.push_back(0.011);
      fitValueError.push_back(0.015);
      fitValueError.push_back(0.024);
      fitValueError.push_back(0.578);
      fitValueError.push_back(0.311);   
    }

    if(plot == 3){      
      fitValue.push_back(0.957);
      fitValue.push_back(0.877);
      fitValue.push_back(0.971);
      fitValue.push_back(0.766);
      fitValue.push_back(0.606);
      fitValue.push_back(0.249);
      
      fitValueError.push_back(0.010);
      fitValueError.push_back(0.016);
      fitValueError.push_back(0.025);
      fitValueError.push_back(0.024);
      fitValueError.push_back(0.957);
      fitValueError.push_back(0.199);   
    }


    if(plot == 4){      
      fitValue.push_back(0.994);
      fitValue.push_back(0.960);
      fitValue.push_back(0.977);
      fitValue.push_back(0.971);
      fitValue.push_back(0.872);
      fitValue.push_back(0.951);
      
      fitValueError.push_back(0.004);
      fitValueError.push_back(0.007);
      fitValueError.push_back(0.009);
      fitValueError.push_back(0.014);
      fitValueError.push_back(0.032);
      fitValueError.push_back(0.042);   
    }

    if(plot == 5){      
      fitValue.push_back(0.971);
      fitValue.push_back(0.959);
      fitValue.push_back(1.074);
      fitValue.push_back(1.027);
      fitValue.push_back(0.455);
      fitValue.push_back(0.687);
      
      fitValueError.push_back(0.011);
      fitValueError.push_back(0.025);
      fitValueError.push_back(0.028);
      fitValueError.push_back(0.042);
      fitValueError.push_back(0.14);
      fitValueError.push_back(0.399);   
    }


  for( int point =0 ; point<6 ; point++){

    //  int point = 0;
  
  std::pair<float, float> currentBin = etaBins[point];
  float etaMean = (currentBin.first + currentBin.second) / 2.;
  
  std::cout<<"Point number  " << point << std::endl;
  std::cout<<"x   " << etaMean << std::endl;
  std::cout<<"y   "<<  fitValue.at(point)<<"  #pm  "<<fitValueError.at(point)<< std::endl;
  
  gr_ratio_vs_eta->SetPoint(point, etaMean, fitValue.at(point) );
  gr_ratio_vs_eta->SetPointError(point, etaRange.at(point) , fitValueError.at(point));
  
  }

  double xmin = 0.;
  double xmax = 5.2;


  TLine* line_one = new TLine(xmin, 1., xmax, 1.);                                                                                                           
  TLine* line_plus = new TLine(xmin, 1.05, xmax, 1.05);                                                                                                 
  TLine* line_minus = new TLine(xmin, 0.95, xmax, 0.95);                                                                                                
  line_one -> SetLineColor(kRed);                                                                                                                            
  line_plus -> SetLineStyle(2);                                                                                                                              
  line_minus -> SetLineStyle(2);                                             

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);                                                                                                         

  double ymin = 0.7;

  //  if(plot == 0) ymin = 0.7;
  //  if(plot == 1) ymin = 0.7;
  //  if(plot == 2) ymin = 0.35;
  //  if(plot == 3) ymin = 0.2;
  //  if(plot == 0) ymin = 0.7;

  TH2D *h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, xmin, xmax, 10, ymin, 1.1);
  gStyle->SetOptStat(kFALSE);
  h2_axes_lo_resp->SetXTitle("#eta (jet)");
  h2_axes_lo_resp->SetYTitle("Data / MC");
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(1.5);                                                                                               
  h2_axes_lo_resp->Draw();
  
  gr_ratio_vs_eta->SetMarkerStyle(20);                                                                                                            
  gr_ratio_vs_eta->SetMarkerSize(1.5);                                                                                                            
  gr_ratio_vs_eta->SetMarkerColor(kBlue - 6);                                                                                                     

  gr_ratio_vs_eta->SetLineColor(kBlue - 6);                                                                                                        
  gr_ratio_vs_eta->Draw("P same");                        
  line_one->Draw("same");                                                                                                                                   
  line_plus->Draw("same");                                                                                                                                  
  line_minus->Draw("same");                                                                                                                                 
  


   if(plot==0)   c1->SaveAs("Ratio_vs_eta_BAL_NoRes_NoExtrap.png");     
   if(plot==1)   c1->SaveAs("Ratio_vs_eta_MPF_NoRes_NoExtrap.png");     
   if(plot==2)   c1->SaveAs("Ratio_vs_eta_BAL_NoRes_SiExtrap.png");     
   if(plot==3)   c1->SaveAs("Ratio_vs_eta_MPF_NoRes_SiExtrap.png");     
   if(plot==4)   c1->SaveAs("Ratio_vs_eta_MPF_SiRes_NoExtrap.png");     
   if(plot==5)   c1->SaveAs("Ratio_vs_eta_MPF_SiRes_SiExtrap.png");     


  }

}






