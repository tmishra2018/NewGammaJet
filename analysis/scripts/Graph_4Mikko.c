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


int main(int argc, char* argv[]) {
  
  bool make_DataMC = true;
  bool make_Pythia = true;

  if(make_DataMC){

  TFile file1("../tuples/Data/PhotonJet_SinglePhoton_GoldenJson_v2_SiPrescale_alphacut030_PFlowAK4chs.root" );
  TFile file2("../tuples/GJET_MC/PhotonJet_GJet_AlphaCut030_PFlowAK4chs.root");

  TGraphErrors* gr_MPF_Data = (TGraphErrors*)file1.Get("analysis/new_extrapolation/extrap_resp_mpf_eta0013_graph");
  TGraphErrors* gr_BAL_Data = (TGraphErrors*)file1.Get("analysis/new_extrapolation/extrap_resp_balancing_eta0013_graph");

  TGraphErrors* gr_MPF_MC = (TGraphErrors*)file2.Get("analysis/new_extrapolation/extrap_resp_mpf_eta0013_graph");
  TGraphErrors* gr_BAL_MC = (TGraphErrors*)file2.Get("analysis/new_extrapolation/extrap_resp_balancing_eta0013_graph");

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
                                                                                                                                                                           
  TGraphErrors* gr_BAL_ratio = 0;                                                                                                                                         
  TGraphErrors* gr_MPF_ratio = 0;                                                                                                                                         
  
  TH2* h2_axes_lo_resp = NULL;                                                                                                                                             
                                                                                                                                                                           
  TLine* line_one = new TLine(0, 1., 0.30, 1.);                                                                                                                         
  TLine* line_plus_resp = new TLine(0, 1.05, 0.30, 1.05);                                                                                                               
  TLine* line_minus_resp = new TLine(0, 0.95, 0.30, 0.95);                                                                                                              
                                                                                                                                                                           
  pad_lo->cd();         

 //  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 0, 5, 10, 0.86, 1.14);                                                                                            
  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 0, 0.30, 10, 0.90, 1.10);                                                                                            
                                                                                                                                                                           
  h2_axes_lo_resp->SetXTitle("#alpha");     //"|#eta (jet)|");                                                                                                               
  h2_axes_lo_resp->SetYTitle("Data / MC");                                                                                                                                 
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);                                                                                                                        
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.55);                                                                                                                       
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);                                                                                                                        
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();                                                                                                                            
  //h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.);                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);                                                                                                                        
  h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.1);                                                                                                                         
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.1);                                                                                                                         
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
  gr_BAL_ratio = get_graphRatio(gr_BAL_Data, gr_BAL_MC);                                                                                                                                
  gr_BAL_ratio->SetName("BAL_ratio");                                                                                                                                
  gr_BAL_ratio->SetMarkerStyle(24);                                                                                                                                       
  gr_BAL_ratio->SetMarkerSize(1.5);                                                                                                                                       
  gr_BAL_ratio->SetMarkerColor(kRed);

  gr_MPF_ratio = get_graphRatio(gr_MPF_Data, gr_MPF_MC);                                                                                                                                
  gr_MPF_ratio->SetName("MPF_ratio");                                                                                                                                
  gr_MPF_ratio->SetMarkerStyle(20);                                                                                                                                       
  gr_MPF_ratio->SetMarkerSize(1.5);                                                                                                                                       
  gr_MPF_ratio->SetMarkerColor(kViolet -6);
                                                                                                                                                             
  line_plus_resp->Draw("same");                                                                                                                                            
  line_minus_resp->Draw("same");                                                                                                                                           
                                                                                                                                                                           
  //  errors->Draw("e3 same");                                                                                                                                             
  //ratioFit->Draw("same");                                                                                                                                                
                                                                                                                                                                           
  gr_BAL_ratio->Draw("Psame");                                                                                                                                            
  gr_MPF_ratio->Draw("Psame");                                                                                                                                            
                                          
  TLegend* legend2 = new TLegend(0.40, 0.88, 0.20, 0.78);                                                                                                                   
  legend2->SetTextFont(42);                                                                                                                                                 
  legend2->SetBorderSize(0);                                                                                                                                                
  legend2->SetFillColor(kWhite);                                                                                                                                            
  legend2->SetFillStyle(0);                                                                                                                                                 
  legend2->SetTextSize(0.1);                                                                                                                                              
  legend2->AddEntry(gr_BAL_ratio, "Balancing", "P");                                                                                                                                     
  //  legend2->AddEntry(gr_MPF_ratio, "MPF", "P");                                                                                                                            

  TLegend* legend3 = new TLegend(0.60, 0.88, 0.40, 0.78);                                                                                                                   
  legend3->SetTextFont(42);                                                                                                                                                 
  legend3->SetBorderSize(0);                                                                                                                                                
  legend3->SetFillColor(kWhite);                                                                                                                                            
  legend3->SetFillStyle(0);                                                                                                                                                 
  legend3->SetTextSize(0.1);                                                                                                                                              
  // legend3->AddEntry(gr_BAL_ratio, "Balancing", "P");                                                                                                                                     
    legend3->AddEntry(gr_MPF_ratio, "MPF", "P");                                                                                                                            


  legend3->Draw("same");   
  legend2->Draw("same");   

                                                                                                                                 
  gPad->RedrawAxis();                                                                                                                                                      
                                                                                                                                                                           
  pad_hi->cd();                                                                                                                                                            
                                                                                                                                                                           
  TH2D* h2_axes = new TH2D("axes_again", "", 10, 0, 0.30, 10, 0.75, 1.055);                                                                                               
  h2_axes->SetYTitle("Jet p_{T} response");                                                                                                                                
  //h2_axes->SetYTitle("< p_{T}^{jet} / p_{T}^{#gamma} >");                                                                                                                
  h2_axes->GetXaxis()->SetTitleOffset(1.1);                                                                                                                                
  h2_axes->GetYaxis()->SetTitleOffset(1.2);                                                                                                                                
  h2_axes->GetYaxis()->SetTitleSize(0.045);                                                                                                                                
  //h2_axes->GetXaxis()->SetMoreLogLabels();                                                                                                                               
  //h2_axes->GetXaxis()->SetNoExponent();                                                                                                                                  
                                                                                                                                                                           
                                                                                                                                                                           
  h2_axes->Draw(); 

  TLegend* legend = new TLegend(0.40, 0.25, 0.20, 0.05);                                                                                                                   
  legend->SetTextFont(42);                                                                                                                                                 
  legend->SetBorderSize(0);                                                                                                                                                
  legend->SetFillColor(kWhite);                                                                                                                                            
  legend->SetFillStyle(0);                                                                                                                                                 
  legend->SetTextSize(0.036);                                                                                                                                              
  //  legend->SetHeader("3.0 <|#eta|< 5.0");                                                                                                                               
  legend->AddEntry(gr_BAL_Data, "Balancing (Data)", "P");                                                                                                                                     
  legend->AddEntry(gr_BAL_MC, "Balancing (MC)", "P");                                                                                                                            
  legend->AddEntry(gr_MPF_Data, "MPF (Data)", "P");                              
  legend->AddEntry(gr_MPF_MC, "MPF (MC)", "P");

  legend->Draw("same");                                                                                                                                                    
                                                                                                                                                                           
  gr_BAL_MC->SetMarkerStyle(25);                                                                                                                                                  
  gr_BAL_MC->SetMarkerSize(1.5);                                                                                                                                                  
  gr_BAL_MC->SetMarkerColor(kOrange +1);                                                                                                                                           
  gr_BAL_MC->SetLineColor(kOrange +1);                                                                                                                                             
  gr_BAL_MC->Draw("Psame");                                                                                                                                                       

  gr_MPF_MC->SetMarkerStyle(21);                                                                                                                                                  
  gr_MPF_MC->SetMarkerSize(1.5);                                                                                                                                                  
  gr_MPF_MC->SetMarkerColor(kBlue);                                                                                                                                           
  gr_MPF_MC->SetLineColor(kBlue);                                                                                                                                             
  gr_MPF_MC->Draw("Psame");                                                                                                                                                       
                                                                                                                                                                           
  gr_BAL_Data->SetMarkerStyle(24);                                                                                                                                                
  gr_BAL_Data->SetMarkerSize(1.5);
  gr_BAL_Data->SetMarkerColor(kRed);  
  gr_BAL_Data->SetLineColor(kRed);
  gr_BAL_Data->Draw("Psame");                                                                                                                                                     

  gr_MPF_Data->SetMarkerStyle(20);                                                                                                                                                
  gr_MPF_Data->SetMarkerSize(1.5);
  gr_MPF_Data->SetMarkerColor(kViolet -6);                                                                                                                                         
  gr_MPF_Data->SetLineColor(kViolet -6);                                                                                                                                         
  gr_MPF_Data->Draw("Psame");                                                                                                                                                     
                                                                                                                                                                           
  gPad->RedrawAxis();                                                                                                                                                      
                                                                                                                                                                           
  c1->SaveAs("Comparison_DataMC.png");

  }

  if(make_Pythia){


  TFile file1("../tuples/GJET_MC/PhotonJet_GJet_AlphaCut030_PFlowAK4chs.root");
  TFile file2("PhotonJet_G_alpha03_PFlowAK5chs.root" );

  TGraphErrors* gr_MPF_MC = (TGraphErrors*)file1.Get("analysis/new_extrapolation/extrap_resp_mpf_eta0013_graph");
  TGraphErrors* gr_BAL_MC = (TGraphErrors*)file1.Get("analysis/new_extrapolation/extrap_resp_balancing_eta0013_graph");


  TGraphErrors* gr_MPF_MC_Pythia6 = (TGraphErrors*)file2.Get("analysis/new_extrapolation/extrap_resp_mpf_eta013_graph");
  TGraphErrors* gr_BAL_MC_Pythia6 = (TGraphErrors*)file2.Get("analysis/new_extrapolation/extrap_resp_balancing_eta013_graph");


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
                                                                                                                                                                           
  TGraphErrors* gr_BAL_ratio = 0;                                                                                                                                         
  TGraphErrors* gr_MPF_ratio = 0;                                                                                                                                         
  
  TH2* h2_axes_lo_resp = NULL;                                                                                                                                             
                                                                                                                                                                           
  TLine* line_one = new TLine(0, 1., 0.30, 1.);                                                                                                                         
  TLine* line_plus_resp = new TLine(0, 1.05, 0.30, 1.05);                                                                                                               
  TLine* line_minus_resp = new TLine(0, 0.95, 0.30, 0.95);                                                                                                              
                                                                                                                                                                           
  pad_lo->cd();         

 //  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 0, 5, 10, 0.86, 1.14);                                                                                            
  h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 10, 0, 0.30, 10, 0.90, 1.10);                                                                                            
                                                                                                                                                                           
  h2_axes_lo_resp->SetXTitle("#alpha");     //"|#eta (jet)|");                                                                                                               
  h2_axes_lo_resp->SetYTitle("MC (P.8) / MC (P.6)");                                                                                                                                 
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);                                                                                                                        
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.55);                                                                                                                       
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);                                                                                                                        
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();                                                                                                                            
  //h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.);                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);                                                                                                                        
  h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);                                                                                                                         
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.1);                                                                                                                         
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.1);                                                                                                                         
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
  gr_BAL_ratio = get_graphRatio(gr_BAL_MC, gr_BAL_MC_Pythia6);                                                                                                                                
  gr_BAL_ratio->SetName("BAL_ratio");                                                                                                                                
  gr_BAL_ratio->SetMarkerStyle(24);                                                                                                                                       
  gr_BAL_ratio->SetMarkerSize(1.5);                                                                                                                                       
  gr_BAL_ratio->SetMarkerColor(kRed);

  gr_MPF_ratio = get_graphRatio(gr_MPF_MC, gr_MPF_MC_Pythia6);                                                                                                                                
  gr_MPF_ratio->SetName("MPF_ratio");                                                                                                                                
  gr_MPF_ratio->SetMarkerStyle(20);                                                                                                                                       
  gr_MPF_ratio->SetMarkerSize(1.5);                                                                                                                                       
  gr_MPF_ratio->SetMarkerColor(kViolet -6);
                                                                                                                                                             
  line_plus_resp->Draw("same");                                                                                                                                            
  line_minus_resp->Draw("same");                                                                                                                                           
                                                                                                                                                                           
  //  errors->Draw("e3 same");                                                                                                                                             
  //ratioFit->Draw("same");                                                                                                                                                
                                                                                                                                                                           
  gr_BAL_ratio->Draw("Psame");                                                                                                                                            
  gr_MPF_ratio->Draw("Psame");                                                                                                                                            
                                          
  TLegend* legend2 = new TLegend(0.40, 0.88, 0.20, 0.78);                                                                                                                   
  legend2->SetTextFont(42);                                                                                                                                                 
  legend2->SetBorderSize(0);                                                                                                                                                
  legend2->SetFillColor(kWhite);                                                                                                                                            
  legend2->SetFillStyle(0);                                                                                                                                                 
  legend2->SetTextSize(0.1);                                                                                                                                              
  legend2->AddEntry(gr_BAL_ratio, "Balancing", "P");                                                                                                                                     
  //  legend2->AddEntry(gr_MPF_ratio, "MPF", "P");                                                                                                                            

  TLegend* legend3 = new TLegend(0.60, 0.88, 0.40, 0.78);                                                                                                                   
  legend3->SetTextFont(42);                                                                                                                                                 
  legend3->SetBorderSize(0);                                                                                                                                                
  legend3->SetFillColor(kWhite);                                                                                                                                            
  legend3->SetFillStyle(0);                                                                                                                                                 
  legend3->SetTextSize(0.1);                                                                                                                                              
  // legend3->AddEntry(gr_BAL_ratio, "Balancing", "P");                                                                                                                                     
  legend3->AddEntry(gr_MPF_ratio, "MPF", "P");                                                                                                                            

  legend3->Draw("same");   
  legend2->Draw("same");   

                                                                                                                                 
  gPad->RedrawAxis();                                                                                                                                                      
                                                                                                                                                                           
  pad_hi->cd();                                                                                                                                                            
                                                                                                                                                                           
  TH2D* h2_axes = new TH2D("axes_again", "", 10, 0, 0.30, 10, 0.75, 1.055);                                                                                               
  h2_axes->SetYTitle("Jet p_{T} response");                                                                                                                                
  //h2_axes->SetYTitle("< p_{T}^{jet} / p_{T}^{#gamma} >");                                                                                                                
  h2_axes->GetXaxis()->SetTitleOffset(1.1);                                                                                                                                
  h2_axes->GetYaxis()->SetTitleOffset(1.2);                                                                                                                                
  h2_axes->GetYaxis()->SetTitleSize(0.045);                                                                                                                                
  //h2_axes->GetXaxis()->SetMoreLogLabels();                                                                                                                               
  //h2_axes->GetXaxis()->SetNoExponent();                                                                                                                                  
                                                                                                                                                                           
                                                                                                                                                                           
  h2_axes->Draw(); 

  TLegend* legend = new TLegend(0.40, 0.25, 0.20, 0.05);                                                                                                                   
  legend->SetTextFont(42);                                                                                                                                                 
  legend->SetBorderSize(0);                                                                                                                                                
  legend->SetFillColor(kWhite);                                                                                                                                            
  legend->SetFillStyle(0);                                                                                                                                                 
  legend->SetTextSize(0.036);                                                                                                                                              
  //  legend->SetHeader("3.0 <|#eta|< 5.0");                                                                                                                               
  legend->AddEntry(gr_BAL_MC, "Balancing (Pythia8)", "P");                                                                                                                                     
  legend->AddEntry(gr_BAL_MC_Pythia6, "Balancing (Pythia6)", "P");                                                                                                                            
  legend->AddEntry(gr_MPF_MC, "MPF (Pythia8)", "P");                              
  legend->AddEntry(gr_MPF_MC_Pythia6, "MPF (Pythia6)", "P");

  legend->Draw("same");                                                                                                                                                    
                                                                                                                                                                           
  gr_BAL_MC->SetMarkerStyle(25);                                                                                                                                                  
  gr_BAL_MC->SetMarkerSize(1.5);                                                                                                                                                  
  gr_BAL_MC->SetMarkerColor(kOrange +1);                                                                                                                                           
  gr_BAL_MC->SetLineColor(kOrange +1);                                                                                                                                             
  gr_BAL_MC->Draw("Psame");                                                                                                                                                       

  gr_MPF_MC->SetMarkerStyle(21);                                                                                                                                                  
  gr_MPF_MC->SetMarkerSize(1.5);                                                                                                                                                  
  gr_MPF_MC->SetMarkerColor(kBlue);                                                                                                                                           
  gr_MPF_MC->SetLineColor(kBlue);                                                                                                                                             
  gr_MPF_MC->Draw("Psame");                                                                                                                                                       
                                                                                                                                                                           
  gr_BAL_MC_Pythia6->SetMarkerStyle(24);                                                                                                                                                
  gr_BAL_MC_Pythia6->SetMarkerSize(1.5);
  gr_BAL_MC_Pythia6->SetMarkerColor(kRed);  
  gr_BAL_MC_Pythia6->SetLineColor(kRed);
  gr_BAL_MC_Pythia6->Draw("Psame");                                                                                                                                                     

  gr_MPF_MC_Pythia6->SetMarkerStyle(20);                                                                                                                                                
  gr_MPF_MC_Pythia6->SetMarkerSize(1.5);
  gr_MPF_MC_Pythia6->SetMarkerColor(kViolet -6);                                                                                                                                         
  gr_MPF_MC_Pythia6->SetLineColor(kViolet -6);                                                                                                                                         
  gr_MPF_MC_Pythia6->Draw("Psame");                                                                                                                                                     
                                                                                                                                                                           
  gPad->RedrawAxis();                                                                                                                                                      
                                                                                                                                                                           
  c1->SaveAs("prova.png");

  }






}






