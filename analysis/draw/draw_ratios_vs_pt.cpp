#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TParameter.h>
#include <TVirtualFitter.h>
#include <TColor.h>
#include <TLatex.h>
#include "drawBase.h"
#include "fitTools.h"
#include "ptBinning.h"
#include "etaBinning.h"
#include "extrapBinning.h"
#include "fineetaBinning.h"

#include <boost/algorithm/string.hpp>

bool ELIF_ = false;
bool FIXM_ = false;
bool OUTPUT_GRAPHS = false;
bool RAW = false;

#define LIGHT_RED TColor::GetColor(0xcf, 0xa0, 0xa1)
#define BALANCING TColor::GetColor(217, 91, 67)
#define MPF TColor::GetColor(192, 41, 66)

void setGraphStyle(TGraphErrors* graph, int markerStyle, int markerColor, int markerSize = 1) {
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetMarkerSize(markerSize);
}

bool SAVE_EPS = false;
bool SAVE_PDF = true;
bool SAVE_PNG = false;
bool SAVE_ROOT = false;
void saveCanvas(TCanvas* canvas, const std::string& name) {
  if (SAVE_EPS)
    canvas->SaveAs((name + ".eps").c_str());
  if (SAVE_PNG)
    canvas->SaveAs((name + ".png").c_str());
  if (SAVE_ROOT)
    canvas->SaveAs((name + ".root").c_str());
  if (SAVE_PDF)
    canvas->SaveAs((name + ".pdf").c_str());
}


void draw_vs_pt_plots(const std::string& resp_reso, const std::string& etaRegion, const std::string& etaRegion_str, const std::string& FIT_RMS, drawBase* db, bool rawJets, const std::string& alphaCut, TFile* outputFile  , double & sf_data, double & sf_mc,  double & sf_dataerr,  double & sf_mcerr, bool isfinebinning, bool isresolution);
void draw_scalefactorVsEta( drawBase* db, double sfMC[], double sfDATA[], double sfMCerr[], double sfDATAerr[], int Npoint, TFile* outputFile);

int main(int argc, char* argv[]) {

  if (argc != 6 && argc != 7) {
    std::cout << "USAGE: ./draw_ratios_vs_pt [DATA_dataset] [mcSignal_dataset] [mcBG_dataset] [recoType] [jetAlgo] [flags=\"\"]" << std::endl;
    exit(53);
  }

  std::string data_dataset(argv[1]);
  std::string mc_photonjet(argv[2]);
  std::string mc_QCD(argv[3]);
  std::string recoType(argv[4]);
  std::string jetAlgo(argv[5]);
  std::string flags = "";
  if (argc == 7) {
    std::string flags_str(argv[6]);
    flags = flags_str;
  }

  jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
  recoType = (recoType == "pf") ? "PFlow" : "PUPPI";
  std::string postFix = recoType + jetAlgo;
  if(recoType == "PFlow") postFix += "chs";

  drawBase* db = new drawBase("PhotonJet", recoType, jetAlgo, OUTPUT_GRAPHS);
  db->set_flags(flags);
  db->set_isCMSArticle(false);

  TString dataFileName;
  if (flags.length() > 0) {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str()/*, flags.c_str()*/);
  } else {
    dataFileName = TString::Format("PhotonJet_%s_%s.root", data_dataset.c_str(), postFix.c_str());
  }

  TFile* dataFile = TFile::Open(dataFileName);
  std::cout << "Opened data file '" << dataFileName << "'." << std::endl;

  db->add_dataFile(dataFile, data_dataset);

  TString mc1FileName;
  if (flags.length() > 0) {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str()/*, flags.c_str()*/);
  } else {
    mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str());
  }
  TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
  std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;

  db->add_mcFile(mcPhotonJetFile, mc_photonjet, "#gamma+jet MC", 46);

  if (mc_QCD != "") {
    TString mc2FileName;
    if (flags.length() > 0) {
      mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc_QCD.c_str(), postFix.c_str()/*, flags.c_str()*/);
    } else {
      mc2FileName = TString::Format("PhotonJet_%s_%s.root", mc_QCD.c_str(), postFix.c_str());
    }
    TFile* mcQCDFile = TFile::Open(mc2FileName);
    std::cout << "Opened mc file '" << mc2FileName << "'." << std::endl;

    if (mc_QCD != mc_photonjet) {
      db->add_mcFile(mcQCDFile, mc_QCD, "QCD MC", 38);
    }
  }

  // MC should already be normalized to a lumi of 1 pb-1
    TParameter<float>* lumi = static_cast<TParameter<float>*>(dataFile->Get("analysis/luminosity"));
    db->set_lumi(lumi->GetVal());
    db->set_lumiNormalization();

  double alpha_cut = static_cast<TParameter<double>*>(dataFile->Get("analysis/alpha_cut"))->GetVal();
  std::stringstream ss;
  ss << ((int) (alpha_cut * 100));
  std::string alphaCut = ss.str();

  std::string fit_rms = "RMS985";//"RMS985";
  std::string outputDir = "PhotonJetPlots_" + db->get_fullSuffix() + "_vs_pt";
  // std::string outputDir = "Plot_vs_pt";
  db->set_outputdir(outputDir);

  

  TFile* output = TFile::Open(std::string(outputDir + "/plots.root").c_str(), "recreate");
  TFile* output_raw = TFile::Open(std::string(outputDir + "/plots_raw.root").c_str(), "recreate");
  
  
  EtaBinning etaBinning;
  size_t s = etaBinning.size();
  double SF_data[s] = {0.};
  double SF_MC[s]= {0.};
  double SF_dataerr[s]= {0.};
  double SF_MCerr[s]= {0.};
  
  double SF_dataspe =  0 ;
  double SF_MCspe =  0 ;
  double SF_dataerrspe =  0 ;
  double SF_MCerrspe =  0 ;
  
  
  
  //scale factor fine eta binning
  
  fineEtaBinning fineetaBinning;
  size_t sfine = fineetaBinning.size();
  double SF_datafine[sfine] = {0.};
  double SF_MCfine[sfine]= {0.};
  double SF_dataerrfine[sfine]= {0.};
  double SF_MCerrfine[sfine]= {0.};
  


  for (size_t i = 0; i < s/*-1*/; i++) { //fixing  // previously excluding the last eta bin ?? why ? 
    std::string etaBin = etaBinning.getBinName(i);
    std::string etaBinTitle = etaBinning.getBinTitle(i);

    draw_vs_pt_plots("response",   etaBin, etaBinTitle, fit_rms, db, false, alphaCut, output, SF_data[i],SF_MC[i], SF_dataerr[i], SF_MCerr[i], false, false); //bool for raw study
    draw_vs_pt_plots("resolution", etaBin, etaBinTitle, fit_rms, db, false, alphaCut, output, SF_data[i],SF_MC[i], SF_dataerr[i], SF_MCerr[i], false,true);
  }
  
  
  

  //special case
  std::string etaBinTitle = "|#eta| #leq 1.3";
  draw_vs_pt_plots("response",   "eta0013", etaBinTitle, fit_rms, db, false, alphaCut, output, SF_dataspe, SF_MCspe, SF_dataerrspe, SF_MCerrspe, false, false);
  draw_vs_pt_plots("resolution", "eta0013", etaBinTitle, fit_rms, db, false, alphaCut, output, SF_dataspe, SF_MCspe, SF_dataerrspe, SF_MCerrspe, false, true);
  draw_scalefactorVsEta(  db, SF_MC, SF_data, SF_MCerr, SF_dataerr, s,output);
  
  
  for (size_t i = 0; i < sfine; i++) { //fixing 
    std::string etaBin = fineetaBinning.getBinName(i);
    std::string etaBinTitle = fineetaBinning.getBinTitle(i);

    draw_vs_pt_plots("response",   etaBin, etaBinTitle, fit_rms, db, false, alphaCut, output, SF_datafine[i],SF_MCfine[i], SF_dataerrfine[i], SF_MCerrfine[i],true, false); //bool for raw study
    draw_vs_pt_plots("resolution", etaBin, etaBinTitle, fit_rms, db, false, alphaCut, output, SF_datafine[i],SF_MCfine[i], SF_dataerrfine[i], SF_MCerrfine[i],true, true);
  }
    draw_scalefactorVsEta(  db, SF_MCfine, SF_datafine, SF_MCerrfine, SF_dataerrfine, sfine,output);
  
  output->Close();
  output_raw->Close();

  delete output;
  delete output_raw;

  return 0;
}

void drawGraphs(TGraphErrors* data, TGraphErrors* mc, double xMin, double xMax, const std::string& methodName, drawBase* db, const std::string& etaRegion, const std::string& legendTitle, const std::string& prefix, const std::string& suffix, bool rawJets) {

  std::string name_base = db->get_outputdir();
  name_base += "/" + prefix;

  if (rawJets)
    name_base += "_RAW";

  if (etaRegion != "")
    name_base += "_" + etaRegion;

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 800);
  c1->cd();

  // Data / MC comparison
  TPad* pad_hi = new TPad("pad_hi", "", 0., 0.33, 0.99, 0.99);
  pad_hi->Draw();
  pad_hi->SetLogx();
  pad_hi->SetLeftMargin(0.12);
  pad_hi->SetBottomMargin(0.015);

  // Data/MC ratio
  TPad* pad_lo = new TPad("pad_lo", "", 0., 0., 0.99, 0.33);
  pad_lo->Draw();
  pad_lo->SetLogx();
  pad_lo->SetLeftMargin(0.12);
  pad_lo->SetTopMargin(1.);
  pad_lo->SetBottomMargin(0.3);

  pad_lo->cd();

  bool isResponse = prefix == "response";

  double lowPadY1 = isResponse ? 0.86 : 0.6;
  double lowPadY2 = isResponse ? 1.14 : 1.4;

  TH2D* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 100, xMin, xMax, 100, lowPadY1, lowPadY2);

  h2_axes_lo_resp->SetXTitle("Photon p_{T} [GeV/c]");
  h2_axes_lo_resp->SetYTitle("Data / MC");
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.2);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.70);
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();
  h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.085);
  h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.09);
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
  h2_axes_lo_resp->GetYaxis()->SetNdivisions(7,true);
  h2_axes_lo_resp->GetYaxis()->SetRangeUser(0.80,1.20);
  h2_axes_lo_resp->Draw();

  float factor = isResponse ? 0.05 : 0.2;

  TLine* line_one = new TLine(xMin, 1., xMax, 1.);
  TLine* line_plus_resp = new TLine(xMin, 1 + factor, xMax, 1 + factor);
  TLine* line_minus_resp = new TLine(xMin, 1 - factor, xMax, 1 - factor);

  line_plus_resp->SetLineWidth(2);
  line_plus_resp->SetLineStyle(2);

  line_minus_resp->SetLineWidth(2);
  line_minus_resp->SetLineStyle(2);

  TGraphErrors* gr_resp_ratio = db->get_graphRatio(data, mc);
  gr_resp_ratio->SetName("response_ratio");
  gr_resp_ratio->SetMarkerStyle(20);
  gr_resp_ratio->SetMarkerSize(1.0);
  gr_resp_ratio->SetMarkerColor(BALANCING);
  gr_resp_ratio->SetLineColor(BALANCING);

  TF1* ratioFit = new TF1("ratioFit", "pol0", xMin, xMax);
  ratioFit->SetParameter(0, 1.);
  ratioFit->SetLineColor(TColor::GetColor("#C02942"));
  ratioFit->SetLineWidth(1.0);
  gr_resp_ratio->Fit(ratioFit, "RQNF EX0");

  TH1D* errors = new TH1D("errors", "errors", xMax - xMin, xMin, xMax);
  //std::cout<<" value of the fit for debugging : "<<ratioFit->GetParameter(0)<<" error : "<<ratioFit->GetParError(0)<<std::endl;
  if(ratioFit->GetParError(0) != 0. ) (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors, 0.68);
  errors->SetStats(false);
  //errors->SetFillColor(TColor::GetColor("#556270"));
  errors->SetFillColor(TColor::GetColor("#ECD078"));
  errors->SetFillStyle(1001);

  double fitValue = ratioFit->GetParameter(0);
  double fitError = ratioFit->GetParError(0);

  //  TString str = TString::Format("Fit: [0] =\\num{%.4f \\pm %.4f} & \\num{%.4f \\pm %.4f}", fitValue, fitError, 1. / fitValue, 1. / fitValue * fitError / fitValue);
  TString str = TString::Format("Fit: [0] = %.4f #pm %.4f", fitValue, fitError);
  std::cout << str << std::endl;

  float height = 0.81 - 0.77;
  float labelYPos = isResponse ? 0.77 : 0.37;

  TPaveText* fitlabel = new TPaveText(0.43, labelYPos, 0.78, labelYPos + height, "brNDC");
  fitlabel->SetTextSize(0.08);
  fitlabel->SetFillColor(0);
  fitlabel->SetTextFont(42);
  //TString fitLabelText = TString::Format("#font[42]{Fit: %.4f #pm %.4f + (%.2e #pm %.2e)x}", fitValue, fitError, ratioFit->GetParameter(1), ratioFit->GetParError(1));
  TString fitLabelText = TString::Format("Fit: %.4f #pm %.4f", fitValue, fitError);
  fitlabel->AddText(fitLabelText);

  gr_resp_ratio->Draw("P same");
  errors->Draw("e3 same");
  line_one->Draw("same");
  line_plus_resp->Draw("same");
  line_minus_resp->Draw("same");
  ratioFit->Draw("same");
  fitlabel->Draw("same");
  gr_resp_ratio->Draw("P same");

  gPad->RedrawAxis();
  pad_hi->cd();

  std::string yTitle;
  if (prefix == "response") {
    yTitle = "Jet p_{T} response";
  } else {
    yTitle = "Jet p_{T} resolution";
  }

  if (rawJets) {
    yTitle += " (raw jets)";
  }

  float ymin, ymax;
  if (prefix == "response") {
    ymin = 0.7;
    ymax = 1.3;
  } else {
    ymin = 0.;
    ymax = 0.3;
  }

  TH2D* h2_axes = new TH2D("axes_again", "", 10, xMin, xMax, 10, ymin, ymax);
  h2_axes->SetYTitle(yTitle.c_str());
  h2_axes->GetXaxis()->SetTitleOffset(1.1);
  h2_axes->GetYaxis()->SetTitleOffset(1.2);
  h2_axes->GetYaxis()->SetTitleSize(0.045);
  h2_axes->GetXaxis()->SetLabelSize(0.);

  h2_axes->Draw();

  //float labelTextSize = 0.035;
  TPaveText* label_algo = db->get_labelAlgo(2);

  float legendY1 = 0.15;
  float legendY2 = 0.38;

  if (prefix == "resolution") {
    legendY1 = 0.60;
    legendY2 = 0.83;
  }

  TLegend* legend = new TLegend(0.50, legendY1, 0.90, legendY2, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);

  TString label = TString::Format("%s (data)", methodName.c_str());
  legend->AddEntry(data, label, "P");

  label = TString::Format("%s (MC)", methodName.c_str());
  legend->AddEntry(mc, label, "P");
  legend->Draw("same");

  float cmsTextSize = 0.043;
  TPaveText* label_cms = db->get_labelCMS(1);
  label_cms->SetTextSize(cmsTextSize);

  TPaveText* label_sqrt = db->get_labelSqrt(1);

  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");

  mc->SetMarkerStyle(24);
  mc->SetMarkerSize(1.0);
  mc->SetMarkerColor(MPF);
  mc->SetLineColor(MPF);
  mc->Draw("Psame");

  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.0);
  data->SetMarkerColor(MPF);
  data->SetLineColor(MPF);

  data->Draw("Psame");

  gPad->RedrawAxis();

  name_base += "_" + suffix;
  saveCanvas(c1, name_base);

  delete c1;
  delete h2_axes_lo_resp;
  delete line_one;
  delete line_plus_resp;
  delete line_minus_resp;
  delete gr_resp_ratio;
  delete ratioFit;
  delete errors;
  delete fitlabel;
  delete h2_axes;
  delete label_algo;
  delete legend;
  delete label_cms;
  delete label_sqrt;
}




void drawGraphsratio(TGraphErrors* ratio, double xMin, double xMax, const std::string& methodName, drawBase* db, const std::string& etaRegion, const std::string& legendTitle, const std::string& prefix, const std::string& suffix, bool is_reso) {

  std::string name_base = db->get_outputdir();
  name_base += "/" + prefix;


  if (etaRegion != "")
    name_base += "_" + etaRegion;

   TCanvas* c1;
  if(is_reso)c1 = new TCanvas("c1", "c1", 600, 600);
  if(!is_reso)c1 = new TCanvas("c1", "c1", 600, 800);
  c1->cd();

  // Data / MC comparison
  TPad* pad_hi ;
  
  if(is_reso) pad_hi = new TPad("pad_hi", "", 0., 0., 0.99, 0.99);
  if(!is_reso)pad_hi = new TPad("pad_hi", "", 0., 0.33, 0.99, 0.99);
  pad_hi->Draw();
  pad_hi->SetLogx();
  pad_hi->SetLeftMargin(0.12);
  pad_hi->cd();
  bool isResponse = prefix == "response";

  double lowPadY1 = isResponse ? 0.86 : 0.6;
  double lowPadY2 = isResponse ? 1.14 : 1.4;

  
  

  float factor = isResponse ? 0.05 : 0.2;


   
  ratio->SetName("response_ratio");
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(1.0);
  ratio->SetMarkerColor(BALANCING);
  ratio->SetLineColor(BALANCING);
  std::cout<<" N point for JER SF "<<ratio->GetN()<<std::endl;
  TF1* ratioFit = new TF1("ratioFit", "pol0", xMin, xMax);
  ratioFit->SetParameter(0, 1.);
  ratioFit->SetLineColor(TColor::GetColor("#C02942"));
  ratioFit->SetLineWidth(1.0);
  ratio->Fit(ratioFit, "RQNF EX0");

  TH1D* errors = new TH1D("errors", "errors", xMax - xMin, xMin, xMax);
  //std::cout<<" value of the fit for debugging : "<<ratioFit->GetParameter(0)<<" error : "<<ratioFit->GetParError(0)<<std::endl;
  if(ratioFit->GetParError(0) != 0. ) (TVirtualFitter::GetFitter())->GetConfidenceIntervals(errors, 0.68);
  errors->SetStats(false);
  //errors->SetFillColor(TColor::GetColor("#556270"));
  errors->SetFillColor(TColor::GetColor("#ECD078"));
  errors->SetFillStyle(1001);

  double fitValue = ratioFit->GetParameter(0);
  double fitError = ratioFit->GetParError(0);

  //  TString str = TString::Format("Fit: [0] =\\num{%.4f \\pm %.4f} & \\num{%.4f \\pm %.4f}", fitValue, fitError, 1. / fitValue, 1. / fitValue * fitError / fitValue);
  TString str = TString::Format("Fit: [0] = %.4f #pm %.4f", fitValue, fitError);
  std::cout << str << std::endl;

  float height = 0.81 - 0.77;
  float labelYPos = isResponse ? 0.77 : 0.37;

  TPaveText* fitlabel = new TPaveText(0.43, labelYPos, 0.78, labelYPos + height, "brNDC");
  fitlabel->SetTextSize(0.08);
  fitlabel->SetFillColor(0);
  fitlabel->SetTextFont(42);
  //TString fitLabelText = TString::Format("#font[42]{Fit: %.4f #pm %.4f + (%.2e #pm %.2e)x}", fitValue, fitError, ratioFit->GetParameter(1), ratioFit->GetParError(1));
  TString fitLabelText = TString::Format("Fit: %.4f #pm %.4f", fitValue, fitError);
  fitlabel->AddText(fitLabelText);
  
  
  TH2D* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 100, xMin, xMax, 100, 0., 2.);

  h2_axes_lo_resp->SetXTitle("Photon p_{T} [GeV/c]");
  h2_axes_lo_resp->SetYTitle("Data / MC ");
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.8);
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
 // h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.08);
 // h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.05);
//  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
  h2_axes_lo_resp->GetYaxis()->SetNdivisions(13,true);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(1.2);
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.045);
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();  

  
  h2_axes_lo_resp->Draw();
  ratio->Draw("P same");
  errors->Draw("e3 same");
  ratioFit->Draw("same");
  fitlabel->Draw("same");
  ratio->Draw("P same");

  gPad->RedrawAxis();

  

 

  //float labelTextSize = 0.035;
  TPaveText* label_algo = db->get_labelAlgo(2);

  float legendY1 = 0.15;
  float legendY2 = 0.38;

  if (prefix == "resolution") {
    legendY1 = 0.60;
    legendY2 = 0.83;
  }

  TLegend* legend = new TLegend(0.50, legendY1, 0.90, legendY2, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);

  TString label = TString::Format("%s (data)", methodName.c_str());

  label = TString::Format("%s (MC)", methodName.c_str());

  float cmsTextSize = 0.043;
  TPaveText* label_cms = db->get_labelCMS(1);
  label_cms->SetTextSize(cmsTextSize);

  TPaveText* label_sqrt = db->get_labelSqrt(1);
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");


  

  gPad->RedrawAxis();

  name_base += "_" + suffix;
  name_base += "JERSF_vspt_";
  saveCanvas(c1, name_base);

  delete c1;
  delete h2_axes_lo_resp;
  delete ratioFit;
  delete errors;
  delete fitlabel;
  delete label_algo;
  delete legend;
  delete label_cms;
  delete label_sqrt;
}


void drawGraphsextrap(TGraphErrors* ratio, double xMin, double xMax, const std::string& methodName, drawBase* db, const std::string& etaRegion, const std::string& legendTitle, const std::string& prefix, const std::string& suffix) {

  std::string name_base = db->get_outputdir();
  name_base += "/" + prefix;


  if (etaRegion != "")
    name_base += "_" + etaRegion;

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  // Data / MC comparison
  TPad* pad_hi = new TPad("pad_hi", "", 0., 0., 0.99, 0.99);
  pad_hi->Draw();
  pad_hi->SetLogx();
  pad_hi->SetLeftMargin(0.12);
  pad_hi->cd();
  bool isResponse = prefix == "response";

  double lowPadY1 = isResponse ? 0.86 : 0.6;
  double lowPadY2 = isResponse ? 1.14 : 1.4;

  
  

  float factor = isResponse ? 0.05 : 0.2;


   
  ratio->SetName("response_ratio");
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(1.0);
  ratio->SetMarkerColor(BALANCING);
  ratio->SetLineColor(BALANCING);
  std::cout<<" N point for JER SF "<<ratio->GetN()<<std::endl;
  

  

  

  //  TString str = TString::Format("Fit: [0] =\\num{%.4f \\pm %.4f} & \\num{%.4f \\pm %.4f}", fitValue, fitError, 1. / fitValue, 1. / fitValue * fitError / fitValue);



  float height = 0.81 - 0.77;
  float labelYPos = isResponse ? 0.77 : 0.37;

  
  
  
  TH2D* h2_axes_lo_resp = new TH2D("axes_lo_resp", "", 100, xMin, xMax, 100, 0., 2.);

  h2_axes_lo_resp->SetXTitle("Photon p_{T} [GeV/c]");
  h2_axes_lo_resp->SetYTitle("Kfsr ");
  h2_axes_lo_resp->GetXaxis()->SetTitleOffset(1.);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(0.8);
  h2_axes_lo_resp->GetXaxis()->SetTickLength(0.06);
  h2_axes_lo_resp->GetXaxis()->SetMoreLogLabels();
 // h2_axes_lo_resp->GetXaxis()->SetLabelSize(0.08);
 // h2_axes_lo_resp->GetYaxis()->SetLabelSize(0.07);
  h2_axes_lo_resp->GetXaxis()->SetTitleSize(0.05);
//  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.08);
  h2_axes_lo_resp->GetYaxis()->SetNdivisions(13,true);
  h2_axes_lo_resp->GetYaxis()->SetTitleOffset(1.2);
  h2_axes_lo_resp->GetYaxis()->SetTitleSize(0.045);
  h2_axes_lo_resp->GetXaxis()->SetNoExponent();  

  
  h2_axes_lo_resp->Draw();
  ratio->Draw("P same");

  ratio->Draw("P same");

  gPad->RedrawAxis();

  

 

  //float labelTextSize = 0.035;
  TPaveText* label_algo = db->get_labelAlgo(2);

  float legendY1 = 0.15;
  float legendY2 = 0.38;

  if (prefix == "resolution") {
    legendY1 = 0.60;
    legendY2 = 0.83;
  }

  TLegend* legend = new TLegend(0.50, legendY1, 0.90, legendY2, legendTitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->SetTextFont(42);

  TString label = TString::Format("%s (data)", methodName.c_str());

  label = TString::Format("%s (MC)", methodName.c_str());

  float cmsTextSize = 0.043;
  TPaveText* label_cms = db->get_labelCMS(1);
  label_cms->SetTextSize(cmsTextSize);

  TPaveText* label_sqrt = db->get_labelSqrt(1);
  legend->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  label_algo->Draw("same");


  

  gPad->RedrawAxis();

  name_base += "_" + suffix;
  name_base += "KFSR_vspt";
  saveCanvas(c1, name_base);

  delete c1;
  delete h2_axes_lo_resp;


  delete label_algo;
  delete legend;
  delete label_cms;
  delete label_sqrt;
}




void draw_vs_pt_plots(const std::string& resp_reso, const std::string& etaRegion, const std::string& etaRegion_str, const std::string& FIT_RMS, drawBase* db, bool rawJets, const std::string& alphaCut, TFile* outputFile, double & sf_data, double & sf_mc,  double & sf_dataerr,  double & sf_mcerr, bool isfinebinning, bool isresolution) {

  //  std::string output = db->get_outputdir();
  //  TFile * outputFile = TFile::Open(TString::Format("%s/plots.root", output.c_str()).Data(), "update");
  
  std::string fullEtaRegion;
  if (etaRegion == "eta0013") fullEtaRegion = "eta00_13";
  else if (etaRegion == "eta0008") fullEtaRegion = "eta00_08";
  else if (etaRegion == "eta0813") fullEtaRegion = "eta08_13";
  else if (etaRegion == "eta1319") fullEtaRegion = "eta13_19";
  else if (etaRegion == "eta1925") fullEtaRegion = "eta19_25";
  else if (etaRegion == "eta2530") fullEtaRegion = "eta25_30";
  else if (etaRegion == "eta3032") fullEtaRegion = "eta30_32";
  else if (etaRegion == "eta3252") fullEtaRegion = "eta32_52";
  else if (etaRegion == "eta3252") fullEtaRegion = "eta32_52";
  
  else if (etaRegion == "eta0003") fullEtaRegion = "eta00_03";
  else if (etaRegion == "eta0305") fullEtaRegion = "eta03_05";
  else if (etaRegion == "eta0508") fullEtaRegion = "eta05_08";
  else if (etaRegion == "eta0810") fullEtaRegion = "eta08_10";
  else if (etaRegion == "eta1013") fullEtaRegion = "eta10_13";
  else if (etaRegion == "eta1315") fullEtaRegion = "eta13_15";
  else if (etaRegion == "eta1517") fullEtaRegion = "eta15_17";
  else if (etaRegion == "eta1719") fullEtaRegion = "eta17_19";
  else if (etaRegion == "eta1922") fullEtaRegion = "eta19_22";
  else if (etaRegion == "eta2223") fullEtaRegion = "eta22_23";
  else if (etaRegion == "eta2325") fullEtaRegion = "eta23_25";
  else if (etaRegion == "eta2526") fullEtaRegion = "eta25_26";
  else if (etaRegion == "eta2629") fullEtaRegion = "eta26_29";
  else if (etaRegion == "eta2930") fullEtaRegion = "eta29_30";
  else if (etaRegion == "eta3031") fullEtaRegion = "eta30_31";
  else if (etaRegion == "eta3135") fullEtaRegion = "eta31_35";
  else if (etaRegion == "eta3538") fullEtaRegion = "eta35_38";
  else if (etaRegion == "eta3852") fullEtaRegion = "eta38_52";
  
  else if (etaRegion == "eta0005") fullEtaRegion = "eta00_05";
  else if (etaRegion == "eta0508") fullEtaRegion = "eta05_08";
  else if (etaRegion == "eta0811") fullEtaRegion = "eta08_11";
  else if (etaRegion == "eta1113") fullEtaRegion = "eta11_13";
  else if (etaRegion == "eta1317") fullEtaRegion = "eta13_17";
  else if (etaRegion == "eta1719") fullEtaRegion = "eta17_19";
  else if (etaRegion == "eta1920") fullEtaRegion = "eta19_20";
  else if (etaRegion == "eta2023") fullEtaRegion = "eta20_23";
  else if (etaRegion == "eta2325") fullEtaRegion = "eta23_25";
  else if (etaRegion == "eta2529") fullEtaRegion = "eta25_29";
  else if (etaRegion == "eta2930") fullEtaRegion = "eta29_30";
  else if (etaRegion == "eta3031") fullEtaRegion = "eta30_31";
  else if (etaRegion == "eta3152") fullEtaRegion = "eta31_52";

  else if (etaRegion == "eta2951") fullEtaRegion = "eta29_51";  
  else if (etaRegion == "eta3051") fullEtaRegion = "eta30_51";
  else if (etaRegion == "eta2329") fullEtaRegion = "eta23_29";
  else if (etaRegion == "eta1720") fullEtaRegion = "eta17_20";  
  else fullEtaRegion = "eta_unknown";

  if (resp_reso != "response" && resp_reso != "resolution") {
    std::cout << "Only 'Response' and 'Resolution' supported. Exiting." << std::endl;
    exit(776);
  }

  
  
  
  std::string file_noextrap_name = "PhotonJetGraphs_" + db->get_fullSuffix() + ".root";
  std::string file_extrap_name = "PhotonJetExtrapGraphs_" + db->get_fullSuffix();
  
  if(isfinebinning){
  
    //file_noextrap_name = "PhotonJetGraphs_fineetabin_" + db->get_fullSuffix() + ".root";
    file_extrap_name = "PhotonJetExtrapGraphs_fineetabin_" + db->get_fullSuffix();
  
  }
  if (etaRegion != "") file_extrap_name += "_" + etaRegion;
  if (rawJets) {
    file_extrap_name += "RAW";
  }
  file_extrap_name += "_" + FIT_RMS + ".root";

  //now open graph files and plot them on on top of the other:
  TFile* file_noextrap = TFile::Open(file_noextrap_name.c_str(), "read");
  if (file_noextrap != 0) std::cout << "-> Opened file: '" << file_noextrap_name << "'" << std::endl;
  else {
    std::cout << "Didn't find file '" << file_noextrap_name << "'. Skipping." << std::endl;
    return;
  }
  TFile* file_extrap = TFile::Open(file_extrap_name.c_str(), "read");
  if (file_extrap != 0) std::cout << "-> Opened file: '" << file_extrap_name << "'" << std::endl;
  else {
    std::cout << "Didn't find file '" << file_extrap_name << "'. Skipping." << std::endl;
    return;
  }

  PtBinning ptBinning;

  std::vector<std::pair<float, float> > ptPhot_binning = ptBinning.getBinning();
  float xMin = ptPhot_binning[0].first;
  float xMax = ptPhot_binning[ptPhot_binning.size() - 1].second;

  int markerSize = 1.25;

  std::string intrName = "gr_intr";
  if (resp_reso == "response"){
    intrName += "Resp";
    sf_data = 1;
    
    sf_mc = 1;
    sf_dataerr = 1;
    sf_mcerr = 1;
  }else
    intrName += "Reso";

  intrName += "_vs_pt";

  std::string prefix = (resp_reso == "response") ? "resp" : "resolution";

  std::string responseBALANCING_name      = prefix + "_balancing";
  std::string responseBALANCINGcount_name = prefix + "_balancing";
  
  if(isfinebinning){
   responseBALANCING_name +="_fine_bining";
   responseBALANCINGcount_name+="_fine_bining";
   }
  
  if (rawJets) {
    responseBALANCING_name += "_raw";
  }
  responseBALANCING_name += "_" + etaRegion + "_data_vs_pt";
  responseBALANCINGcount_name+= "_" + etaRegion + "_RawNEvents_data_vs_pt";
  
  TH1F* H_count_vs_pt   = NULL;
  TH1F* H_countMC_vs_pt = NULL;
  
  std::cout << responseBALANCING_name<<" " << responseBALANCINGcount_name << std::endl;
  if(resp_reso == "response"){
  H_count_vs_pt =(TH1F*)file_noextrap->Get(responseBALANCINGcount_name.c_str());
  H_count_vs_pt->SetName(TString::Format("%s_PtBalchs_DATA_a%s_%s_RawNEvents_data_vs_pt", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  }

  TGraphErrors* gr_responseBALANCING_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCING_name.c_str());
  gr_responseBALANCING_vs_pt->SetMarkerStyle(21);
  gr_responseBALANCING_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCING_vs_pt->SetMarkerColor(kGray + 2);
  gr_responseBALANCING_vs_pt->SetName(TString::Format("%s_PtBalchs_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  std::string responseBALANCINGMC_name = responseBALANCING_name;
  std::string responseBALANCINGcountMC_name = responseBALANCINGcount_name;
  boost::replace_all(responseBALANCINGMC_name, "data", "mc");
  boost::replace_all(responseBALANCINGcountMC_name, "data", "mc");
  
  if(resp_reso == "response"){
  H_countMC_vs_pt =(TH1F*)file_noextrap->Get(responseBALANCINGcountMC_name.c_str());
  H_countMC_vs_pt->SetName(TString::Format("%s_PtBalchs_MC_a%s_%s_RawNEvents_data_vs_pt", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  }

  TGraphErrors* gr_responseBALANCINGMC_vs_pt = (TGraphErrors*)file_noextrap->Get(responseBALANCINGMC_name.c_str());
  gr_responseBALANCINGMC_vs_pt->SetMarkerStyle(25);
  gr_responseBALANCINGMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseBALANCINGMC_vs_pt->SetMarkerColor(kGray + 2);
  gr_responseBALANCINGMC_vs_pt->SetName(TString::Format("%s_PtBalchs_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  TGraphErrors* gr_responseMPF_vs_pt = NULL;
  TGraphErrors* gr_responseMPFMC_vs_pt = NULL;

  std::string responseMPF_name = prefix + "_mpf";
  if(isfinebinning) responseMPF_name +="_fine_bining";
  if (rawJets) {
    responseMPF_name += "_raw";
  }
  responseMPF_name += "_" + etaRegion + "_data_vs_pt";
  
  boost::replace_all(responseBALANCING_name, "balancing", "mpf");
  boost::replace_all(responseBALANCINGMC_name, "balancing", "mpf");
  
  TH1F* H_count_vs_pt_mpf   = NULL;
  TH1F* H_countMC_vs_pt_mpf = NULL;
  
  if(resp_reso == "response"){
  H_count_vs_pt_mpf =(TH1F*)file_noextrap->Get(responseBALANCINGcount_name.c_str());
  H_count_vs_pt_mpf->SetName(TString::Format("%s_MPFchs_DATA_a%s_%s_RawNEvents_data_vs_pt", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  
  H_countMC_vs_pt_mpf =(TH1F*)file_noextrap->Get(responseBALANCINGcountMC_name.c_str());
  H_countMC_vs_pt_mpf->SetName(TString::Format("%s_MPFchs_MC_a%s_%s_RawNEvents_data_vs_pt", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  }

  gr_responseMPF_vs_pt = (TGraphErrors*)file_noextrap->Get(responseMPF_name.c_str());
  gr_responseMPF_vs_pt->SetMarkerStyle(20);
  gr_responseMPF_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPF_vs_pt->SetMarkerColor(38);
  gr_responseMPF_vs_pt->SetName(TString::Format("%s_MPFchs_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  std::string responseMPFMC_name = responseMPF_name;
  boost::replace_all(responseMPFMC_name, "data", "mc");

  gr_responseMPFMC_vs_pt = (TGraphErrors*)file_noextrap->Get(responseMPFMC_name.c_str());
  gr_responseMPFMC_vs_pt->SetMarkerStyle(24);
  gr_responseMPFMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPFMC_vs_pt->SetMarkerColor(38);
  gr_responseMPFMC_vs_pt->SetName(TString::Format("%s_MPFchs_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  std::string resp_reso_short = (resp_reso == "response") ? "Resp" : "Reso";

  std::string responseEXTRAP_name = "gr_DATA" + resp_reso_short + "_vs_pt";
  
  
  
  TGraphErrors* gr_responseEXTRAP_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAP_name.c_str());
  gr_responseEXTRAP_vs_pt->SetMarkerStyle(22);
  gr_responseEXTRAP_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAP_vs_pt->SetMarkerColor(46);
  gr_responseEXTRAP_vs_pt->SetName(TString::Format("%s_PtBalchs_extrap_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  //
  //  gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  //  gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove second point also

  std::string responseEXTRAPMC_name = "gr_extrap" + resp_reso_short + "_vs_pt";
    
  TGraphErrors* gr_responseKFSR_vs_pt = new TGraphErrors(0);
  if(isfinebinning && isresolution){
  gr_responseKFSR_vs_pt = (TGraphErrors*)file_extrap->Get("gr_DATAkfr_vs_pt");
  gr_responseKFSR_vs_pt->SetMarkerStyle(22);
  gr_responseKFSR_vs_pt->SetMarkerSize(markerSize);
  gr_responseKFSR_vs_pt->SetMarkerColor(46);
  gr_responseKFSR_vs_pt->SetName(TString::Format("%s_PtBalchs_KFSR_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  gr_responseKFSR_vs_pt->SetTitle(TString::Format("%s_PtBalchs_KFSR_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  }
  
  TGraphErrors* gr_MC_responseKFSR_vs_pt = new TGraphErrors(0);
  if(isfinebinning && isresolution){
  gr_MC_responseKFSR_vs_pt = (TGraphErrors*)file_extrap->Get("gr_extrapkfr_vs_pt");
  gr_MC_responseKFSR_vs_pt->SetMarkerStyle(22);
  gr_MC_responseKFSR_vs_pt->SetMarkerSize(markerSize);
  gr_MC_responseKFSR_vs_pt->SetMarkerColor(46);
  gr_MC_responseKFSR_vs_pt->SetName(TString::Format("%s_PtBalchs_KFSR_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  gr_MC_responseKFSR_vs_pt->SetTitle(TString::Format("%s_PtBalchs_KFSR_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  }
  
  TGraphErrors* gr_EXTRAP_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAP_name.c_str());
  gr_EXTRAP_vs_pt->SetMarkerStyle(22);
  gr_EXTRAP_vs_pt->SetMarkerSize(markerSize);
  gr_EXTRAP_vs_pt->SetMarkerColor(46);
  gr_EXTRAP_vs_pt->SetName(TString::Format("%s_PtBalchs_KFSR_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  
  //
  //  gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  //  gr_responseEXTRAP_vs_pt->RemovePoint(0); //remove second point also

  //std::string responseEXTRAPMC_name = "gr_extrap" + resp_reso_short + "_vs_pt";

   /* TF1* reso_mc_Fit = new TF1("reso_mc", "[0]", 40., 3000.);
    reso_mc_Fit->SetParameter(1,0.);
    reso_mc_Fit->SetLineColor(46);
    reso_mc_Fit->SetLineWidth(2);*/


  TGraphErrors* gr_responseEXTRAPMC_vs_pt = (TGraphErrors*)file_extrap->Get(responseEXTRAPMC_name.c_str());
  gr_responseEXTRAPMC_vs_pt->SetMarkerStyle(26);
  gr_responseEXTRAPMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseEXTRAPMC_vs_pt->SetMarkerColor(46);
  gr_responseEXTRAPMC_vs_pt->SetName(TString::Format("%s_PtBalchs_extrap_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
 // if (resp_reso != "response"){gr_responseEXTRAPMC_vs_pt->Fit(reso_mc_Fit,"RQ");}
  
    
   
  //  gr_responseEXTRAPMC_vs_pt->RemovePoint(0); //remove first point (cant extrapolate at such low pt)
  //  gr_responseEXTRAPMC_vs_pt->RemovePoint(0); //remove second point also

  TGraphErrors* gr_responseMPFExtrap_vs_pt = (TGraphErrors*) file_extrap->Get(std::string("gr_DATA" + resp_reso_short + "MPF_vs_pt").c_str());
  gr_responseMPFExtrap_vs_pt->SetMarkerStyle(23);
  gr_responseMPFExtrap_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPFExtrap_vs_pt->SetMarkerColor(49);
  gr_responseMPFExtrap_vs_pt->SetName(TString::Format("%s_MPFchs_extrap_DATA_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  //  gr_responseMPFExtrap_vs_pt->RemovePoint(0);
  //  gr_responseMPFExtrap_vs_pt->RemovePoint(0);

  TGraphErrors* gr_responseMPFExtrapMC_vs_pt = (TGraphErrors*) file_extrap->Get(std::string("gr_extrap" + resp_reso_short + "MPF_vs_pt").c_str());
  gr_responseMPFExtrapMC_vs_pt->SetMarkerStyle(32);
  gr_responseMPFExtrapMC_vs_pt->SetMarkerSize(markerSize);
  gr_responseMPFExtrapMC_vs_pt->SetMarkerColor(49);
  gr_responseMPFExtrapMC_vs_pt->SetName(TString::Format("%s_MPFchs_extrap_MC_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  //  gr_responseMPFExtrapMC_vs_pt->RemovePoint(0);
  //  gr_responseMPFExtrapMC_vs_pt->RemovePoint(0);

  float ymin, ymax;
  if (resp_reso == "response") {
    ymin =  0.7;
    ymax = 1.4;
  } else {
    ymin = 0.;
    ymax = 0.3;
  }

  std::string plotVarName = (resp_reso == "response") ? "Response" : "Resolution";
  std::string legendTrue = "True " + plotVarName;

  TH2D* axes = new TH2D("axes", "", 10, xMin, xMax, 10, ymin, ymax);
  axes->SetXTitle("Photon p_{T} [GeV]");
  std::string yTitle = plotVarName;
  if (resp_reso == "response") {
    yTitle = "Jet p_{T} response";
  } else {
    yTitle = "Jet p_{T} resolution";
  }

  if (rawJets) {
    yTitle += " (raw jets)";
  }

  axes->SetYTitle(yTitle.c_str());
  axes->GetXaxis()->SetMoreLogLabels();
  axes->GetXaxis()->SetNoExponent();

  std::string legendTitle = "  " + etaRegion_str;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  // new stuff -- dataMC ratio
  TGraphErrors* gr_dataMC_BALANCING = fitTools::get_graphRatio(gr_responseBALANCING_vs_pt, gr_responseBALANCINGMC_vs_pt);
  gr_dataMC_BALANCING->SetMarkerStyle(21);
  gr_dataMC_BALANCING->SetMarkerSize(markerSize);
  gr_dataMC_BALANCING->SetMarkerColor(kGray + 2);
  gr_dataMC_BALANCING->SetName(TString::Format("%s_PtBalchs_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));

  TGraphErrors* gr_dataMC_MPF = fitTools::get_graphRatio(gr_responseMPF_vs_pt, gr_responseMPFMC_vs_pt);
  gr_dataMC_MPF->SetMarkerStyle(20);
  gr_dataMC_MPF->SetMarkerSize(markerSize);
  gr_dataMC_MPF->SetMarkerColor(38);
  gr_dataMC_MPF->SetName(TString::Format("%s_MPFchs_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  
  
  TF1* reso_data_Fit = new TF1("reso", "[0]", 40., 3000.);
  reso_data_Fit->SetParameter(1,0.);
  reso_data_Fit->SetLineColor(46);
  reso_data_Fit->SetLineWidth(2);

 
  TGraphErrors* gr_dataMC_EXTRAP = fitTools::get_graphRatio(gr_responseEXTRAP_vs_pt, gr_responseEXTRAPMC_vs_pt);
  gr_dataMC_EXTRAP->SetMarkerStyle(22);
  gr_dataMC_EXTRAP->SetMarkerSize(markerSize);
  gr_dataMC_EXTRAP->SetMarkerColor(46);
  gr_dataMC_EXTRAP->SetName(TString::Format("%s_PtBalchs_extrap_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  if (resp_reso != "response"){ gr_dataMC_EXTRAP->Fit(reso_data_Fit,"RQ");}
  
  // JER (scale factor)
  if (resp_reso != "response"){ sf_data = reso_data_Fit->GetParameter(0);
   sf_mc   = 1;
   sf_dataerr = reso_data_Fit->GetParError(0);
   sf_mcerr   = 1;}
   
   
   std::cout<<" SF  : " <<sf_data<<" ERRORS : "<< sf_dataerr<<std::endl;
   //std::cout<<" SF MC : " <<sf_mc<<" ERRORS : "<< sf_mcerr<<std::endl; 

  TGraphErrors* gr_dataMC_MPF_EXTRAP = fitTools::get_graphRatio(gr_responseMPFExtrap_vs_pt, gr_responseMPFExtrapMC_vs_pt);
  gr_dataMC_MPF_EXTRAP->SetMarkerStyle(23);
  gr_dataMC_MPF_EXTRAP->SetMarkerSize(markerSize);
  gr_dataMC_MPF_EXTRAP->SetMarkerColor(49);
  gr_dataMC_MPF_EXTRAP->SetName(TString::Format("%s_MPFchs_extrap_a%s_%s", prefix.c_str(), alphaCut.c_str(), fullEtaRegion.c_str()));
  
  if (outputFile) {
    outputFile->cd();
    gr_dataMC_BALANCING->Write();
    gr_dataMC_MPF->Write();
    gr_dataMC_EXTRAP->Write();
    gr_dataMC_MPF_EXTRAP->Write();

    gr_responseBALANCING_vs_pt->Write();
    gr_responseBALANCINGMC_vs_pt->Write();
    gr_responseMPF_vs_pt->Write();
    gr_responseMPFMC_vs_pt->Write();
    gr_responseEXTRAP_vs_pt->Write();
    gr_responseEXTRAPMC_vs_pt->Write();
    gr_responseMPFExtrap_vs_pt->Write();
    gr_responseMPFExtrapMC_vs_pt->Write();
    if(resp_reso == "response"){
    H_countMC_vs_pt_mpf->Write();
    H_count_vs_pt_mpf->Write();
    H_count_vs_pt->Write();
    H_countMC_vs_pt->Write();
    
    }
    if(isfinebinning && isresolution){
      gr_responseKFSR_vs_pt    -> Write();
      gr_MC_responseKFSR_vs_pt -> Write();
    }
  }
  
  drawGraphs(gr_responseBALANCING_vs_pt, gr_responseBALANCINGMC_vs_pt, xMin, xMax, "p_{T} Balance", db, etaRegion, legendTitle, resp_reso, "balancing", rawJets);
  drawGraphs(gr_responseMPF_vs_pt, gr_responseMPFMC_vs_pt, xMin, xMax, "MPF", db, etaRegion, legendTitle, resp_reso, "mpf", rawJets);
  drawGraphs(gr_responseEXTRAP_vs_pt, gr_responseEXTRAPMC_vs_pt, xMin, xMax, "p_{T} Balance extrap.", db, etaRegion, legendTitle, resp_reso, "balancing_extrap", rawJets);
  drawGraphs(gr_responseMPFExtrap_vs_pt, gr_responseMPFExtrapMC_vs_pt, xMin, xMax, "MPF extrap.", db, etaRegion, legendTitle, resp_reso, "mpf_extrap", rawJets);
  //  drawGraphs(gr_responseMPFExtrap_vs_pt, gr_responseMPFExtrapMC_vs_pt, xMin, xMax, "MPF extrap.", db, etaRegion, legendTitle, resp_reso, "mpf_extrap", rawJets);
  drawGraphsratio(gr_dataMC_EXTRAP, xMin, xMax, "p_{T} Balance #sigma extrap.", db, etaRegion, legendTitle, resp_reso, "balancing_reso_extrap",true);
 if(isresolution) drawGraphsratio(gr_dataMC_BALANCING, xMin, xMax, "p_{T} Balance #sigma extrap.", db, etaRegion, legendTitle, resp_reso, "balancing_reso",true); 
 if(isfinebinning && isresolution){ drawGraphsextrap(gr_responseKFSR_vs_pt, xMin, xMax, "Kfsr", db, etaRegion, legendTitle, resp_reso, "balancing_Kfsr_data");}
 if(isfinebinning && isresolution){ drawGraphsextrap(gr_MC_responseKFSR_vs_pt, xMin, xMax, "Kfsr", db, etaRegion, legendTitle, resp_reso, "balancing_Kfsr_MC");}
  file_noextrap->Close();
  file_extrap->Close();
}

void draw_scalefactorVsEta( drawBase* db, double  sfMC[], double sfDATA[], double sfMCerr[], double sfDATAerr[], int Npoint, TFile* outputFile)
{
      std::string name_base = db->get_outputdir();
      std::string graph_name ;
      name_base += "/";
    
    double Scalefactor [Npoint]    = {0.}; 
    double Scalefactorerr [Npoint] = {0.};
    double x[Npoint] = {0.};
    
    double Scalefactor_incl [Npoint]    = {0.}; 
    double Scalefactorerr_incl [Npoint] = {0.};
    
      fineEtaBinning fineetaBinning;
      EtaBinning     etaBinning;
    
    if(Npoint <= 7 ){
    name_base += "wide_eta_binning_";
    graph_name += "wide_eta_binning_"; 
   /* x[0] = 0.4;
    x[1] = 1.05;
    x[2] = 1.6;
    x[3] = 2.2;
    x[4] = 2.7;
    x[5] = 3.1;*/
    for(size_t i = 0 ; i < etaBinning.size() ; ++i){
     
     x[i] = etaBinning.getBinValue(i).first + ( etaBinning.getBinValue(i).second - etaBinning.getBinValue(i).first)/2. ;
     std::cout<<" wide eta bin "<<x[i]<<std::endl;
    }
    
    }else{
    name_base += "fine_eta_binning_";
    graph_name += "fine_eta_binning_";
    for(size_t i = 0 ; i < fineetaBinning.size() ; ++i){
     
     x[i] = fineetaBinning.getBinValue(i).first + ( fineetaBinning.getBinValue(i).second - fineetaBinning.getBinValue(i).first)/2. ;
     std::cout<<" fine eta bin "<<x[i]<<std::endl;
    }
    
    
    }
   // x[6] = 4.2;
    
    
    double x_err[Npoint] = {0.};
    for(int i =0 ; i < Npoint ; i++){
    
      Scalefactor [i]    = sfDATA[i]  ;
      Scalefactorerr [i] = sfDATAerr[i];
      
       std::cout<<" eta SF "<<Scalefactor [i]<<" +- "<<Scalefactorerr [i]<<std::endl;      
      
    }
    
   TLatex* box = new TLatex();
  box->SetNDC();
  box->SetTextSize(0.04);
  box->SetTextAlign(31);
  box->SetTextFont(42);
  box->SetLineWidth(2);
 
  TLatex* boxcms = new TLatex();
  boxcms->SetNDC();
  boxcms->SetTextFont(61);
  boxcms->SetLineWidth(2);
  boxcms->SetTextSize(0.05);
  
  TLatex* boxprel = new TLatex();
  
  boxprel->SetNDC();
  boxprel->SetTextAlign(13);
  boxprel->SetTextFont(52);
  boxprel->SetTextSize(0.04);
  boxprel->SetLineWidth(2);
    
  graph_name += "Scale_factor_res_vs_ETA";
  name_base  +="Scale_factor_res_vs_ETA.pdf";
  TGraphErrors* gr_sf = new TGraphErrors(Npoint, x, Scalefactor , x_err, Scalefactorerr);
    
  gr_sf->SetTitle(graph_name.c_str());
  gr_sf->SetLineWidth(2);
  gr_sf->SetLineColor(4);
  gr_sf->SetMarkerSize(1);
  gr_sf->SetMarkerStyle(34);
  gr_sf->SetMarkerColor(kBlue);
  
  gr_sf->SetName(graph_name.c_str());
  gr_sf->GetXaxis()->SetTitle("#eta");
  gr_sf->GetYaxis()->SetTitle("Resolution Scale Factor");
  gr_sf->GetYaxis()->SetRangeUser(0.8,2.);
  
  if (outputFile) {
    outputFile->cd();
    gr_sf->Write();
    
    }
  
  float cmsTextSize = 0.043;
  TPaveText* label_cms = db->get_labelCMS(1);
  label_cms->SetTextSize(cmsTextSize);

  TPaveText* label_sqrt = db->get_labelSqrt(1);

  
  

  TCanvas *screen_1 = new TCanvas("screen1", "screen1", 800, 800);
  gr_sf->Draw("APE");
  //box->DrawLatex(0.95,0.91,"8.6 pb^{-1} (13 TeV)");
  //boxcms->DrawLatex(0.12,0.91,"CMS");
  //boxprel->DrawLatex(0.24,0.935,"Preliminary");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  screen_1->SaveAs(name_base.c_str());
  delete screen_1;
}
