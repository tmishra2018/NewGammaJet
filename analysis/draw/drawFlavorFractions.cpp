#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include <sstream>

#include <TColor.h>
#include <TStyle.h>
#include <THStack.h>
#include <TPad.h>
#include <TAttFill.h>
#include <TLegend.h>

#include "TParameter.h"
#include "TError.h"
#include "drawBase.h"
#include "fitTools.h"
#include "etaBinning.h"
#include "ptBinning.h"

#include <boost/algorithm/string.hpp>

int main(int argc, char* argv[]) {

	if (argc != 4) {
		std::cout << "USAGE: ./drawPhotonJet [mc_SIGNAL_dataset] [recoType] [jetAlgo]" << std::endl;
		exit(23);
	}

	std::string mc_photonjet(argv[1]);
	std::string recoType(argv[2]);
	std::string jetAlgo(argv[3]);
	recoType = (recoType == "pf") ? "PFlow" : "PUPPI";
	jetAlgo = (jetAlgo == "ak4") ? "AK4" : "AK8";
	std::string postFix = recoType + jetAlgo;
	if(recoType == "PFlow") postFix += "chs";


	TString mc1FileName;
	mc1FileName = TString::Format("PhotonJet_%s_%s.root", mc_photonjet.c_str(), postFix.c_str()); 
	TFile* mcPhotonJetFile = TFile::Open(mc1FileName);
	if (mcPhotonJetFile) {
		std::cout << "Opened mc file '" << mc1FileName << "'." << std::endl;
	}

	TString outputDir = TString::Format("PhotonJetPlots_%s_%s_FlavorFractions", mc_photonjet.c_str(), postFix.c_str());
	mkdir(outputDir, 0777);

	TFile * inputfile;
	TString type;  
	EtaBinning etaBinning;
	size_t etaBinningSize = etaBinning.size() + 1;

	double alpha_cut = static_cast<TParameter<double>*>(mcPhotonJetFile->Get("analysis/alpha_cut"))->GetVal();
	std::stringstream ss;
	ss << ((int) (alpha_cut * 100));
	std::string alphaCut = ss.str();


	inputfile = mcPhotonJetFile; 
	type = "mc"; 

	PtBinning ptBinning;
	std::vector<std::pair<float, float> > ptBins = ptBinning.getBinning();
	size_t ptBinningSize = ptBinning.size();
	std::pair<float, float> currentBin;

	Float_t ptphot_bins[ptBinningSize+1];
	TH1F *hudFrac[etaBinningSize];
	TH1F *hsFrac[etaBinningSize];
	TH1F *hcFrac[etaBinningSize];
	TH1F *hbFrac[etaBinningSize];
	TH1F *hgluFrac[etaBinningSize];
	TH1F *hundefFrac[etaBinningSize];

	for(size_t j = 0; j<ptBinningSize; j++) {
		currentBin = ptBinning.getBinValue(j);
		ptphot_bins[j] = currentBin.first;
		if(j == ptBinningSize-1) ptphot_bins[j+1] = currentBin.second;
	}

	TString histoName;
	TString histoName_tmp;
	TH1F *h_tmp;

	float udfraction = -999.;
	float sfraction = -999.;
	float cfraction = -999.;
	float bfraction = -999.;
	float glufraction  = -999.;
	float undeffraction = -999.;
	float totalevts = -999;


	for (size_t i = 0; i < etaBinningSize-1; i++) {

		histoName = TString::Format("flavorFraction_ud_a%s_%s",  alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hudFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
		histoName = TString::Format("flavorFraction_s_a%s_%s", alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hsFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
		histoName = TString::Format("flavorFraction_c_a%s_%s", alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hcFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
		histoName = TString::Format("flavorFraction_b_a%s_%s", alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hbFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
		histoName = TString::Format("flavorFraction_glu_a%s_%s", alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hgluFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);
		histoName = TString::Format("flavorFraction_Undefined_a%s_%s", alphaCut.c_str(), etaBinning.getBinName(i).c_str());
		hundefFrac[i] = new TH1F(histoName,histoName,ptBinningSize,ptphot_bins);

		for(size_t j = 0; j<ptBinningSize; j++) {
                histoName_tmp=TString::Format("fSumEntries_%s", etaBinning.getBinName(i).c_str());
			currentBin = ptBinning.getBinValue(j);
			histoName_tmp.Append(TString::Format("_ptPhot_%i",int(currentBin.first)));
			histoName_tmp.Append(TString::Format("_%i", int(currentBin.second)));
			h_tmp=(TH1F*)inputfile->Get("analysis/flavorcomposition/"+histoName_tmp);

			totalevts=h_tmp->GetEntries();
			if(totalevts>0.) {
				udfraction = (h_tmp->GetBinContent(h_tmp->FindBin(1)) + h_tmp->GetBinContent(h_tmp->FindBin(2)))/totalevts;
				sfraction = (h_tmp->GetBinContent(h_tmp->FindBin(3)))/totalevts;
				cfraction = (h_tmp->GetBinContent(h_tmp->FindBin(4)))/totalevts;
				bfraction = (h_tmp->GetBinContent(h_tmp->FindBin(5)))/totalevts;
				glufraction = (h_tmp->GetBinContent(h_tmp->FindBin(21)))/totalevts;
				undeffraction = (h_tmp->GetBinContent(h_tmp->FindBin(0)))/totalevts;
			}
			hudFrac[i]->SetBinContent(j+1,udfraction);
			hsFrac[i]->SetBinContent(j+1,sfraction);
			hcFrac[i]->SetBinContent(j+1,cfraction);
			hbFrac[i]->SetBinContent(j+1,bfraction);
			hgluFrac[i]->SetBinContent(j+1,glufraction);
			hundefFrac[i]->SetBinContent(j+1,undeffraction);

		}//loop on pt bins


 hudFrac[i]->SetMarkerStyle(20);
 hsFrac[i]->SetMarkerStyle(21);
 hcFrac[i]->SetMarkerStyle(22);
 hbFrac[i]->SetMarkerStyle(23);
 hgluFrac[i]->SetMarkerStyle(33);
 hundefFrac[i]->SetMarkerStyle(34);

	}      

//	TCanvas *c = new TCanvas("c","c",600,600);    
//	c->SaveAs(outputDir+"/"+type+"_"+stackName+".pdf");
//	c->Destructor();


TFile out(outputDir+"/"+type+"_alpha"+alphaCut+"_FlavorComposition.root","recreate");
out.cd();
for (size_t i = 0; i < etaBinningSize-1; i++) {
hudFrac[i]->Write();
hsFrac[i]->Write();
hcFrac[i]->Write();
hbFrac[i]->Write();
hgluFrac[i]->Write();
hundefFrac[i]->Write();
}
out.Close();
}


