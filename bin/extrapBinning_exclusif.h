#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

#include "ptBinning.h"

class ExtrapBinning {
  public:
    ExtrapBinning() {}

    void initialize(PtBinning ptBinning, const std::string& recoType) {
      
      /*size_t s = ptBinning.size();
      for (size_t i = 0; i < s; i++) {

	//        std::pair<float, float> bin = ptBinning.getBinValue(i);
        float minPt = 0.;
        float ptStep = 0.;

	/*
	//old binning
        // In percent 
          if (bin.first <= 80.) {
            minPt = 8.;
            ptStep = 3.;
          } else if (bin.first <= 350) {
            minPt = 5.;
            ptStep = 2.5;
          } else {
            minPt = 2.;
            ptStep = 2.;
          }

        // In percent -- federico
	  minPt = 0.0;
	  ptStep = 3.0;

        minPt /= 100.;
        ptStep /= 100.;
        float maxPt = minPt + ptStep;

	mMapping.push_back(std::make_pair(minPt, maxPt));

      }*///end ciclo su ptBin
      mMapping.push_back(std::make_pair(0.,0.05));
      mMapping.push_back(std::make_pair(0.05,0.1));
      mMapping.push_back(std::make_pair(0.1,0.15));
      mMapping.push_back(std::make_pair(0.15,0.2));
      mMapping.push_back(std::make_pair(0.2,0.25));
      mMapping.push_back(std::make_pair(0.25,0.3));
    } //end initialize

    int getBin(float ptPhoton, float ptSecondJet, int ptBin) const {


      
     // std::pair<float, float> mapping = mMapping[ptBin]; // mPtBinning.getBinValue(ptBin);

      double alpha = ptSecondJet / ptPhoton;
     
      
      size_t extrapBin = 0 ;//(size_t) floor((alpha - mapping.first) / (mapping.second - mapping.first));
      if(alpha < 0.3 && alpha >=0.25) extrapBin = 5 ;
      if(alpha < 0.25 && alpha >=0.2) extrapBin = 4 ;
      if(alpha < 0.20 && alpha >=0.15) extrapBin = 3 ;
      if(alpha < 0.15 && alpha >=0.1) extrapBin = 2 ;
      if(alpha < 0.1 && alpha >=0.05) extrapBin = 1 ;
      if(alpha < 0.05)   extrapBin = 0 ;
      
      
      return /*(extrapBin >= size()) ? -1 :*/ extrapBin;
    }

    size_t size() const {
      return mSize;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mMapping[bin];
    }

  private:
    PtBinning mPtBinning;
    static const int mSize = 6;

    std::vector<std::pair<float, float> > mMapping; // first is minPt, second is maxPt
};
