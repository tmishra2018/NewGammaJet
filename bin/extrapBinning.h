#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

#include "ptBinning.h"

class ExtrapBinning {
  public:
    ExtrapBinning() {}

    void initialize() {
      
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
     // mMapping.push_back(std::make_pair(0.,0.));
      mMapping.push_back(std::make_pair(0.,0.1));
      mMapping.push_back(std::make_pair(0.,0.15));
      mMapping.push_back(std::make_pair(0.,0.2));
      mMapping.push_back(std::make_pair(0.,0.25));
      mMapping.push_back(std::make_pair(0.,0.3));
      mMapping.push_back(std::make_pair(0.,0.5));
      mMapping.push_back(std::make_pair(0.,0.7));
      mMapping.push_back(std::make_pair(0.,1.0));	
    } //end initialize

    int getBin(float ptPhoton, float ptSecondJet, int ptBin) const {


      
     // std::pair<float, float> mapping = mMapping[ptBin]; // mPtBinning.getBinValue(ptBin);

      double alpha = ptSecondJet / ptPhoton;
     
      
      size_t extrapBin = 0 ;//(size_t) floor((alpha - mapping.first) / (mapping.second - mapping.first));
      if(alpha < 1.0) extrapBin = 7 ;
      if(alpha < 0.7) extrapBin = 6 ;
      if(alpha < 0.5) extrapBin = 5 ;
      if(alpha < 0.3) extrapBin = 4 ;
      if(alpha < 0.25) extrapBin = 3 ;
      if(alpha < 0.20) extrapBin = 2 ;
      if(alpha < 0.15) extrapBin = 1 ;
      if(alpha < 0.1) extrapBin = 0 ;
     // if(alpha ==0)   extrapBin = 0 ;
      
      
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
    static const int mSize = 8;

    std::vector<std::pair<float, float> > mMapping; // first is minPt, second is maxPt
};
