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
      
      if( ptSecondJet == 0.  ) alpha = 2.3  / ptPhoton ;
     
      
      size_t extrapBin = -1 ;//(size_t) floor((alpha - mapping.first) / (mapping.second - mapping.first));
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

    std::vector<std::pair<float, float> > getBinning(){
      return mMapping;
    }
  private:
    PtBinning mPtBinning;
    static const int mSize = 6;

    std::vector<std::pair<float, float> > mMapping; // first is minPt, second is maxPt
};
