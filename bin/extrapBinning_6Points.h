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

      /*
      size_t s = ptBinning.size();
      for (size_t i = 0; i < s; i++) {

        std::pair<float, float> bin = ptBinning.getBinValue(i);
        float minPt = 0.;
        float ptStep = 0.;

		//old binning
        // In percent 
        if (recoType == "PFlow" || recoType == "JPT") {
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
        } else {
          if (bin.first <= 60.) {
            minPt = 8.;
            ptStep = 1.5;
          } else if (bin.first <= 350.) {
            minPt = 6.;
            ptStep = 1.5;
          } else {
            minPt = 2.;
            ptStep = 2.;
          }
	}

        // In percent -- federico
        if (recoType == "PFlow" || recoType == "JPT") {
	  minPt = 2.0;
	  ptStep = 3.5;
        } else {
          if (bin.first <= 60.) {
            minPt = 8.;
            ptStep = 1.5;
          } else if (bin.first <= 350.) {
            minPt = 6.;
            ptStep = 1.5;
          } else {
            minPt = 2.;
            ptStep = 2.;
          }
	}
	

        minPt /= 100.;
        ptStep /= 100.;
        float maxPt = minPt + ptStep;

	mMapping.push_back(std::make_pair(minPt, maxPt));

      }//end ciclo su ptBin
	*/

	//Mikko's Test
      mMapping.push_back(std::make_pair(0. , 0.05 ));
      mMapping.push_back(std::make_pair(0.05 , 0.1 ));
      mMapping.push_back(std::make_pair(0.1 , 0.15 ));
      mMapping.push_back(std::make_pair(0.15 , 0.20 ));
      mMapping.push_back(std::make_pair(0.20 , 0.25 ));
      mMapping.push_back(std::make_pair(0.25 , 0.30 ));


    } //end initialize

    int getBin(float ptPhoton, float ptSecondJet, int ptBin) const {
      
      //      std::pair<float, float> mapping = mMapping[ptBin]; // mPtBinning.getBinValue(ptBin);
      
      double alpha = ptSecondJet / ptPhoton;
      //      size_t extrapBin = (size_t) floor((alpha - mapping.first) / (mapping.second - mapping.first));
      
      //      return (extrapBin >= size()) ? -1 : extrapBin;
      
      std::vector<std::pair<float, float> >::const_iterator it = mMapping.begin();                                                                                           
      for (; it != mMapping.end(); ++it) {                                                                                                                                   
	std::pair<float, float> bin = *it;                                                                                                                                  
	if (alpha > bin.first && alpha <= bin.second) {                                                                                                                           
	  return it - mMapping.begin();                                                                                                                                      
	}                                
      }
      return -1;
	}
    
    size_t size() const {
      return mSize;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mMapping[bin];
    }

  private:
    PtBinning mPtBinning;
    // static const int mSize = 10;
    //federico
    static const int mSize = 6;

    std::vector<std::pair<float, float> > mMapping; // first is minPt, second is maxPt
};
