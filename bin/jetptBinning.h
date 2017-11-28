#pragma once

#include <cmath>
#include <vector>
#include <utility>

class JetPtBinning {
  public:
    JetPtBinning() {
      fillPtBins();
    }

    int getPtBin(float pt) {
      std::vector<std::pair<float, float> >::const_iterator it = mPtBins.begin();
      for (; it != mPtBins.end(); ++it) {
        std::pair<float, float> bin = *it;
        if (pt >= bin.first && pt < bin.second) {
          return it - mPtBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mPtBins.size();
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mPtBins[bin];
    }

    std::vector<std::pair<float, float> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<float, float> >(mPtBins.begin(), mPtBins.begin() + n);
    }

    std::vector<std::pair<float, float> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<float, float> >(mPtBins.begin() + from, mPtBins.begin() + to);
    }

  private:
    std::vector<std::pair<float, float> > mPtBins;

    void fillPtBins() { //2016

    for (double ptforbins = 15; ptforbins < 50. - 0.25; ptforbins += 0.5)
	    mPtBins.push_back(std::make_pair(ptforbins,ptforbins+0.5));
    
    for (double ptforbins = 50.; ptforbins < 300. - 0.5; ptforbins += 1.)
	    mPtBins.push_back(std::make_pair(ptforbins,ptforbins+1.));
 
    for (double ptforbins = 300.; ptforbins < 1000. - 1.; ptforbins += 5.)
        mPtBins.push_back(std::make_pair(ptforbins,ptforbins+5.));
    
    for (double ptforbins = 1000.; ptforbins < 3500 - 1.; ptforbins += 10.)
        mPtBins.push_back(std::make_pair(ptforbins,ptforbins+10.));
    
/*for (const auto& p : mPtBins)
{
  std::cout <<"bin = "<< p.first << ", " << p.second << std::endl;
  } 
*/
 
    } 
};
