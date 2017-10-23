#pragma once

#include <cmath>
#include <vector>
#include <utility>

class PhotonPtBinning {
  public:
    PhotonPtBinning() {
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

    void fillPtBins() { 
         
      mPtBins.push_back(std::make_pair(40., 45.));
      mPtBins.push_back(std::make_pair(45., 50.));
      mPtBins.push_back(std::make_pair(50., 55.));
      mPtBins.push_back(std::make_pair(55., 60.));
      mPtBins.push_back(std::make_pair(60., 70.));
      mPtBins.push_back(std::make_pair(70., 85.));
      mPtBins.push_back(std::make_pair(85., 95.));
      mPtBins.push_back(std::make_pair(95., 105.));
      mPtBins.push_back(std::make_pair(105., 115.));
      mPtBins.push_back(std::make_pair(115., 130.));
      mPtBins.push_back(std::make_pair(130., 150.));
      mPtBins.push_back(std::make_pair(150., 175.));
      mPtBins.push_back(std::make_pair(175., 200.));
      mPtBins.push_back(std::make_pair(200., 230.));
      mPtBins.push_back(std::make_pair(230., 265.));
      mPtBins.push_back(std::make_pair(265., 300.));
      mPtBins.push_back(std::make_pair(300., 350.));
      mPtBins.push_back(std::make_pair(350., 400.));
      mPtBins.push_back(std::make_pair(400., 450.));
      mPtBins.push_back(std::make_pair(450., 500.));
      mPtBins.push_back(std::make_pair(500., 600.));
      mPtBins.push_back(std::make_pair(600., 700.));
      mPtBins.push_back(std::make_pair(700., 850.));
      mPtBins.push_back(std::make_pair(850., 1000.));
      mPtBins.push_back(std::make_pair(1000.,2000.));
      mPtBins.push_back(std::make_pair(2000.,3000.));
      
      
    } 
};
