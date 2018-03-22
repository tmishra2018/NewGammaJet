#pragma once

#include <cmath>
#include <vector>
#include <utility>

class PhotonPtBinning {
  public:
    PhotonPtBinning() {
      fillPtBins();
    }

    int getPtBin(double pt) {
      std::vector<std::pair<double, double> >::const_iterator it = mPtBins.begin();
      for (; it != mPtBins.end(); ++it) {
        std::pair<double, double> bin = *it;
        if (pt >= bin.first && pt < bin.second) {
          return it - mPtBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mPtBins.size();
    }

    std::pair<double, double> getBinValue(int bin) const {
      return mPtBins[bin];
    }

    std::vector<std::pair<double, double> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<double, double> >(mPtBins.begin(), mPtBins.begin() + n);
    }

    std::vector<std::pair<double, double> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<double, double> >(mPtBins.begin() + from, mPtBins.begin() + to);
    }

  private:
    std::vector<std::pair<double, double> > mPtBins;

    void fillPtBins() { 

//    for (double ptforbins = 40; ptforbins < 400. - 2.5; ptforbins += 5.)
    for (double ptforbins = 60; ptforbins < 400. - 2.5; ptforbins += 5.)
            mPtBins.push_back(std::make_pair(ptforbins,ptforbins+5.));
    for (double ptforbins = 400; ptforbins < 1000. - 5.; ptforbins += 10.)
            mPtBins.push_back(std::make_pair(ptforbins,ptforbins+10.));
    for (double ptforbins = 1000; ptforbins < 2950. - 25.; ptforbins += 50.)
            mPtBins.push_back(std::make_pair(ptforbins,ptforbins+50.));
//fill last bin by hand to be sure it's exactly 3000
            mPtBins.push_back(std::make_pair(2950.,3000.00001));

/*
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
   */   
      
    } 
};
