#pragma once

#include <cmath>
#include <vector>
#include <utility>

class RunBinning {
  public:
    RunBinning() {
      fillRunBins();
    }

    int getRunBin(int n) {
      std::vector<std::pair<int, int> >::const_iterator it = mRunBins.begin();
      for (; it != mRunBins.end(); ++it) {
        std::pair<int, int> bin = *it;
        if (n >= bin.first && n <= bin.second) {
          return it - mRunBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mRunBins.size();
    }

    std::pair<int, int> getBinValue(int bin) const {
      return mRunBins[bin];
    }

    std::vector<std::pair<int, int> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<int, int> >(mRunBins.begin(), mRunBins.begin() + n);
    }

    std::vector<std::pair<int, int> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<int, int> >(mRunBins.begin() + from, mRunBins.begin() + to);
    }

  private:
    std::vector<std::pair<int, int> > mRunBins;

    void fillRunBins() {
      // bins reader: >= && <
      
     /*
     2016 collisions
      mRunBins.push_back(std::make_pair(272007, 274250)); // 210.285 pb-1
      mRunBins.push_back(std::make_pair(274251, 274421)); // 223.401 pb-1
      mRunBins.push_back(std::make_pair(274422, 274999)); // 246.562 pb-1
      mRunBins.push_back(std::make_pair(275000, 275309)); // 229.871 pb-1
      mRunBins.push_back(std::make_pair(275310, 275376)); // 205.806 pb-1
      mRunBins.push_back(std::make_pair(275377, 275657)); // 234.228 pb-1
      mRunBins.push_back(std::make_pair(275658, 276243)); // 208.717 pb-1
      mRunBins.push_back(std::make_pair(276244, 276315)); // 246.805 pb-1
      mRunBins.push_back(std::make_pair(276316, 276525)); // 287.615 pb-1
      mRunBins.push_back(std::make_pair(276526, 276587));
      mRunBins.push_back(std::make_pair(276588, 276811));
      mRunBins.push_back(std::make_pair(276812, 276831));
      mRunBins.push_back(std::make_pair(276832, 277072));
      mRunBins.push_back(std::make_pair(277073, 277166));
      mRunBins.push_back(std::make_pair(277167, 277772));
      mRunBins.push_back(std::make_pair(277773, 278308));
      mRunBins.push_back(std::make_pair(278309, 278509));
      mRunBins.push_back(std::make_pair(278510, 278802));
      mRunBins.push_back(std::make_pair(278803, 278820));
      mRunBins.push_back(std::make_pair(278821, 279694));
      mRunBins.push_back(std::make_pair(279695, 279767));
      mRunBins.push_back(std::make_pair(279768, 279931));
      mRunBins.push_back(std::make_pair(279932, 280191));
      mRunBins.push_back(std::make_pair(280192, 280384));
      mRunBins.push_back(std::make_pair(280385, 280919));
      mRunBins.push_back(std::make_pair(280920, 282037));
      mRunBins.push_back(std::make_pair(282038, 282814));
      mRunBins.push_back(std::make_pair(282815, 283270));
      mRunBins.push_back(std::make_pair(283271, 283408));
      mRunBins.push_back(std::make_pair(283409, 283830));
      mRunBins.push_back(std::make_pair(283831, 284044)); // 2016 Collisions
     */
     
      mRunBins.push_back(std::make_pair(297046,  297219)); // almost 1/fb in each bin
      mRunBins.push_back(std::make_pair(297220, 297435)); // 
      mRunBins.push_back(std::make_pair(297436, 297563)); // 
      mRunBins.push_back(std::make_pair(297564, 299061)); // RunB
      mRunBins.push_back(std::make_pair(299062, 299368)); //Run C starting from now on 
      mRunBins.push_back(std::make_pair(299369, 300122)); // 
      mRunBins.push_back(std::make_pair(300123, 300237));
      mRunBins.push_back(std::make_pair(300238, 300461));
      mRunBins.push_back(std::make_pair(300462, 300576));
      mRunBins.push_back(std::make_pair(300577, 300780));
      mRunBins.push_back(std::make_pair(300781, 301298)); 
      mRunBins.push_back(std::make_pair(301299, 301461)); 
      mRunBins.push_back(std::make_pair(301462, 301941)); 
      mRunBins.push_back(std::make_pair(301942, 302030));
      mRunBins.push_back(std::make_pair(302031, 302277));
      mRunBins.push_back(std::make_pair(302278, 302448));
      mRunBins.push_back(std::make_pair(302449, 302572));
      mRunBins.push_back(std::make_pair(302573, 302663));
      mRunBins.push_back(std::make_pair(302664, 303824));
      mRunBins.push_back(std::make_pair(303825, 304062));
      mRunBins.push_back(std::make_pair(304063, 304144));
      mRunBins.push_back(std::make_pair(304145, 304204));
      mRunBins.push_back(std::make_pair(304205, 304366));
      mRunBins.push_back(std::make_pair(304367, 304508));
      mRunBins.push_back(std::make_pair(304509, 304655));
      mRunBins.push_back(std::make_pair(304656, 304738));
      mRunBins.push_back(std::make_pair(304739, 304797));
      mRunBins.push_back(std::make_pair(304798, 305040));
      mRunBins.push_back(std::make_pair(305041, 305188));
      mRunBins.push_back(std::make_pair(305189, 305237));
      mRunBins.push_back(std::make_pair(305238, 305312));
      mRunBins.push_back(std::make_pair(305313, 305366));
      mRunBins.push_back(std::make_pair(305367, 305406));
      mRunBins.push_back(std::make_pair(305407, 305590));     
      mRunBins.push_back(std::make_pair(305591, 305814));
      mRunBins.push_back(std::make_pair(305815, 305862));
      mRunBins.push_back(std::make_pair(305863, 306091));
      mRunBins.push_back(std::make_pair(306092, 306135));
      mRunBins.push_back(std::make_pair(306136, 306155));
      mRunBins.push_back(std::make_pair(306156, 306155));
      mRunBins.push_back(std::make_pair(306156, 306462));
      // 
       // 2017 Collisions
     
      // end of pp collision 2015
    }

};
