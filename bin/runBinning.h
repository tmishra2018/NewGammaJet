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
     
     
     /* mRunBins.push_back(std::make_pair(297046,  297223)); // almost 1/fb in each bin
      mRunBins.push_back(std::make_pair(297223, 297484)); // 
      mRunBins.push_back(std::make_pair(297485, 299602)); // 
      mRunBins.push_back(std::make_pair(297603, 299064)); // RunB
      mRunBins.push_back(std::make_pair(299065, 299367)); //Run C starting from now on 
      mRunBins.push_back(std::make_pair(299368, 300122)); // 
      mRunBins.push_back(std::make_pair(300123, 300463));
      mRunBins.push_back(std::make_pair(300464, 300575));
      mRunBins.push_back(std::make_pair(300576, 300784));
      mRunBins.push_back(std::make_pair(300785, 301297));
      mRunBins.push_back(std::make_pair(301298, 301460)); 
      mRunBins.push_back(std::make_pair(301461, 301958)); 
      mRunBins.push_back(std::make_pair(301959, 302029)); 
      mRunBins.push_back(std::make_pair(302030, 302343));*/ // 
       // 2017 Collisions
     
      // end of pp collision 2015
    }

};
