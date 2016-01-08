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
      mRunBins.push_back(std::make_pair(256630, 257613)); // 210.285 pb-1
      mRunBins.push_back(std::make_pair(257614, 257969)); // 223.401 pb-1
      mRunBins.push_back(std::make_pair(258129, 258177)); // 246.562 pb-1
      mRunBins.push_back(std::make_pair(258211, 258448)); // 229.871 pb-1
      mRunBins.push_back(std::make_pair(258655, 258713)); // 205.806 pb-1
      mRunBins.push_back(std::make_pair(258714, 259685)); // 234.228 pb-1
      mRunBins.push_back(std::make_pair(259686, 259891)); // 208.717 pb-1
      mRunBins.push_back(std::make_pair(260373, 260532)); // 246.805 pb-1
      mRunBins.push_back(std::make_pair(260533, 260627)); // 287.615 pb-1
      // end of pp collision 2015
    }

};
