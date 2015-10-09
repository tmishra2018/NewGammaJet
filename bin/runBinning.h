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
        if (n >= bin.first && n < bin.second) {
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
      // bins read: >= && <
      mRunBins.push_back(std::make_pair(256630, 256678));
      mRunBins.push_back(std::make_pair(256728, 256730));
      mRunBins.push_back(std::make_pair(256734, 256844));
      mRunBins.push_back(std::make_pair(256866, 257600));
    }

};
