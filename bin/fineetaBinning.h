#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <utility>

struct fineEtaBin {
  std::pair<float, float> bin;
  std::string name;
  std::string title;
};

class fineEtaBinning {
  public:
    fineEtaBinning() {
      fillEtaBins();
    }

    int getBin(float eta) const {
      eta = fabs(eta);
      std::vector<fineEtaBin>::const_iterator it = mEtaBins.begin();
      for (; it != mEtaBins.end(); ++it) {
        fineEtaBin bin = *it;
        if (eta >= bin.bin.first && eta < bin.bin.second) {
          return it - mEtaBins.begin();
        }
      }

      return -1;
    }

    std::string getBinName(int bin) const {
      return mEtaBins[bin].name;
    }

    std::string getBinTitle(int bin) const {
      return mEtaBins[bin].title;
    }

    std::pair<float, float> getBinValue(int bin) const {
      return mEtaBins[bin].bin;
    }


    size_t size() const {
      return mEtaBins.size();
    }
    
    std::vector< fineEtaBin > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<fineEtaBin >(mEtaBins.begin(), mEtaBins.begin() + n);
    }
    void Is_Jer_computation(bool is_Jer ){
      if(is_Jer){
        mJER_bining = true ; 
      }else{
        mJER_bining = false ;
      }
    
    
    }
    

  private:
    std::vector<fineEtaBin> mEtaBins;
    bool mJER_bining = false;
    void fillEtaBins() {

      // official fine bin 
      fineEtaBin bin;

      
      /*
      bin.bin = std::make_pair(0., 0.522);
      bin.name = "eta0005";
      bin.title = "|#eta| < 0.52";  
      mEtaBins.push_back(bin);
      
      
      bin.bin = std::make_pair(0.522,0.783);
      bin.name = "eta0508";
      bin.title = "0.5 #leq |#eta| < 0.8";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(0.783,1.131);
      bin.name = "eta0811";
      bin.title = "0.8 #leq |#eta| < 1.1";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.131,1.305);
      bin.name = "eta1113";
      bin.title = "1.1 #leq |#eta| < 1.3";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.305,1.740);
      bin.name = "eta1317";
      bin.title = "1.3 #leq |#eta| < 1.7";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.740,2.043);
      bin.name = "eta1720";
      bin.title = "1.7 #leq |#eta| < 2.0";
      mEtaBins.push_back(bin);
          
      bin.bin = std::make_pair(2.043,2.322);
      bin.name = "eta2023";
      bin.title = "2.0 #leq |#eta| < 2.3";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.322,2.853);
      bin.name = "eta2329";
      bin.title = "2.3 #leq |#eta| < 2.9";
      mEtaBins.push_back(bin);
      
      
      
      bin.bin = std::make_pair(2.853,5.191);
      bin.name = "eta2951";
      bin.title = "2.9 #leq |#eta| < 5.1";
      mEtaBins.push_back(bin);
      */
     
     // if(mJER_bining){
     
      bin.bin = std::make_pair(0., 0.261);
      bin.name = "eta0003";
      bin.title = "|#eta| < 0.26";  
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(0.261,0.522);
      bin.name = "eta0305";
      bin.title = "0.3 #leq |#eta| < 0.5";  
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(0.522,0.783);
      bin.name = "eta0508";
      bin.title = "0.5 #leq |#eta| < 0.8";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(0.783,1.044);
      bin.name = "eta0810";
      bin.title = "0.7 #leq |#eta| < 1.0";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.044,1.305);
      bin.name = "eta1013";
      bin.title = "1.0 #leq |#eta| < 1.3";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.305,1.479);
      bin.name = "eta1315";
      bin.title = "1.3 #leq |#eta| < 1.5";
      mEtaBins.push_back(bin);

      bin.bin = std::make_pair(1.479,1.653);
      bin.name = "eta1517";
      bin.title = "1.5 #leq |#eta| < 1.7";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(1.653,1.930);
      bin.name = "eta1719";
      bin.title = "1.7 #leq |#eta| < 1.9";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(1.930,2.172);
      bin.name = "eta1922";
      bin.title = "1.9 #leq |#eta| < 2.1";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.172,2.322);
      bin.name = "eta2223";
      bin.title = "2.1 #leq |#eta| < 2.3";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.322,2.500);
      bin.name = "eta2325";
      bin.title = "2.3 #leq |#eta| < 2.5";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.500,2.650);
      bin.name = "eta2526";
      bin.title = "2.5 #leq |#eta| < 2.7";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.650,2.853);
      bin.name = "eta2629";
      bin.title = "2.7 #leq |#eta| < 2.9";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.853,2.964);
      bin.name = "eta2930";
      bin.title = "2.9 #leq |#eta| < 3.0";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(2.964,3.139);
      bin.name = "eta3031";
      bin.title = "3.0 #leq |#eta| < 3.1";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(3.139,3.489);
      bin.name = "eta3135";
      bin.title = "3.1 #leq |#eta| < 3.5";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(3.489,3.839);
      bin.name = "eta3538";
      bin.title = "3.5 #leq |#eta| < 3.8";
      mEtaBins.push_back(bin);
      
      bin.bin = std::make_pair(3.839,5.191);
      bin.name = "eta3852";
      bin.title = "3.8 #leq |#eta| < 5.1";
      mEtaBins.push_back(bin);
      
      
      
      

    }
};
