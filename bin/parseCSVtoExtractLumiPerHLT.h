#ifndef __parseCSVtoExtractLumiPerHLT_C__
#define __parseCSVtoExtractLumiPerHLT_C__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <utility>

using namespace std;
map<int, map<int, double> > m_LUMI;
map<int, map<int, double> > m_30;
map<int, map<int, double> > m_50;
map<int, map<int, double> > m_75;
map<int, map<int, double> > m_90;
map<int, map<int, double> > m_120;
map<int, map<int, double> > m_175;

double getLUMItot(int run, int ls) {

  return m_LUMI[run][ls];
}

double getPrescaleperHLT(int run, int ls, int HLT) {

  if(HLT==1)return m_30[run][ls];
  if(HLT==2)return m_50[run][ls];
  if(HLT==3)return m_75[run][ls];
  if(HLT==4)return m_90[run][ls];
  if(HLT==5)return m_120[run][ls];
  if(HLT==6)return m_175[run][ls];
  return -999;
}

int parsePrescalefromJson(string filename="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_NoHLT.csv") {

  string filename30="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_30.csv";
  string filename50="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_50.csv";
  string filename75="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_75.csv";
  string filename90="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_90.csv";
  string filename120="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_120.csv";
  string filename175="/afs/cern.ch/work/h/hlattaud/private/Ploting_area/CMSSW_8_0_25/src/JetMETCorrections/GammaJetFilter/bin/LUMdataHLT2016_runBCDEFGH_HLT_photon_175.csv";
  cout << "Opening " << filename << "..."<<filename30<< "..."<<filename50<< "..."<<filename75<< "..."<<filename90<< "..."<<filename120<< "..."<<filename175<<std::endl;

  string line,line30,line50,line75,line90,line120,line175;
  ifstream file(filename);
  ifstream file30(filename30);
  ifstream file50(filename50);
  ifstream file75(filename75);
  ifstream file90(filename90);
  ifstream file120(filename120);
  ifstream file175(filename175);

  if (file.is_open() && file30.is_open() && file50.is_open() && file75.is_open() && file90.is_open() && file120.is_open() && file175.is_open()){
    cout << "ok" << endl;

    //loop over lines in file
    while ( getline(file,line) ){

      string str, run_str, ls_str;
      int delim_pos;
      double LUMIrecorded  = -1;
      
      if ( line.at(0) != '#' ){

        //loop over strings in line
        for (int string_num=0; (delim_pos = line.find(",")) != -1; string_num++){

          str = line.substr(0, delim_pos);
          line.erase(0, delim_pos + 1);

          if (string_num == 0)  //first string holds run number
            run_str = str.substr(0, str.find(":"));

          else if (string_num == 1) //second string has ls
            ls_str = str.substr(0, str.find(":"));

          else if (string_num == 6) //seventth string has recorded lumi
            LUMIrecorded = stod( str );
        }
        
       
      //   cout<<"test arguments for string run nb : "<<run_str << " lumi section : "<<ls_str<<" recorded LUMI : "<<LUMIrecorded<<endl;
       
        int run = stoi( run_str );
        int ls = stoi( ls_str );
     //   cout<<"test arguments for stoi run nb : "<<run << " lumi section : "<<ls<<" PU : "<<PU<<endl;
        m_LUMI[run][ls] = LUMIrecorded;
      }
    }
    
    file.close();
    
    
    while ( getline(file30,line30) && getline(file50,line50) && getline(file75,line75) && getline(file90,line90) && getline(file120,line120) && getline(file175,line175) ){

      string  str30, run_str30, ls_str30, str50, run_str50, ls_str50, str75, run_str75, ls_str75, str90, run_str90, ls_str90, str120, run_str120, ls_str120, str175, run_str175, ls_str175;
      int delim_pos;
      int delim_pos50;
      int delim_pos75;
      int delim_pos90;
      int delim_pos120;
      int delim_pos175;
      double prescale30  = -1;
      double prescale50  = -1;
      double prescale75  = -1;
      double prescale90  = -1;
      double prescale120 = -1;
      double prescale175 = -1;
      if ( line30.at(0) != '#'){

        //loop over strings in line
        for (int string_num=0; (delim_pos = line30.find(",")) != -1; string_num++){

          delim_pos50 = line50.find(",");
          delim_pos75 = line75.find(",");
          delim_pos90 = line90.find(",");
          delim_pos120 = line120.find(",");
          delim_pos175 = line175.find(",");
          
          str30 = line30.substr(0, delim_pos);
          str50 = line50.substr(0, delim_pos50);
          str75 = line75.substr(0, delim_pos75);
          str90 = line90.substr(0, delim_pos90);
          str120 = line120.substr(0, delim_pos120);
          str175 = line175.substr(0, delim_pos175);
          line30.erase(0, delim_pos + 1);
          line50.erase(0, delim_pos50 + 1);
          line75.erase(0, delim_pos75 + 1);
          line90.erase(0, delim_pos90 + 1);
          line120.erase(0, delim_pos120 + 1);
          line175.erase(0, delim_pos175 + 1);

          if (string_num == 0)  //first string holds run number
            run_str30 = str30.substr(0, str30.find(":"));

          else if (string_num == 1) //second string has ls
            ls_str30 = str30.substr(0, str30.find(":"));

          else if (string_num == 5){ //sixth string has delivred LUMI
          
          //std::cout<<" test prescale "<<str30<<std::endl;
            prescale30 = stod( str30 );
                  //    std::cout<<" test prescale "<<str50<<std::endl;
            prescale50 = stod( str50 );
                 //     std::cout<<" test prescale "<<str75<<std::endl;
            prescale75 = stod( str75 );
                 //     std::cout<<" test prescale "<<str90<<std::endl;
            prescale90 = stod( str90 );
                //      std::cout<<" test prescale "<<str120<<std::endl;
            prescale120 = stod( str120 );
                //      std::cout<<" test prescale "<<str175<<std::endl;
            prescale175 = stod( str175 );
            
            }
        
       }
       
        int run = stoi( run_str30 );
        int ls = stoi( ls_str30 );

        m_30[run][ls] = prescale30 / getLUMItot(run,ls);
        m_50[run][ls] = prescale50 / getLUMItot(run,ls);
        m_75[run][ls] = prescale75 / getLUMItot(run,ls);
        m_90[run][ls] = prescale90 / getLUMItot(run,ls);
        m_120[run][ls] = prescale120 / getLUMItot(run,ls);
        m_175[run][ls] = prescale175 / getLUMItot(run,ls);
       
      }
    }
    
   file30.close();
   file90.close();
   file50.close();
   file75.close();
   file120.close();
   file175.close(); 
    
  }
  else
    cout << "Unable to open files" << endl;

  return 0;
}

#endif //__parseCSVtoExtractLumiPerHLT_C__
