#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime
from ROOT import *

parser = argparse.ArgumentParser(description='Generate PU profile')
parser.add_argument('--Cert_json',type=str,dest="Cert_json",default=1, help=" cert Json file of the analysis")
parser.add_argument('--whichrun',type=str,dest="whichrun",default=1, help="Run era")
args = parser.parse_args()
Cert_json = args.Cert_json #path to the processed lumi JSON file
Runsuffix = args.whichrun # which run
print ( Cert_json)
if not os.path.exists(Cert_json):
  print("Error: %s not found" % lumi_json)
  exit(1)
  
print("\tDownloading the latest pileup file")
cmd_1="rm pileup_latest.txt"
cmd_2="cp /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt ."
print(cmd_2)

os.system(cmd_1)
os.system(cmd_2)
pileup_latest = "pileup_latest.txt"
  
    
HLT_list=['HLT_Photon33_v*',
     'HLT_Photon50_R9Id90_HE10_IsoM_v*',
     'HLT_Photon75_R9Id90_HE10_IsoM_v*',
     'HLT_Photon90_R9Id90_HE10_IsoM_v*',
     'HLT_Photon120_R9Id90_HE10_IsoM_v*',
     'HLT_Photon165_R9Id90_HE10_IsoM_v*',
'HLT_Photon200_v*']
HLT_prefix_list=['33','50','75','90','120','165', '200']

today = datetime.date.today()
today.strftime('%d-%m-%Y')

outputdirectory="PUfiles_"+Runsuffix+"_"+str(today)

cmd="mkdir -p "+outputdirectory
os.system(cmd)

print("\tRunning brilcalc for each HLT")
i_HLT=0
for iHLT in HLT_list:
	
	YourOutput="PUdataHLTphoton"+HLT_prefix_list[i_HLT]+Runsuffix
	print("\tRunning brilcalc for"+HLT_prefix_list[i_HLT])
	cmd1="brilcalc lumi -i "+Cert_json+" --hltpath "+iHLT+" --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json --byls --minBiasXsec 69200 -o "+outputdirectory+"/"+YourOutput+".csv"
	print(cmd1)
	os.system(cmd1)
	print("\tRunning pileupReCalc_HLTpaths.py for"+HLT_prefix_list[i_HLT])
	cmd2="pileupReCalc_HLTpaths.py -i "+outputdirectory+"/"+YourOutput+".csv --inputLumiJSON pileup_latest.txt -o "+outputdirectory+"/"+YourOutput+".txt --runperiod Run2"
	print(cmd2)
	os.system(cmd2)
	outputROOT = "pu_truth_data2017UL_100bins_HLTphoton"+HLT_prefix_list[i_HLT]+Runsuffix+".root"
	print("\tRunning pileupCalc "+HLT_prefix_list[i_HLT])
	cmd3= "pileupCalc.py -i "+Cert_json+" --inputLumiJSON "+outputdirectory+"/"+YourOutput+".txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 "+outputROOT
	print(cmd3)
	os.system(cmd3)
	print ("Create %s" % outputROOT)
	i_HLT+=1








 


