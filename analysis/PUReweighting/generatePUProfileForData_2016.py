#! /bin/env python

import argparse, os, tempfile, shutil, sys
from subprocess import call, PIPE, STDOUT, Popen

parser = argparse.ArgumentParser(description='Generate PU profile')
parser.add_argument('--Cert_json',type=str,dest="Cert_json",default=1, help=" cert Json file of the analysis")
parser.add_argument('--pileup_latestHLT',type=str,dest="pileup_latestHLT",default=1, help=" pileup json file per HLT")
parser.add_argument('--whichHLT',type=str,dest="whichHLT",default=1, help="Photon HLT in use")
parser.add_argument('--whichrun',type=str,dest="whichrun",default=1, help="Run era")
#parser.add_argument('pileup_latest', nargs=1)
args = parser.parse_args()
Cert_json = args.Cert_json
Pileuplatest = args.pileup_latestHLT
HLTsuffix = args.whichHLT
Runsuffix = args.whichrun
#pileup_latest = args.pileup_latest[0]

print ( Cert_json)
#print ( pileup_latest)

if not os.path.exists(Cert_json):
  print("Error: %s not found" % lumi_json)
  exit(1)
#if not os.path.exists(pileup_latest):
#  print("Error: %s not found" % pileup_latest)
#  exit(1)

outputROOT = "pu_truth_data2016_100bins_HLTphoton"+HLTsuffix+Runsuffix+".root"

pileup_latest = "/afs/cern.ch/user/h/hlattaud/private/tutorial/CMSSW_8_0_24_patch1/src/JetMETCorrections/GammaJetFilter/PUReweightingHLT/"+Pileuplatest
 
print("\tRunning pileupCalc...")
cmd= "pileupCalc.py -i "+Cert_json+" --inputLumiJSON "+pileup_latest+" --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 "+outputROOT
print(cmd)

os.system(cmd)
print ("Create %s" % outputROOT)

