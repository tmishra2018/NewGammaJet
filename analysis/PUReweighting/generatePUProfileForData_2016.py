#! /bin/env python

import argparse, os, tempfile, shutil, sys
from subprocess import call, PIPE, STDOUT, Popen

parser = argparse.ArgumentParser(description='Generate PU profile')
parser.add_argument('Cert_json', nargs=1)
#parser.add_argument('pileup_latest', nargs=1)
args = parser.parse_args()

Cert_json = args.Cert_json[0]
#pileup_latest = args.pileup_latest[0]

print ( Cert_json)
#print ( pileup_latest)

if not os.path.exists(Cert_json):
  print("Error: %s not found" % lumi_json)
  exit(1)
#if not os.path.exists(pileup_latest):
#  print("Error: %s not found" % pileup_latest)
#  exit(1)

outputROOT = "pu_truth_data2016_100bins.root"

pileup_latest = "/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt"
 
print("\tRunning pileupCalc...")
cmd= "pileupCalc.py -i "+Cert_json+" --inputLumiJSON "+pileup_latest+" --calcMode true --minBiasXsec 71300 --maxPileupBin 100 --numPileupBins 100 "+outputROOT
print(cmd)

os.system(cmd)
print ("Create %s" % outputROOT)

