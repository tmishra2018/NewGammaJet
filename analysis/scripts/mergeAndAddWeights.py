#!/usr/bin/env python

import os
import sys
import argparse
import math
from ROOT import *
from array import array
import datetime

usage = "usage:  python mergeAndAddWeights.py --inputList list_of_files_to_merge.txt -o ./"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
        help="input list of files to be merged",
	    )
parser.add_argument("-o", "--outputDir", type=str, dest="outputDir", default="./",
        help="output directory",
	    )
parser.add_argument("--xsec", "--xsec", type=float, dest="Xsec", default="0",
        help="output directory",
	    )

args = parser.parse_args()
print args 

inputList = args.inputList
outputDir = args.outputDir
xsec=args.Xsec
###################
#read input file
ins = open(args.inputList,"r")
files = " "

for line in ins:
  file = line.strip()
  files += str(" "+file)
 # pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
 # fullname = #os.path.split(pathT2)[1]
 # name = fullname.split("_")
  
today = datetime.date.today()
today.strftime('%d-%m-%Y')

filename_out = outputDir+"/PhotonJet_2ndLevel_.root"#+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+name[5]+str(today)+".root" 
os.system("hadd -f "+filename_out+"  "+files )

inputFile = TFile(filename_out,"UPDATE")
h_sumOfWeights = inputFile.Get("h_sumW")
sumOfWeights = h_sumOfWeights.Integral()

##### update total luminosity ######
lumi = inputFile.Get("totallumi")
lumi.SetVal(1) # in /pb

##### update tree with weight for total normalization #####
analysis_tree = inputFile.Get("rootTupleTree/tree")
Putree = inputFile.Get("puvariable")

#analysis_tree.GetEntry(0)
#xsec = analysis_tree.crossSection
 
print xsec

evtWeightTot = array("f", [0.] )
b_evtWeightTot = analysis_tree.Branch("evtWeightTot", evtWeightTot,"evtWeightTot/F")
#import pdb; pdb.set_trace()
for event in analysis_tree:
  evtWeightTot[0] = (xsec / sumOfWeights)
  b_evtWeightTot.Fill()

inputFile.cd("rootTupleTree")
analysis_tree.Write("",TObject.kOverwrite)
lumi.Write()
inputFile.cd("../")
evtWeightTotpu = array("f", [0.] )
b_evtWeightTotpu = Putree.Branch("evtWeightTot", evtWeightTotpu,"evtWeightTotpu/F")
#import pdb; pdb.set_trace()
for event in Putree:
  evtWeightTotpu[0] = (xsec / sumOfWeights)
  b_evtWeightTotpu.Fill()

Putree.Write("",TObject.kOverwrite)




