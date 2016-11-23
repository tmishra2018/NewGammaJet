#!/usr/bin/env python

import os
import sys
import argparse
import math
from ROOT import *
from array import array
import datetime

usage = "usage:  python mergeData.py --inputList list_of_files_to_merge.txt -o ./  --lumi_tot 1000. [in /pb]"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
        help="input list of files to be merged",
	    )
parser.add_argument("-o", "--outputDir", type=str, dest="outputDir", default="./",
        help="output directory",
	    )
parser.add_argument("--lumi_tot", type=float, dest="lumi_tot", default=1,
        help="total luminosity",
	    )

args = parser.parse_args()
print args 
inputList = args.inputList
outputDir = args.outputDir
#xsec = args.xsec
lumi_tot = args.lumi_tot
###################
#read input file
ins = open(args.inputList,"r")
files1 = " "
files2 = " "

i = 0



for line in ins:
  if i < 550:
    file = line.strip()
    files1 += str(" "+file)
    pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
    fullname = os.path.split(pathT2)[1]
    name = fullname.split("_")
    i= i+1
  if i >= 550:
    file = line.strip()
    files2 += str(" "+file)
    pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
    fullname = os.path.split(pathT2)[1]
    name = fullname.split("_")
    i= i+1

print i
print files1
print files2

today = datetime.date.today()
today.strftime('%d-%m-%Y')

pwd = os.environ['PWD']

if i < 550:
  filename_out = outputDir+"/PhotonJet_2ndLevel_"+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+".root" 
  print filename_out
  os.system("hadd -f "+filename_out+"  "+files1 )
else:
  filename_out = outputDir+"/PhotonJet_2ndLevel_"+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+".root" 
  filename_out1 = outputDir+"/PhotonJet_2ndLevel_"+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+"_Part1.root" 
  filename_out2 = outputDir+"/PhotonJet_2ndLevel_"+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+"_Part2.root" 
  print filename_out1
  print filename_out2
  print filename_out
  os.system("hadd "+filename_out1+" "+files1 )
  os.system("hadd "+filename_out2+"  "+files2 )
  os.system("hadd "+filename_out+"  "+filename_out1+" "+filename_out2 )

#  command1 = "hadd "+filename_out1+" "+files1 
#  command2 = "hadd "+filename_out2+" "+files2
#  command3 = "hadd "+filename_out+"  "+filename_out1+" "+filename_out2
#  logfile = "logfile_Merging.log"
#  outputname = "submit_Merging.src"
#  outputfile = open(outputname,'w')
#  outputfile.write('#!/bin/bash\n')
#  outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
#  outputfile.write('cd '+pwd+'\n')
#  outputfile.write('eval `scramv1 runtime -sh`\n')
#  outputfile.write(command1+"\n")
#  outputfile.write(command2+"\n")
#  outputfile.write(command3+"\n")
#  print outputname 
#  print logfile 
#  os.system("bsub -q cmslong -o "+logfile+" source "+pwd+"/"+outputname)

##### update total luminosity ######
inputFile = TFile(filename_out,"UPDATE")
lumi = inputFile.Get("totallumi")
lumi.SetVal(lumi_tot) # in /pb
inputFile.cd("rootTupleTree")
lumi.Write()

