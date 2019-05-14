#!/usr/bin/env python

import os
import sys
import argparse
import math
from ROOT import *
from array import array
import datetime

usage = "usage:  python mergeData.py --inputList list_of_files_to_merge.txt -o ./  --lumi_tot 1000. [in /pb] --run run era"

parser = argparse.ArgumentParser(description='Process options.')

parser.add_argument("-c", "--cleaning", type=bool, dest="cleaning", default="False",
        help="remove corrupted files from list or not",
	    )
parser.add_argument("-i", "--inputList", type=str, dest="inputList", default="",
        help="input list of files to be merged",
	    )
parser.add_argument("-o", "--outputDir", type=str, dest="outputDir", default="./",
        help="output directory",
	    )
parser.add_argument("--lumi_tot", type=float, dest="lumi_tot", default=1,
        help="total luminosity",
	    )
parser.add_argument("--run", type=str, dest="run", default=1,
        help="run era",
	    )

args = parser.parse_args()
print args 
inputList = args.inputList
outputDir = args.outputDir
#xsec = args.xsec
lumi_tot = args.lumi_tot
run = args.run
cleaning = args.cleaning
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
   # pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
   # fullname = os.path.split(pathT2)[1]
   # name = fullname.split("_")
    i= i+1
  if i >= 550:
    file = line.strip()
    files2 += str(" "+file)
   # pathT2 = file.split("dcap://cmsrm-se01.roma1.infn.it/")[1]
   # fullname = os.path.split(pathT2)[1]
   # name = fullname.split("_")
    i= i+1

print i
print files1
print files2

today = datetime.date.today()
today.strftime('%d-%m-%Y')

pwd = os.environ['PWD']

bad_files = []

ignored_cleaning_runs = ['ABC', 'ABCD']

if i < 550:
  filename_out = outputDir+"/PhotonJet_2ndLevel_DATA_RUN_"+run+"_"+str(today)+".root" #+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+".root" 
  print filename_out
  exitcode = os.system("hadd -f "+filename_out+"  "+files1 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
  while exitcode > 0 and cleaning and not (run in ignored_cleaning_runs):
    print 'Error while merging, removing a file...'
    bad_file = os.popen("cat tmp_py_merge_files_from_step2_{}.out".format(inputList)).read()[:-1].split('\n')[-1].split(':')[1][1:]
    print bad_file
    bad_files.append(bad_file)
    if bad_file in files1:
      files1 = files1.replace(bad_file,'')
    exitcode = os.system("hadd -f "+filename_out+"  "+files1 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
else:
  filename_out = outputDir+"/PhotonJet_2ndLevel_DATA_RUN_"+run+"_"+str(today)+".root"#+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+".root" 
  filename_out1 = outputDir+"/PhotonJet_2ndLevel_DATA_RUN_"+run+"_"+str(today)+"_Part1.root"#+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+"_Part1.root" 
  filename_out2 = outputDir+"/PhotonJet_2ndLevel_DATA_RUN_"+run+"_"+str(today)+"_Part2.root"#+name[0]+"_"+name[1]+"_"+name[2]+"_"+name[3]+"_"+name[4]+"_"+str(today)+"_Part2.root" 
  print filename_out1
  print filename_out2
  print filename_out
  exitcode = os.system("hadd -f "+filename_out1+" "+files1 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
  while exitcode > 0 and cleaning and not (run in ignored_cleaning_runs):
    print 'Error while merging, removing a file...'
    bad_file = os.popen("cat tmp_py_merge_files_from_step2_{}.out".format(inputList)).read()[:-1].split('\n')[-1].split(':')[1][1:]
    print bad_file
    bad_files.append(bad_file)
    if bad_file in files1:
      files1 = files1.replace(bad_file,'')
    exitcode = os.system("hadd -f "+filename_out1+" "+files1 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
  exitcode = os.system("hadd -f "+filename_out2+" "+files2 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
  while exitcode > 0 and cleaning and not (run in ignored_cleaning_runs):
    print 'Error while merging, removing a file...'
    bad_file = os.popen("cat tmp_py_merge_files_from_step2_{}.out".format(inputList)).read()[:-1].split('\n')[-1].split(':')[1][1:]
    print bad_file
    bad_files.append(bad_file)
    if bad_file in files2:
      files2 = files2.replace(bad_file,'')
    exitcode = os.system("hadd -f "+filename_out2+" "+files2 +" > tmp_py_merge_files_from_step2_{}.out".format(inputList))
  os.system("hadd -f "+filename_out+"  "+filename_out1+" "+filename_out2 )

all_files = os.popen("cat "+args.inputList).read()[:-1].split('\n')
good_files = [f for f in all_files if f not in bad_files]

if cleaning and not (run in ignored_cleaning_runs):
  os.system("rm "+args.inputList)
  for files in good_files:
    os.system("echo "+files+" >> "+args.inputList)

os.system("rm tmp_py_merge_files_from_step2_{}.out".format(inputList))
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

print lumi.GetVal()
lumi.Write()
analysis_tree = inputFile.Get("rootTupleTree/tree")

evtWeightTotA = array("f", [0.] )
b_evtWeightTot = analysis_tree.Branch("evtWeightTotA", evtWeightTotA,"evtWeightTotA/F")

PassGenmatching = array("f", [0.] )
b_PassGenmatching = analysis_tree.Branch("PassGenmatching", PassGenmatching,"PassGenmatching/D")

for event in analysis_tree:
  evtWeightTotA[0] = 0
  PassGenmatching[0] = 0
  b_evtWeightTot.Fill()
  b_PassGenmatching.Fill()
inputFile.cd("rootTupleTree")
analysis_tree.Write("",TObject.kOverwrite)


