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
all_files = open(args.inputList,"r").readlines()
all_files = [path[:-1] for path in all_files]
splitted_files = []
splitting_files_by = 500

if len(all_files) <= splitting_files_by:
  splitted_files = [all_files]
else:
  for k in range(int(1.*len(all_files)/splitting_files_by)):
    splitted_files.append(all_files[int(k*splitting_files_by):int((k+1)*splitting_files_by)])
  if not int((k+1)*splitting_files_by) == len(all_files)-1:
    k+=1
  splitted_files.append(all_files[int(k*splitting_files_by):])

print "Will merge {} files, splitting by groups of {} so will use {} part{}.".format(
  len(all_files),
  splitting_files_by,
  len(splitted_files),
  "s" if len(splitted_files) > 1 else ""
)

today = datetime.date.today()
today.strftime('%d-%m-%Y')

pwd = os.environ['PWD']

bad_files = []

filename_out = outputDir+"/PhotonJet_2ndLevel_DATA_RUN_"+run+"_"+str(today)+".root"
filename_out_parts = []

for k in range(len(splitted_files)):
  splitted_files[k] = ' '.join(splitted_files[k])

for k in range(len(splitted_files)):
  part_num = k+1
  print "Merging part {}.".format(part_num)
  filename_out_part = "{}/PhotonJet_2ndLevel_DATA_RUN_{}_{}_Part{}.root".format(
    outputDir,
    run,
    str(today),
    part_num)
  command = " ".join(
    ["hadd -f {}".format(filename_out_part), splitted_files[k], " > tmp_py_merge_files_from_step2_{}.out".format(inputList)])
  exitcode = os.system(command)
  bad_file = os.popen("cat tmp_py_merge_files_from_step2_{}.out".format(inputList)).read()[:-1].split('\n')[-1].split(':')[1][1:]
  while exitcode > 0 and cleaning and bad_file in splitted_files[k]:
    print 'Error while merging, removing a file...'
    print bad_file
    bad_files.append(bad_file)
    splitted_files[k] = splitted_files[k].replace(bad_file,'')
    command = " ".join(
      ["hadd -f {}".format(filename_out_part), splitted_files[k], " > tmp_py_merge_files_from_step2_{}.out".format(inputList)])
    exitcode = os.system(command)
    bad_file = os.popen("cat tmp_py_merge_files_from_step2_{}.out".format(inputList)).read()[:-1].split('\n')[-1].split(':')[1][1:]
  filename_out_parts.append(filename_out_part)
if k > 0:
  print "Merging parts together..."
  os.system(
    " ".join(["hadd -f {}".format(filename_out)]+filename_out_parts)
  )
  for file_part in filename_out_parts:
    os.system("rm -f {} &".format(file_part))
else :
  os.system("mv {} {}".format(filename_out_parts[0], filename_out))

all_files = os.popen("cat "+args.inputList).read()[:-1].split('\n')
good_files = [f for f in all_files if f not in bad_files]

if cleaning :
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


