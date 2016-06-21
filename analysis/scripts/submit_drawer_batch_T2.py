#! /usr/bin/env python

import os
import sys
import optparse
import re

usage = "usage: To be run from directory containing the root file:  python scripts/submit_batch_T2.py -i 2ndLevel.root -v TagName"

parser = optparse.OptionParser("submitAllGJetsID.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
                  help='run in batch in queue specified as option (default -q cmslong)', 
                  default='cmslong',
                  metavar="QUEUE")

parser.add_option("-d",  type=str, dest="DataSample", default="",
                  help="template of cmssw file to execute",
                  )

parser.add_option("-s",  type=str, dest="MCSample", default="",
                  help="template of cmssw file to execute",
                  )

parser.add_option('--all',      
                  action='store_true',
                  dest='AllAlphaCut',      
                  help='draw All AlphaCut',
                  default=False)

(opt, args) = parser.parse_args()
################################################

os.system("mkdir -p batch")
pwd = os.environ['PWD']

for j in range (0,2):
  if j == 0:
    jetType ="pf"
  if j == 1:
    jetType ="puppi"
  print jetType
  if opt.AllAlphaCut==True:
    for i in range (0,4):
      if i == 0:
        alphacut = 0.1
      if i == 1:
        alphacut = 0.15
      if i == 2:
        alphacut = 0.2
      if i == 3:
        alphacut = 0.3
      os.system("mkdir -p AlphaCut0"+str(int(alphacut*100) ) )
      dataSample = opt.DataSample+"_alphacut0"+str(int(alphacut*100))
      mcSample   = opt.MCSample+"_alphacut0"+str(int(alphacut*100))
      command = "../../../draw/drawPhotonJet_2bkg "+dataSample+" "+mcSample+" "+mcSample+" "+jetType+" ak4 LUMI"
      command2 = "../../../draw/drawPhotonJetExtrap --type "+jetType+" --algo ak4 "+dataSample+" "+mcSample+" "+mcSample
      command3 = "../../../draw/draw_ratios_vs_pt "+dataSample+" "+mcSample+" "+mcSample+" "+jetType+" ak4"
      print "submit "+command
      print "submit "+command2
      print "submit "+command3
      print ""
      logfile = "logfile_drawer"+"_"+jetType+"_alphacut0"+str( int(alphacut*100))+".log"
      outputname = "batch/submit_drawer_"+jetType+"_alphacut0"+str( int(alphacut*100) )+".src"
      outputfile = open(outputname,'w')
      outputfile.write('#!/bin/bash\n')
      outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
      outputfile.write('cd '+pwd+'/AlphaCut0'+str(int(alphacut*100))+' \n')
      outputfile.write('eval `scramv1 runtime -sh`\n')
      outputfile.write(command+"\n")
      outputfile.write(command2+"\n")
      outputfile.write(command3+"\n")
      print outputname 
      print logfile 
      os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
  else:
    alphacut = 0.3
    os.system("mkdir -p AlphaCut0"+str(int(alphacut*100) ) )
    dataSample = opt.DataSample+"_alphacut0"+str(int(alphacut*100))
    mcSample   = opt.MCSample+"_alphacut0"+str(int(alphacut*100))
    command = "../../../draw/drawPhotonJet_2bkg "+dataSample+" "+mcSample+" "+mcSample+" "+jetType+" ak4 LUMI"
    command2 = "../../../draw/drawPhotonJetExtrap --type "+jetType+" --algo ak4 "+dataSample+" "+mcSample+" "+mcSample
    command3 = "../../../draw/draw_ratios_vs_pt "+dataSample+" "+mcSample+" "+mcSample+" "+jetType+" ak4"
    print "submit "+command
    print "submit "+command2
    print "submit "+command3
    print ""
    logfile = "logfile_drawer"+"_"+jetType+"_alphacut0"+str( int(alphacut*100))+".log"
    outputname = "batch/submit_drawer_"+jetType+"_alphacut0"+str( int(alphacut*100) )+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
    outputfile.write('cd '+pwd+'/AlphaCut0'+str(int(alphacut*100))+' \n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write(command+"\n")
    outputfile.write(command2+"\n")
    outputfile.write(command3+"\n")
    print outputname 
    print logfile 
    os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
    
