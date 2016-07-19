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

parser.add_option("-i", "--inputList", dest="inputlist",
                  help="list of samples to be analyzed",
                  )

parser.add_option("-v",  type=str, dest="tagname", default="TAG",
                  help="template of cmssw file to execute",
                  )

parser.add_option('--all',      
                  action='store_true',
                  dest='AllAlphaCut',      
                  help='draw All AlphaCut',
                  default=False)

parser.add_option('--mc',      
                  action='store_true',
                  dest='IsMC',      
                  help='run on MonteCarlo',
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
      if opt.IsMC:
        command = "gammaJetFinalizer -i ../"+opt.inputlist+" -d "+opt.tagname+"_alphacut0"+str( int(alphacut*100) )+" --type "+jetType+" --algo ak4 --alpha "+str(alphacut)+" --mc"
      else:
        command = "gammaJetFinalizer -i ../"+opt.inputlist+" -d "+opt.tagname+"_alphacut0"+str( int(alphacut*100) )+" --type "+jetType+" --algo ak4 --alpha "+str(alphacut)
      print "submit "+command
      print ""
      logfile = "logfile_"+opt.tagname+"_"+jetType+"_alphacut0"+str( int(alphacut*100))+".log"
      outputname = "batch/submit_"+opt.tagname+"_"+jetType+"_alphacut0"+str( int(alphacut*100) )+".src"
      outputfile = open(outputname,'w')
      outputfile.write('#!/bin/bash\n')
      outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
      outputfile.write('cd '+pwd+'/AlphaCut0'+str(int(alphacut*100))+' \n')
      outputfile.write('eval `scramv1 runtime -sh`\n')
      outputfile.write(command+"\n")
      print outputname 
      print logfile 
      os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
  else:
    alphacut = 0.3
    os.system("mkdir -p AlphaCut0"+str(int(alphacut*100) ) )
    if opt.IsMC:
      command = "gammaJetFinalizer -i ../"+opt.inputlist+" -d "+opt.tagname+"_alphacut0"+str( int(alphacut*100) )+" --type "+jetType+" --algo ak4 --alpha "+str(alphacut)+" --mc"
    else:
      command = "gammaJetFinalizer -i ../"+opt.inputlist+" -d "+opt.tagname+"_alphacut0"+str( int(alphacut*100) )+" --type "+jetType+" --algo ak4 --alpha "+str(alphacut)
    print "submit "+command
    print ""
    logfile = "logfile_"+opt.tagname+"_"+jetType+"_alphacut0"+str( int(alphacut*100))+".log"
    outputname = "batch/submit_"+opt.tagname+"_"+jetType+"_alphacut0"+str( int(alphacut*100) )+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
    outputfile.write('cd '+pwd+'/AlphaCut0'+str(int(alphacut*100))+' \n')
    outputfile.write('eval `scramv1 runtime -sh`\n')
    outputfile.write(command+"\n")
    print outputname 
    print logfile 
    os.system("bsub -q "+opt.queue+" -o batch/"+logfile+" source "+pwd+"/"+outputname)
