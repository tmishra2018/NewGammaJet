#! /usr/bin/env python

import os
import sys
import optparse
import re

usage = "usage: To be run from directory containing the root file:  python scripts/submit_batch_T2.py -i list.txt"

parser = optparse.OptionParser("submitAllGJetsID.py")
parser.add_option('-q', '--queue',       action='store',     dest='queue',       
                  help='run in batch in queue specified as option (default -q cmslong)', 
                  default='cmslong',
                  metavar="QUEUE")

parser.add_option("-i", "--inputList", dest="inputlist",
                  help="list of samples to be analyzed",
                  )

(opt, args) = parser.parse_args()
################################################

os.system("mkdir -p batch")
pwd = os.environ['PWD']
  
command = "../../scripts/mergeAndAddWeights.py -i "+opt.inputlist+" -o ."
print "submit "+command
print ""
logfile = "logfile_Merging.log"
outputname = "batch/submit_Merging.src"
outputfile = open(outputname,'w')
outputfile.write('#!/bin/bash\n')
outputfile.write('export SCRAM_ARCH=slc6_amd64_gcc530\n')
outputfile.write('cd '+pwd+' \n')
outputfile.write('eval `scramv1 runtime -sh`\n')
outputfile.write(command+"\n")
print outputname 
print logfile 
os.system("bsub -q cmslong -o batch/"+logfile+" source "+pwd+"/"+outputname)
