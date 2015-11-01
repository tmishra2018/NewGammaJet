#! /bin/env python

import argparse, os, tempfile, shutil, sys
from subprocess import call, PIPE, STDOUT, Popen

# JSON files for Winter13 rereco. 
jsonFiles = [
    ["Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt", [256630, 258750]],
 ]

print(jsonFiles)

# Trigger list for 2015
triggers = {
    tuple([256630, 258750]) : ["HLT_Photon30_R9Id90_HE10_IsoM_v.*",
                        "HLT_Photon50_R9Id90_HE10_IsoM_v.*",
                        "HLT_Photon75_R9Id90_HE10_IsoM_v.*",
                        "HLT_Photon90_R9Id90_HE10_IsoM_v.*",
                        "HLT_Photon120_R9Id90_HE10_IsoM_v.*",
                        "HLT_Photon165_R9Id90_HE10_IsoM_v.*"]
    }


#def parallelize(cmd):
#  print("Launching parallel for %r" % cmd)
#  p_cmd = ["parallel", "-j", "20", "-u"]
#  p_cmd.extend(cmd)
#  p = Popen(p_cmd, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
#  out = p.communicate(input = "\n".join(validJSONFiles))[0]
#  print(out)


parser = argparse.ArgumentParser(description='Generate PU profile for each prescaled triggers')
#parser.add_argument('crab_json', nargs=1)
parser.add_argument('lumi_json', nargs=1)
args = parser.parse_args()

lumi_json = args.lumi_json[0]

if not os.path.exists(lumi_json):
  print("Error: %s not found" % lumi_json)
  exit(1)

#tempFolder = tempfile.mkdtemp(dir="scratch")
#print(tempFolder)

for trigger_runs, triggerslist in triggers.items():
  # Get list of JSON for this run range
  validJSONFiles = [];
  for json in jsonFiles:
    json_name = json[0]
    json_runs = json[1]

    print(json_name)
    print(json_runs)

    if (json_runs[0] >= trigger_runs[0] and json_runs[0] <= trigger_runs[1]) or (json_runs[1] >= trigger_runs[0] and json_runs[1] <= trigger_runs[1]):
      validJSONFiles.append(json[0])

  for trigger in triggerslist:
    print("Generate PU profile for '%s'" % trigger)

    tempFilename = "%s_%d_%d" % (trigger.rstrip("_*"), trigger_runs[0], trigger_runs[1])
#    outputCSV = os.path.join(tempFolder, "{.}_%s.csv" % tempFilename)
#    outputJSON = os.path.join(tempFolder, "{.}_%s_JSON.txt" % tempFilename)
    outputROOT = "pu_truth_data_photon_2015_true_%s_50bins.root" % tempFilename

#    print("\tRunning lumiCalc2...")

#    parallelize(["pixelLumiCalc.py", "--begin", str(trigger_runs[0]), "--end", str(trigger_runs[1]), "lumibyls", "-i", "{}", "--hltpath", trigger, "-o", outputCSV])
    #parallelize(["echo", "lumibyls", "-i", "{}", "--hltpath", trigger, "-o", outputCSV])

#    print("\tRunning pileupReCalc_HLTpaths.py...")

 #   parallelize(["pileupReCalc_HLTpaths.py", "-i", outputCSV, "--inputLumiJSON", lumi_json, "-o", outputJSON])

    cmd= "pileupCalc.py -i ../../json/"+json_name+" --inputLumiJSON "+lumi_json+" --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 "+outputROOT

    os.system(cmd)

    print("\tRunning pileupCalc...")
 

  
#    parallelize(["pileupCalc.py", "-i", json, "--inputLumiJSON", lumi_json, "--calcMode", "true", "--minBiasXsec", "69000", "--maxPileupBin", "50", "--numPileupBins", "50", outputROOT])

    print

 # Merge everything
toMerge = {}
for trigger_runs, triggerslist in triggers.items():
  for json in jsonFiles:
    json_name = os.path.splitext(json[0])[0]

    for trigger in triggerslist:
      root_file = "pu_truth_data_photon_2012_true_%s_%s_%d_%d_50bins.root" % (json_name, trigger.rstrip("_*"), trigger_runs[0], trigger_runs[1])
      if os.path.exists(root_file):
        toMerge.setdefault(trigger.rstrip("_*"), []).append(root_file)

for trigger, files in toMerge.items():
  output_file = "pu_truth_data_photon_2012_true_%s_50bins.root" % trigger
  input_files = " ".join(files)
  os.system("hadd %s %s" % (output_file, input_files))

#shutil.rmtree(tempFolder)
