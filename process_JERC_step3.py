import os
import subprocess
from functools import partial
import multiprocessing as mp

from optparse import OptionParser
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-d", "--dates", dest = "dates", default=None,
                  help = "Give the dates of the Step2 outputs, by default it takes the most recent. Format is -d YYY-MM-DD,YYY-MM-DD,...")

(options,args) = parser.parse_args()

samples = {
    'B' : 'Run2017B-09Aug2019_UL2017-v1',
    'C' : 'Run2017C-09Aug2019_UL2017-v1',
    'D' : 'Run2017D-09Aug2019_UL2017-v1',
    'E' : 'Run2017E-09Aug2019_UL2017-v1',
    'F' : 'Run2017F-09Aug2019_UL2017-v1',
    'MC1': 'GJets_HT-40To100_RunIISummer19MiniAOD-106X',
    'MC2': 'GJets_HT-100To200_RunIISummer19MiniAOD-106X',
    'MC3': 'GJets_HT-200To400_RunIISummer19MiniAOD-106X',
    'MC4': 'GJets_HT-400To600_RunIISummer19MiniAOD-106X',
    'MC5': 'GJets_HT-600ToInf_RunIISummer19MiniAOD-106X',
}

merged_runs = ["BC", "DE", "BCDEF"]

lumis_per_pb = {
    'Run2017B-09Aug2019_UL2017-v1' : 4823,
    'Run2017C-09Aug2019_UL2017-v1' : 9664,
    'Run2017D-09Aug2019_UL2017-v1' : 4252,
    'Run2017E-09Aug2019_UL2017-v1' : 9278,
    'Run2017F-09Aug2019_UL2017-v1' : 13540,
}
# recomputed with brilcalc :
lumis_per_pb = {
    'Run2017B-09Aug2019_UL2017-v1' : 4793.961426839,
    'Run2017C-09Aug2019_UL2017-v1' : 9631.214820913,
    'Run2017D-09Aug2019_UL2017-v1' : 4247.682053046,
    'Run2017E-09Aug2019_UL2017-v1' : 9313.642401775,
    'Run2017F-09Aug2019_UL2017-v1' : 13539.378417564,
}

xsec_pb = {
    'GJets_HT-40To100_RunIISummer19MiniAOD-106X' : 18700.0,
    'GJets_HT-100To200_RunIISummer19MiniAOD-106X' : 8640.0,
    'GJets_HT-200To400_RunIISummer19MiniAOD-106X' : 2185.0,
    'GJets_HT-400To600_RunIISummer19MiniAOD-106X' : 259.9,
    'GJets_HT-600ToInf_RunIISummer19MiniAOD-106X': 85.31,
}

L1Offset_approaches = ["ComplexL1", "SimpleL1"]
JECs = ['wo_L2Res', 'only_L2Res', 'L2L3Res']
JERs = ['JER']

Step1_outputs_crab_dir = "/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/CMSSW_10_6_3/src/CMSDIJET/DijetRootTreeMaker/prod/crab/CERN_crab/prod_2020-03-24/"
Step2_outputs_base_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/lists_2017UL/"#+dirlist+'/'+datetime.date.today().isoformat()+'/'
Step3_outputs_base_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2017UL/"

##############################################################################
###################### Do not change the following lines #####################
##############################################################################

Step3_dir = "$CMSSW_BASE/src/JetMETCorrections/GammaJetFilter/"
Step3_scripts_dir = Step3_dir+"/analysis/scripts/"
Step3_PU_dir = Step3_dir+"/analysis/PUReweighting/"

def print_log_started_step(run, JERC, step):
    command = "echo 'run "+run+' '+JERC+" "+step+" started' >> "+Step3_dir+'/tmp-process_logs.log'
    os.system(command)

def get_most_recent(directory):
    return os.popen('ls -t {}'.format(directory)).read().split('\n')[0]

def find_files_from_step2(run, JERC):
    print_log_started_step(run, JERC, 'find files')

    Step2_files_dir = "/".join([Step2_outputs_base_dir, samples[run]])
    grep_path = "{}".format(JERC)

    if options.dates is None:
        Step2_files_dir = "/".join([Step2_files_dir, get_most_recent(Step2_files_dir)])
    else:
        grep_path += " | grep -e "+" -e ".join(options.dates.split(","))

    list_file = "Step2_files_{}_{}.txt".format(run,JERC)

    command = "cd {Step3_scripts_dir} && find {Step2_files_dir} -type f -iname \*reduced_skim.root | grep {grep_path} > {list_file}".format(
        Step3_scripts_dir = Step3_scripts_dir,
        Step2_files_dir = Step2_files_dir,
        grep_path = grep_path,
        list_file = list_file,
    )
    os.system(command)
        
def merge_files_from_step2(run, JERC, cleaning=True):
    print_log_started_step(run, JERC, 'merge files')

    output_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2017UL/"+JERC+'/'
    command_mkdir = "mkdir -p {}".format(output_dir)

    command_cd = 'cd {}'.format(Step3_scripts_dir)

    input_list = "Step2_files_{}_{}.txt".format(run,JERC)

    if samples[run] in lumis_per_pb.keys():
        script = "mergeData.py"
        extras = "--lumi_tot {} --run {}".format(lumis_per_pb[samples[run]], run)
    elif samples[run] in xsec_pb.keys():
        script = "mergeAndAddWeights.py"
        extras = "--xsec {} --JEC {}".format(xsec_pb[samples[run]], "_".join([samples[run], JERC]))

    command = " && ".join([
        command_mkdir, command_cd,
        "python {script} -i {input_list} -o {output_dir} {extras}".format(
            script = script,
            input_list = input_list,
            output_dir = output_dir,
            extras = extras,
        )
    ])
    print command
    os.system(command)

def make_pileup_data(run,JERC):
    print_log_started_step(run, JERC, 'make pileup')
    cmds = []
    cmds.append('cd {}'.format(Step3_PU_dir))

    if samples[run] in lumis_per_pb.keys(): # data
        cmds.append("python ComputePU_perHLT.py --Cert_json {Cert_json} --whichrun {whichrun}".format(
            Cert_json = "{}/crab_{}/results/processedLumis.json".format(Step1_outputs_crab_dir, samples[run]),
            whichrun = "_".join([run, JERC]),
        )
    )
    command = " && ".join(cmds)
    print command
    os.system(command)

def make_pileup_data_merged(run,JERC):
    print_log_started_step(run, JERC, 'make pileup')
    cmds = []
    cmds.append('cd {}'.format(Step3_PU_dir))

    processedLumis = ''
    for letter in run :
        processedLumis += os.popen('cat '+Step1_outputs_crab_dir+'/crab_'+samples[letter]+'/results/processedLumis.json').read()[:-1]
    processedLumis = processedLumis.replace('}{',', ')
    cmds.append("echo '{}' > ./tmp_processedLumis_{}_{}.json".format(processedLumis, run, JERC))
    cmds.append(
        "python ComputePU_perHLT.py --Cert_json {Cert_json} --whichrun {whichrun}".format(
            Cert_json = "./tmp_processedLumis_{}_{}.json".format(run, JERC),
            whichrun = "_".join([run, JERC])
        )
    )

    command = " && ".join(cmds)
    print command
    #os.system(command)

def make_pileup_MC(JERC):
    run = "MC"
    print_log_started_step(run, JERC, 'make pileup')
    files_list = "files_PUreweight_2017UL_{}_{}".format(run, JERC)
    cmds = []
    cmds.append('cd {}'.format(Step3_PU_dir))
    cmds.append("rm -f {list}".format(list = files_list))
    for mc_sample in xsec_pb.keys():
        cmds.append(
            "echo {merged_file} >> {list}.list".format(
                merged_file = get_most_recent(
                    "{}/{}/*_MC_{}*".format(Step3_outputs_base_dir, JERC, "_".join([mc_sample, JERC]))
                    ),
                list = files_list,
            )
        )
    cmds.append(
        "./generate_mc_pileup.exe {list}".format(list = files_list)
    )
    command = " && ".join(cmds)
    print command
    os.system(command)

def Produce_Combination_File_and_plots(run,JERC):
    print_log_started_step(run, JERC, 'make plots')

    os.system("cd {} ; ln -sf computed_mc_files_PUreweight_2017UL_MC_{}_pu_truth_100bins.root computed_mc_files_2017UL_MC_for_{}_{}_pu_truth_100bins.root".format(
        Step3_PU_dir,
        JERC,
        run,
        JERC,
    )
    )

    output = JERC
    IOV = run+'_'+JERC
    Pu_profile = run+'_'+JERC
    if len(run) == 1:
        script = "Produce_Combination_File_and_plots-MC_list.py"
        Input_data = get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_DATA_RUN_'+run+'_*')
    else:
        script ="Produce_Combination_File_and_plots-DATA_and_MC_list.py"
        Input_data = [get_most_recent("{}/{}/PhotonJet_2ndLevel_DATA_RUN_merged_*_for_{}_*".format(Step3_outputs_base_dir,JERC,mrun))]
        for r in run[1:]:
            Input_data.append(get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_DATA_RUN_'+run+'_*'))
        os.system("cd {} && rm -f input_data_{}_{}.list".format(Step3_dir, run, JERC))
        for file in Input_data:
            os.system("cd {} && echo '{}' >> input_data_{}_{}.list".format(Step3_dir, file, run, JERC))
        Input_data = '{}/input_data_{}_{}.list'.format(Step3_dir, run, JERC)
    Input_mc = []
    for mc_sample in [m for m in samples.keys() if 'MC' in m]:
        Input_mc.append(get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_MC_{}*{}*'.format(samples[mc_sample], JERC)))
    os.system("cd {} && rm -f input_mc_{}_{}.list".format(Step3_dir, run, JERC))
    for file in Input_mc:
        os.system("cd {} && echo '{}' >> input_mc_{}_{}.list".format(Step3_dir, file, run, JERC))
    command = 'cd '+Step3_dir
    command += ' ; python {}'.format(script)
    command += ' --output='+output
    command += ' --IOV='+IOV
    command += ' --Pu_profile='+Pu_profile
    command += ' --Input_data='+Input_data
    command += ' --Input_mc='+'{}/input_mc_{}_{}.list'.format(Step3_dir, run, JERC)
    print command
    #os.system(command)

def get_copy_for_merge(mrun, JERC):
    run_to_copy = mrun[0]
    cmds = []
    cmds.append("cd {}/{}".format(Step3_outputs_base_dir,JERC))
    cmds.append(
        "cp {} {}".format(
            get_most_recent("PhotonJet_2ndLevel_DATA_RUN_{}_*".format(run_to_copy)),
            get_most_recent(
                "PhotonJet_2ndLevel_DATA_RUN_{}_*".format(run_to_copy)
            ).replace(run_to_copy, "merged_{}_for_{}".format(run_to_copy, mrun))
        )
    )
    command = " && ".join(cmds)
    print(command)

def set_lumi_for_merge(mrun, JERC):
    lumi = 0
    for run in mrun:
        lumi += lumis_per_pb[samples[run]]
    file = get_most_recent("{}/{}/PhotonJet_2ndLevel_DATA_RUN_*_for_merged_{}_*".format(Step3_outputs_base_dir,JERC,mrun))
    set_lumi(file, lumi)

def set_lumi(file, lumi):
    from ROOT import TFile
    tfile = TFile(file, "UPDATE")
    _lumi = tfile.Get("totallumi")
    _lumi.SetVal(lumi)
    _lumi.Write()
    tfile.Close()

##############################################################################
##############################################################################
##############################################################################

CMSSW_BASE = os.popen('echo $CMSSW_BASE').read()[:-1]
Step3_BASE_JEC = os.popen('echo /afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/$CMSSW_VERSION').read()[:-1]
Step3_BASE_JER = os.popen('echo /afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/JER_CMSSWs/$CMSSW_VERSION').read()[:-1]

if CMSSW_BASE == Step3_BASE_JEC:
    JERCs = JECs
    print "Processing JECs:"
    print JERCs
elif CMSSW_BASE == Step3_BASE_JER:
    JERCs = JERs
    print "Processing JERs:"
    print JERCs
else:
    print "Ouch, the directories are not as they are supposed to be. Did you followed the git installation instructions correctly?"

os.system('rm -f '+Step3_dir+'/tmp-process_logs.log')
os.system('touch '+Step3_dir+'/tmp-process_logs.log')

run_JERCs_MC = []
run_JERCs_data = []
MC_pileup_JERCs = set()

for L1Offset_approach in L1Offset_approaches:
    for JERC in JERCs:
        L1Offset_JERC = '_'.join([L1Offset_approach, JERC])
        for run in samples.keys():
            if samples[run] in lumis_per_pb.keys():
                run_JERCs_data.append((run, L1Offset_JERC))
            elif samples[run] in xsec_pb.keys():
                run_JERCs_MC.append((run, L1Offset_JERC))

run_JERCs = run_JERCs_MC + run_JERCs_data

print run_JERCs

os.system("sleep 10s")
for run_JERC in run_JERCs_MC:
    if run_JERC in run_JERCs:
        run, JERC = run_JERC
        find_files_from_step2(run, JERC)
        merge_files_from_step2(run, JERC)
        MC_pileup_JERCs.add(JERC)

for MC_pileup_JERC in MC_pileup_JERCs:
    make_pileup_MC(MC_pileup_JERC)

for run_JERC in run_JERCs_data:
    if run_JERC in run_JERCs:
        run, JERC = run_JERC

        find_files_from_step2(run, JERC)
        merge_files_from_step2(run, JERC)

        make_pileup_data(run,JERC)

        Produce_Combination_File_and_plots(run,JERC)

for merged_run in merged_runs:
    done_jercs = {}
    for mrun in merged_run:
        done_jercs[mrun] = set()
        for run_JERC in run_JERCs:
            run, JERC = run_JERC
            if mrun == run:
                done_jercs[mrun].add(JERC)
    todo_jercs = set([jerc for jerc in done_jercs[mrun]])
    for mrun in merged_run:
        todo_jercs = todo_jercs.intersection(done_jercs[mrun])

    for JERC in todo_jercs:
        get_copy_for_merge(mrun, JERC)
        set_lumi_for_merge(mrun, JERC)

        make_pileup_data_merged(mrun,JERC)

        Produce_Combination_File_and_plots(mrun,JERC)

os.system('rm -f '+Step3_dir+'/tmp-process_logs.log')
