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
    'A' : 'Run2018A-17Sep2018-v2',
    'B' : 'Run2018B-17Sep2018-v1',
    'C' : 'Run2018C-17Sep2018-v1',
    'D' : 'Run2018D-PromptReco-v2',
    'MC': 'GJet_Pt-15To6000_RunIIAutumn18MiniAOD-102X',
}

merged_runs = ["ABC", "ABCD"]

lumis_per_pb = {}
# recomputed with brilcalc :
lumis_per_pb = {
    'Run2018A-17Sep2018-v2' : 13654.355526985,
    'Run2018B-17Sep2018-v1' : 7057.825158567,
    'Run2018C-17Sep2018-v1' : 6894.770971269,
    'Run2018D-PromptReco-v2' : 31066.589629726,
}

xsec_pb = {
    'GJet_Pt-15To6000_RunIIAutumn18MiniAOD-102X' : 283000.0,
}

JECs = ['wo_L2Res', 'only_L2Res', 'L2L3Res']
JERs = ['JER']

Step1_outputs_crab_dir = "/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/CMSSW_10_2_5/src/CMSDIJET/DijetRootTreeMaker/prod/crab/CERN_crab/old_crab_jobs/2019-02-28/"
Step2_outputs_base_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/lists_2018/"#+dirlist+'/'+datetime.date.today().isoformat()+'/'
Step3_outputs_base_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2018/"

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

    output_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2018/"+JERC+'/'
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
    os.system(command)

def make_pileup_MC(JERC):
    run = "MC"
    print_log_started_step(run, JERC, 'make pileup')
    files_list = "files_PUreweight_2018_{}_{}".format(run, JERC)
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

    os.system("cd {} ; ln -sf computed_mc_files_PUreweight_2018_MC_{}_pu_truth_100bins.root computed_mc_files_2018_MC_for_{}_{}_pu_truth_100bins.root".format(
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
        Input_data = [get_most_recent("{}/{}/PhotonJet_2ndLevel_DATA_RUN_merged_{}_for_{}_*".format(Step3_outputs_base_dir,JERC,run[0],run))]
        for r in run[1:]:
            Input_data.append(get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_DATA_RUN_'+r+'_*'))
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
    os.system(command)

def get_copy_for_merge(mrun, JERC):
    run_to_copy = mrun[0]
    command = "cp {} {}".format(
        get_most_recent("{}/{}/PhotonJet_2ndLevel_DATA_RUN_{}_*".format(Step3_outputs_base_dir,JERC,run_to_copy)),
        get_most_recent(
            "{}/{}/PhotonJet_2ndLevel_DATA_RUN_{}_*".format(Step3_outputs_base_dir,JERC,run_to_copy)
        ).replace("_RUN_{}_".format(run_to_copy), "_RUN_merged_{}_for_{}_".format(run_to_copy, mrun))
    )
    print(command)
    os.system(command)

def set_lumi_for_merge(mrun, JERC):
    lumi = 0
    for run in mrun:
        lumi += lumis_per_pb[samples[run]]
    file = get_most_recent("{}/{}/PhotonJet_2ndLevel_DATA_RUN_merged_*_for_{}_*".format(Step3_outputs_base_dir,JERC,mrun))
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

for JERC in JERCs:
    for run in samples.keys():
        if samples[run] in lumis_per_pb.keys():
            run_JERCs_data.append((run, JERC))
        elif samples[run] in xsec_pb.keys():
            run_JERCs_MC.append((run, JERC))

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
        get_copy_for_merge(merged_run, JERC)
        set_lumi_for_merge(merged_run, JERC)

        make_pileup_data_merged(merged_run,JERC)

        Produce_Combination_File_and_plots(merged_run,JERC)

os.system('rm -f '+Step3_dir+'/tmp-process_logs.log')
