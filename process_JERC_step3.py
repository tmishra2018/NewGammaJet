import os
import subprocess
from functools import partial
import multiprocessing as mp

samples = {
    'A' : 'Run2018A-17Sep2018-v2',
    'B' : 'Run2018B-17Sep2018-v1',
    'C' : 'Run2018C-17Sep2018-v1',
    'D' : 'Run2018D-PromptReco-v2',
    'MC': 'GJet_Pt-15To6000_RunIIAutumn18MiniAOD-102X',
}

lumis_or_xsec_pb = {
    'A' : 13654.355526985,
    'B' : 7057.825158567,
    'C' : 6894.770971269,
    'D' : 31066.589629726,
    'MC': 283000.0,
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

def wait_for(runs, JERC, step, time='1m'):
    waiting = 0
    for run in runs:
        command = 'if cat '+Step3_dir+'/tmp-process_logs.log | grep -q '+"'run "+run+' '+JERC+" "+step+"'"
        command += " ; then echo '0' ; else echo '1' ; fi"
        waiting += int(os.popen(command).read()[:-1])
    if not waiting == 0:
        command = 'echo '+"'waiting "+time+" for runs "
        for run in runs:
            command += run+' '
        command += JERC+" "+step+"'"
        os.system(command)
        os.system('sleep '+time)
        wait_for(runs, JERC, step, time=time)

def get_most_recent(directory):
    return os.popen('ls -t {}'.format(directory)).read().split('\n')[0]

def find_files_from_step2(run, JERC):
    print_log_started_step(run, JERC, 'find files')
    list_file = "Step2_files_{}_{}.txt".format(run,JERC)
    command = 'cd '+Step3_scripts_dir
    if run in samples.keys():
        command +=' ; find '
        directory = Step2_outputs_base_dir+'/'
        directory += samples[run]+'/'
        directory += get_most_recent(directory)+'/'
        command += directory
        command += ' -type f '
        command += " -iname \*{}\*reduced_skim.root ".format(JERC)
        command += " > "+list_file
    if run in ['ABC', 'ABCD']:
        wait_for([letter for letter in run], JERC, 'merged OK',time='10m')
        command += ' ; rm -f '+list_file
        for letter in run :
            # find_files_from_step2(letter, JERC)
            command += ' ; cat '+"Step2_files_{}_{}.txt".format(letter,JERC)
            command += " >> "+list_file
    os.system(command)
        
def merge_files_from_step2(run, JERC, cleaning=True):
    print_log_started_step(run, JERC, 'merge files')
    command = 'cd '+Step3_scripts_dir
    output_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2018/"+JERC+'/'
    os.system("mkdir -p {}".format(output_dir))
    if run in lumis_or_xsec_pb.keys():
        lumi_or_xsec_pb = lumis_or_xsec_pb[run]
    elif run in ['ABC','ABCD']:
        lumi_or_xsec_pb = 0
        for letter in run:
            lumi_or_xsec_pb += lumis_or_xsec_pb[letter]
    if run == 'MC':
        command += ' ; python mergeAndAddWeights.py -c {} '.format(cleaning)
    else:
        command += ' ; python mergeData.py -c {} '.format(cleaning)
    command += " -i Step2_files_{}_{}.txt".format(run,JERC)
    command += " -o "+output_dir
    if run == 'MC':
        command += " --xsec "
    else:
        command += " --lumi_tot "
    command += str(lumi_or_xsec_pb)
    if not run == 'MC':
        command += " --run "+run
        command += " && echo 'run "+run+' '+JERC+" merged OK' >> "+Step3_dir+'/tmp-process_logs.log'
    print command
    os.system(command)

def make_pileup(run,JERC):
    print_log_started_step(run, JERC, 'make pileup')
    command = 'cd '+Step3_PU_dir
    if run == 'MC':
        command += ' ; echo '+get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/*_MC_*')+' > files_2018_MC_'+JERC+'.list'
        command += ' ; ./generate_mc_pileup.exe files_2018_MC_'+JERC
    if not run in samples.keys():
        processedLumis = ''
        for letter in run :
            processedLumis += os.popen('cat '+Step1_outputs_crab_dir+'/crab_'+samples[letter]+'/results/processedLumis.json').read()[:-1]
        processedLumis = processedLumis.replace('}{',', ')
        command += " ; echo '"+processedLumis+"' > ./tmp_processedLumis_"+run+'_'+JERC+'.json'
    if not run == 'MC':
        command += ' ; python ComputePU_perHLT.py'
        if run in samples.keys():
            command += ' --Cert_json '+Step1_outputs_crab_dir+'/crab_'+samples[run]+'/results/processedLumis.json'
        else:
            command += ' --Cert_json ./tmp_processedLumis_'+run+'_'+JERC+'.json'
        command += ' --whichrun '+run+'_'+JERC
    if run == 'MC':
        command += " && echo 'run MC "+JERC+" PU OK' >> "+Step3_dir+'/tmp-process_logs.log'
    print command
    os.system(command)

def Produce_Combination_File_and_plots(run,JERC):
    print_log_started_step(run, JERC, 'make plots')
    if not run == 'MC':
        wait_for(['MC'], JERC, 'PU OK')
        os.system("cd {} ; ln -sf computed_mc_files_2018_MC_{}_pu_truth_100bins.root computed_mc_files_2018_MC_for_{}_{}_pu_truth_100bins.root".format(
            Step3_PU_dir,
            JERC,
            run,
            JERC
            )
        )
    output = JERC
    IOV = run+'_'+JERC
    Pu_profile = run+'_'+JERC
    Input_data = get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_DATA_RUN_'+run+'_*')
    Input_mc = get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/'+'PhotonJet_2ndLevel_MC_*')
    command = 'cd '+Step3_dir
    command += ' ; python Produce_Combination_File_and_plots.py'
    command += ' --output='+output
    command += ' --IOV='+IOV
    command += ' --Pu_profile='+Pu_profile
    command += ' --Input_data='+Input_data
    command += ' --Input_mc='+Input_mc
    os.system(command)

def multithreadmap(process,X,ncores=8, **kwargs):
    """
    multithreading map of a function, default on 8 cpu cores.
    """
    func = partial(process, **kwargs)
    p=mp.Pool(ncores)
    Xout = p.map(func,X)
    p.terminate()
    return(Xout)

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


def process(run_JERC):
    run, JERC = run_JERC
    find_files_from_step2(run, JERC)
    merge_files_from_step2(run, JERC)
    
    make_pileup(run,JERC)
    
    if not run == 'MC' :
        Produce_Combination_File_and_plots(run,JERC)

os.system('rm -f '+Step3_dir+'/tmp-process_logs.log')
os.system('touch '+Step3_dir+'/tmp-process_logs.log')


done_run_JERCs = []
# for run in ['A','B', 'MC']:
#     for JERC in JERCs:
#         done_run_JERCs.append((run,JERC))
# for run in ['C']:
#     for JERC in ['wo_L2Res', 'L2L3Res', 'JER']:
#         done_run_JERCs.append((run,JERC))
# for run in ['D']:
#     for JERC in ['wo_L2Res', 'L2L3Res', 'only_L2Res']:
#         done_run_JERCs.append((run,JERC))

run_JERCs = []
sorted_runs_for_multithreadmap = ['MC', 'C', 'B']
sorted_runs_for_multithreadmap = [run for run in sorted_runs_for_multithreadmap if run in samples.keys()]
sorted_runs_for_multithreadmap += [key for key in samples.keys() if not key in sorted_runs_for_multithreadmap]
sorted_runs_for_multithreadmap += ['ABC', 'ABCD']
for run in sorted_runs_for_multithreadmap:
    for JERC in JERCs:
        if not (run,JERC) in done_run_JERCs:
                run_JERCs.append((run,JERC))
print run_JERCs
os.system("sleep 10s")
multithreadmap(process, run_JERCs)
os.system('rm -f '+Step3_dir+'/tmp-process_logs.log')
