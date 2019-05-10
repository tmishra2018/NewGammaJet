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

def get_most_recent(directory):
    return os.popen('ls -t {}'.format(directory)).read().split('\n')[0]

def find_files_from_step2(run, JERC):
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
        command += ' ; rm '+list_file
        for letter in run :
            # find_files_from_step2(letter, JERC)
            command += ' ; cat '+"Step2_files_{}_{}.txt".format(letter,JERC)
            command += " >> "+list_file
    os.system(command)
        
def merge_files_from_step2(run, JERC, cleaning=True):
    command = 'cd '+Step3_scripts_dir
    output_dir = "/eos/user/${USER:0:1}/$USER/JEC-task/Step3_outputs/2018/"+JERC+'/'
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
    print command
    os.system(command)

def make_pileup(run,JERC):
    command = 'cd '+Step3_PU_dir
    if run == 'MC':
        command += ' ; echo '+get_most_recent(Step3_outputs_base_dir+'/'+JERC+'/*_MC_*')+' > files_2018_MC_'+JERC+'.list'
        command += ' ; ./generate_mc_pileup.exe files_2018_MC_'+JERC
    else:
        command += ' ; python ComputePU_perHLT.py'
        command += ' --Cert_json '+Step1_outputs_crab_dir+'/crab_'+samples[run]+'/results/processedLumis.json'
        command += ' --whichrun '+run+'_'+JERC
    print command
    os.system(command)

def Produce_Combination_File_and_plots(run,JERC):
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

steps = [['MC'], ['A', 'B', 'C'], ['ABC', 'D'],['ABCD']]
for step in steps:
    run_JERCs = []
    for run in step:
        for JERC in JERCs:
            run_JERCs.append((run,JERC))
    multithreadmap(process, run_JERCs)
