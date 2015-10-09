from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'Output_Run2//New_GJet_Pt_15_6000_v2_20150727_0223/workdir'
config.General.requestName = 'GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8__RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v3__MINIAODSIM'
config.section_('JobType')
config.JobType.psetName = 'Output_Run2//New_GJet_Pt_15_6000_v2_20150727_0223/cfg/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8_cmssw.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8__RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v3__MINIAODSIM.root']
config.section_('Data')
config.Data.inputDataset = '/GJet_Pt-15To6000_TuneCUETP8M1-Flat_13TeV_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v3/MINIAODSIM'
config.Data.unitsPerJob = 10
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fpreiato/JetMET2015/MC/Spring15/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
