from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'Output_Run2//SinglePhoton_GoldenJson_ResClosure_v4_20150814_1952/workdir'
config.General.requestName = 'SinglePhoton__Run2015B-PromptReco-v1__MINIAOD'
config.section_('JobType')
config.JobType.psetName = 'Output_Run2//SinglePhoton_GoldenJson_ResClosure_v4_20150814_1952/cfg/SinglePhoton_cmssw.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['SinglePhoton__Run2015B-PromptReco-v1__MINIAOD.root']
config.section_('Data')
config.Data.inputDataset = '/SinglePhoton/Run2015B-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
#config.Data.lumiMask = '/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_4_3/src/JetMETCorrections/GammaJetFilter/analysis/2ndLevel/jsonFile_245155.json'
config.Data.lumiMask = '/cmshome/fpreiato/GammaJet/CMSSW_7_4_5/src/JetMETCorrections/GammaJetFilter/json/Cert_246908-251883_13TeV_PromptReco_Collisions15_JSON_v2.json'
#config.Data.lumiMask = '/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_4_3/src/JetMETCorrections/GammaJetFilter/json/jsonFile_246908.json'
#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#config.Data.runRange = '251244-251252' # '193093-194075'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/fpreiato/JetMET2015/data/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
