import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

#process.load("Configuration/StandardSequences/GeometryDB_cff")
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1)) # all events
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000)) # 5000 events
###################################### Run on AOD instead of MiniAOD? ########
runOnAOD = False #Federico
###################################### Run on RECO instead of MiniAOD? ########
runOnRECO = False
if runOnRECO: runOnAOD=True


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#    'root://xrootd.unl.edu//store/data/Run2015A/Commissioning/AOD/PromptReco-v1/000/246/865/00000/E69D98C7-150B-E511-9118-02163E014206.root' #AOD
      #'root://xrootd.unl.edu//store/data/Run2015A/Commissioning/RECO/PromptReco-v1/000/246/865/00000/3296E66B-160B-E511-AD26-02163E013395.root' #RECO
#    'file:Express_AOD_1.root' # AOD
   'file:../tuples/Data/SinglePhoton_Run2015D_file3.root'  #Express_miniAOD.root' # miniAOD
#      THISINPUTFILE
      )
    )

#process.options = cms.untracked.PSet(
#   allowUnscheduled = cms.untracked.bool(True)
#)
## Production Info                                                             
#process.configurationMetadata = cms.untracked.PSet(    
#      annotation = cms.untracked.string('Express_miniAOD nevts:-1'),    
#      name = cms.untracked.string('Applications'),                 
#      version = cms.untracked.string('$Revision: 1.19 $')      
#      )      
#
#############
# Output definition
#### save miniAOD ####
if runOnAOD:
  process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
       compressionAlgorithm = cms.untracked.string('LZMA'),
       compressionLevel = cms.untracked.int32(4),
       dataset = cms.untracked.PSet(
       dataTier = cms.untracked.string(''),
       filterName = cms.untracked.string('')
        ),
      dropMetaData = cms.untracked.string('ALL'),
      eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
      fastCloning = cms.untracked.bool(False),
      fileName = cms.untracked.string('Express_miniAOD_Giulia.root'),
      outputCommands = process.MINIAODEventContent.outputCommands,
      overrideInputFileSplitLevels = cms.untracked.bool(True)
      )

### RUN MINIAOD SEQUENCE
if runOnAOD:
  from FWCore.ParameterSet.Utilities import convertToUnscheduled
  process=convertToUnscheduled(process)
  process.load('Configuration.StandardSequences.PAT_cff')
  from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
  process = miniAOD_customizeAllData(process)

if runOnRECO:
  ### RUN PFCLUSTERJETS
  process.load("RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff")
  process.load("RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi")
  process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
  process.pfClusterRefsForJetsHCAL = cms.EDProducer("PFClusterRefCandidateProducer",
      src          = cms.InputTag('particleFlowClusterHCAL'),
      particleType = cms.string('pi+')
      )
  process.pfClusterRefsForJetsECAL = cms.EDProducer("PFClusterRefCandidateProducer",
      src          = cms.InputTag('particleFlowClusterECAL'),
      particleType = cms.string('pi+')
      )
  process.pfClusterRefsForJetsHF = cms.EDProducer("PFClusterRefCandidateProducer",
      src          = cms.InputTag('particleFlowClusterHF'),
      particleType = cms.string('pi+')
      )
  process.pfClusterRefsForJetsHO = cms.EDProducer("PFClusterRefCandidateProducer",
      src          = cms.InputTag('particleFlowClusterHO'),
      particleType = cms.string('pi+')
					  )
  process.pfClusterRefsForJets = cms.EDProducer("PFClusterRefCandidateMerger",
      src = cms.VInputTag("pfClusterRefsForJetsHCAL", "pfClusterRefsForJetsECAL", "pfClusterRefsForJetsHF", "pfClusterRefsForJetsHO")
      )
  process.load("RecoJets.JetProducers.ak4PFClusterJets_cfi")
  process.pfClusterRefsForJets_step = cms.Sequence(
      process.particleFlowRecHitECAL*
      process.particleFlowRecHitHBHE*
      process.particleFlowRecHitHF*
      process.particleFlowRecHitHO*
      process.particleFlowClusterECALUncorrected*
      process.particleFlowClusterECAL*
      process.particleFlowClusterHBHE*
      process.particleFlowClusterHCAL*
      process.particleFlowClusterHF*
      process.particleFlowClusterHO*
      process.pfClusterRefsForJetsHCAL*
      process.pfClusterRefsForJetsECAL*
      process.pfClusterRefsForJetsHF*
      process.pfClusterRefsForJetsHO*
      process.pfClusterRefsForJets*
      process.ak4PFClusterJets
      )
  ### RUN PFCALOJETS
  ############ need the following setup when running in <CMSSW_7_5_X:
  # git cms-addpkg RecoParticleFlow/PFProducer
  # git cherry-pick af5c1ba33e88b3be627c262eb93d678f9f70e729
  process.hltParticleFlowBlock = cms.EDProducer("PFBlockProducer",
      debug = cms.untracked.bool(False),
      verbose = cms.untracked.bool(False),
      elementImporters = cms.VPSet(
	cms.PSet(
	  source = cms.InputTag("particleFlowClusterECAL"),
	  #source = cms.InputTag("particleFlowClusterECALUncorrected"), #we use uncorrected
	  importerName = cms.string('GenericClusterImporter')
	  ),
	cms.PSet(
	  source = cms.InputTag("particleFlowClusterHCAL"),
	  importerName = cms.string('GenericClusterImporter')
	  ),
	cms.PSet(
	  source = cms.InputTag("particleFlowClusterHO"),
	  importerName = cms.string('GenericClusterImporter')
	  ),
	cms.PSet(
	  source = cms.InputTag("particleFlowClusterHF"),
	  importerName = cms.string('GenericClusterImporter')
	  )
	),
      linkDefinitions = cms.VPSet(
	cms.PSet(
	  linkType = cms.string('ECAL:HCAL'),
	  useKDTree = cms.bool(False),
	  #linkerName = cms.string('ECALAndHCALLinker')
	  linkerName = cms.string('ECALAndHCALCaloJetLinker') #new ECal and HCal Linker for PFCaloJets
	  ),
	cms.PSet(
	  linkType = cms.string('HCAL:HO'),
	  useKDTree = cms.bool(False),
	  linkerName = cms.string('HCALAndHOLinker')
	  ),
	cms.PSet(
	  linkType = cms.string('HFEM:HFHAD'),
	  useKDTree = cms.bool(False),
	  linkerName = cms.string('HFEMAndHFHADLinker')
	  ),
	cms.PSet(
	  linkType = cms.string('ECAL:ECAL'),
	  useKDTree = cms.bool(False),
	  linkerName = cms.string('ECALAndECALLinker')
	  )
	)
      )
  from RecoParticleFlow.PFProducer.particleFlow_cfi import particleFlowTmp
  process.hltParticleFlow = particleFlowTmp.clone(
      GedPhotonValueMap = cms.InputTag(""),
      useEGammaFilters = cms.bool(False),
      useEGammaElectrons = cms.bool(False), 
      useEGammaSupercluster = cms.bool(False),
      rejectTracks_Step45 = cms.bool(False),
      usePFNuclearInteractions = cms.bool(False),  
      blocks = cms.InputTag("hltParticleFlowBlock"), 
      egammaElectrons = cms.InputTag(""),
      useVerticesForNeutral = cms.bool(False),
      PFEGammaCandidates = cms.InputTag(""),
      useProtectionsForJetMET = cms.bool(False),
      usePFConversions = cms.bool(False),
      rejectTracks_Bad = cms.bool(False),
      muons = cms.InputTag(""),
      postMuonCleaning = cms.bool(False),
      usePFSCEleCalib = cms.bool(False)
      )
  from RecoJets.JetProducers.PFJetParameters_cfi import *
  process.PFCaloJetParameters = PFJetParameters.clone(
      src = cms.InputTag('hltParticleFlow')
      )
  from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
  process.ak4PFCaloJets = cms.EDProducer(
      "FastjetJetProducer",
      process.PFCaloJetParameters,
      AnomalousCellParameters,
      jetAlgorithm = cms.string("AntiKt"),
      rParam       = cms.double(0.4)
      )
  process.pfClusterRefsForJets_step += process.hltParticleFlowBlock
  process.pfClusterRefsForJets_step += process.hltParticleFlow
  process.pfClusterRefsForJets_step += process.ak4PFCaloJets
  #######
  process.MINIAODoutput.outputCommands.append("keep *_ak4PFCaloJets_*_*")
  process.MINIAODoutput.outputCommands.append("keep *_ak4PFClusterJets_*_*")


# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, THISGLOBALTAG, '')
process.GlobalTag.globaltag = cms.string("74X_dataRun2_Prompt_v2")
#process.GlobalTag.globaltag = cms.string("GR_P_V56")
#process.GlobalTag = GlobalTag(process.GlobalTag, "GR_P_V56", '')

#from FWCore.ParameterSet.VarParsing import VarParsing
#options = VarParsing()
#
#options.register ('datasetName',
#    '',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "The dataset currently processed. A folder named 'datasetName' must exists")
#
#options.register ('globalTag',
#    '',
#    VarParsing.multiplicity.singleton,
#    VarParsing.varType.string,
#    "The globaltag to use")
#
#options.parseArguments()
#if len(options.globalTag) == 0:
#  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")
#
#process.GlobalTag.globaltag = cms.string("PHYS14_25_V2::All")
#process.GlobalTag.globaltag = cms.string("%s::All" % options.globalTag) ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)

###############

#fullPath = os.path.join(os.getcwd(), options.datasetName)

##these statements have to be here

if runOnAOD:
#  process.endjob_step = cms.EndPath(process.endOfProcess)
  process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)



##################
process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi") # federico -> add for the photon ID/isolation

# Do some CHS stuff
process.ak4PFchsL1Fastjet  = process.ak4PFL1Fastjet.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL2Relative = process.ak4PFL2Relative.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL3Absolute = process.ak4PFL3Absolute.clone(algorithm = 'AK4PFchs')
process.ak4PFchsResidual   = process.ak4PFResidual.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL1FastL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak4PFchsL1Fastjet', 'ak4PFchsL2Relative','ak4PFchsL3Absolute')
    )
process.ak4PFchsL1FastL2L3Residual = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak4PFchsL1Fastjet', 'ak4PFchsL2Relative','ak4PFchsL3Absolute','ak4PFchsResidual')
    )


##-------------------- User analyzer  --------------------------------

if runOnAOD:
  calo_collection='ak4CaloJets'
else:
  calo_collection=''

if runOnRECO:
  cluster_collection='ak4PFClusterJets'
  pfcalo_collection='ak4PFCaloJets'
else:
  cluster_collection=''
  pfcalo_collection=''

##### User analyzer ######
process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("slimmedPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    #json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    #csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
     json = cms.string( "file:lumiSummary.json"),
     csv = cms.string( "file:lumibyls.csv"),
   
     filterData = cms.untracked.bool(False),

    runOnNonCHS   = cms.untracked.bool(False),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK4    = cms.untracked.bool(True),
    runOnPFAK8    = cms.untracked.bool(False),

    runOnCaloAK4  = cms.untracked.bool(False),
    runOnCaloAK8  = cms.untracked.bool(False),

    runOnPFClusterAK4  = cms.untracked.bool(False),

    # federico -> ValueMap names from the producer upstream
    #    full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"), # from rel73 ok in photon class
     phoChargedIsolation           = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
     phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
     phoPhotonIsolation             = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
     prescales = cms.InputTag("patTrigger"),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),
#    correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3"), #federico
    #if you want to correct data also with residual
    #correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3Residual"),

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(True)


    )



########## Path ##########
process.p = cms.Path()
if runOnRECO:
     process.p += process.pfClusterRefsForJets_step
process.p +=process.photonIDValueMapProducer  # federico -> add process for isolation
process.p += process.gammaJet
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.TFileService = cms.Service("TFileService",
#     fileName = cms.string(THISROOTFILE) # run on crab
     fileName = cms.string("output_singleFile_Data.root") # run in local
     )

#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
