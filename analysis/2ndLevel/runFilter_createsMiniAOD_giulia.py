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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
      #'/store/user/sbrochet/Photon/Photon_Run2012A-22Jan2013_24Apr13-v1/37e3bf2409397e623ffd52beab84a202/patTuple_PF2PAT_92_1_v1Q.root' 
      #'/store/user/sbrochet/SinglePhotonParked/SinglePhoton_Run2012D-22Jan2013_16May13-v1/d8690a0448f852b4d70656ff31f27990/patTuple_PF2PAT_981_4_yaZ.root'
      #"file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_4_3/src/JetMETCorrections/GammaJetFilter/analysis/2ndLevel/miniAOD-data_test.root"
      #'file:/cmshome/fpreiato/CMSSW_7_4_3/test/341F1ECF-D700-E511-81F9-02163E014295.root'
      THISINPUTFILE
      )
    )
#process.options = cms.untracked.PSet(
#    allowUnscheduled = cms.untracked.bool(True)
#)
## Production Info                                                             
#process.configurationMetadata = cms.untracked.PSet(    
#      annotation = cms.untracked.string('miniAOD-prod nevts:-1'),    
#      name = cms.untracked.string('Applications'),                 
#      version = cms.untracked.string('$Revision: 1.19 $')      
#      )      
#
#############
# Output definition

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
    fileName = cms.untracked.string('miniAOD-prod_PAT.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True)
    )

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, THISGLOBALTAG, '')
#process.GlobalTag.globaltag = cms.string("MCRUN2_74_V8")

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
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)


############
##   MINIAOD
############
#
##do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)
#process.load('Configuration.StandardSequences.PAT_cff')
#
## customisation of the process.
#
## Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 
#
##call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#process = miniAOD_customizeAllData(process)
#
## End of customisation functions
#
##################
process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")


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


process.gammaJet = cms.EDFilter('GammaJetFilter',
    isMC = cms.untracked.bool(False),
    photons = cms.untracked.InputTag("slimmedPhotons"),
    firstJetPtCut = cms.untracked.bool(False),

    #json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
    #csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
    json = cms.string( "file:lumiSummary.json"),
    csv = cms.string( "file:lumibyls.csv"),
    #filterData = cms.untracked.bool(True),
    filterData = cms.untracked.bool(False),

    runOnNonCHS   = cms.untracked.bool(False),
    runOnCHS      = cms.untracked.bool(True),

    runOnPFAK4    = cms.untracked.bool(True),
    runOnPFAK8    = cms.untracked.bool(False),

    runOnCaloAK4  = cms.untracked.bool(False),
    runOnCaloAK8  = cms.untracked.bool(False),

    # JEC
    doJetCorrection = cms.untracked.bool(True),
    correctJecFromRaw = cms.untracked.bool(True),
    correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3"),
    #if you want to correct data also with residual
    #correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3Residual"),

    # MET
    redoTypeIMETCorrection = cms.untracked.bool(True)
    )

process.p = cms.Path(process.gammaJet)
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.TFileService = cms.Service("TFileService",
     #fileName = cms.string("output.root")
     fileName = cms.string(THISROOTFILE)
     )

#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
