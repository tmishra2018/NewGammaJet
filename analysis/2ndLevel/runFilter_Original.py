import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration/StandardSequences/GeometryDB_cff")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
    #'/store/user/sbrochet/Photon/Photon_Run2012A-22Jan2013_24Apr13-v1/37e3bf2409397e623ffd52beab84a202/patTuple_PF2PAT_92_1_v1Q.root' 
    #      '/store/user/sbrochet/SinglePhotonParked/SinglePhoton_Run2012D-22Jan2013_16May13-v1/d8690a0448f852b4d70656ff31f27990/patTuple_PF2PAT_981_4_yaZ.root'
    'file:/cmshome/fpreiato/GammaJet/CMSSW_7_3_2/test/Commissioning2015_Data_1.root'
    )
                            )

from FWCore.ParameterSet.VarParsing import VarParsing


options = VarParsing()

options.register ('datasetName',
    '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "The dataset currently processed. A folder named 'datasetName' must exists")

options.register ('globalTag',
                  '',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

process.GlobalTag.globaltag = cms.string("%s::All" % options.globalTag) ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#process.GlobalTag.globaltag = cms.string("PHYS14_25_V2::All" % options.globalTag) ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)

fullPath = os.path.join(os.getcwd(), options.datasetName)

process.gammaJet = cms.EDFilter('GammaJetFilter',
                                isMC = cms.untracked.bool(False),
                                #    photons = cms.untracked.InputTag("selectedPatPhotons"),
                                photons = cms.untracked.InputTag("slimmedPhotons"),
                                firstJetPtCut = cms.untracked.bool(False),
                                
                                json = cms.string(os.path.join(fullPath, "lumiSummary.json")),
                                csv = cms.string(os.path.join(fullPath, "lumibyls.csv")),
                                filterData = cms.untracked.bool(True),
                                
                                runOnNonCHS   = cms.untracked.bool(False),
                                runOnCHS      = cms.untracked.bool(True),
                                
                                runOnPFAK4    = cms.untracked.bool(True),
                                runOnPFAK8    = cms.untracked.bool(True),
                                
                                runOnCaloAK4  = cms.untracked.bool(False),
                                runOnCaloAK8  = cms.untracked.bool(False),
                                
                                # JEC
                                doJetCorrection = cms.untracked.bool(False),
                                correctJecFromRaw = cms.untracked.bool(True),
                                #correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3"),
                                correctorLabel = cms.untracked.string("ak4PFchsL1FastL2L3Residual"),
                                
                                # MET
                                redoTypeIMETCorrection = cms.untracked.bool(False)
                                )

process.p = cms.Path(process.gammaJet)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output_data.root")
                                   )

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
