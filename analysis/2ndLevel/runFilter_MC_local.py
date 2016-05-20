import FWCore.ParameterSet.Config as cms
import os
from CondCore.CondDB.CondDB_cfi  import *

process = cms.Process("GAMMAJET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000 #1000

#--- import of standard configurations
#process.load("Configuration/StandardSequences/GeometryDB_cff")
#process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag('slimmedJets')    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel       = cms.string('QGL_AK4PFchs')        # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

process.GlobalTag.globaltag = cms.string("80X_mcRun2_asymptotic_2016_v3") # run in local

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

process.chs = cms.EDFilter('CandPtrSelector', src = cms.InputTag('packedPFCandidates'), cut = cms.string('fromPV'))

# Do some CHS stuff
process.ak4PFchsL1Fastjet  = process.ak4PFL1Fastjet.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL2Relative = process.ak4PFL2Relative.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL3Absolute = process.ak4PFL3Absolute.clone(algorithm = 'AK4PFchs')
process.ak4PFchsResidual   = process.ak4PFResidual.clone(algorithm = 'AK4PFchs')
process.ak4PFchsL1FastL2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak4PFchsL1Fastjet', 'ak4PFchsL2Relative','ak4PFchsL3Absolute')
    )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) ) #run over all events
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) ) # run only on # events

from FWCore.ParameterSet.VarParsing import VarParsing
#readFiles = cms.untracked.vstring(
#    )

#readFiles.extend( [
#  'file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_3_2/test/test_file_MINIAOD_for_JEC2015.root',
#  ])

process.source = cms.Source (
    "PoolSource", 
    fileNames = cms.untracked.vstring(
        #'file:/cmshome/gdimperi/GammaJet/JetCorrections/CMSSW_7_3_2/test/test_file_MINIAOD_for_JEC2015.root'
        'file:../tuples/GJET_Pythia/GJet_Pythia_80X_file1.root'
        #'file:../tuples/QCD_Pythia/QCD_Pythia_80X_file1.root'
      )
    )

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('crossSection',
				  '',
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.float,
				  "Dataset cross section")

options.parseArguments()

crossSection = float(options.crossSection) if isinstance(options.crossSection, float) and float(options.crossSection) != 0 else 1000

print("Running on sample with:")
print("\tCross-section: %f" % crossSection)

process.gammaJet = cms.EDFilter('GammaJetFilter',
                                isMC = cms.untracked.bool(True),
                                firstJetPtCut = cms.untracked.bool(False),
                                
                                crossSection = cms.double(crossSection),
                                                               
                                full5x5SigmaIEtaIEtaMap   = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
                                phoChargedIsolation           = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
                                phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
                                phoPhotonIsolation             = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"), 
                                prescalesTag = cms.InputTag("patTrigger"),  
                                triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),  
                                generatorTag = cms.InputTag("generator"),  
                                vertexTag = cms.InputTag("offlineSlimmedPrimaryVertices"),  
                                photonsTag = cms.InputTag("slimmedPhotons"),
                                jetsTag = cms.InputTag("slimmedJets"),
                                jetsAK8Tag = cms.InputTag("slimmedJetsAK8"),
                                metTag = cms.InputTag("slimmedMETs"),
                                electronsTag = cms.InputTag("slimmedElectrons"),
                                muonsTag = cms.InputTag("slimmedMuons"),
                                rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
                                PUInfoTag = cms.InputTag("slimmedAddPileupInfo"),
                                pfCands = cms.InputTag("packedPFCandidates"),                                                               

                                runOnNonCHS   = cms.untracked.bool(False),
                                runOnCHS      = cms.untracked.bool(True),
                                
                                runOnPFAK4    = cms.untracked.bool(True),
                                runOnPFAK8    = cms.untracked.bool(False),
                                
                                # MET
                                redoTypeIMETCorrection = cms.untracked.bool(True),
                                doFootprintMETCorrection = cms.untracked.bool(True),                             

                                # JEC
                                doJetCorrection = cms.untracked.bool(True),
                                correctJecFromRaw = cms.untracked.bool(True),
                                                                                               
                                L1corr_MC = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV1_MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt'),
                                L2corr_MC = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV1_MC/Spring16_25nsV1_MC_L2Relative_AK4PFchs.txt'),
                                L3corr_MC = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV1_MC/Spring16_25nsV1_MC_L3Absolute_AK4PFchs.txt'),
                                L1RCcorr_MC =cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV1_MC/Spring16_25nsV1_MC_L1FastJet_AK4PFchs.txt') #L1RC
                               
               )

process.p = cms.Path(
    process.chs *
    process.photonIDValueMapProducer *
    process.QGTagger *
    process.gammaJet)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("delete_me.root"),
    SelectEvents = cms.untracked.PSet(
      SelectEvents = cms.vstring('p')
      )
    )

process.out.outputCommands = cms.untracked.vstring('keep *',
#    'drop *_selectedPatJets*_*_*',
#    'drop *_selectedPatPhotons*_*_*',
#    'keep *_selectedPatJets*_genJets_*',
#    'keep *_selectedPatJets*_caloTowers_*',
#    # Drop CHS
#    'drop *_*chs*_*_*'
)

process.TFileService = cms.Service("TFileService",
# federico
#          fileName = cms.string("output_singleFile_QCD.root") # run in local
          fileName = cms.string("output_singleFile_GJet.root") # run in local
    )

#process.out.fileName = 'patTuple_cleaned.root'
# set True if you want long output
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

#process.outpath = cms.EndPath(process.out)
