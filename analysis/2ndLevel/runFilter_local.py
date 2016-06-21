import FWCore.ParameterSet.Config as cms
import os
from CondCore.CondDB.CondDB_cfi  import *

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
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets          = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel       = cms.string('QGL_AK4PFchs') #see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion

process.GlobalTag.globaltag = cms.string("80X_dataRun2_Prompt_v8")

process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000))

process.source = cms.Source(
  "PoolSource",
  fileNames = cms.untracked.vstring(
    'file:../tuples/SinglePhoton_file1_80X.root'
    )
  )

process.gammaJet = cms.EDFilter('GammaJetFilter',
                                isMC = cms.untracked.bool(False),
                                firstJetPtCut = cms.untracked.bool(False),

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
                                jetsPUPPITag = cms.InputTag("slimmedJetsPuppi"),
                                metTag = cms.InputTag("slimmedMETs"),
                                metPUPPITag = cms.InputTag("slimmedMETsPuppi"),
                                electronsTag = cms.InputTag("slimmedElectrons"),
                                muonsTag = cms.InputTag("slimmedMuons"),
                                rhoTag = cms.InputTag("fixedGridRhoFastjetAll"),
                                PUInfoTag = cms.InputTag("slimmedAddPileupInfo"),
                                pfCands = cms.InputTag("packedPFCandidates"),                                                               

                                runOnPFAK4    = cms.untracked.bool(True),
                                runOnPFAK8    = cms.untracked.bool(False),
                                runOnPUPPIAK4    = cms.untracked.bool(True),
                                
                                # MET
                                redoTypeIMETCorrection = cms.untracked.bool(True),
                                doFootprintMETCorrection = cms.untracked.bool(False),
                                
                                # JEC
                                doJetCorrection = cms.untracked.bool(True),
                                correctJecFromRaw = cms.untracked.bool(True),
                                applyL2Res = cms.untracked.bool(True),
                                applyL2L3Res = cms.untracked.bool(False),
                                
                                L1corr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L1FastJet_AK4PFchs.txt'),
                                L2corr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2Relative_AK4PFchs.txt'),
                                L3corr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L3Absolute_AK4PFchs.txt'),
                                L1RCcorr_DATA      = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L1RC_AK4PFchs.txt'),
                                L2Rescorr_DATA     = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2Residual_AK4PFchs.txt'),
                                L2L3Rescorr_DATA = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2L3Residual_AK4PFchs.txt'),

                                L1PUPPIcorr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L1FastJet_AK4PFPuppi.txt'),
                                L2PUPPIcorr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2Relative_AK4PFPuppi.txt'),
                                L3PUPPIcorr_DATA           = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L3Absolute_AK4PFPuppi.txt'),
                                L1RCPUPPIcorr_DATA      = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L1FastJet_AK4PFPuppi.txt'),
                                L2ResPUPPIcorr_DATA     = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2Residual_AK4PFPuppi.txt'),
                                L2L3ResPUPPIcorr_DATA = cms.FileInPath('JetMETCorrections/GammaJetFilter/data/Spring16_25nsV4_DATA/Spring16_25nsV4_DATA_L2L3Residual_AK4PFPuppi.txt')

                                )



########## Path ##########
process.p = cms.Path(
  process.photonIDValueMapProducer *
  process.QGTagger *
  process.gammaJet)

#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output_singleFile_Data.root")
                                   )

# set True if you want long output
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
