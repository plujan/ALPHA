import FWCore.ParameterSet.Config as cms

process = cms.Process('GenParticleListDrawer')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
#process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:/uscms/home/tucker/nobackup/store/mc/Summer11/ZprimeSSMToMuMu_M-2250_TuneZ2_7TeV-pythia6/AODSIM/PU_S4_START42_V11-v1/0000/7032435D-1A93-E011-94D5-0017A4770008.root'))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('/store/user/lbenato/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_PRIVATE-MC/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/160515_095125/0000/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_1.root') )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.ParticleListDrawer = cms.EDAnalyzer('ParticleListDrawer',
                                            maxEventsToPrint = cms.untracked.int32(100),
                                            src = cms.InputTag('prunedGenParticles'),
                                            printOnlyHardInteraction = cms.untracked.bool(False),
                                            useMessageLogger = cms.untracked.bool(False)
                                            )
process.p = cms.Path(process.ParticleListDrawer)
