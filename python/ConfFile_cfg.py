import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("ALPHA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/lustre/cmswork/pazzini/CMSSW_8_0_4/src/2455D4FC-A5F0-E511-88AC-0025905A48C0.root'
    )
)

#output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("deleteme.root"),
    closeFileFast = cms.untracked.bool(True)
)

#electrons first attempt
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.cleanedMuons = cms.EDProducer("PATMuonCleanerBySegments",
    src = cms.InputTag("slimmedMuons"),#("calibratedMuons"),
    preselection = cms.string("track.isNonnull"),
    passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
    fractionOfSharedSegments = cms.double(0.499)
)

process.ntuple = cms.EDAnalyzer('Ntuple',
    pileupSet = cms.PSet(
        pileup = cms.InputTag("slimmedAddPileupInfo"),
        dataFileName = cms.string('%s/src/Analysis/ALPHA/data/Prod6.root' % os.environ['CMSSW_BASE']),
        mcFileName = cms.string('%s/src/Analysis/ALPHA/data/MC_True.root' % os.environ['CMSSW_BASE']),
        dataName = cms.string('pileup'),
        mcName = cms.string('S10'),
    ),
    electronSet = cms.PSet(
        electrons = cms.InputTag("slimmedElectrons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
        electron1id = cms.int32(1), # 0: veto, 1: loose, 2: medium, 3: tight, 4: HEEP
        electron2id = cms.int32(1),
        electron1iso = cms.int32(1), # 0: veto, 1: standard
        electron2iso = cms.int32(1),
        electron1pt = cms.double(20.),
        electron2pt = cms.double(10.),
    ),
    muonSet = cms.PSet(
        muons = cms.InputTag("cleanedMuons"),#("slimmedMuons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        muonIdFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonID_Z_RunCD_Reco74X_Dec1.root' % os.environ['CMSSW_BASE']),
        muonIsoFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonIso_Z_RunCD_Reco74X_Dec1.root' % os.environ['CMSSW_BASE']),
        muonHighptFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonHighPt_Z_RunCD_Reco74X_Dec17.root' % os.environ['CMSSW_BASE']),
        muonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/SingleMuonTrigger_Z_RunCD_Reco74X_Dec1.root' % os.environ['CMSSW_BASE']),
        doubleMuonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root' % os.environ['CMSSW_BASE']),#obsolete
        muon1id = cms.int32(3), # 0: no selections, 1: loose, 2: medium, 3: tight, 4: soft, 5: high pt
        muon2id = cms.int32(3),
        muon1iso = cms.int32(2), # 0: no selections, 1: loose (0.25), 2: tight (0.15) (pfIso in cone 0.4)
        muon2iso = cms.int32(2),
        muon1pt = cms.double(20.),
        muon2pt = cms.double(10.),
    ),
    photonSet = cms.PSet(
        photons = cms.InputTag("slimmedPhotons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        photon1id = cms.int32(1), # 0: veto, 1: loose, 2: medium, 3: tight
        photon2id = cms.int32(1),
        photon1iso = cms.int32(0),
        photon2iso = cms.int32(0),
        photon1pt = cms.double(20.),
        photon2pt = cms.double(10.),
    ),
    jetSet = cms.PSet(
        jets = cms.InputTag("slimmedJetsAK8"), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(1), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet1pt = cms.double(20.),
        jet2pt = cms.double(20.),
        btag = cms.string("combinedSecondaryVertexBJetTags"),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag("slimmedMETs"),
        metRecoil = cms.bool(True),
        metRecoilMC = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsMC_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
        metRecoilData = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsData_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
    ),
    writeNElectrons = cms.int32(0),
    writeNMuons = cms.int32(0),
    writeNLeptons = cms.int32(2),
    writeNJets = cms.int32(2),
    verbose  = cms.bool(True),
)


#process.output = cms.OutputModule("PoolOutputModule",
#  splitLevel = cms.untracked.int32(0),
#  fileName = cms.untracked.string('deleteme.root'),
#  dataset = cms.untracked.PSet(
#    filterName = cms.untracked.string(''),
#    dataTier = cms.untracked.string('')
#  )
#)

process.seq = cms.Sequence(process.egmGsfElectronIDSequence * process.cleanedMuons * process.ntuple)
#process.seq = cms.Sequence(process.cleanedMuons * process.ntuple)
process.p = cms.Path(process.seq)
