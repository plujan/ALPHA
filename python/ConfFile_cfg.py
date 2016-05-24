import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("ALPHA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/lustre/cmswork/pazzini/CMSSW_8_0_4/src/2455D4FC-A5F0-E511-88AC-0025905A48C0.root'
        '/store/user/lbenato/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_PRIVATE-MC/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/160515_095125/0000/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_1.root'
    )
)

#output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("deleteme.root"),
    closeFileFast = cms.untracked.bool(True)
)

#electrons upstream modules
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']
for ele_idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,ele_idmod,setupVIDElectronSelection)

#photons upstream modules
switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)
ph_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']
for ph_idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,ph_idmod,setupVIDPhotonSelection)

#muons upstream modules
process.cleanedMuons = cms.EDProducer("PATMuonCleanerBySegments",
    src = cms.InputTag("slimmedMuons"),#("calibratedMuons"),
    preselection = cms.string("track.isNonnull"),
    passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
    fractionOfSharedSegments = cms.double(0.499)
)

process.ntuple = cms.EDAnalyzer('Ntuple',
    genSet = cms.PSet(
        genparticles = cms.InputTag("packedGenParticles"),
    ),
    pileupSet = cms.PSet(
        pileup = cms.InputTag("slimmedAddPileupInfo"),
        dataFileName = cms.string('%s/src/Analysis/ALPHA/data/Prod6.root' % os.environ['CMSSW_BASE']),
        mcFileName = cms.string('%s/src/Analysis/ALPHA/data/MC_True.root' % os.environ['CMSSW_BASE']),
        dataName = cms.string('pileup'),
        mcName = cms.string('S10'),
    ),
    triggerSet = cms.PSet(
        trigger = cms.InputTag("TriggerResults"),
        paths = cms.vstring('HLT_Mu45_eta2p1_v', 'HLT_Mu50_v2', 'HLT_IsoMu20_v', 'HLT_Mu27_TkMu8_v', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v', 'HLT_Ele105_CaloIdVT_GsfTrkIdT_v', 'HLT_Ele23_WPLoose_Gsf_v', 'HLT_Ele27_WPLoose_Gsf_v', 'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 'HLT_DoubleEle33_CaloIdL_v'),
    ),
    electronSet = cms.PSet(
        electrons = cms.InputTag("slimmedElectrons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
        eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
        eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
        eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
        eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
        eleMVANonTrigMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
        eleMVANonTrigTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
        eleMVATrigMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
        eleMVATrigTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
        eleVetoIdFileName = cms.string('%s/src/Analysis/ALPHA/data/CutBasedID_VetoWP_76X_18Feb.txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleLooseIdFileName = cms.string('%s/src/Analysis/ALPHA/data/CutBasedID_LooseWP_76X_18Feb.txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/CutBasedID_MediumWP_76X_18Feb.txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleTightIdFileName = cms.string('%s/src/Analysis/ALPHA/data/CutBasedID_TightWP_76X_18Feb.txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleMVATrigMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/ScaleFactor_GsfElectronToRECO_passingTrigWP90.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        eleMVATrigTightIdFileName = cms.string('%s/src/Analysis/ALPHA/data/ScaleFactor_GsfElectronToRECO_passingTrigWP80.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        eleRecoEffFileName = cms.string('%s/src/Analysis/ALPHA/data/eleRECO.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        electron1id = cms.int32(1), # 0: veto, 1: loose, 2: medium, 3: tight, 4: HEEP, 5: MVA medium nonTrig, 6: MVA tight nonTrig, 7: MVA medium Trig, 8: MVA tight Trig
        electron2id = cms.int32(1),
        #electron1iso = cms.int32(1), # 0: veto, 1: standard
        #electron2iso = cms.int32(1),
        electron1pt = cms.double(20.),
        electron2pt = cms.double(10.),
    ),
    muonSet = cms.PSet(
        muons = cms.InputTag("cleanedMuons"),#("slimmedMuons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        muonIdFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonID_Z_RunCD_Reco76X_Feb15.root' % os.environ['CMSSW_BASE']),
        muonIsoFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonIso_Z_RunCD_Reco76X_Feb15.root' % os.environ['CMSSW_BASE']),
        muonHighptFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonHighPt_Z_RunCD_Reco74X_Dec17.root' % os.environ['CMSSW_BASE']),
        muonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/SingleMuonTrigger_Z_RunCD_Reco76X_Feb15.root' % os.environ['CMSSW_BASE']),
        doubleMuonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root' % os.environ['CMSSW_BASE']),#obsolete
        muon1id = cms.int32(3), # 0: tracker id, 1: loose, 2: medium, 3: tight, 4: high pt
        muon2id = cms.int32(3),
        muon1iso = cms.int32(2), # 0: no selections, 1: loose (0.25), 2: tight (0.15) (pfIso in cone 0.4)
        muon2iso = cms.int32(2),
        muon1pt = cms.double(20.),
        muon2pt = cms.double(10.),
    ),
    photonSet = cms.PSet(
        photons = cms.InputTag("slimmedPhotons"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        phoLooseIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
        phoMediumIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
        phoTightIdMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
        phoMVANonTrigMediumIdMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
        phoLooseIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Loosenumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Mediumnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoTightIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Tightnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoMVANonTrigMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/MVAnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        photonid = cms.int32(1), # 1: loose, 2: medium, 3: tight, 4:MVA NonTrig medium
        photonpt = cms.double(20.),
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
    writeNPhotons = cms.int32(1),
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

process.seq = cms.Sequence(process.egmGsfElectronIDSequence * process.egmPhotonIDSequence * process.cleanedMuons * process.ntuple)
process.p = cms.Path(process.seq)
