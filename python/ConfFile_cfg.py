import FWCore.ParameterSet.Config as cms
import os

isMC = True

process = cms.Process("ALPHA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

# input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/lustre/cmswork/zucchett/CMSSW_8_0_5/src/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
        '/store/user/lbenato/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_PRIVATE-MC/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/160515_095125/0000/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_1.root'
    )
)

#output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("deleteme.root"),
    closeFileFast = cms.untracked.bool(True)
)

#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
import FWCore.PythonUtilities.LumiList as LumiList
if not isMC:
    process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()

# Trigger filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.HLTFilter = cms.EDFilter("HLTHighLevel",
    TriggerResultsTag = cms.InputTag("TriggerResults", "", "HLT"),
    HLTPaths = cms.vstring(
        'HLT_Mu45_eta2p1_v*',
        'HLT_Mu50_v*',
        'HLT_IsoMu20_v*',
        'HLT_Mu27_TkMu8_v*',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*',
        'HLT_Ele105_CaloIdVT_GsfTrkIdT_v*',
        'HLT_Ele23_WPLoose_Gsf_v*',
        'HLT_Ele27_WPLoose_Gsf_v*',
        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*',
        'HLT_DoubleEle33_CaloIdL_v*'
    ),
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),    # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False)    # throw exception on unknown path names
)

# Primary vertex
import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24), 
    maxd0 = cms.double(2) 
)


#-----------------------#
#        OBJECTS        #
#-----------------------#

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


#-----------------------#
#        NTUPLE         #
#-----------------------#

process.ntuple = cms.EDAnalyzer('Ntuple',
    genSet = cms.PSet(
        genProduct = cms.InputTag("generator"),
        lheProduct = cms.InputTag("externalLHEProducer"),
        genParticles = cms.InputTag("prunedGenParticles"),
        pdgId = cms.vint32(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 1000022, 9100000, 9000001, 9000002, 9100012, 9100022, 9900032, 1023),
    ),
    pileupSet = cms.PSet(
        pileup = cms.InputTag("slimmedAddPileupInfo"),
        dataFileName = cms.string('%s/src/Analysis/ALPHA/data/Prod6.root' % os.environ['CMSSW_BASE']),
        mcFileName = cms.string('%s/src/Analysis/ALPHA/data/PU_MC.root' % os.environ['CMSSW_BASE']),
        dataName = cms.string('pileup'),
        mcName = cms.string('2016_25ns_SpringMC_PUScenarioV1'),
    ),
    triggerSet = cms.PSet(
        trigger = cms.InputTag("TriggerResults", "", "HLT"),
        paths = cms.vstring('HLT_Mu45_eta2p1_v', 'HLT_Mu50_v', 'HLT_IsoMu20_v', 'HLT_Mu27_TkMu8_v', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v', 'HLT_Ele105_CaloIdVT_GsfTrkIdT_v', 'HLT_Ele23_WPLoose_Gsf_v', 'HLT_Ele27_WPLoose_Gsf_v', 'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 'HLT_DoubleEle33_CaloIdL_v'),
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
    tauSet = cms.PSet(
        taus = cms.InputTag("slimmedTaus"),
        vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
        taupt = cms.double(20.),
        taueta = cms.double(9999.),
        tauIdByDecayMode = cms.int32(1),# 0: not set, 1: old, 2: new
        tauIdByDeltaBetaIso = cms.int32(0),# 0: not set, 1: loose, 2: medium, 3: tight
        tauIdByMVAIso = cms.int32(0),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight, 5: V tight
        tauIdByMuonRejection = cms.int32(0),# 0: not set, 1: loose, 2: tight
        tauIdByElectronRejection = cms.int32(0),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight
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
        jets = cms.InputTag("slimmedJets"), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(0), # 0: no selection, 1: loose, 2: tight, 3: tight lepton veto
        jet1pt = cms.double(20.),
        jet2pt = cms.double(20.),
        jeteta = cms.double(2.4),
        btag = cms.string("combinedSecondaryVertexBJetTags"),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag("slimmedMETs"),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsMC_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
        metRecoilData = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsData_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
    ),
    writeNElectrons = cms.int32(0),
    writeNMuons = cms.int32(0),
    writeNLeptons = cms.int32(2),
    writeNTaus = cms.int32(0),
    writeNPhotons = cms.int32(1),
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

process.seq = cms.Sequence(
    process.HLTFilter *
    process.primaryVertexFilter *
    process.egmGsfElectronIDSequence *
    process.egmPhotonIDSequence *
    process.cleanedMuons *
    process.ntuple
)

process.p = cms.Path(process.seq)
