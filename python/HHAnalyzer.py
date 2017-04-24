#FIXME JSON file here
# New electron SF
# Check muon SF
# TO EXECUTE: 
# cmsRun python/HHAnalyzer.py maxEvents=10

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os
import sys

options = VarParsing ('analysis')
options.register ('tCut', 0, options.multiplicity.singleton, options.varType.int,
                  "Trigger cut on/off")
options.register ('Data', False, options.multiplicity.singleton, options.varType.bool,
                  "Data on/off")
options.parseArguments()
#Un-comment to change output name - default is 'output' (useful for debug)
#options.outputFile = 'test'

# Determine sample name for MC stitching
sample = (options.inputFiles[0]).split('/')[-1].replace('.txt', '') if len(options.inputFiles) > 0 else ''
if sample=='list': sample = (options.inputFiles[0]).split('/')[-3]

process = cms.Process('ALPHA')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# input
# default: if no filelist from command line, run on specified samples
if len(options.inputFiles) == 0:
    #SIGNAL - full list
    filelist = ['file:miniAOD-prod_PAT.root'] 
# production: read externally provided filelist
else:
    filelist = open(options.inputFiles[0], 'r').readlines()
process.source = cms.Source ('PoolSource', fileNames = cms.untracked.vstring(filelist) )

#output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root' if not options.outputFile else options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)

# Determine whether we are running on data or MC
isData = ('/store/data/' in process.source.fileNames[0] or options.Data==1)
isCustom = ('GluGluToAToZhToLLBB' in process.source.fileNames[0])
isReHLT = ('_reHLT_' in process.source.fileNames[0])
isPythiaLO = (True if (sample=='' or sample=='') else False)
isReRecoBCD = ('Run2016B-23Sep' in sample or 'Run2016C-23Sep' in sample or 'Run2016D-23Sep' in sample or options.Data==1)
isReRecoEF = ('Run2016E-23Sep' in sample or 'Run2016F-23Sep' in sample)
isReRecoG = ('Run2016G-23Sep' in sample)
isReRecoH = ('Run2016H-PromptReco' in sample)
isPromptReco = (('PromptReco' in sample) and (not isReRecoH))
print 'Running on', ('data' if isData else 'MC'),', sample is', sample
if isReHLT: print '-> re-HLT sample'
if isPythiaLO: print '-> Pythia LO sample'
if isReRecoBCD:  JECstring = "Summer16_23Sep2016BCDV3"
elif isReRecoEF: JECstring = "Summer16_23Sep2016EFV3"
elif isReRecoG:  JECstring = "Summer16_23Sep2016GV3"
elif isReRecoH:  JECstring = "Summer16_23Sep2016HV3"
elif isPromptReco: JECstring = "Spring16_25nsV6"
else: JECstring = "Summer16_23Sep2016BCDV3" # default value needed
JECstring_MC = "Summer16_23Sep2016V3" # unique default option for MC
JERstring = 'Spring16_25nsV10' # okay for moriond
print ' JEC:', JECstring, JECstring_MC
print ' JER:', JERstring
print ' maxEvents:', options.maxEvents

# Print trigger cut status
if not isData: print ' Trigger cut is', ('off' if options.tCut == 0 else 'on')

#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
import FWCore.PythonUtilities.LumiList as LumiList
if isData:
    #process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange() #12.9
     process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange() #36.26 (B,C,D,E,F,G reReco + H v2,v3 promptReco)

process.counter = cms.EDAnalyzer('CounterAnalyzer',
    lheProduct = cms.InputTag('externalLHEProducer' if not isCustom else 'source'),
    pythiaLOSample = cms.bool(True if isPythiaLO else False)
)

# Trigger filter
#import HLTrigger.HLTfilters.hltHighLevel_cfi
#triggerTag = 'HLT2' if isReHLT else 'HLT'
#process.HLTFilter = cms.EDFilter('HLTHighLevel',
 #   TriggerResultsTag = cms.InputTag('TriggerResults', '', triggerTag),
  #  HLTPaths = cms.vstring(
   #     'HLT_QuadJet45_TripleBTagCSV_p087_v*',
   #     'HLT_DoubleJet90_Double30_TripleBTagCSV_p087_v*',
   # ),
   # eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
   # andOr = cms.bool(True),    # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
   # throw = cms.bool(False)    # throw exception on unknown path names
#)

# METFilters
#process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
#process.BadPFMuonFilter.muons = cms.InputTag('slimmedMuons')
#process.BadPFMuonFilter.PFCandidates = cms.InputTag('packedPFCandidates')
#process.load('RecoMET.METFilters.BadPFMuonSummer16Filter_cfi')
#process.BadPFMuonSummer16Filter.muons = cms.InputTag("slimmedMuons")
#process.BadPFMuonSummer16Filter.PFCandidates = cms.InputTag("packedPFCandidates")

#process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
#process.BadChargedCandidateFilter.muons = cms.InputTag('slimmedMuons')
#process.BadChargedCandidateFilter.PFCandidates = cms.InputTag('packedPFCandidates')
#process.load('RecoMET.METFilters.BadChargedCandidateSummer16Filter_cfi')
#process.BadChargedCandidateSummer16Filter.muons = cms.InputTag('slimmedMuons')
#process.BadChargedCandidateSummer16Filter.PFCandidates = cms.InputTag('packedPFCandidates')

#MET corrections and uncertainties
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#if isData:
 #   jecFile = cms.string('{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_Uncertainty_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring))
#else:
#    jecFile = cms.string('{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_Uncertainty_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring))
#runMetCorAndUncFromMiniAOD(process,
                            #metType="PF",
                            #correctionLevel=["T1","Smear"],
                            #computeUncertainties=True,
                            #produceIntermediateCorrections=False,
                            #addToPatDefaultSequence=False,
 #                           isData=isData,
                            #onMiniAOD=True,
                            #reapplyJEC=reapplyJEC,
                            #reclusterJets=reclusterJets,
                            #jetSelection=jetSelection,
                            #recoMetFromPFCs=recoMetFromPFCs,
                            #autoJetCleaning=jetCleaning,
                            #manualJetConfig=manualJetConfig,
                            #jetFlavor=jetFlavor,
                            #jetCorLabelUpToL3=jetCorLabelL3,
                            #jetCorLabelL3Res=jetCorLabelRes,
                            #jecUnFile=jecFile,
                            #CHS=CHS,
                            #postfix=postfix,
  #                         )

if isData:
    filterString = "RECO"
else:
    filterString = "PAT"

# Primary vertex
#import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
#process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
#    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
#    minimumNDOF = cms.uint32(4) ,
#    maxAbsZ = cms.double(24), 
#    maxd0 = cms.double(2) 
#)


#-----------------------#
#        OBJECTS        #
#-----------------------#

#GLOBALTAG
process.load('Configuration.StandardSequences.Services_cff')#Lisa
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.trackerGeometryDB_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
GT = ''
#if isData:          GT = '80X_dataRun2_2016SeptRepro_v7'
GT = '90X_upgrade2023_realistic_v1'   #'80X_mcRun2_asymptotic_2016_TrancheIV_v8'
process.GlobalTag = GlobalTag(process.GlobalTag, GT)
print 'GlobalTag', GT

#electron/photon regression modules
#from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
#process = regressionWeights(process)
#process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

#electrons upstream modules
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
#ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
#                  'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                  'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff']

#for ele_idmod in ele_id_modules:
#    setupAllVIDIdsInModule(process,ele_idmod,setupVIDElectronSelection)

#muons upstream modules
#process.cleanedMuons = cms.EDProducer('PATMuonCleanerBySegments',
 #   src = cms.InputTag('slimmedMuons'),#('calibratedMuons'),
  #  preselection = cms.string('track.isNonnull'),
   # passthrough = cms.string('isGlobalMuon && numberOfMatches >= 2'),
   # fractionOfSharedSegments = cms.double(0.499)
#)

# Jet corrector https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrOnTheFly
#process.load('JetMETCorrections.Configuration.JetCorrectors_cff')



#### new
process.load('CommonTools.PileupAlgos.Puppi_cff')
process.puppi.candName           = cms.InputTag('packedPFCandidates')
process.puppi.vertexName         = cms.InputTag('offlineSlimmedPrimaryVertices')
process.puppi.useExistingWeights = cms.bool(True)
process.puppi.clonePackedCands   = cms.bool(True)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsPuppi'),
   jetCorrections = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], ''),
   pfCandidates = cms.InputTag('puppi'),
   btagDiscriminators = ["pfCombinedSecondaryVertexV2BJetTags", "pfCombinedInclusiveSecondaryVertexV2BJetTags"]
)
##############################


#https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer#ECAL_scale_and_resolution_correc
#process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#calibratedPatElectrons = cms.EDProducer('CalibratedPatElectronProducerRun2',
 #                                       # input collections
  #                                      electrons = cms.InputTag('slimmedElectrons'),
  #                                      gbrForestName = cms.string('gedelectron_p4combination_25ns'),
                                        # data or MC corrections
                                        # if isMC is false, data corrections are applied
  #                                      isMC = cms.bool(False) if isData else cms.bool(True),
   #                                     # set to True to get special 'fake' smearing for synchronization. Use JUST in case of synchronization
    #                                    isSynchronization = cms.bool(False),
     #                                   correctionFile = cms.string('80X_DCS05July_plus_Golden22')
      #                                  )

#-----------------------#
#        NTUPLE         #
#-----------------------#

process.ntuple = cms.EDAnalyzer('HHAnalyzer',
    genSet = cms.PSet(
        genProduct = cms.InputTag('generator'),
        lheProduct = cms.InputTag('externalLHEProducer' if not isCustom else 'source'),
        genParticles = cms.InputTag('prunedGenParticles'),
        pdgId = cms.vint32(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 36, 39, 1000022, 9100000, 9000001, 9000002, 9100012, 9100022, 9900032, 1023),
        samplesDYJetsToLL = cms.vstring(),
        samplesZJetsToNuNu = cms.vstring(),
        samplesWJetsToLNu = cms.vstring(),
        samplesDir = cms.string('%s/src/Analysis/ALPHA/data/Stitch/' % os.environ['CMSSW_BASE']),
        sample = cms.string( sample ),
        ewkFile = cms.string('%s/src/Analysis/ALPHA/data/scalefactors_v4.root' % os.environ['CMSSW_BASE']),
        applyEWK = cms.bool(True if sample.startswith('DYJets') or sample.startswith('WJets') else False),
        applyTopPtReweigth = cms.bool(True if sample.startswith('TT_') else False),
        pythiaLOSample = cms.bool(True if isPythiaLO else False)
    ),

  

    jetSet = cms.PSet(
        #https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Jet_Energy_Corrections_in_Run2
        jets = cms.InputTag('slimmedJetsPuppi'), # 'slimmedJetsPuppi'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(1), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet1pt = cms.double(20.),
        jet2pt = cms.double(20.),
        jeteta = cms.double(2.8),
        addQGdiscriminator = cms.bool(False), #True
        recalibrateJets = cms.bool(False),#True
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        smearJets = cms.bool(False),#True
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),#'offlineSlimmedPrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),
        jecShift = cms.int32(0),  # -1: down, 0: nominal, 1: up
        jecUncertaintyDATA = cms.string('{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_Uncertainty_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring) ),
        jecUncertaintyMC = cms.string('{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_Uncertainty_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC) ),
        jecCorrectorDATA = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L1FastJet_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2Relative_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L3Absolute_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2L3Residual_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
        ),
        jecCorrectorMC = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L1FastJet_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L2Relative_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L3Absolute_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
        ),
        massCorrectorDATA = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2Relative_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L3Absolute_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2L3Residual_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
        ),
        massCorrectorMC = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L2Relative_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L3Absolute_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
        ),
        massCorrectorPuppi = cms.string('%s/src/Analysis/ALPHA/data/puppiCorrSummer16.root' % os.environ['CMSSW_BASE']),
        reshapeBTag = cms.bool(False), #True
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('{0}/src/Analysis/ALPHA/data/CSVv2_Moriond17_B_H.csv'.format(os.environ['CMSSW_BASE'])),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        #met = cms.InputTag("patPFMetT1Smear"),#
        #metRecoil = cms.bool(False),
        #metRecoilMC = cms.string('{0}/src/Analysis/ALPHA/data/recoilfit_gjetsMC_Zu1_pf_v5.root'.format(os.environ['CMSSW_BASE'])),
        #metRecoilData = cms.string('{0}/src/Analysis/ALPHA/data/recoilfit_gjetsData_Zu1_pf_v5.root'.format(os.environ['CMSSW_BASE'])),
      #  metTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/MET_trigger_eff_data_SingleMuRunBH.root' % os.environ['CMSSW_BASE']),
        jerShift = cms.int32(0),  # -1: down, 0: nominal, 1: up
        jerNameRes = cms.string('{0}/src/Analysis/ALPHA/data/JER/{1}_MC_PtResolution_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JERstring)),
        jerNameSf = cms.string('{0}/src/Analysis/ALPHA/data/JER/{1}_MC_SF_AK4PFchs.txt'.format(os.environ['CMSSW_BASE'], JERstring)),
    
    ),
    fatJetSet = cms.PSet(
        jets = cms.InputTag('slimmedJetsAK8'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(1), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet1pt = cms.double(150.),
        jet2pt = cms.double(150.),
        jeteta = cms.double(2.5),
        addQGdiscriminator = cms.bool(True),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(True),
        recalibratePuppiMass = cms.bool(True),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        smearJets = cms.bool(True),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),       
        jecShift = cms.int32(0),  # -1: down, 0: nominal, 1: up      
        jecUncertaintyDATA = cms.string('{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_Uncertainty_AK8PFPuppi.txt'.format(os.environ['CMSSW_BASE'], JECstring)),
        jecUncertaintyMC = cms.string('{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_Uncertainty_AK8PFPuppi.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC)),
        jecCorrectorDATA = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L1FastJet_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2Relative_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L3Absolute_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2L3Residual_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
        ),
        jecCorrectorMC = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L1FastJet_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L2Relative_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L3Absolute_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
        ),
        massCorrectorDATA = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2Relative_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L3Absolute_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
            '{0}/src/Analysis/ALPHA/data/{1}_DATA/{1}_DATA_L2L3Residual_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring),
        ),
        massCorrectorMC = cms.vstring(
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L2Relative_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
            '{0}/src/Analysis/ALPHA/data/{1}_MC/{1}_MC_L3Absolute_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JECstring_MC),
        ),
        massCorrectorPuppi = cms.string('%s/src/Analysis/ALPHA/data/puppiCorrSummer16.root' % os.environ['CMSSW_BASE']),
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('{0}/src/Analysis/ALPHA/data/CSVv2_Moriond17_B_H.csv'.format(os.environ['CMSSW_BASE'])),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETs','','ALPHA'),#("patPFMetT1Smear"),#
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string(''),
        metRecoilData = cms.string(''),
        metTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/MET_trigger_eff_data_SingleMuRunBH.root' % os.environ['CMSSW_BASE']),
        jerShift = cms.int32(0),  # -1: down, 0: nominal, 1: up
        jerNameRes = cms.string('{0}/src/Analysis/ALPHA/data/JER/{1}_MC_PtResolution_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JERstring)),
        jerNameSf = cms.string('{0}/src/Analysis/ALPHA/data/JER/{1}_MC_SF_AK8PFchs.txt'.format(os.environ['CMSSW_BASE'], JERstring)),

    ),
   # writeNElectrons = cms.int32(0),
   # writeNMuons = cms.int32(0),
    writeNLeptons = cms.int32(2),
    writeNJets = cms.int32(0),
    writeNFatJets = cms.int32(1),
    histFile = cms.string('{0}/src/Analysis/ALPHA/data/HistList_HH.dat'.format(os.environ['CMSSW_BASE'])),
    verbose  = cms.bool(False),
)

process.seq = cms.Sequence(
        process.counter *

        #process.BadPFMuonFilter * process.BadPFMuonSummer16Filter *
        #process.BadChargedCandidateFilter * process.BadChargedCandidateSummer16Filter *        
        #process.fullPatMetSequence *
       # process.primaryVertexFilter *
        ##process.regressionApplication * #debug
        #process.egmGsfElectronIDSequence *
        #process.calibratedPatElectrons *
        #process.cleanedMuons *
        #process.ak4PFL2L3ResidualCorrectorChain *
        #process.QGTagger *
        process.ntuple
)

process.p = cms.Path(process.seq)
