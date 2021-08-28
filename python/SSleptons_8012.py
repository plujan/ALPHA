import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import os

options = VarParsing ('analysis')
options.parseArguments()

# Determine sample name for MC stitching
sample = (options.inputFiles[0]).split('/')[-1].replace('.txt', '') if len(options.inputFiles) > 0 else ''
if sample=='list': sample = (options.inputFiles[0]).split('/')[-3]

# try...
#options.Pset('SkipEvent = cms.untracked.vstring('ProductNotFound')')

# isMC = True

process = cms.Process('ALPHA')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.threshold = 'ERROR'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000) )

# input
# default: if no filelist from command line, run on specified samples

if len(options.inputFiles) == 0:
   process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
#        'file:/lustre/cmswork/zucchett/CMSSW_8_0_5/src/00F0B3DC-211B-E611-A6A0-001E67248A39.root'
#        '/store/user/lbenato/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_PRIVATE-MC/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_10000ev/160515_095125/0000/BulkGraviton_ZZ_ZlepZhad_narrow_M800_13TeV-madgraph_MINIAOD_1.root'
# Wjets
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/50000/92B964DB-E023-E611-8FDF-0025901D4C94.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/50000/2079DB4D-E123-E611-88F2-0025901D4D54.root',
# DYjets
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/7443D5EA-3E2A-E611-AA3F-0CC47A4D7634.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/76942F47-3F2A-E611-BA56-0025905A60BE.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/DC246147-3F2A-E611-B7F0-0025905B857A.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/901D3229-3F2A-E611-AE6D-0CC47A78A426.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/425C3E4A-3F2A-E611-87AD-0025905A6104.root',
#
# WZ 3lept 
#(from WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1.txt file list)
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/56BA9F82-4B2B-E611-9F68-C4346BC7AAE0.root',
#(from WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_list.txt  file list)
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/52FB58B4-D81B-E611-B678-0CC47A57CD00.root',
#
# ttLL
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/70001/E01231E5-CE1C-E611-9EB8-002590E2F8B6.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/70001/6A622162-CE1C-E611-B245-0CC47A13D216.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/TTTo2L2Nu_13TeV-powheg/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/70000/546F0BC8-CD1C-E611-BF83-002590494E34.root',
#
# ttHadr
'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1/20000/0CC5C039-0632-E611-B597-0025905A48C0.root',
'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1/20000/40A77157-FE31-E611-9D84-0025905A60EE.root',
#
# double Muons, RunD
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016D/DoubleMuon/MINIAOD/PromptReco-v2/000/276/544/00000/58DF1812-024A-E611-9127-02163E01414B.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016D/DoubleMuon/MINIAOD/PromptReco-v2/000/276/542/00000/40641D94-E649-E611-8C75-02163E0146B5.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016D/DoubleMuon/MINIAOD/PromptReco-v2/000/276/543/00000/1410D0F3-EF49-E611-9EC6-02163E01374E.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016D/DoubleMuon/MINIAOD/PromptReco-v2/000/276/542/00000/3AC8848D-E649-E611-BC06-02163E0129A8.root',
# WpWp ------------------
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/DA134B54-8032-E611-AC58-047D7B881D10.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/5C7C555A-8032-E611-A812-0CC47A7E018E.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/2CFAEA8F-8032-E611-8CB5-047D7BD6DF00.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/3EA59203-8032-E611-8E24-0CC47A4D76BE.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/90000/C46CE521-8032-E611-B72D-0025901D4844.root',
# WmWm ------------------
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/00000/1C2655B1-5A38-E611-A3C0-C4346BC8C638.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/00000/5E5047B7-5A38-E611-8785-001E67F336A4.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/20000/16DA9C99-A438-E611-A10B-6CC2173BBA30.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/20000/36047DC6-3F3C-E611-AF07-003048CBA444.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/20000/3808F29E-A438-E611-B947-00266CFEFDE0.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v2/20000/E6FBDAE6-B338-E611-8EF3-C4346BC8C310.root',
# 
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/834/00000/5A8E357E-E83D-E611-818D-02163E01449E.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/834/00000/6C74613E-E13D-E611-A542-02163E013596.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/833/00000/6EA228A0-D83D-E611-B79A-02163E01349F.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/275/601/00000/4423F253-7B3A-E611-8707-02163E013706.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/276/097/00000/BE923E81-0741-E611-9BD7-02163E011A1C.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016C/DoubleMuon/MINIAOD/PromptReco-v2/000/276/097/00000/863B98EE-1041-E611-AF4D-02163E0129BD.root',
#
#double Muons, RunH
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/281/639/00000/9C6748A8-9885-E611-B2F5-02163E014305.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/281/639/00000/22163C41-B985-E611-81FD-02163E013406.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v2/000/281/680/00000/16DED305-9C85-E611-A36F-02163E0135C6.root',
#
# SingleMuons_3.root
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/E4DA02F6-1621-E611-8376-02163E014656.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/84C81BE6-1621-E611-BC41-02163E011EB0.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/0E9B9CE5-1621-E611-AF8E-02163E012474.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/D67EF0E5-1621-E611-B823-02163E011DE7.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/88C138FE-1621-E611-8B04-02163E01349C.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/BEC7F2E9-1621-E611-95E2-02163E014409.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/F21552E5-1621-E611-B733-02163E01441A.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/725/00000/541B0BFC-1621-E611-B383-02163E01220D.root',
#
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/AA551594-941B-E611-98A0-00304865C426.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/2AE73891-941B-E611-9BDB-0CC47A6C0682.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/9C180B94-941B-E611-9F0F-002590E2DD10.root',
#'dcap://t2-srm-02.lnl.infn.it/pnfs/lnl.infn.it/data/cms//store/mc/RunIISpring16MiniAODv2/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/1EE14398-941B-E611-A2C6-003048F5B69E.root',
    )
)
# production: read externally provided filelist
else:
    filelist = open(options.inputFiles[0], 'r').readlines()
    process.source = cms.Source ('PoolSource', fileNames = cms.untracked.vstring(filelist) )

#output
#process.TFileService = cms.Service('TFileService',
#    fileName = cms.string('SingleMuons_3.root'),
#    closeFileFast = cms.untracked.bool(True)
#)

#output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root' if len(options.outputFile) == 0 else options.outputFile),
    closeFileFast = cms.untracked.bool(True)
)

# Determine whether we are running on data or MC
isData = ('/store/data/' in process.source.fileNames[0])
#isCustom = ('GluGluToAToZhToLLBB' in process.source.fileNames[0])
#isCustom = ('WpWpJJ_13TeV-powheg' in process.source.fileNames[0])
# ***needed for WpWp, WmWm and ST_t-channel productions, 
#isCustom = True
isCustom = False
isDibosonInclusive = (True if (sample=='WW_TuneCUETP8M1_13TeV-pythia8_v0-v1' or sample=='WZ_TuneCUETP8M1_13TeV-pythia8_v0-v1' or sample=='ZZ_TuneCUETP8M1_13TeV-pythia8_v0-v1') else False)
print 'Running on', ('data' if isData else 'MC')

#-----------------------#
#        FILTERS        #
#-----------------------#

# JSON filter
import FWCore.PythonUtilities.LumiList as LumiList
# if not isMC:
if isData:
#  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-273450_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-274443_13TeV_PromptReco_Collisions16_JSON.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#   process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-276811_13TeV_PromptReco_Collisions16_JSON.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-279588_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
#   process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()
  process.source.lumisToProcess = LumiList.LumiList(filename = '%s/src/Analysis/ALPHA/data/JSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt' % os.environ['CMSSW_BASE']).getVLuminosityBlockRange()

process.counter = cms.EDAnalyzer('CounterAnalyzer',
    lheProduct = cms.InputTag('externalLHEProducer' if not isCustom else 'source'),
    pythiaLOSample = cms.bool(True if isDibosonInclusive else False),
)

# Trigger filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.HLTFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag = cms.InputTag('TriggerResults', '', 'HLT'),
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

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag('slimmedMuons')
process.BadPFMuonFilter.PFCandidates = cms.InputTag('packedPFCandidates')

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag('slimmedMuons')
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag('packedPFCandidates')

process.METFilter = cms.EDFilter('HLTHighLevel',
    TriggerResultsTag = cms.InputTag('TriggerResults', '', 'RECO'),
    HLTPaths = cms.vstring(
        'Flag_HBHENoiseFilter',
        'Flag_HBHENoiseIsoFilter',
        'Flag_EcalDeadCellTriggerPrimitiveFilter',
        'Flag_goodVertices',
        'Flag_eeBadScFilter',
        'Flag_globalTightHalo2016Filter',
    ),
    eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
    andOr = cms.bool(True),    # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
    throw = cms.bool(False)    # throw exception on unknown path names
)

# Primary vertex
import RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi
process.primaryVertexFilter = cms.EDFilter('GoodVertexFilter',
    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24), 
    maxd0 = cms.double(2) 
)


#-----------------------#
#        OBJECTS        #
#-----------------------#


#electrons upstream modules
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
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
process.cleanedMuons = cms.EDProducer('PATMuonCleanerBySegments',
    src = cms.InputTag('slimmedMuons'),#('calibratedMuons'),
    preselection = cms.string('track.isNonnull'),
    passthrough = cms.string('isGlobalMuon && numberOfMatches >= 2'),
    fractionOfSharedSegments = cms.double(0.499)
)

# Jet corrector https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrOnTheFly
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')


#quark gluon likelihood upstream modules
qgDatabaseVersion = 'v2b' # check https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion
from CondCore.DBCommon.CondDBSetup_cfi import *
QGPoolDBESSource = cms.ESSource('PoolDBESSource',
      CondDBSetup,
      toGet = cms.VPSet(),
      connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'),
)
for type in ['AK4PFchs','AK4PFchs_antib']:
    QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
        record = cms.string('QGLikelihoodRcd'),
        tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
        label  = cms.untracked.string('QGL_'+type)
    )))
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag('slimmedJets') # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs') # Other options: see https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGDataBaseVersion


process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
#process.selectedElectrons = cms.EDFilter('PATElectronSelector', 
                                         #src = cms.InputTag('slimmedElectrons'), 
                                         #cut = cms.string('pt > 5 && abs(eta)<2.5') 
                                         #)
calibratedPatElectrons = cms.EDProducer('CalibratedPatElectronProducerRun2',
                                        # input collections
                                        electrons = cms.InputTag('slimmedElectrons'),
                                        gbrForestName = cms.string('gedelectron_p4combination_25ns'),
                                        # data or MC corrections
                                        # if isMC is false, data corrections are applied
                                        isMC = cms.bool(False) if isData else cms.bool(True),
                                        # set to True to get special 'fake' smearing for synchronization. Use JUST in case of synchronization
                                        isSynchronization = cms.bool(False),
                                        correctionFile = cms.string('80X_DCS05July_plus_Golden22')
                                        )


#-----------------------#
#        NTUPLE         #
#-----------------------#

process.ntuple = cms.EDAnalyzer('SSleptons',
    genSet = cms.PSet(
        genProduct = cms.InputTag('generator'),
  #     lheProduct = cms.InputTag('externalLHEProducer'),
	lheProduct = cms.InputTag('externalLHEProducer' if not isCustom else 'source'),
        genParticles = cms.InputTag('prunedGenParticles'),
 # this is the list of particles to be stored in GenPVect (see code in ../plugins/SSlepton.cc	)
  #      pdgId = cms.vint32(1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 1000022, 9100000, 9000001, 9000002, 9100012, 9100022, 9900032, 1023),
        pdgId = cms.vint32( 11, 12, 13, 14, 15, 16, 21, 23, 24, 25, 211, 311, 321, 411,421,431, 443, 4122, 511, 521, 531, 5122 ),
        samplesDYJetsToLL = cms.vstring(
            'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v2',
            'DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
            'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
            'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
            'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
            'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
        ),
        samplesZJetsToNuNu = cms.vstring(),
        samplesWJetsToLNu = cms.vstring(),
        samplesDir = cms.string('%s/src/Analysis/ALPHA/data/Stitch/' % os.environ['CMSSW_BASE']),
        sample = cms.string( sample ),
        ewkFile = cms.string('%s/src/Analysis/ALPHA/data/scalefactors_v4.root' % os.environ['CMSSW_BASE']),
	applyEWK = cms.bool(True if sample.startswith('DYJets') or sample.startswith('WJets') else False),
        pythiaLOSample = cms.bool(True if isDibosonInclusive else False),
    ),
    pileupSet = cms.PSet(
        pileup = cms.InputTag('slimmedAddPileupInfo'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
#        dataFileName = cms.string('%s/src/Analysis/ALPHA/data/PU_69000.root' % os.environ['CMSSW_BASE']),
	dataFileName = cms.string('%s/src/Analysis/ALPHA/data/PU_71300.root' % os.environ['CMSSW_BASE']),
#        dataFileName = cms.string('%s/src/Analysis/ALPHA/data/Prod6.root' % os.environ['CMSSW_BASE']),
        mcFileName = cms.string('%s/src/Analysis/ALPHA/data/PU_MC.root' % os.environ['CMSSW_BASE']),
        dataName = cms.string('pileup'),
        mcName = cms.string('2016_25ns_SpringMC_PUScenarioV1'),
    ),
    triggerSet = cms.PSet(
        trigger = cms.InputTag('TriggerResults', '', 'HLT'),
        paths = cms.vstring('HLT_Mu45_eta2p1_v', 'HLT_Mu50_v', 'HLT_IsoMu20_v', 'HLT_Mu27_TkMu8_v', 'HLT_Mu30_TkMu11_v', 'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v', 'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v', 'HLT_Ele105_CaloIdVT_GsfTrkIdT_v', 'HLT_Ele23_WPLoose_Gsf_v', 'HLT_Ele27_WPLoose_Gsf_v', 'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v', 'HLT_DoubleEle33_CaloIdL_v'),
    ),
    electronSet = cms.PSet(
        electrons = cms.InputTag('slimmedElectrons'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        eleVetoIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto'),
        eleLooseIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose'),
        eleMediumIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium'),
        eleTightIdMap = cms.InputTag('egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight'),
        eleHEEPIdMap = cms.InputTag('egmGsfElectronIDs:heepElectronID-HEEPV60'),
        eleMVANonTrigMediumIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
        eleMVANonTrigTightIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
        eleMVATrigMediumIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90'),
        eleMVATrigTightIdMap = cms.InputTag('egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80'),
        eleSingleTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/SingleEleTriggerEff.root' % os.environ['CMSSW_BASE']),
        eleVetoIdFileName = cms.string('%s/src/Analysis/ALPHA/data/runB_p2_passingVeto_egammaEffi_txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleLooseIdFileName = cms.string('%s/src/Analysis/ALPHA/data/runB_p2_passingLoose_egammaEffi_txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/runB_p2_passingMedium_egammaEffi_txt_SF2D.root' % os.environ['CMSSW_BASE']),
        eleTightIdFileName = cms.string('%s/src/Analysis/ALPHA/data/runB_p2_passingTight_egammaEffi_txt_SF2D.root' % os.environ['CMSSW_BASE']),
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
        muons = cms.InputTag('cleanedMuons'),#('slimmedMuons'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        muonTrkFileName = cms.string('%s/src/Analysis/ALPHA/data/TrkEff.root' % os.environ['CMSSW_BASE']),
        muonIdFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonID_Z_RunBCD_prompt80X_7p65.root' % os.environ['CMSSW_BASE']),
        muonIsoFileName = cms.string('%s/src/Analysis/ALPHA/data/MuonIso_Z_RunBCD_prompt80X_7p65.root' % os.environ['CMSSW_BASE']),
        muonTrkHighptFileName = cms.string('%s/src/Analysis/ALPHA/data/trackHighPtID_effSF_80X.root' % os.environ['CMSSW_BASE']),
        muonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root' % os.environ['CMSSW_BASE']),
        doubleMuonTriggerFileName = cms.string('%s/src/Analysis/ALPHA/data/MuHLTEfficiencies_Run_2012ABCD_53X_DR03-2.root' % os.environ['CMSSW_BASE']),#obsolete
        muon1id = cms.int32(3), # 0: tracker id, 1: loose, 2: medium, 3: tight, 4: high pt
        muon2id = cms.int32(3),
        muon1iso = cms.int32(-1), # 0: tk high pT iso, 1: loose (0.25), 2: tight (0.15) (pfIso in cone 0.4), 3: no selection
        muon2iso = cms.int32(-1),
        muon1pt = cms.double(20.),
        muon2pt = cms.double(5.),
        useTuneP = cms.bool(True),
        doRochester = cms.bool(False),
    ),
    tauSet = cms.PSet(
        taus = cms.InputTag('slimmedTaus'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        taupt = cms.double(20.),
        taueta = cms.double(2.3),
        tauIdByDecayMode = cms.int32(1),# 0: not set, 1: old, 2: new
        tauIdByDeltaBetaIso = cms.int32(3),# 0: not set, 1: loose, 2: medium, 3: tight
        tauIdByMVAIso = cms.int32(4),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight, 5: V tight
        tauIdByMuonRejection = cms.int32(2),# 0: not set, 1: loose, 2: tight
        tauIdByElectronRejection = cms.int32(4),# 0: not set, 1: V loose, 2: loose, 3: medium, 4: tight
    ),
    photonSet = cms.PSet(
        photons = cms.InputTag('slimmedPhotons'),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        phoLooseIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose'),
        phoMediumIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium'),
        phoTightIdMap = cms.InputTag('egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight'),
        phoMVANonTrigMediumIdMap = cms.InputTag('egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90'),
        phoLooseIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Loosenumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Mediumnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoTightIdFileName = cms.string('%s/src/Analysis/ALPHA/data/Tightnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        phoMVANonTrigMediumIdFileName = cms.string('%s/src/Analysis/ALPHA/data/MVAnumbers.txt.egamma_SF2D.root' % os.environ['CMSSW_BASE']),
        photonid = cms.int32(1), # 1: loose, 2: medium, 3: tight, 4:MVA NonTrig medium
        photonpt = cms.double(20.),
    ),
    jetSet = cms.PSet(
        jets = cms.InputTag('slimmedJets'),#('slimmedJetsAK8'), #selectedPatJetsAK8PFCHSPrunedPacked
        jetid = cms.int32(1), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet1pt = cms.double(30.),
        jet2pt = cms.double(20.),
        jeteta = cms.double(4.5),
        addQGdiscriminator = cms.bool(True),
        recalibrateJets = cms.bool(True),
        recalibrateMass = cms.bool(False),
        recalibratePuppiMass = cms.bool(False),
        vertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        rho = cms.InputTag('fixedGridRhoFastjetAll'),        
        jecUncertaintyDATA = cms.string('%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt' % os.environ['CMSSW_BASE']),
        jecUncertaintyMC = cms.string('%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt' % os.environ['CMSSW_BASE']),
        jecCorrectorDATA = cms.vstring(
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
        ),
        jecCorrectorMC = cms.vstring(
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
        ),
        massCorrectorDATA = cms.vstring(
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_DATA/Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
        ),
        massCorrectorMC = cms.vstring(
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
            '%s/src/Analysis/ALPHA/data/Spring16_25nsV6_MC/Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt' % os.environ['CMSSW_BASE'],
        ),
        massCorrectorPuppi = cms.string('%s/src/Analysis/ALPHA/data/puppiJecCorr.root' % os.environ['CMSSW_BASE']),
        reshapeBTag = cms.bool(True),
        btag = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
        btagDB = cms.string('%s/src/Analysis/ALPHA/data/CSVv2.csv' % os.environ['CMSSW_BASE']),
        jet1btag = cms.int32(0), # 0: no selection, 1: loose, 2: medium, 3: tight
        jet2btag = cms.int32(0),
        met = cms.InputTag('slimmedMETs'),
        metRecoil = cms.bool(False),
        metRecoilMC = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsMC_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
        metRecoilData = cms.string('%s/src/Analysis/ALPHA/data/recoilfit_gjetsData_Zu1_pf_v5.root' % os.environ['CMSSW_BASE']),
    ),
    writeNElectrons = cms.int32(4),
    writeNMuons = cms.int32(4),
    writeNLeptons = cms.int32(0),
    writeNTaus = cms.int32(2),
    writeNPhotons = cms.int32(1),
    writeNJets = cms.int32(8),
 #   verbose  = cms.bool(True),
    verbose  = cms.bool(False),
)


#process.output = cms.OutputModule('PoolOutputModule',
#  splitLevel = cms.untracked.int32(0),
#  fileName = cms.untracked.string('deleteme.root'),
#  dataset = cms.untracked.PSet(
#    filterName = cms.untracked.string(''),
#    dataTier = cms.untracked.string('')
#  )
#)



if isData:
    process.seq = cms.Sequence(
        process.counter *
     #   process.HLTFilter *
	process.METFilter *
        process.BadPFMuonFilter *
        process.BadChargedCandidateFilter *
        process.primaryVertexFilter *
        process.egmGsfElectronIDSequence *
	process.calibratedPatElectrons *
        process.egmPhotonIDSequence *
        process.cleanedMuons *
     #   process.ak4PFL2L3ResidualCorrectorChain *
        process.QGTagger *
        process.ntuple
    )
else:
    process.seq = cms.Sequence(
        process.counter *
        process.BadPFMuonFilter *
        process.BadChargedCandidateFilter *	
        process.primaryVertexFilter *
        process.egmGsfElectronIDSequence *
	process.calibratedPatElectrons *
        process.egmPhotonIDSequence *
        process.cleanedMuons *
     #   process.ak4PFL2L3ResidualCorrectorChain *
        process.QGTagger *
        process.ntuple
    )

process.p = cms.Path(process.seq)
