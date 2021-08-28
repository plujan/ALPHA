#!/usr/bin/env python

samplelists = [ 
####### DATA
    'SingleMuonRun2016B-23Sep2016-v3',
    'SingleMuonRun2016C-23Sep2016-v1',
    'SingleMuonRun2016D-23Sep2016-v1',
    'SingleMuonRun2016E-23Sep2016-v1',
    'SingleMuonRun2016F-23Sep2016-v1',
    'SingleMuonRun2016G-23Sep2016-v1',
    'SingleMuonRun2016H-PromptReco-v3',
    
    'SingleElectronRun2016B-23Sep2016-v3',
    'SingleElectronRun2016C-23Sep2016-v1',
    'SingleElectronRun2016D-23Sep2016-v1',
    'SingleElectronRun2016E-23Sep2016-v1',
    'SingleElectronRun2016F-23Sep2016-v1',
    'SingleElectronRun2016G-23Sep2016-v1',
    'SingleElectronRun2016H-PromptReco-v3',

    'DoubleMuonRun2016B-23Sep2016-v3',
    'DoubleMuonRun2016C-23Sep2016-v1',
    'DoubleMuonRun2016D-23Sep2016-v1',
    'DoubleMuonRun2016E-23Sep2016-v1',
    'DoubleMuonRun2016F-23Sep2016-v1',
    'DoubleMuonRun2016G-23Sep2016-v1',
    'DoubleMuonRun2016H-PromptReco-v3',
    
######## MC
    #### QCD
    'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v2',
    'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
    'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    
    # W->lnu+jets
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v3',
    'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v2',
    'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',

    #### Z->ll+jets
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',


    #### TTbar and single top
    'TT_TuneCUETP8M1_13TeV-powheg-pythia8-v2',
    'TTToSemiLeptonic_13TeV-powheg_v0_ext1-v1',
    'TTTo2L2Nu_13TeV-powheg_v0_ext1-v1',

    'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_v0-v1',
    'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v1',
    'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v1',
    'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_v0-v1',
    'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1_v0-v1',

    #### VV
    'VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v0-v1',
    'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8_v0-v1',
    'ZZ_TuneCUETP8M1_13TeV-pythia8_v0-v1',
    'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v0-v1',
    'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8_v0-v2',
    'WpWpJJ_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v1',
    'WmWmJJ_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v1'

]
