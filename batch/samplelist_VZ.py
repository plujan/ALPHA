#!/usr/bin/env python

samplelists = [ 
####### DATA
    # 'SingleMuonRun2016B-23Sep2016-v1',
    # 'SingleMuonRun2016B-23Sep2016-v3',
    # 'SingleMuonRun2016C-23Sep2016-v1',
    # 'SingleMuonRun2016D-23Sep2016-v1',
    # 'SingleMuonRun2016E-23Sep2016-v1',
    # 'SingleMuonRun2016F-23Sep2016-v1',
    # 'SingleMuonRun2016G-23Sep2016-v1',
    # 'SingleMuonRun2016H-PromptReco-v2',
    # 'SingleMuonRun2016H-PromptReco-v3',

    'SingleMuonRun2016B-03Feb2017_ver1-v1',
    'SingleMuonRun2016B-03Feb2017_ver2-v2',
    'SingleMuonRun2016C-03Feb2017-v1',
    'SingleMuonRun2016D-03Feb2017-v1',
    'SingleMuonRun2016E-03Feb2017-v1',
    'SingleMuonRun2016F-03Feb2017-v1',
    'SingleMuonRun2016G-03Feb2017-v1',
    'SingleMuonRun2016H-03Feb2017_ver3-v1',
    'SingleMuonRun2016H-03Feb2017_ver2-v1',
    
    # 'SingleElectronRun2016B-23Sep2016-v2',
    # 'SingleElectronRun2016B-23Sep2016-v3',
    # 'SingleElectronRun2016C-23Sep2016-v1',
    # 'SingleElectronRun2016D-23Sep2016-v1',
    # 'SingleElectronRun2016E-23Sep2016-v1',
    # 'SingleElectronRun2016F-23Sep2016-v1',
    # 'SingleElectronRun2016G-23Sep2016-v1',
    # 'SingleElectronRun2016H-PromptReco-v2',
    # 'SingleElectronRun2016H-PromptReco-v3',

    'SingleElectronRun2016B-03Feb2017_ver1-v1',
    'SingleElectronRun2016B-03Feb2017_ver2-v2',
    'SingleElectronRun2016C-03Feb2017-v1',
    'SingleElectronRun2016D-03Feb2017-v1',
    'SingleElectronRun2016E-03Feb2017-v1',
    'SingleElectronRun2016F-03Feb2017-v1',
    'SingleElectronRun2016G-03Feb2017-v1',
    'SingleElectronRun2016H-03Feb2017_ver2-v1',
    'SingleElectronRun2016H-03Feb2017_ver3-v1',

######## MC
    #### Z(ll)+jets
    'DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext5-v1',
    'DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext5-v1',
    'DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext2-v1',
    'DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1',

    #### TTbar
    'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
    'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
  
    #### ST
    'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1-v1',
    'ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-v1',
    'ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1-v1',
    'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1-v1',                              # <-----
    'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1-v1',                          # <-----

    #### VV
    'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8-v1',
    'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8-v1',

    #### BulGravZZ
    'BulkGravToZZToZlepZhad_narrow_M-600_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph-v1',  
    'BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1400_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-2500_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph-v1',
    'BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph-v1',
    
    #### WprimeWZ
    'WprimeToWZToWhadZlep_narrow_M-600_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-800_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-1000_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-1200_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-1400_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-1600_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-1800_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-2000_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-2500_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-3000_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-3500_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-4000_13TeV-madgraph-v1',
    'WprimeToWZToWhadZlep_narrow_M-4500_13TeV-madgraph-v1',
  
]
