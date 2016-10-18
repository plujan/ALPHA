#!/usr/bin/env python

samplelists = [ 
####### DATA
  #### Run2016B
  #'SingleMuonRun2016B-PromptReco-v1', 
  #'SingleMuonRun2016B-PromptReco-v2',
  #'SingleElectronRun2016B-PromptReco-v1',
  #'SingleElectronRun2016B-PromptReco-v2',

  #### Run2016C
  #'SingleMuonRun2016C-PromptReco-v2', 
  #'SingleElectronRun2016C-PromptReco-v2',

  #### Run2016D
  #'SingleMuonRun2016D-PromptReco-v2', 
  #'SingleElectronRun2016D-PromptReco-v2',

  #### Run2016E
  #'SingleMuonRun2016E-PromptReco-v2', 
  #'SingleElectronRun2016E-PromptReco-v2',
    
######## MC
  #### Z(ll)+jets
  ##'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1',
  #'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
  #'DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
  #'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
  #'DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
  #'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
  #'DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
  #'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1',
  #'DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1',
  #### TTbar
  #'TTToSemiLeptonic_13TeV-powheg_ext1-v2',
  #'TTTo2L2Nu_13TeV-powheg_ext1-v1',
  #### ST
  #'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1-v1',
  #'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1-v1',
  #'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1-v1',
  #'ST_tW_antitop_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1-v1',
  #'ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1-v1',
  'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_withHLT_ext1-v1',
  #### VV
  #'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8-v1',
  #'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8-v1',
  #'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8-v1',
  #### Signal
  #'BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-1400_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-2500_13TeV-madgraph_reHLT-v2',
  #'BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-600_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-650_13TeV-madgraph_reHLT-v2',
  #'BulkGravToZZToZlepZhad_narrow_M-700_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-750_13TeV-madgraph_reHLT-v1',
  #'BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_reHLT-v1',  
   #### Signal
  #'WprimeToWZToWhadZlep_narrow_M-1000_13TeV-madgraph_reHLT-v3',
  #'WprimeToWZToWhadZlep_narrow_M-1200_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-1400_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-1600_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-1800_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-2000_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-2500_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-3000_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-3500_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-4000_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-4500_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-600_13TeV-madgraph_reHLT-v1',
  #'WprimeToWZToWhadZlep_narrow_M-800_13TeV-madgraph_reHLT-v1',
  
]
