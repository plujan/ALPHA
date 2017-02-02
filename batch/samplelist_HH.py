#!/usr/bin/env python
# Sample lists for HH->4b

samplelists = [ 
####### DATA

    # MORIOND GOLDEN JSON -- 36 fb-1
    'BTagCSVRun2016B-23Sep2016-v2',
    'BTagCSVRun2016B-23Sep2016-v3',
    'BTagCSVRun2016C-23Sep2016-v1',
    'BTagCSVRun2016D-23Sep2016-v1',
    'BTagCSVRun2016E-23Sep2016-v1',
    'BTagCSVRun2016F-23Sep2016-v1',
    'BTagCSVRun2016G-23Sep2016-v1',
    'BTagCSVRun2016H-PromptReco-v2',
    'BTagCSVRun2016H-PromptReco-v3',

    # MORIOND GOLDEN JSON -- 36 fb-1  -- for trigger study
    'SingleMuonRun2016B-23Sep2016-v2',
    'SingleMuonRun2016B-23Sep2016-v3',
    'SingleMuonRun2016C-23Sep2016-v1',
    'SingleMuonRun2016D-23Sep2016-v1',
    'SingleMuonRun2016E-23Sep2016-v1',
    'SingleMuonRun2016F-23Sep2016-v1',
    'SingleMuonRun2016G-23Sep2016-v1',
    'SingleMuonRun2016H-PromptReco-v2',
    'SingleMuonRun2016H-PromptReco-v3',

####### MC
	'GluGluToHHTo4B_node_2_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_3_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_4_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_5_13TeV-madgraph-v1',
	#'GluGluToHHTo4B_node_6_13TeV-madgraph-v1', - not ready
	'GluGluToHHTo4B_node_7_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_8_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_9_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_10_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_11_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_12_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_13_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_box_13TeV-madgraph-v1',
	'GluGluToHHTo4B_node_SM_13TeV-madgraph-v1',
    'VBFHHTo4B_CV_1_C2V_1_C3_1_13TeV-madgraph-v1',
	
	# Background QCD
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
	
    # TOP
    'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-v1',

####### FOR TRIGGER STUDIES
    # Single Top leptonic
    #'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1-v1',
    #'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_withHLT_ext1-v1',

    # W Jets
    #'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v3',
    #'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    #'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    #'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    #'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    #'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
    #'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    #'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v2',
    #'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',

]
