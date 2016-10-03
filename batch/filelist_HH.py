#!/usr/bin/env python
# Filelist for HH->4b

filelists = [ 
####### DATA
        'BTagCSVRun2016D-PromptReco-v2',

####### MC
	'GluGluToHHTo4B_node_2_13TeV-madgraph',
	'GluGluToHHTo4B_node_3_13TeV-madgraph',
	'GluGluToHHTo4B_node_4_13TeV-madgraph',
	'GluGluToHHTo4B_node_5_13TeV-madgraph',
	'GluGluToHHTo4B_node_6_13TeV-madgraph',
	'GluGluToHHTo4B_node_7_13TeV-madgraph',
	'GluGluToHHTo4B_node_8_13TeV-madgraph',
	'GluGluToHHTo4B_node_9_13TeV-madgraph',
	'GluGluToHHTo4B_node_10_13TeV-madgraph',
	'GluGluToHHTo4B_node_11_13TeV-madgraph',
	'GluGluToHHTo4B_node_12_13TeV-madgraph',
	'GluGluToHHTo4B_node_13_13TeV-madgraph',
	'GluGluToHHTo4B_node_box_13TeV-madgraph',
	'GluGluToHHTo4B_node_SM_13TeV-madgraph',
	
	# Background QCD
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v2',
	'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v3',
	'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
	
	# Background QCD bEnriched
	'QCD_bEnriched_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
	'QCD_bEnriched_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_bEnriched_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',

        # Background QCD BGenFilter
        'QCD_HT1000to1500_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT100to200_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT1500to2000_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT2000toInf_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT200to300_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT300to500_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT500to700_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',
        'QCD_HT700to1000_BGenFilter_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0-v1',

        # TOP
        'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_v0-v1',
        'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v2',
        'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v0-v1',

        # di-bosons
        'WZ_TuneCUETP8M1_13TeV-pythia8_v0-v1',
        'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8_v0-v1',

####### FOR TRIGGER STUDIES
        # Single Top leptonic
        'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1_v0-v1',

        # W Jets
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

]
