#! /usr/bin/env python

#voms-proxy-init -voms cms
#

sample = {
    ########## DATA ##########
    'METRun2016B-23Sep2016-v2' : {
        'nevents' : 583427,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016B-23Sep2016-v3' : {
        'nevents' : 35987712,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016C-23Sep2016-v1' : {
        'nevents' : 17381222,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016D-23Sep2016-v1' : {
        'nevents' : 20947429,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016E-23Sep2016-v1' : {
        'nevents' : 22348402,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016F-23Sep2016-v1' : {
        'nevents' : 13319829,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016G-23Sep2016-v1' : {
        'nevents' : 26974131,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
        },
    'METRun2016H-PromptReco-v1' : {
        'nevents' : 13873,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'METRun2016H-PromptReco-v2' : {
        'nevents' : 39496995,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'METRun2016H-PromptReco-v3' : {
        'nevents' : 1129511,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },    

    'SingleMuonRun2016B-23Sep2016-v1' : {
        'nevents' : 2789243,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016B-23Sep2016-v3' : {
        'nevents' : 158145722,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016C-23Sep2016-v1' : {
        'nevents' : 67441308,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016D-23Sep2016-v1' : {
        'nevents' : 98017996,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016E-23Sep2016-v1' : {
        'nevents' : 90984718,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016F-23Sep2016-v1' : {
        'nevents' : 65425108,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016G-23Sep2016-v1' : {
        'nevents' : 149916849,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016H-PromptReco-v2' : {
        'nevents' : 171134793,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016H-PromptReco-v3' : {
        'nevents' : 4393222,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },

    'BTagCSVRun2016B-23Sep2016-v2' : {
        'nevents' : 1972666,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016B-23Sep2016-v3' : {
        'nevents' : 77890616,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016C-23Sep2016-v1' : {
        'nevents' : 30358567,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016D-23Sep2016-v1' : {
        'nevents' : 56527008,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016E-23Sep2016-v1' : {
        'nevents' : 60415444,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016F-23Sep2016-v1' : {
        'nevents' : 37608672,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016G-23Sep2016-v1' : {
        'nevents' : 100834056,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016H-PromptReco-v2' : {
        'nevents' : 64179785,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BTagCSVRun2016H-PromptReco-v3' : {
        'nevents' : 1695733,
        'xsec'    : 1.,
        'matcheff': 1.,
        'kfactor' : 1.,  
    },

    ########## MC ##########
    #DYJetsToNuNu_PtZ
    'DYJetsToNuNu_PtZ-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 21953584,
        'xsec'    : 593.9,#taken from mcm
        'matcheff': 0.349,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 5353639,
        'xsec'    : 170.4,#Check!!! Totally different on mcm!!!!
        'matcheff': 0.38,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 5356674,
        'xsec'    : 170.4,
        'matcheff': 0.38,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 1059634,
        'xsec'    : 6.636, #2.082,
        'matcheff': 0.4,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 1050705,
        'xsec'    : 0.9372, #0.2816,
        'matcheff': 0.41,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 1050592,
        'xsec'    : 0.9372, #0.2816,
        'matcheff': 0.41,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 1022595,
        'xsec'    : 0.1042, #0.02639,
        'matcheff': 0.42,
        'kfactor' : 1.,
    },
    'DYJetsToNuNu_PtZ-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 1024620,
        'xsec'    : 0.1042, #0.02639,
        'matcheff': 0.42,
        'kfactor' : 1.,
    },
    #WJetsToLNu_Pt*
    'WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 10089661,
        'xsec'    : 676.3,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 10088599,
        'xsec'    : 676.3,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 1001250,
        'xsec'    : 23.94,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 1000132,
        'xsec'    : 23.94,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext4-v1' : {
        'nevents' : 10021205,
        'xsec'    : 23.94,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 951713,
        'xsec'    : 3.031,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-400To600_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 988234,
        'xsec'    : 3.031,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-v1' : {
        'nevents' : 989482,
        'xsec'    : 0.4524,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'WJetsToLNu_Pt-600ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_ext1-v1' : {
        'nevents' : 985127,
        'xsec'    : 0.4524 ,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
##TT
    'TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-v1' : {
        'nevents' : 77229341,
        'xsec'    : 730.,#taken from McM: check from ichep!it is different!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
##TTJets backup
#    'TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8_backup-v1' : {
#        'nevents' : 43561608,
#        'xsec'    : 0.,#NO GENSIM PARENTS
#        'matcheff': 1.,
#        'kfactor' : 1.,
#    },
##ST
    'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_ext1-v1' : {
        'nevents' : 6933094,
        'xsec'    : 38.09,#taken from McM: check from ichep!it is different!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4-v1' : {
        'nevents' : 998276,
        'xsec'    : 38.09,#taken from McM: check from ichep!it is different!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_ext1-v1' : {
        'nevents' : 6952830,
        'xsec'    : 38.09,#taken from McM: check from ichep!it is different!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4-v1' : {
        'nevents' : 992024,
        'xsec'    : 38.09,#taken from McM: check from ichep!it is different!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
##WW
    'WW_TuneCUETP8M1_13TeV-pythia8-v1' : {
        'nevents' : 994012,
        'xsec'    : 63.21,#Taken from mcm!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
#    'WW_TuneCUETP8M1_13TeV-pythia8_ext1-v1' : {
#        'nevents' : 6987124,
#        'xsec'    : 63.21,#Taken from mcm!
#        'matcheff': 1.,
#        'kfactor' : 1.,
#    },
##WZ
    'WZ_TuneCUETP8M1_13TeV-pythia8-v1' : {
        'nevents' : 1000000,
        'xsec'    : 22.82,#Taken from mcm!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
#    'WZ_TuneCUETP8M1_13TeV-pythia8_ext1-v1' : {
#        'nevents' : 2995828,
#        'xsec'    : 22.82,#Taken from mcm!
#        'matcheff': 1.,
#        'kfactor' : 1.,
#    },
#ZZ
    'ZZ_TuneCUETP8M1_13TeV-pythia8-v1' : {
        'nevents' : 990064,
        'xsec'    : 10.32,#Taken from mcm!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'ZZ_TuneCUETP8M1_13TeV-pythia8_ext1-v1' : {
        'nevents' : 998034,
        'xsec'    : 10.32,#Taken from mcm!
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    # FIXME - check xs for QCD (taken from ichep)
    'QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 80684349,
        'xsec'    : 1712000.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 18722416,
        'xsec'    : 1712000.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 38857977,
        'xsec'    : 1712000.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 17035891,
        'xsec'    : 347700.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 37502012,
        'xsec'    : 347700.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 18929951,
        'xsec'    : 32100.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v2' : {
        'nevents' : 43341392,
        'xsec'    : 32100.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 15629253,
        'xsec'    : 6831.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 29783527,
        'xsec'    : 6831.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },   
    'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 4767100,
        'xsec'    : 1207.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 10360193,
        'xsec'    : 1207.,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 3970819,
        'xsec'    : 119.9,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 7855883,
        'xsec'    : 119.9,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8-v1' : {
        'nevents' : 1991645,
        'xsec'    : 25.24,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_ext1-v1' : {
        'nevents' : 4047360,
        'xsec'    : 25.24,
        'matcheff': 1.,
        'kfactor' : 1.,
    },


    ### Signal ###
    #ZHadZinv
    'BulkGravToZZToZhadZinv_narrow_M-600_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-800_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-1000_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-1200_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-1400_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-1600_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-1800_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-2000_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-2500_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-3000_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-3500_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-4000_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZhadZinv_narrow_M-4500_13TeV-madgraph-v1' : {
        'nevents' : 100000,
        'xsec'    : 2.*0.6991*0.20,
        'matcheff': 1.,
        'kfactor' : 1.,
    },

    #HH FIXME - bsm xs
    'GluGluToHHTo4B_node_2_13TeV-madgraph-v1' : {
        'nevents' : 299600,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_3_13TeV-madgraph-v1' : {
        'nevents' : 299800,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_4_13TeV-madgraph-v1' : {
        'nevents' : 297000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_5_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_6_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_7_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_8_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_9_13TeV-madgraph-v1' : {
        'nevents' : 298200,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_10_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_11_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_12_13TeV-madgraph-v1' : {
        'nevents' : 299400,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_13_13TeV-madgraph-v1' : {
        'nevents' : 300000,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_box_13TeV-madgraph-v1' : {
        'nevents' : 278600,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'GluGluToHHTo4B_node_SM_13TeV-madgraph-v1' : {
        'nevents' : 299800,
        'xsec'    : 0.03345*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    'VBFHHTo4B_CV_1_C2V_1_C3_1_13TeV-madgraph-v1' : {
        'nevents' : 297969,
        'xsec'    : 1.*0.577*0.577,
        'matcheff': 1.,
        'kfactor' : 1.,
    },
    
}





samples = {
    'data_obs' : {
        'order' : 0,
        'files' : ['METRun2016B-23Sep2016-v2', 'METRun2016B-23Sep2016-v3', 'METRun2016C-23Sep2016-v1', 'METRun2016D-23Sep2016-v1', 'METRun2016E-23Sep2016-v1', 'METRun2016F-23Sep2016-v1', 'METRun2016G-23Sep2016-v1', 'METRun2016H-PromptReco-v1', 'METRun2016H-PromptReco-v2', 'METRun2016H-PromptReco-v3'],
        'fillcolor' : 0,
        'fillstyle' : 1,
        'linecolor' : 1,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : "Data",
        'weight': 1.,
        'plot': True,
    },
    
    # Dummy entry for background sum
    'BkgSum' : {
        'order' : 0,
        'files' : [],
        'fillcolor' : 1,
        'fillstyle' : 3003,
        'linecolor' : 1,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : "MC stat.",
        'weight': 1.,
        'plot': True,
    },
    
    ### BulkG -> Zhad Zinv
    'XZZInv_M600' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-600_13TeV-madgraph-v1'],
        'fillcolor' : 51,
        'fillstyle' : 3005,
        'linecolor' : 51,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 600 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M800' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-800_13TeV-madgraph-v1'],
        'fillcolor' : 52,
        'fillstyle' : 3005,
        'linecolor' : 52,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 800 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M1000' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-1000_13TeV-madgraph-v1'],
        'fillcolor' : 53,
        'fillstyle' : 3005,
        'linecolor' : 53,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 1000 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M1200' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-1200_13TeV-madgraph-v1'],
        'fillcolor' : 55,
        'fillstyle' : 3005,
        'linecolor' : 55,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 1200 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M1400' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-1400_13TeV-madgraph-v1'],
        'fillcolor' : 57,
        'fillstyle' : 3005,
        'linecolor' : 57,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 1400 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M1600' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-1600_13TeV-madgraph-v1'],
        'fillcolor' : 59,
        'fillstyle' : 3005,
        'linecolor' : 59,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 1600 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M1800' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-1800_13TeV-madgraph-v1'],
        'fillcolor' : 61,
        'fillstyle' : 3005,
        'linecolor' : 61,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 1800 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M2000' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-2000_13TeV-madgraph-v1'],
        'fillcolor' : 63,
        'fillstyle' : 3005,
        'linecolor' : 63,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 2000 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M2500' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-2500_13TeV-madgraph-v1'],
        'fillcolor' : 65,
        'fillstyle' : 3005,
        'linecolor' : 65,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 2500 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M3000' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-3000_13TeV-madgraph-v1'],
        'fillcolor' : 67,
        'fillstyle' : 3005,
        'linecolor' : 67,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 3000 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M3500' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-3500_13TeV-madgraph-v1'],
        'fillcolor' : 69,
        'fillstyle' : 3005,
        'linecolor' : 69,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 3500 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M4000' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-4000_13TeV-madgraph-v1'],
        'fillcolor' : 71,
        'fillstyle' : 3005,
        'linecolor' : 71,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 4000 GeV",
        'weight': 1.,
        'plot': True,
    },
    'XZZInv_M4500' : {
        'order' : 1001,
        'files' : ['BulkGravToZZToZhadZinv_narrow_M-4500_13TeV-madgraph-v1'],
        'fillcolor' : 73,
        'fillstyle' : 3005,
        'linecolor' : 73,
        'linewidth' : 3,
        'linestyle' : 1,
        'label' : "m_{X} = 4500 GeV",
        'weight': 1.,
        'plot': True,
    },
    
}
