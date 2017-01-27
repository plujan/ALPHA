#! /usr/bin/env python

#voms-proxy-init -voms cms
#

sample = {
    #Re-Reco
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
    # Signal
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
