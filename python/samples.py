#! /usr/bin/env python

sample = {
    'DoubleMuonRun2016B-PromptReco-v2' : {
        'nevents' : 19821894,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016B-PromptReco-v1' : {
        'nevents' : 2816842,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'DoubleEGRun2016B-PromptReco-v2' : {
        'nevents' : 29217360,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'METRun2016B-PromptReco-v1' : {
        'nevents' : 602650,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SingleMuonRun2016B-PromptReco-v2' : {
        'nevents' : 18201601,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SingleElectronRun2016B-PromptReco-v2' : {
        'nevents' : 33670715,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SingleElectronRun2016B-PromptReco-v1' : {
        'nevents' : 1422819,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'METRun2016B-PromptReco-v2' : {
        'nevents' : 4176245,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'DoubleEGRun2016B-PromptReco-v1' : {
        'nevents' : 5704443,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SinglePhotonRun2016B-PromptReco-v2' : {
        'nevents' : 1000000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'SinglePhotonRun2016B-PromptReco-v1' : {
        'nevents' : 1000000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'DoubleMuonRun2016B-PromptReco-v1' : {
        'nevents' : 4201017,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1' : {
        'nevents' : 49877138,
        'xsec'    : 6025.2,
        'kfactor' : 1.,
    },
    
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1' : {
        'nevents' : 9908534,
        'xsec'    : 61526.7,
        'kfactor' : 1.,
    },
    
    'TTTo2L2Nu_13TeV-powheg_v0_ext1-v1' : {
        'nevents' : 104607105,
        'xsec'    : 92.42,
        'kfactor' : 1.,
    },
    
    # Signal
    'BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-1600_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    'BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_PRIVATE-MC_v0-v1' : {
        'nevents' : 10000,
        'xsec'    : 1.,
        'kfactor' : 1.,
    },
    
}





samples = {
    'data_obs' : {
        'order' : 0,
        'files' : ['SingleMuonRun2016B-PromptReco-v2', 'DoubleEGRun2016B-PromptReco-v2',],
        'fillcolor' : 0,
        'fillstyle' : 1,
        'linecolor' : 1,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : "Data",
        'weight': 1.,
        'plot': True,
    },
    'DYJetsToLL' : {
        'order' : 1,
        'files' : ['DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1', ],
        'fillcolor' : 418,
        'fillstyle' : 1001,
        'linecolor' : 418,
        'linewidth' : 2,
        'linestyle' : 1,
        'label' : "Z(ll) + jets",
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
}
