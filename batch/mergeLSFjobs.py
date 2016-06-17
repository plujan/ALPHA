#!/usr/bin/env python
import os, re
import multiprocessing
import commands
import math, time
import sys
from ROOT import TFile, TH1, TH1F
from Analysis.ALPHA.samples import sample

LUMI        = 2070 # in pb-1

# use the following lists to include/exclude samples to be merged

blacklist = []

whitelist = [
    'SingleMuonRun2016B-PromptReco-v1', 
    'SingleMuonRun2016B-PromptReco-v2',
#    'SingleElectronRun2016B-PromptReco-v1',
#    'SingleElectronRun2016B-PromptReco-v2',
    'METRun2016B-PromptReco-v1',
    'METRun2016B-PromptReco-v2',
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'TTToSemiLeptonic_13TeV-powheg_v0_ext1_red1-v1',
    'TTTo2L2Nu_13TeV-powheg_v0_ext1_red1-v1',
    'VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v0-v1',
     # Signal
    'BulkGravToZZToZlepZhad_narrow_M-600_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-800_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1600_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_PRIVATE-MC_v0-v1',
    'BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_PRIVATE-MC_v0-v1',
    # AZh
    'GluGluToAToZhToLLBB_M225_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M250_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M275_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M300_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M325_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M350_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M400_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M500_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M600_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M700_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M800_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M900_13TeV-amcatnlo',
    'GluGluToAToZhToLLBB_M1000_13TeV-amcatnlo',
]



########## DO NOT TOUCH BELOW THIS POINT ##########

import argparse

parser = argparse.ArgumentParser(description='combine the LSF outputs into one tree')
parser.add_argument('folder', help='the folder containing the LSF output')
args = parser.parse_args()

if not os.path.exists(os.path.expandvars(args.folder)):
    print '--- ERROR ---'
    print '  \''+args.folder+'\' path not found'
    print '  please point to the correct path to the folder containing the LSF output' 
    print 
    exit()

jobs = []

def hadd(name):
    if len(whitelist)>0 and not name in whitelist: return
    if len(blacklist)>0 and name in blacklist: return
    
    os.system('hadd '+name+'.root '+name+'/*/*.root')
    # Add weight
    file = TFile(name+'.root', "UPDATE")
    tree = file.Get("ntuple/tree")
    if 'Run2016' in name: weight = 1.
    else:
        nevents = file.Get("counter/c_nEvents").GetBinContent(1)
        xs = sample[name]['xsec'] * sample[name]['kfactor']
        weight = LUMI * xs / nevents
    tree.SetWeight(weight)
    tree.AutoSave()
    file.Close()
pass


subdirs = [x for x in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, x))]

os.chdir(args.folder)

for s in subdirs:
    print s
    p = multiprocessing.Process(target=hadd, args=(s,))
    jobs.append(p)
    p.start()

os.system('cd ..')
