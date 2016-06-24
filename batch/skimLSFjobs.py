#!/usr/bin/env python
import os, re
import multiprocessing
import commands
import math, time
import sys
from ROOT import TObject, TFile, TH1, TH1F

import argparse

parser = argparse.ArgumentParser(description='skim the LSF outputs into another tree')
parser.add_argument('folder', help='the folder containing the LSF output')
args = parser.parse_args()

if not os.path.exists(os.path.expandvars(args.folder)):
    print '--- ERROR ---'
    print '  \''+args.folder+'\' path not found'
    print '  please point to the correct path to the folder containing the LSF output' 
    print 
    exit()

jobs = []


def skim(name):
    oldFile = TFile(name+'.root', "READ")
    oldTree = oldFile.Get("ntuple/tree")
    
    newFile = TFile("Skim/"+name+'.root', "RECREATE")
    newFile.mkdir("ntuple/")
    newFile.cd("ntuple/")
    newTree = oldTree.CopyTree("V.pt>200 && FatJet1.pt>200")
    newTree.AutoSave()
    
    oldFile.Close()
    newFile.Close()



subdirs = [x for x in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, x))]

os.chdir(args.folder)

os.mkdir('Skim')

for s in subdirs:
#    skim(s)
    p = multiprocessing.Process(target=skim, args=(s,))
    jobs.append(p)
    p.start()

os.system('cd ..')

