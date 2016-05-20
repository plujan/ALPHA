#!/usr/bin/env python
import os, re
import multiprocessing
import commands
import math, time
import sys

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
    os.system('hadd '+name+'.root '+name+'/*/*.root')
pass


subdirs = [x for x in os.listdir(args.folder) if os.path.isdir(os.path.join(args.folder, x))]

os.chdir(args.folder)

for s in subdirs:
    print s
    p = multiprocessing.Process(target=hadd, args=(s,))
    jobs.append(p)
    p.start()

os.system('cd ..')
