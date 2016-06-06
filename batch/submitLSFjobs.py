#!/usr/bin/env python
import os, re
import commands
import math, time
import sys
from Analysis.ALPHA.samples import sample

########## FILELIST ##########

filelists = [ 
    'SingleMuonRun2016B-PromptReco-v1', 
    'SingleMuonRun2016B-PromptReco-v2',
    'DoubleMuonRun2016B-PromptReco-v1',
    'DoubleMuonRun2016B-PromptReco-v2',
    'SingleElectronRun2016B-PromptReco-v1',
    'SingleElectronRun2016B-PromptReco-v2',
    'DoubleEGRun2016B-PromptReco-v1',
    'DoubleEGRun2016B-PromptReco-v2',
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v0_ext1-v1',
    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v0-v1',
    'TTTo2L2Nu_13TeV-powheg_v0_ext1-v1',
#    # Signal
##    'BulkGravToZZToZlepZhad_narrow_M-1000_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-1200_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-1600_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-1800_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-2000_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-3000_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-3500_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-4000_13TeV-madgraph_PRIVATE-MC_v0-v1',
##    'BulkGravToZZToZlepZhad_narrow_M-4500_13TeV-madgraph_PRIVATE-MC_v0-v1',
]

########## DO NOT TOUCH BELOW THIS POINT ##########

########## OPTIONS ##########

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-b', '--base',         action='store', type='string', dest='base',         default='$CMSSW_BASE/src/Analysis/ALPHA/')
parser.add_option('-o', '--output',       action='store', type='string', dest='output',       default='output')
parser.add_option('-c', '--cfg',          action='store', type='string', dest='cfg',          default='python/Diboson.py')
parser.add_option('-q', '--queue',        action='store', type='string', dest='queue',        default='local-cms-short')
parser.add_option('-m', '--maxlsftime',   action='store', type='int',    dest='maxlsftime',   default=5)
parser.add_option('-e', '--eventspersec', action='store', type='int',    dest='eventspersec', default=100)
(options, args) = parser.parse_args()

print
if not os.path.exists(os.path.expandvars(options.base)):
    print '--- ERROR ---'
    print '  \''+options.base+'\' path not found expanding '+options.base
    print '  please point to the correct path to ALPHA/ using option \'-b PATH-TO-ALPHA\'' 
    print 
    exit()

if not os.path.exists(os.path.expandvars(options.base+options.cfg)):
    print '--- ERROR ---'
    print '  \''+options.cfg+'\' file not found in '+options.base+options.cfg
    print '  please point to a valid cfg file using option \'-c CFG-FILENAME\'' 
    print 
    exit()

path = os.getcwd()
if os.path.exists(options.output):
    print '--- ERROR ---'
    print '  \''+options.output+'\' folder already exists!'
    print '  please delete it or use a different name using option \'-o FOLDER-NAME\'' 
    print 
    exit()
os.system('mkdir '+options.output)


########## LOOP ON FILELISTS ##########
for l in filelists:
    if not l in sample:
        print l
        continue
    dir= 'Run2016' if 'Run2016' in l else 'Spring16'
    file=open(os.path.expandvars(options.base+'filelists/'+dir+'/'+l+'.txt'),'r')
    filelist = file.readlines()
    splitting= max(int(float(sample[l]['nevents'])/(options.maxlsftime*3600*options.eventspersec)),1)
    njobs    = int(len(filelist)/splitting)+1
    sublists = [filelist[i:i+njobs] for i in range(0, len(filelist), njobs)]
    print '\nSplitting',l,'in',len(sublists),'chunk(s) of approximately',njobs,'files each'
    lfold = options.output+'/'+l
    os.system('mkdir '+lfold)

    ########## LOOP ON LSF JOB ##########
    for x in range(len(sublists)):
        lsubfold = lfold+'/'+str(x).zfill(4)
        os.system('mkdir '+lsubfold)
        os.chdir(lsubfold)
        splitlist=open('list.txt','w')  
        splitlist.write(''.join(str(x) for x in sublists[x]))
        splitlist.close()
        
        with open('job.sh', 'w') as fout:
            #fout.write('#!/bin/bash\n')
            #fout.write('#BSUB -J '+l+'_'+str(x).zfill(4)+'\n')
            fout.write('echo "PWD:"\n')
            fout.write('pwd\n')
            fout.write('export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n')
            fout.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
            fout.write('echo "environment:"\n')
            fout.write('echo\n')
            fout.write('env > local.env\n')
            fout.write('env\n')
            fout.write('# ulimit -v 3000000 # NO\n')
            fout.write('echo "copying job dir to worker"\n')
            fout.write('eval `scram runtime -sh`\n')
            fout.write('ls\n')
            fout.write('echo "running"\n')
            fout.write('cmsRun '+options.base+options.cfg+' inputFiles='+options.base+lsubfold+'/list.txt\n')
            #fout.write('cmsRun '+options.base+options.cfg+' outputFile_clear outputFile='+outname+'.root inputFiles_clear inputFiles_load='+lsubfold+'/list.txt\n')
            fout.write('exit $?\n') 
            fout.write('echo ""\n')
        os.system('chmod 755 job.sh')
        
        ########## SEND JOB ON LSF QUEUE ##########
        os.system('bsub -q '+options.queue+' -o logs < job.sh')
        #print 'filelist ' + l + ' - job nr ' + str(x).zfill(4) + ' -> submitted'
        os.chdir('../../../')
   
print
print 'CURRENT JOB SUMMARY:'
os.system('bjobs')
print
