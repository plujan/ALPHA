# ALPHA
A Light Post-Heppy Analyzer

## Git instructions
In your working area, first set up the CMSSW release:
```bash
cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src/
cmsenv
git cms-init
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
mkdir Analysis
cd Analysis
```
then, clone the git repository:
```bash
git clone https://github.com/CMS-PD/ALPHA
```
See also the TWiki for developers git instructions: [https://twiki.cern.ch/twiki/bin/view/CMS/ALPHA](https://twiki.cern.ch/twiki/bin/view/CMS/ALPHA)

## Run instructions
Compile the code:
```bash
scram b -j 8
```
and run it:
```bash
cmsRun python/Diboson.py
```
