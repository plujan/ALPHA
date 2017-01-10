# ALPHA 
A Light Post-Heppy Analyzer

# FROZEN TO CMSSW_8_0_12 (JAN-2017)

## Git instructions
# Prerequisites
git account
git environment set

# Git instructions
In your working area, first set up the CMSSW release:
```bash
cmsrel CMSSW_8_0_12
cd CMSSW_8_0_12/src/
cmsenv
git cms-init
```
Packages needed by ALPHA:
Merge the most recent MET filters and EGM smearing, scale, and IDs
```bash
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
git cms-merge-topic -u emanueledimarco:ecal_smear_fix_80X
git cms-addpkg EgammaAnalysis/ElectronTools
cd EgammaAnalysis/ElectronTools/data
git clone -b ICHEP2016_approval_7p65fb https://github.com/emanueledimarco/ScalesSmearings.git
cd $CMSSW_BASE/src
git cms-merge-topic -u ikrav:egm_id_80X_v1 
mkdir Analysis
cd Analysis
```
then, clone the ALPHA git repository:
```bash
git clone https://github.com/CMS-PD/ALPHA
```
and setup the code for KinFitter:
```bash
cd $CMSSW_BASE/src/Analysis/ALPHA
sh setup.sh
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
