# ALPHA
A Light Post-Heppy Analyzer

## Git instructions
# Prerequisites
git account
git environment set

# Git instructions
In your working area, first set up the CMSSW release:
```bash
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src/
cmsenv
git cms-init
```
Packages needed by ALPHA:
Merge the most recent MET filters and EGM smearing, scale, and IDs
```bash
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
### for future use, still not blessed by EGAMMA ###git cms-merge-topic -u rafaellopesdesa:Regression80XEgammaAnalysis
git cms-merge-topic -u shervin86:Moriond2017_JEC_energyScales
cd EgammaAnalysis/ElectronTools/data
# download the txt files with the corrections
git clone git@github.com:ECALELFS/ScalesSmearings.git
### EgammaID
cd $CMSSW_BASE/src
git cms-merge-topic -u ikrav:egm_id_80X_v2
### Egamma HEEP
git cms-merge-topic -u Sam-Harper:HEEPV70VID_8010_ReducedCheckout # HEEPV70VID_8010_ReducedCheckout
git cms-merge-topic -u ikrav:egm_id_80X_v3
git cms-merge-topic -u Sam-Harper:PackedCandNoPuppi
mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/ #we need this for the mva weights which runs in VID regardless if you need it or not
git clone git@github.com:cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data #we need this for the mva weights which runs in VID regardless if you need it or not
### PhotonID
cd $CMSSW_BASE/src
git cms-merge-topic -u ikrav:egm_id_80X_v3_photons

cd $CMSSW_BASE/src
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

## Run instructions
Compile the code:
```bash
scram b -j 8
```
and run it:
```bash
cmsRun python/Diboson.py
```

