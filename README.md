# ALPHA
A Light Post-Heppy Analyzer

## Git instructions
In your working area, first set up the CMSSW release:
```cmsrel CMSSW_8_0_4
cmsenv
mkdir Analysis
cd Analysis
```
then, clone the git repository:
```git clone https://github.com/CMS-PD/ALPHA
```
my knowledge of git ends here.

## Run instructions
Compile the code:
```scram b
```
and run it:
```cmsRun python/ConfFile_cfg.py
```
