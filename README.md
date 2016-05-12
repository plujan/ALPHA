# ALPHA
A Light Post-Heppy Analyzer

## Git instructions
In your working area, first set up the CMSSW release:
```bash
cmsrel CMSSW_8_0_4
cd CMSSW_8_0_4/src/
cmsenv
```
then, clone the git repository:
```bash
git clone https://github.com/CMS-PD/ALPHA
```
my knowledge of git ends here, but more will be added in the future!

## Run instructions
Compile the code:
```bash
scram b
```
and run it:
```bash
cmsRun python/ConfFile_cfg.py
```
