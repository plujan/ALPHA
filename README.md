# ALPHA
A Light Post-Heppy Analyzer

## Git instructions
# Prerequisites
git account
git environment set

# Git instructions
In your working area, first set up the CMSSW release:
```bash
cmsrel CMSSW_8_0_26_patch1
cd CMSSW_8_0_26_patch1/src
cmsenv
git cms-init
```

Packages needed by ALPHA:
```bash
### MET + FILTERS
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes?rev=19#MET_Recipes
git cms-merge-topic cms-met:METRecipe_8020 -u
git cms-merge-topic cms-met:METRecipe_80X_part2 -u

### ELECTRON REGRESSION + SMEARER
# https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression?rev=20#Consistent_EGMSmearer
git cms-merge-topic -u rafaellopesdesa:EgammaAnalysis80_EGMSmearer_Moriond17_23Jan 
cd EgammaAnalysis/ElectronTools/data
git clone https://github.com/ECALELFS/ScalesSmearings.git
cd $CMSSW_BASE/src
# https://twiki.cern.ch/twiki/bin/view/CMS/EGMRegression?rev=20#Application_of_regression_outsid
git cms-merge-topic -u rafaellopesdesa:RegressionCheckNegEnergy 

### ELECTRON ID (CUT-BASED AND MVA)
# https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=41#Recipe_for_regular_users_for_8_0
git cms-merge-topic -u ikrav:egm_id_80X_v2

### ELECTRON ID (HEEP)
# https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2?rev=27#Instructions_to_checkout_HEEPV70
git cms-merge-topic -u Sam-Harper:HEEPV70VID_8010_ReducedCheckout 
git cms-merge-topic -u Sam-Harper:PackedCandNoPuppi 
mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/ 
git clone git@github.com:cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data 

### PHOTON ID (CUT-BASED AND MVA)
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2?rev=22#Recipes_for_regular_users
git cms-merge-topic -u ikrav:egm_id_80X_v3_photons
cd $CMSSW_BASE/external
cd slc6_amd64_gcc530/ # you may have a different architecture, this is just one example from lxplus
git clone https://github.com/ikrav/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_BASE/src
```

then, clone the ALPHA git repository:
```bash
### CLONE THE ALPHA REPO
cd $CMSSW_BASE/src
mkdir Analysis
cd Analysis
git clone -b REMINIAOD --single-branch https://github.com/CMS-PD/ALPHA
sh setup.sh
```

and compile the code:
```bash
### COMPILE
cd $CMSSW_BASE/src
scram b -j 8
```