#!/bin/bash

#mkdir $CMSSW_BASE/src/PhysicsTools/
cp -r KinFitter/ $CMSSW_BASE/src/PhysicsTools/

cp /lustre/cmswork/pazzini/ALPHA/ElectronIDValueMapProducer_cfi.py $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/python/ElectronIDValueMapProducer_cfi.py

cp /lustre/cmswork/pazzini/ALPHA/ElectronMVAValueMapProducer_cfi.py $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/python/ElectronMVAValueMapProducer_cfi.py

cp /lustre/cmswork/pazzini/ALPHA/PhotonIDValueMapProducer_cfi.py $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/python/PhotonIDValueMapProducer_cfi.py

cp /lustre/cmswork/pazzini/ALPHA/PhotonMVAValueMapProducer_cfi.py $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/python/PhotonMVAValueMapProducer_cfi.py