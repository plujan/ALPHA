#!/bin/bash

#mkdir $CMSSW_BASE/src/PhysicsTools/
cp -r KinFitter/ $CMSSW_BASE/src/PhysicsTools/

cp $CMSSW_BASE/src/Analysis/ALPHA/data/ElectronIDValueMapProducer $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/python/ElectronIDValueMapProducer_cfi.py

cp $CMSSW_BASE/src/Analysis/ALPHA/data/ElectronMVAValueMapProducer $CMSSW_BASE/src/RecoEgamma/ElectronIdentification/python/ElectronMVAValueMapProducer_cfi.py

cp $CMSSW_BASE/src/Analysis/ALPHA/data/PhotonIDValueMapProducer $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/python/PhotonIDValueMapProducer_cfi.py

cp $CMSSW_BASE/src/Analysis/ALPHA/data/PhotonMVAValueMapProducer $CMSSW_BASE/src/RecoEgamma/PhotonIdentification/python/PhotonMVAValueMapProducer_cfi.py

