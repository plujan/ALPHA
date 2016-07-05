#!/usr/bin/env python
import os, re
import multiprocessing
import commands
import math, time
import sys
from ROOT import TObject, TFile, TH1, TH1F
from array import array

import argparse

parser = argparse.ArgumentParser(description='skim the LSF outputs into another tree')
parser.add_argument('folder', help='the folder containing the LSF output')
args = parser.parse_args()

if not os.path.exists(os.path.expandvars(args.folder)):
    print '--- ERROR ---'
    print '  \''+args.folder+'\' path not found'
    print '  please point to the correct path to the folder containing the LSF output' 
    print 
    exit()

jobs = []


def skim(name):
    
    isMC = array('b', [0]) 
    isZtoEE = array('b', [0]) 
    isZtoMM = array('b', [0]) 
    isTtoEM = array('b', [0]) 
    isWtoEN = array('b', [0]) 
    isWtoMN = array('b', [0]) 
    isZtoNN = array('b', [0]) 
    isMerged = array('b', [0]) 
    isResolved = array('b', [0]) 
    
    EventNumber = array('l', [-1]) 
    RunNumber = array('l', [-1]) 
    LumiNumber = array('l', [-1]) 

    EventWeight = array('f', [-1.0]) 
    StitchWeight = array('f', [-1.0]) 
    ZewkWeight = array('f', [-1.0]) 
    WewkWeight = array('f', [-1.0]) 
    PUWeight = array('f', [-1.0]) 
    TriggerWeight = array('f', [-1.0]) 
    LeptonWeight = array('f', [-1.0]) 

    FacWeightUp = array('f', [-1.0]) 
    FacWeightDown = array('f', [-1.0]) 
    RenWeightUp = array('f', [-1.0]) 
    RenWeightDown = array('f', [-1.0]) 
    ScaleWeightUp = array('f', [-1.0]) 
    ScaleWeightDown = array('f', [-1.0]) 

    nPV = array('d', [-1]) 
    nElectrons = array('d', [-1]) 
    nMuons = array('d', [-1]) 
    nTaus = array('d', [-1]) 
    nPhotons = array('d', [-1]) 
    nJets = array('d', [-1]) 
    nFatJets = array('d', [-1]) 
    nBTagJets = array('d', [-1]) 

    MaxJetBTag = array('f', [-1.0]) 
    MaxFatJetBTag = array('f', [-1.0]) 
    MinJetMetDPhi = array('f', [-1.0]) 

    CosThetaStar = array('f', [-1.0]) 
    CosTheta1 = array('f', [-1.0]) 
    CosTheta2 = array('f', [-1.0]) 
    Phi = array('f', [-1.0]) 
    Phi1 = array('f', [-1.0]) 
    AngularLD = array('f', [-1.0]) 

    Lepton1_isMuon = array('b', [0]) 
    Lepton1_isElectron = array('b', [0]) 
    Lepton1_isLoose = array('b', [0]) 
    Lepton1_isHighPt = array('b', [0]) 
    Lepton1_isTrackerHighPt = array('b', [0]) 
    Lepton1_isTight = array('b', [0]) 
    Lepton1_pt = array('f', [-1.0]) 
    Lepton1_trkIso = array('f', [-1.0]) 

    Lepton2_isMuon = array('b', [0]) 
    Lepton2_isElectron = array('b', [0]) 
    Lepton2_isLoose = array('b', [0]) 
    Lepton2_isHighPt = array('b', [0]) 
    Lepton2_isTrackerHighPt = array('b', [0]) 
    Lepton2_isTight = array('b', [0]) 
    Lepton2_pt = array('f', [-1.0]) 
    Lepton2_trkIso = array('f', [-1.0]) 

    MEt_pt = array('f', [-1.0]) 

    V_pt = array('f', [-1.0]) 
    V_dPhi = array('f', [-1.0]) 
    V_mass = array('f', [-1.0]) 
    V_tmass = array('f', [-1.0]) 

    X_pt = array('f', [-1.0]) 
    X_dPhi = array('f', [-1.0]) 
    X_mass = array('f', [-1.0]) 
    X_tmass = array('f', [-1.0]) 

    FatJet1_isTight = array('b', [0]) 

    FatJet1_pt = array('f', [-1.0]) 
    FatJet1_prunedMass = array('f', [-1.0]) 
    FatJet1_softdropMass = array('f', [-1.0]) 
    FatJet1_softdropPuppiMass = array('f', [-1.0]) 
    FatJet1_prunedMassCorr = array('f', [-1.0]) 
    FatJet1_softdropMassCorr = array('f', [-1.0]) 
    FatJet1_softdropPuppiMassCorr = array('f', [-1.0]) 
    FatJet1_tau21 = array('f', [-1.0]) 
    FatJet1_CSV1 = array('f', [-1.0]) 
    FatJet1_CSV2 = array('f', [-1.0]) 

    oldFile = TFile(name, "READ")
    oldTree = oldFile.Get("ntuple/tree")

    oldTree.SetBranchAddress("isMC", isMC)
    oldTree.SetBranchAddress("EventNumber", EventNumber)
    oldTree.SetBranchAddress("LumiNumber", LumiNumber)
    oldTree.SetBranchAddress("RunNumber", RunNumber)
    oldTree.SetBranchAddress("EventWeight", EventWeight)
    oldTree.SetBranchAddress("FacWeightUp", FacWeightUp)
    oldTree.SetBranchAddress("FacWeightDown", FacWeightDown)
    oldTree.SetBranchAddress("RenWeightUp", RenWeightUp)
    oldTree.SetBranchAddress("RenWeightDown", RenWeightDown)
    oldTree.SetBranchAddress("ScaleWeightUp", ScaleWeightUp)
    oldTree.SetBranchAddress("ScaleWeightDown", ScaleWeightDown)
    oldTree.SetBranchAddress("StitchWeight", StitchWeight)
    oldTree.SetBranchAddress("ZewkWeight", ZewkWeight)
    oldTree.SetBranchAddress("WewkWeight", WewkWeight)
    oldTree.SetBranchAddress("PUWeight", PUWeight)
    oldTree.SetBranchAddress("TriggerWeight", TriggerWeight)
    oldTree.SetBranchAddress("LeptonWeight", LeptonWeight)
    
    # Analysis variables
    oldTree.SetBranchAddress("isZtoEE", isZtoEE)
    oldTree.SetBranchAddress("isZtoMM", isZtoMM)
    oldTree.SetBranchAddress("isTtoEM", isTtoEM)
    oldTree.SetBranchAddress("isWtoEN", isWtoEN)
    oldTree.SetBranchAddress("isWtoMN", isWtoMN)
    oldTree.SetBranchAddress("isZtoNN", isZtoNN)
    oldTree.SetBranchAddress("isMerged", isMerged)
    oldTree.SetBranchAddress("isResolved", isResolved)
    
    # Objects
    oldTree.SetBranchAddress("nPV", nPV)
    oldTree.SetBranchAddress("nElectrons", nElectrons)
    oldTree.SetBranchAddress("nMuons", nMuons)
    oldTree.SetBranchAddress("nTaus", nTaus)
    oldTree.SetBranchAddress("nPhotons", nPhotons)
    oldTree.SetBranchAddress("nJets", nJets)
    oldTree.SetBranchAddress("nFatJets", nFatJets)
    oldTree.SetBranchAddress("nBTagJets", nBTagJets)
    
    oldTree.SetBranchAddress("MaxJetBTag", MaxJetBTag)
    oldTree.SetBranchAddress("MaxFatJetBTag", MaxFatJetBTag)
    oldTree.SetBranchAddress("MinJetMetDPhi", MinJetMetDPhi)
    # Angular variables
    oldTree.SetBranchAddress("CosThetaStar", CosThetaStar)
    oldTree.SetBranchAddress("CosTheta1", CosTheta1)
    oldTree.SetBranchAddress("CosTheta2", CosTheta2)
    oldTree.SetBranchAddress("Phi", Phi)
    oldTree.SetBranchAddress("Phi1", Phi1)
      
    # Lepton1
    oldTree.SetBranchAddress("Lepton1_isMuon", Lepton1_isMuon)
    oldTree.SetBranchAddress("Lepton1_isElectron", Lepton1_isElectron)
    oldTree.SetBranchAddress("Lepton1_isLoose", Lepton1_isLoose)
    oldTree.SetBranchAddress("Lepton1_isHighPt", Lepton1_isHighPt)
    oldTree.SetBranchAddress("Lepton1_isTrackerHighPt", Lepton1_isTrackerHighPt)
    oldTree.SetBranchAddress("Lepton1_isTight", Lepton1_isTight)
    oldTree.SetBranchAddress("Lepton1_isMuon", Lepton1_isMuon)
    oldTree.SetBranchAddress("Lepton1_pt", Lepton1_pt)
    oldTree.SetBranchAddress("Lepton1_trkIso", Lepton1_trkIso)

    # Lepton2
    oldTree.SetBranchAddress("Lepton2_isMuon", Lepton2_isMuon)
    oldTree.SetBranchAddress("Lepton2_isElectron", Lepton2_isElectron)
    oldTree.SetBranchAddress("Lepton2_isLoose", Lepton2_isLoose)
    oldTree.SetBranchAddress("Lepton2_isHighPt", Lepton2_isHighPt)
    oldTree.SetBranchAddress("Lepton2_isTrackerHighPt", Lepton2_isTrackerHighPt)
    oldTree.SetBranchAddress("Lepton2_isTight", Lepton2_isTight)
    oldTree.SetBranchAddress("Lepton2_isMuon", Lepton2_isMuon)
    oldTree.SetBranchAddress("Lepton2_pt", Lepton2_pt)
    oldTree.SetBranchAddress("Lepton2_trkIso", Lepton2_trkIso)

    # MET        
    oldTree.SetBranchAddress("MEt_pt", MEt_pt)

    # V        
    oldTree.SetBranchAddress("V_pt", V_pt)
    oldTree.SetBranchAddress("V_dPhi", V_dPhi)
    oldTree.SetBranchAddress("V_mass", V_mass)
    oldTree.SetBranchAddress("V_tmass", V_tmass)

    # X        
    oldTree.SetBranchAddress("X_pt", X_pt)
    oldTree.SetBranchAddress("X_dPhi", X_dPhi)
    oldTree.SetBranchAddress("X_mass", X_mass)
    oldTree.SetBranchAddress("X_tmass", X_tmass)

    # FatJet1
    oldTree.SetBranchAddress("FatJet1_isTight", FatJet1_isTight)
    oldTree.SetBranchAddress("FatJet1_pt", FatJet1_pt)
    oldTree.SetBranchAddress("FatJet1_prunedMass", FatJet1_prunedMass)
    oldTree.SetBranchAddress("FatJet1_softdropMass", FatJet1_softdropMass)
    oldTree.SetBranchAddress("FatJet1_softdropPuppiMass", FatJet1_softdropPuppiMass)
    oldTree.SetBranchAddress("FatJet1_prunedMassCorr", FatJet1_prunedMassCorr)
    oldTree.SetBranchAddress("FatJet1_softdropMassCorr", FatJet1_softdropMassCorr)
    oldTree.SetBranchAddress("FatJet1_softdropPuppiMassCorr", FatJet1_softdropPuppiMassCorr)
    oldTree.SetBranchAddress("FatJet1_tau21", FatJet1_tau21)
    oldTree.SetBranchAddress("FatJet1_CSV1", FatJet1_CSV1)
    oldTree.SetBranchAddress("FatJet1_CSV2", FatJet1_CSV2)
    

    oldTree.SetBranchStatus("*",0)

    oldTree.SetBranchStatus("isMC",1)
    oldTree.SetBranchStatus("EventNumber",1)
    oldTree.SetBranchStatus("LumiNumber",1)
    oldTree.SetBranchStatus("RunNumber",1)
    oldTree.SetBranchStatus("EventWeight",1)
    oldTree.SetBranchStatus("FacWeightUp",1)
    oldTree.SetBranchStatus("FacWeightDown",1)
    oldTree.SetBranchStatus("RenWeightUp",1)
    oldTree.SetBranchStatus("RenWeightDown",1)
    oldTree.SetBranchStatus("ScaleWeightUp",1)
    oldTree.SetBranchStatus("ScaleWeightDown",1)
    oldTree.SetBranchStatus("StitchWeight",1)
    oldTree.SetBranchStatus("ZewkWeight",1)
    oldTree.SetBranchStatus("WewkWeight",1)
    oldTree.SetBranchStatus("PUWeight",1)
    oldTree.SetBranchStatus("TriggerWeight",1)
    oldTree.SetBranchStatus("LeptonWeight",1)
    
    # Analysis variables
    oldTree.SetBranchStatus("isZtoEE",1)
    oldTree.SetBranchStatus("isZtoMM",1)
    oldTree.SetBranchStatus("isTtoEM",1)
    oldTree.SetBranchStatus("isWtoEN",1)
    oldTree.SetBranchStatus("isWtoMN",1)
    oldTree.SetBranchStatus("isZtoNN",1)
    oldTree.SetBranchStatus("isMerged",1)
    oldTree.SetBranchStatus("isResolved",1)
    
    # Objects
    oldTree.SetBranchStatus("nPV",1)
    oldTree.SetBranchStatus("nElectrons",1)
    oldTree.SetBranchStatus("nMuons",1)
    oldTree.SetBranchStatus("nTaus",1)
    oldTree.SetBranchStatus("nPhotons",1)
    oldTree.SetBranchStatus("nJets",1)
    oldTree.SetBranchStatus("nFatJets",1)
    oldTree.SetBranchStatus("nBTagJets",1)
    
    oldTree.SetBranchStatus("MaxJetBTag",1)
    oldTree.SetBranchStatus("MaxFatJetBTag",1)
    oldTree.SetBranchStatus("MinJetMetDPhi",1)
    oldTree.SetBranchStatus("Chi2",1)
    # Angular variables
    oldTree.SetBranchStatus("CosThetaStar",1)
    oldTree.SetBranchStatus("CosTheta1",1)
    oldTree.SetBranchStatus("CosTheta2",1)
    oldTree.SetBranchStatus("Phi",1)
    oldTree.SetBranchStatus("Phi1",1)
      
    # Lepton1
    oldTree.SetBranchStatus("Lepton1_isMuon",1)
    oldTree.SetBranchStatus("Lepton1_isElectron",1)
    oldTree.SetBranchStatus("Lepton1_isLoose",1)
    oldTree.SetBranchStatus("Lepton1_isHighPt",1)
    oldTree.SetBranchStatus("Lepton1_isTrackerHighPt",1)
    oldTree.SetBranchStatus("Lepton1_isTight",1)
    oldTree.SetBranchStatus("Lepton1_isMuon",1)
    oldTree.SetBranchStatus("Lepton1_pt",1)
    oldTree.SetBranchStatus("Lepton1_trkIso",1)

    # Lepton2
    oldTree.SetBranchStatus("Lepton2_isMuon",1)
    oldTree.SetBranchStatus("Lepton2_isElectron",1)
    oldTree.SetBranchStatus("Lepton2_isLoose",1)
    oldTree.SetBranchStatus("Lepton2_isHighPt",1)
    oldTree.SetBranchStatus("Lepton2_isTrackerHighPt",1)
    oldTree.SetBranchStatus("Lepton2_isTight",1)
    oldTree.SetBranchStatus("Lepton2_isMuon",1)
    oldTree.SetBranchStatus("Lepton2_pt",1)
    oldTree.SetBranchStatus("Lepton2_trkIso",1)

    # MET        
    oldTree.SetBranchStatus("MEt_pt",1)

    # V        
    oldTree.SetBranchStatus("V_pt",1)
    oldTree.SetBranchStatus("V_dPhi",1)
    oldTree.SetBranchStatus("V_mass",1)
    oldTree.SetBranchStatus("V_tmass",1)

    # X        
    oldTree.SetBranchStatus("X_pt",1)
    oldTree.SetBranchStatus("X_dPhi",1)
    oldTree.SetBranchStatus("X_mass",1)
    oldTree.SetBranchStatus("X_tmass",1)

    # FatJet1
    oldTree.SetBranchStatus("FatJet1_isTight",1)
    oldTree.SetBranchStatus("FatJet1_pt",1)
    oldTree.SetBranchStatus("FatJet1_prunedMass",1)
    oldTree.SetBranchStatus("FatJet1_softdropMass",1)
    oldTree.SetBranchStatus("FatJet1_softdropPuppiMass",1)
    oldTree.SetBranchStatus("FatJet1_prunedMassCorr",1)
    oldTree.SetBranchStatus("FatJet1_softdropMassCorr",1)
    oldTree.SetBranchStatus("FatJet1_softdropPuppiMassCorr",1)
    oldTree.SetBranchStatus("FatJet1_tau21",1)
    oldTree.SetBranchStatus("FatJet1_CSV1",1)
    oldTree.SetBranchStatus("FatJet1_CSV2",1)

    
    newFile = TFile("Pruned/"+name, "RECREATE")
    newFile.mkdir("ntuple/")
    newFile.cd("ntuple/")
    
    newTree = oldTree.CopyTree("V_pt>200 && FatJet1_pt>170")
    newTree.AutoSave()
    
    oldFile.Close()
    newFile.Close()


print args.folder
subfiles = [x for x in os.listdir(args.folder) if os.path.isfile(os.path.join(args.folder, x))]

os.chdir(args.folder)

os.mkdir('Pruned')

for s in subfiles:
    print s
#    skim(s)
    p = multiprocessing.Process(target=skim, args=(s,))
    jobs.append(p)
    p.start()

os.system('cd ..')

