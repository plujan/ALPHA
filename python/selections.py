#! /usr/bin/env python

selection = {
    "triggerMET" : "(isMC?1:(HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v||HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v||HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v||HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v))",
    "triggerIsoEle" : "(isMC?1:HLT_Ele23_WPLoose_Gsf_v)",
    "triggerIsoMuo" : "(isMC?1:(HLT_IsoMu20_v||HLT_IsoTkMu20_v))",
    "triggerEle" : "(isMC?1:(HLT_Ele105_CaloIdVT_GsfTrkIdT_v))",
    "triggerMuo" : "(isMC?1:(HLT_Mu45_eta2p1_v||HLT_Mu50_v))",
    # Leptons
    "singleEle" : "isWtoEN && Lepton1.pt>135 && Lepton1.tightId && V.dPhi<2 && X.dPhi>2 && MEt.pt>80", #  && nTaus==0
    "singleMuo" : "isWtoMN && Lepton1.pt>55 && Lepton1.highptId && Lepton1.trkIso<0.1 && V.dPhi<2 && X.dPhi>2", # && nTaus==0
    "singleIsoLep" : "((isWtoEN && Lepton1.isElectron && Lepton1.looseId && Lepton2.isMuon && Lepton2.highptId) || (isWtoMN && Lepton2.isElectron && Lepton2.looseId && Lepton1.isMuon && Lepton1.highptId))
    "doubleEle" : "isZtoEE && Lepton1.pt>135 && Lepton2.pt>35 && Lepton1.looseId && Lepton2.looseId", # && Lepton1.highptId && Lepton2.highptId && Lepton1.miniIso<0.1 && Lepton2.miniIso<0.1",
    "doubleMuo" : "isZtoMM && (Lepton1.highptId || Lepton2.highptId) && Lepton1.pt>55 && Lepton2.pt>20 && Lepton1.trkIso<0.1 && Lepton2.trkIso<0.1",
    "doubleIsoEle" : "isZtoEE && Lepton1.pt>25 && Lepton2.pt>10 && Lepton1.looseId && Lepton2.looseId",
    "doubleIsoMuo" : "isZtoMM && (Lepton1.highptId || Lepton2.highptId) && Lepton1.pt>25 && Lepton2.pt>10 && Lepton1.trkIso<0.1 && Lepton2.trkIso<0.1",
    "noLeptons" : "MEt.pt>200 && nMuons==0 && nElectrons==0 && nPhotons==0 && FatJet1.tightId && MinJetMetDPhi>0.5 && X.dPhi>2", #nTaus==0 && 
    # V
    "Boost" : "V.pt>200 && FatJet1.pt>200",
    "Zcut" : "V.mass>70 && V.mass<110",
    "Topcut" : "MaxFatJetBTag>0.460",
    "TopVetocut" : "MaxFatJetBTag<0.460",
    "SRcut" : "(FatJet1.prunedMassCorr>105 && FatJet1.prunedMassCorr<135)",
    "SBcut" : "((FatJet1.prunedMassCorr>30 && FatJet1.prunedMassCorr<65) || (FatJet1.prunedMassCorr>135))", # && FatJet1.prunedMassCorr<300
    "VRcut" : "(FatJet1.prunedMassCorr>65 && FatJet1.prunedMassCorr<105)",
    "LSBcut" : "FatJet1.prunedMassCorr>30 && FatJet1.prunedMassCorr<65",
    "HSBcut" : "FatJet1.prunedMassCorr>135 && FatJet1.prunedMassCorr<300",
    "HPcut" : "FatJet1.tau21<0.55",
    "LPcut" : "FatJet1.tau21>0.55",
    "1Btag" : "((FatJet1.CSV1>0.460 && FatJet1.CSV2<0.460) || (FatJet1.CSV1<0.460 && FatJet1.CSV2>0.460))",
    "2Btag" : "(FatJet1.CSV1>0.460 && FatJet1.CSV2>0.460)",
    
        #------------------------------#
    #----------    VZ    ----------#
    #------------------------------#
    # 2 electrons
    "XVZeePre"  : "triggerEle && doubleEle && Boost && Zcut",
    "XZVeeInc"  : "triggerEle && doubleEle && Boost && Zcut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZVeelpSB" : "triggerEle && doubleEle && Boost && Zcut && LPcut && SBcut",
    "XZVeehpSB" : "triggerEle && doubleEle && Boost && Zcut && HPcut && SBcut",
    #
    "XZVeelpSR" : "triggerEle && doubleEle && Boost && Zcut && LPcut && SRcut",
    "XZVeehpSR" : "triggerEle && doubleEle && Boost && Zcut && HPcut && SRcut",
    # 2 muons
    "XZVmmPre"  : "triggerMuo && doubleMuo && Boost && Zcut",
    "XZVmmInc"  : "triggerMuo && doubleMuo && Boost && Zcut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZVmmlpSB" : "triggerMuo && doubleMuo && Boost && Zcut && LPcut && SBcut",
    "XZVmmhpSB" : "triggerMuo && doubleMuo && Boost && Zcut && HPcut && SBcut",
    #
    "XZVmmlpSR" : "triggerMuo && doubleMuo && Boost && Zcut && LPcut && SRcut",
    "XZVmmhpSR" : "triggerMuo && doubleMuo && Boost && Zcut && HPcut && SRcut",
    # 0 leptons
    "XZVnnPre"  : "triggerMET && noLeptons && Boost",
    "XZVnnInc"  : "triggerMET && noLeptons && Boost && TopVetocut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZVnnlpSB" : "triggerMET && noLeptons && Boost && TopVetocut && LPcut && SBcut",
    "XZVnnhpSB" : "triggerMET && noLeptons && Boost && TopVetocut && HPcut && SBcut",
    #
    "XZVnnlpSR" : "triggerMET && noLeptons && Boost && TopVetocut && LPcut && SRcut",
    "XZVnnhpSR" : "triggerMET && noLeptons && Boost && TopVetocut && HPcut && SRcut",
    # 2 lepton, 1 additional btag
    "XZVlnTR"   : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut",
    "XZVlnlpTR" : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut && LPcut && (SBcut || VRcut || SRcut)",
    "XZVlnhpTR" : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut && HPcut && (SBcut || VRcut || SRcut)",
    # 0 lepton, 1 additional btag
    "XZVnnbTR"  : "triggerMET && noLeptons && Topcut && LPcut && (SBcut || VRcut || SRcut)",
    "XZVnnbbTR" : "triggerMET && noLeptons && Topcut && HPcut && (SBcut || VRcut || SRcut)",
    
    #------------------------------#
    #----------    VH    ----------#
    #------------------------------#
    # 2 electrons
    "XZheePre"  : "triggerEle && doubleEle && Boost && Zcut",
    "XZheeInc"  : "triggerEle && doubleEle && Boost && Zcut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZheebSB"  : "triggerEle && doubleEle && Boost && Zcut && 1Btag && SBcut",
    "XZheebbSB" : "triggerEle && doubleEle && Boost && Zcut && 2Btag && SBcut",
    #
    "XZheebSR"  : "triggerEle && doubleEle && Boost && Zcut && 1Btag && SRcut",
    "XZheebbSR" : "triggerEle && doubleEle && Boost && Zcut && 2Btag && SRcut",
    # 2 muons
    "XZhmmPre"  : "triggerMuo && doubleMuo && Boost && Zcut",
    "XZhmmInc"  : "triggerMuo && doubleMuo && Boost && Zcut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZhmmbSB"  : "triggerMuo && doubleMuo && Boost && Zcut && 1Btag && SBcut",
    "XZhmmbbSB" : "triggerMuo && doubleMuo && Boost && Zcut && 2Btag && SBcut",
    #
    "XZhmmbSR"  : "triggerMuo && doubleMuo && Boost && Zcut && 1Btag && SRcut",
    "XZhmmbbSR" : "triggerMuo && doubleMuo && Boost && Zcut && 2Btag && SRcut",
    # 1 electron
    "XWhenPre"  : "triggerEle && singleEle && Boost",
    "XWhenInc"  : "triggerEle && singleEle && Boost && TopVetocut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XWhenbSB"  : "triggerEle && singleEle && Boost && TopVetocut && 1Btag && SBcut",
    "XWhenbbSB" : "triggerEle && singleEle && Boost && TopVetocut && 2Btag && SBcut",
    #
    "XWhenbSR"  : "triggerEle && singleEle && Boost && TopVetocut && 1Btag && SRcut",
    "XWhenbbSR" : "triggerEle && singleEle && Boost && TopVetocut && 2Btag && SRcut",
    # 1 muon
    "XWhmnPre"  : "triggerMuo && singleMuo && Boost",
    "XWhmnInc"  : "triggerMuo && singleMuo && Boost && TopVetocut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XWhmnbSB"  : "triggerMuo && singleMuo && Boost && TopVetocut && 1Btag && SBcut",
    "XWhmnbbSB" : "triggerMuo && singleMuo && Boost && TopVetocut && 2Btag && SBcut",
    #
    "XWhmnbSR"  : "triggerMuo && singleMuo && Boost && TopVetocut && 1Btag && SRcut",
    "XWhmnbbSR" : "triggerMuo && singleMuo && Boost && TopVetocut && 2Btag && SRcut",
    # 0 leptons
    "XZhnnPre"  : "triggerMET && noLeptons && Boost",
    "XZhnnInc"  : "triggerMET && noLeptons && Boost && TopVetocut && (FatJet1.prunedMassCorr<65 || FatJet1.prunedMassCorr>135)",
    #
    "XZhnnbSB"  : "triggerMET && noLeptons && Boost && TopVetocut && 1Btag && SBcut",
    "XZhnnbbSB" : "triggerMET && noLeptons && Boost && TopVetocut && 2Btag && SBcut",
    #
    "XZhnnbSR"  : "triggerMET && noLeptons && Boost && TopVetocut && 1Btag && SRcut",
    "XZhnnbbSR" : "triggerMET && noLeptons && Boost && TopVetocut && 2Btag && SRcut",
    # 1 lepton, 1 additional btag
    "XWhlnTR"   : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut",
    "XWhlnbTR"  : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut && 1Btag && (SBcut || VRcut || SRcut)",
    "XWhlnbbTR" : "((triggerEle && singleEle) || (triggerMuo && singleMuo)) && Boost && Topcut && 2Btag && (SBcut || VRcut || SRcut)",
    # 0 lepton, 1 additional btag
    "XZhnnbTR"  : "triggerMET && noLeptons && Boost && Topcut && 1Btag && (SBcut || VRcut || SRcut)",
    "XZhnnbbTR" : "triggerMET && noLeptons && Boost && Topcut && 2Btag && (SBcut || VRcut || SRcut)",
    
    #------------------------------#
    #----------    AZh   ----------#
    #------------------------------#
    "AZheeInc"    : "triggerIsoEle && doubleIsoEle",
    "AZhmmInc"    : "triggerIsoMuo && doubleIsoMuo",
    "AZheejjSB" : "triggerIsoEle && doubleIsoEle && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && Jet1.CSV<0.460 && Jet2.CSV<0.460 && MEt.pt<60",
    "AZhmmjjSB" : "triggerIsoMuo && doubleIsoMuo && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && Jet1.CSV<0.460 && Jet2.CSV<0.460 && MEt.pt<60",
    "AZheebjSB" : "triggerIsoEle && doubleIsoEle && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && ((Jet1.CSV>0.935 && Jet2.CSV<0.460) || (Jet1.CSV<0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "AZhmmbjSB" : "triggerIsoMuo && doubleIsoMuo && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && ((Jet1.CSV>0.935 && Jet2.CSV<0.460) || (Jet1.CSV<0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "AZheebbSB" : "triggerIsoEle && doubleIsoEle && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && ((Jet1.CSV>0.935 && Jet2.CSV>0.460) || (Jet1.CSV>0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "AZhmmbbSB" : "triggerIsoMuo && doubleIsoMuo && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass<90 || H.mass>140) && ((Jet1.CSV>0.935 && Jet2.CSV>0.460) || (Jet1.CSV>0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "TopSB" : "triggerIsoMuo && doubleIsoMuo && nJets>=2 && (V.mass<70 || V.mass>110) && MEt.pt>60", #Jet1.CSV>0.935 && Jet2.CSV>0.935 && 
    "AZheebbSR" : "triggerIsoEle && doubleIsoEle && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass>90 && H.mass<140) && ((Jet1.CSV>0.935 && Jet2.CSV>0.460) || (Jet1.CSV>0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "AZhmmbbSR" : "triggerIsoMuo && doubleIsoMuo && nJets>=2 && V.mass>70 && V.mass<110 && (H.mass>90 && H.mass<140) && ((Jet1.CSV>0.935 && Jet2.CSV>0.460) || (Jet1.CSV>0.460 && Jet2.CSV>0.935)) && MEt.pt<60",
    "AZhllbbSR": "(AZhmmbbSR) || (AZheebbSR)",
}

