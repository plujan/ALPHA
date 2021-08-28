
#! /usr/bin/env python

import os, sys, getopt, multiprocessing
import copy, math
from array import array
from ROOT import gROOT, gSystem, gStyle, gRandom
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph, TGaxis
from ROOT import TStyle, TCanvas, TPad, TLegend, TLatex, TText

# Import PDF library and PDF diagonalizer
gSystem.Load('PDFs/HWWLVJRooPdfs_cxx.so')
gSystem.Load('PDFs/PdfDiagonalizer_cc.so')
# if modified or recompiled, copy .cxx and .h to $COMBINE_BASE/src/HiggsAnalysis/CombinedLimit/src/

from ROOT import RooFit, RooRealVar, RooDataHist, RooDataSet, RooAbsData, RooAbsReal, RooAbsPdf, RooPlot, RooBinning, RooCategory, RooSimultaneous, RooArgList, RooArgSet, RooWorkspace, RooMsgService
from ROOT import RooFormulaVar, RooGenericPdf, RooGaussian, RooExponential, RooPolynomial, RooChebychev, RooBreitWigner, RooCBShape, RooExtendPdf, RooAddPdf, RooProdPdf, RooNumConvPdf, RooFFTConvPdf
from ROOT import PdfDiagonalizer, RooAlphaExp, RooErfExpPdf, Roo2ExpPdf, RooAlpha42ExpPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooErfExpTailPdf, RooAlpha4ExpTailPdf, RooExpTail2Pdf, RooAlpha4ExpTail2Pdf, RooErfExpNPdf, RooExpN2Pdf, RooAlpha4ErfExpNPdf, RooAlpha4ExpN2Pdf, RooAlpha4ErfExpTailPdf, RooAlpha4ErfExpNPdf, RooAlpha, RooDoubleCrystalBall, RooErfPow2Pdf

from Analysis.ALPHA.LdrawUtils import *
from Analysis.ALPHA.variables import *
from Analysis.ALPHA.selectionsForAlpha import *
from Analysis.ALPHA.samples import samples
from Lalpha_settings import *

#from parameters import param

import optparse
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-a', '--all', action='store_true', default=False, dest='all')
parser.add_option('-b', '--batch', action='store_true', default=False, dest='batch')
parser.add_option('-c', '--channel', action='store', type='string', dest='channel', default='')
parser.add_option('-d', '--different', action='store_true', default=False, dest='different')
parser.add_option('-e', '--extrapolate', action='store_true', default=False, dest='extrapolate')
parser.add_option('-s', '--scan', action='store_true', default=False, dest='scan')
parser.add_option('-n', '--name', action='store', type='string', default='', dest='signalname')
parser.add_option('-v', '--verbose', action='store_true', default=False, dest='verbose')
parser.add_option("-U", "--unblind", action="store_true", default=False, dest="unblind")
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)

########## SETTINGS ##########

# Silent RooFit
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)
#RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)

#gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetPadTopMargin(0.06)
gStyle.SetPadRightMargin(0.05)

#ALTERNATIVE = options.different
#EXTRAPOLATE = options.extrapolate
#SCAN        = options.scan
#SIGNALNAME  = options.signalname
##NTUPLEDIR   = '/lustre/cmswork/lbenato/VZAnalyses/CMSSW_8_0_12/src/Analysis/ALPHA/ProdDummy_v5/Skim/'#'/lustre/cmsdata/ALPHA/v9/Skim/'
##NTUPLEDIR   = '/lustre/cmswork/lbenato/VZAnalyses/CMSSW_8_0_12/src/Analysis/ALPHA/Prod_v6/Skim/'#'/lustre/cmsdata/ALPHA/v9/Skim/'
#NTUPLEDIR   = '/lustre/cmswork/lbenato/VZAnalyses/CMSSW_8_0_25/src/Analysis/ALPHA/v1/Skim/'
##LUMI        = 12900 # in pb-1
#LUMI        = 36814 # in pb-1
#PLOTDIR     = 'v1/plotsAlphaNorm' if not EXTRAPOLATE else 'v1/plotsAlphaExt'
#RATIO       = 4
#SHOWERR     = True
#BLIND       = False if (EXTRAPOLATE or options.unblind) else True
#VERBOSE     = options.verbose
#LISA        = True
#CHAN        = options.channel


#channelList = ['XVZnnlp', 'XVZnnhp']#['XVZmmlp', 'XVZmmhp', 'XVZeelp', 'XVZeehp']

#if (SIGNALNAME == 'XWZ') or (SIGNALNAME == 'XZZ'):
#    signName  = SIGNALNAME
#else :
#    signName  = 'XVZ' # 'XZZ'

#genPoints = [600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000 ]
##massPoints = [x for x in range(600, 4500+1, 100)] #if not HVTMODEL else genPoints
#massPoints = [x for x in range(600, 3000+1, 200)] #if not HVTMODEL else genPoints
##massPoints = [600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000]


##gROOT.ProcessLine('.L %s/src/Analysis/ALPHA/plugins/Objects.h+' % os.environ['CMSSW_BASE'])
##from ROOT import LeptonType, JetType, FatJetType, MEtType, MEtFullType, CandidateType, LorentzType

#LOWMIN = 30.
#LOWMAX = 50. if EXTRAPOLATE else 65. 

#SIGMIN = 50. if EXTRAPOLATE else 65.
#SIGMAX = 65. if EXTRAPOLATE else 105.

#HMIN = 105.
#HMAX = 135.

#HIGMIN = 135.
#HIGMAX = 300.

#VBINMIN= 30.
#VBINMAX= 300.
#VBINS  = 54

#XBINMIN= 650.#550.
#XBINMAX= 3150.
#XBINS_PLOT  = 25
#XBINS  = 1000

## Binning
#binsJmass = RooBinning(VBINS, VBINMIN, VBINMAX)
#binsXmass = RooBinning(XBINS_PLOT, XBINMIN, XBINMAX)
#binsSignal = RooBinning(XBINS_PLOT*5, XBINMIN, XBINMAX)
#binsSignal.addUniform(XBINS_PLOT*5, XBINMIN, XBINMAX)

## WSF[higP]76X with v2 JEC
##WSF = [1.154, 0.960, ]
##WSFErr = [0.305, 0.090, ]
## WSF[higP]80X with ICHEP JEC
#WSF = [0.88, 1.03, ]
#WSFErr = [0.24, 0.078,]

########## ######## ##########

def alpha_norm(channel):

    nElec = channel.count('e')
    nMuon = channel.count('m')
    nLept = nElec + nMuon
    lowP  = channel.count('lp')
    higP  = channel.count('hp')
    allSB = channel.count('allSB')#L
    nBtag = 0
    bkgSF = 1.
            

    if lowP:
        TAUMIN = 0.40
        TAUMAX = 0.75
    else:
        TAUMIN = 0.00
        TAUMAX = 0.40

    ## Channel-dependent settings
    ## Background function. Semi-working options are: EXP, EXP2, EXPN, EXPTAIL
    if nMuon == 0 and nElec == 0:
        treeName  = 'alpha'
        colorVjet = samples['DYJetsToNuNu_PtZ']['linecolor']
        cut       = 'isZtoNN'
        trigger   = 'HLT_PFMET'
        massVar   = 'X_tmass'
        jetVar    = 'FatJet1_softdropPuppiMassCorr'
        tauVar    = 'FatJet1_puppiTau21'
        binFact   = 1
        # --- mX ---
        fitFunc         = 'EXPTAIL2'#provo la 2
        fitAltFunc      = 'EXPN'
        fitFuncBkg      = 'EXPTAIL' if lowP else 'EXPTAIL'
        # --- mJ ---
        fitFuncVjet     = 'TEST' if lowP else 'ERFEXP'#da fit migliori all'8 novembre: 'EXPGAUS' if lowP else 'ERFEXP'#L#funziona prima del 2 novembre: 'ERFEXP' if lowP else 'ERFEXP'#: EXPGAUS best, ERFEXPGAUS no convergence
        fitAltFuncVjet  = 'ERFEXP' if lowP else 'EXPGAUS'#da fit migliori all'8 novembre: 'ERFEXP' if lowP else 'EXPGAUS'#funziona prima del 2 novembre: 'EXPGAUS' if lowP else 'EXPGAUS'#'POL3' forse per la hp usavo altro...tenta exptail,expgaus 
        fitFuncVV       = 'EXPGAUS' # 'EXPGAUS2' WHEN THE VH SAMPLES WILL BE INCLUDED
        fitFuncTop      = 'ERFEXPGAUS2'# if lowP else 'GAUS3'#L# erfexp-gaus2 converge,erfexpgaus va meglio 'ERFEXPGAUS' if lowP else 'ERFEXPGAUS2'

    '''
    if nBtag == 2:
        btagCut = selection['1Btag']
    elif nBtag == 1:
        btagCut = selection['2Btag']
    else :
        btagCut = '0==0'
    '''

    print 'Opening ntuples in ', NTUPLEDIR
    print '--- Channel', channel, '---'
    print '  number of electrons:', nElec, ' muons:', nMuon, ' b-tags:', nBtag
    if(lowP): print '  low purity :', lowP
    if(higP):print '  high purity:', higP
    if(allSB):print '  all SB:', allSB
    print '  read tree  :', treeName, 'and cut:', cut
    print '  signal     :', signName
    if ALTERNATIVE: print '  using ALTERNATIVE fit functions'
    if EXTRAPOLATE: print '  using VALIDATION mode'
    #print '  fit func mX primary bkg: ', fitFunc
    #print '  fit func mX secondary bkg: ', fitFuncBkg
    print '  fit func mJ Vjet: ',fitFuncVjet
    print '  fit alt func mJ Vjet: ',fitAltFuncVjet
    print '-'*11*2

    #*******************************************************#
    #                                                       #
    #              Variables and selections                 #
    #                                                       #
    #*******************************************************#
            
    # Define all the variables from the trees that will be used in the cuts and fits
    # this steps actually perform a 'projection' of the entire tree on the variables in thei ranges, so be careful once setting the limits
    X_mass                  = RooRealVar( massVar,      'm_{VZ}^{T}' if "XVZnn" in channel else 'm_{VZ}' ,          XBINMIN, XBINMAX, 'GeV')
    J_mass                  = RooRealVar( jetVar,       'jet mass',         VBINMIN, VBINMAX, 'GeV')
    J_tau21                 = RooRealVar( tauVar,       'jet #tau_{21}',    0.,      1.            )
    isZtoNN                 = RooRealVar( 'isZtoNN',    '',                 0,       2             )
    weight                  = RooRealVar( 'EventWeight','',                -1.e9,    1.e9          )
        
    # Define the RooArgSet which will include all the variables defined before
    # there is a maximum of 9 variables in the declaration, so the others need to be added with 'add'
    variables = RooArgSet(X_mass, J_mass, J_tau21, isZtoNN, weight, 'variables')
        
    # set reasonable ranges for J_mass and X_mass
    # these are used in the fit in order to avoid ROOFIT to look in regions very far away from where we are fitting 
    # (honestly, it is not clear to me why it is necessary, but without them the fit often explodes)
    X_mass.setRange('X_reasonable_range', XBINMIN, XBINMAX)
    J_mass.setRange('V_reasonable_range', LOWMIN, HIGMAX)
    J_mass.setRange('V_extended_reasonable_range', 0, HIGMAX)
    
    # Set RooArgSets once for all, see https://root.cern.ch/phpBB3/viewtopic.php?t=11758
    jetMassArg = RooArgSet(J_mass)
    # Define the ranges in fatJetMass - these will be used to define SB and SR
    J_mass.setRange('LSBrange', LOWMIN, LOWMAX)
    J_mass.setRange('HSBrange', HIGMIN, HIGMAX)
    J_mass.setRange('SRrange',  SIGMIN, SIGMAX)
    J_mass.setRange('HRrange',  SIGMAX, HIGMIN)
    
    # Set binning for plots
    J_mass.setBins(VBINS)
    X_mass.setBins(binFact*XBINS)
    
    # Define the selection for the various categories (base + SR / LSBcut / HSBcut )
    baseCut = '{3} && {0}>{1} && {0}<{2}'.format(J_tau21.GetName(), TAUMIN, TAUMAX, cut)
    
    # Cuts
    SRcut  = baseCut + ' && {0}>{1} && {0}<{2}'.format(J_mass.GetName(), SIGMIN, SIGMAX)
    LSBcut = baseCut + ' && {0}>{1} && {0}<{2}'.format(J_mass.GetName(), LOWMIN, LOWMAX)
    HSBcut = baseCut + ' && {0}>{1} && {0}<{2}'.format(J_mass.GetName(), HIGMIN, HIGMAX)
    SBcut  = baseCut + ' && (({0}>{1} && {0}<{2}) || ({0}>{3} && {0}<{4}))'.format(J_mass.GetName(), LOWMIN, LOWMAX, HIGMIN, HIGMAX)
    HRcut  = baseCut + ' && {0}>{1} && {0}<{2}'.format(J_mass.GetName(), SIGMAX, HIGMIN)
            

    #*******************************************************#
    #                                                       #
    #                      Input files                      #
    #                                                       #
    #*******************************************************#
    
    # Import the files using TChains (separately for the bkg 'classes' that we want to describe: here DY and VV+ST+TT)
    treeData = TChain(treeName)
    treeMC   = TChain(treeName)
    treeVjet = TChain(treeName)
    treeVV   = TChain(treeName)
    treeTop  = TChain(treeName)
    
    # Read data
    pd = getPrimaryDataset(trigger)
    if len(pd)==0: raw_input('Warning: Primary Dataset not recognized, continue?')
    for i, s in enumerate(pd): 
        treeData.Add(NTUPLEDIR + s + '.root')
    
    # Read V+jets backgrounds
    #for i, s in enumerate(['WJetsToLNu_HT', 'DYJetsToNuNu_PtZ', 'DYJetsToLL_HT']):
    for i, s in enumerate(['DYJetsToNuNu_PtZ', 'WJetsToLNu_Pt']):
        for j, ss in enumerate(samples[s]['files']): 
            treeVjet.Add(NTUPLEDIR + ss + '.root')
    
    # Read VV backgrounds
    for i, s in enumerate(['VV']):
        for j, ss in enumerate(samples[s]['files']): 
            treeVV.Add(NTUPLEDIR + ss + '.root')
    
    # Read Top backgrounds
    for i, s in enumerate(['TTbar','ST']):
        for j, ss in enumerate(samples[s]['files']): 
            treeTop.Add(NTUPLEDIR + ss + '.root')
        
    # Sum all background MC
    treeMC.Add(treeVjet)
    treeMC.Add(treeVV)
    treeMC.Add(treeTop)

    # create a dataset to host data in sideband (using this dataset we are automatically blind in the SR!)    
    setDataSB = RooDataSet('setDataSB', 'setDataSB', variables, RooFit.Cut(SBcut), RooFit.WeightVar(weight), RooFit.Import(treeData))
    setDataLSB = RooDataSet('setDataLSB', 'setDataLSB', variables, RooFit.Import(setDataSB), RooFit.Cut(LSBcut), RooFit.WeightVar(weight))
    setDataHSB = RooDataSet('setDataHSB', 'setDataHSB', variables, RooFit.Import(setDataSB), RooFit.Cut(HSBcut), RooFit.WeightVar(weight))
    
    # Observed data (WARNING, BLIND!)
    setDataSR = RooDataSet('setDataSR', 'setDataSR', variables, RooFit.Cut(SRcut), RooFit.WeightVar(weight), RooFit.Import(treeData))
    setDataHR = RooDataSet('setDataHR', 'setDataHR', variables, RooFit.Cut(HRcut), RooFit.WeightVar(weight), RooFit.Import(treeData)) # Observed in the VH mass, just for plotting purposes

    setDataSRSB = RooDataSet('setDataSRSB', 'setDataSRSB', variables, RooFit.Cut('('+SRcut+') || ('+SBcut+')'), RooFit.WeightVar(weight), RooFit.Import(treeData))
    
    # same for the bkg datasets from MC, where we just apply the base selections (not blind)
    setVjet = RooDataSet('setVjet', 'setVjet', variables, RooFit.Cut(baseCut), RooFit.WeightVar(weight), RooFit.Import(treeVjet))
    setVjetSB = RooDataSet('setVjetSB', 'setVjetSB', variables, RooFit.Import(setVjet), RooFit.Cut(SBcut), RooFit.WeightVar(weight))
    setVjetSR = RooDataSet('setVjetSR', 'setVjetSR', variables, RooFit.Import(setVjet), RooFit.Cut(SRcut), RooFit.WeightVar(weight))
    setVV = RooDataSet('setVV', 'setVV', variables, RooFit.Cut(baseCut), RooFit.WeightVar(weight), RooFit.Import(treeVV))
    setVVSB = RooDataSet('setVVSB', 'setVVSB', variables, RooFit.Import(setVV), RooFit.Cut(SBcut), RooFit.WeightVar(weight))
    setVVSR = RooDataSet('setVVSR', 'setVVSR', variables, RooFit.Import(setVV), RooFit.Cut(SRcut), RooFit.WeightVar(weight))
    setTop = RooDataSet('setTop', 'setTop', variables, RooFit.Cut(baseCut), RooFit.WeightVar(weight), RooFit.Import(treeTop))
    setTopSB = RooDataSet('setTopSB', 'setTopSB', variables, RooFit.Import(setTop), RooFit.Cut(SBcut), RooFit.WeightVar(weight))
    setTopSR = RooDataSet('setTopSR', 'setTopSR', variables, RooFit.Import(setTop), RooFit.Cut(SRcut), RooFit.WeightVar(weight))
    
    print '  Data events SB: %.2f' % setDataSB.sumEntries()
    print '  V+jets entries: %.2f' % setVjet.sumEntries()
    print '  VV, VH entries: %.2f' % setVV.sumEntries()
    print '  Top,ST entries: %.2f' % setTop.sumEntries()
    print 'IN SIGNAL REGION:::::::::::::::::::::::::'
    print '  V+jets entries: %.2f' % setVjetSR.sumEntries()
    print '  VV, VH entries: %.2f' % setVVSR.sumEntries()
    print '  Top,ST entries: %.2f' % setTopSR.sumEntries()
    
    nVV   = RooRealVar('nVV',  'VV normalization',   setVV.sumEntries(),   0., 1.e7)
    nTop  = RooRealVar('nTop', 'Top normalization',  setTop.sumEntries(),  0., 1.e7)
    nMain = RooRealVar('nMain','Vjet normalization', setVjet.sumEntries(), 0., 1.e7)
    nMain2 = RooRealVar('nMain2','Vjet2 normalization', setVjet.sumEntries(), 0., 1.e7)
    nMainEasy = RooRealVar('nMainEasy','Vjet normalization', setVjet.sumEntries(), 0., 1.e7)
    
    # Apply Top SF #Dopo! #L
    nTop.setVal(nTop.getVal()*WSF[higP])
    nTop.setError(nTop.getVal()*WSFErr[higP])
    
    # Apply VV SF #Dopo! #L
    nVV.setVal(nVV.getVal()*WSF[higP])
    nVV.setError(nVV.getVal()*WSFErr[higP])

    # Define entries
    entryVjet = RooRealVar('entryVjets',  'V+jets normalization', setVjet.sumEntries(), 0., 1.e7)
    entryVV = RooRealVar('entryVV',  'VV normalization', setVV.sumEntries(), 0., 1.e7)
    entryTop = RooRealVar('entryTop',  'Top normalization', setTop.sumEntries(), 0., 1.e7)
    entryVjetSB = RooRealVar('entryVjetsSB',  'V+jets normalization', setVjet.sumEntries(SBcut), 0., 1.e7)
    entryVVSB = RooRealVar('entryVVSB',  'VV normalization', setVV.sumEntries(SBcut), 0., 1.e7)
    entryTopSB = RooRealVar('entryTopSB',  'Top normalization', setTop.sumEntries(SBcut), 0., 1.e7)
    
    entrySB = RooRealVar('entrySB',  'Data SB normalization', setDataSB.sumEntries(SBcut), 0., 1.e7)
    entrySB.setError(math.sqrt(entrySB.getVal()))
    
    entryLSB = RooRealVar('entryLSB',  'Data LSB normalization', setDataSB.sumEntries(LSBcut), 0., 1.e7)
    entryLSB.setError(math.sqrt(entryLSB.getVal()))

    entryHSB = RooRealVar('entryHSB',  'Data HSB normalization', setDataSB.sumEntries(HSBcut), 0., 1.e7)
    entryHSB.setError(math.sqrt(entryHSB.getVal()))

    #Normalization

    ###################################################################################
    #        _   _                                                                    #
    #       | \ | |                          | (_)         | | (_)                    #
    #       |  \| | ___  _ __ _ __ ___   __ _| |_ ___  __ _| |_ _  ___  _ __          #
    #       | . ` |/ _ \| '__| '_ ` _ \ / _` | | / __|/ _` | __| |/ _ \| '_ \         #
    #       | |\  | (_) | |  | | | | | | (_| | | \__ \ (_| | |_| | (_) | | | |        #
    #       |_| \_|\___/|_|  |_| |_| |_|\__,_|_|_|___/\__,_|\__|_|\___/|_| |_|        #
    #                                                                                 #
    ###################################################################################
    # fancy ASCII art thanks to, I guess, Jose
    
    # start by creating the fit models to get the normalization: 
    # * MAIN and SECONDARY bkg are taken from MC by fitting the whole J_mass range
    # * The two PDFs are added together using the relative normalizations of the two bkg from MC
    # * DATA is then fit in the sidebands only using the combined bkg PDF
    # * The results of the fit are then estrapolated in the SR and the integral is evaluated.
    # * This defines the bkg normalization in the SR
    
    #*******************************************************#
    #                                                       #
    #                 V+jets normalization                  #
    #                                                       #
    #*******************************************************#
    
    # Variables for V+jets
    a0Vjet = RooRealVar('a0Vjet', 'width of the erf', -0.2, -5, 5)
    a1Vjet = RooRealVar('a1Vjet', 'width of the erf', 0.2,  -2, 2)
    a2Vjet = RooRealVar('a2Vjet', 'width of the erf', -0.2, -2, 2)
    a3Vjet = RooRealVar('a3Vjet', 'width of the erf', 0.2, -5, 5)
    a4Vjet = RooRealVar('a4Vjet', 'width of the erf', -0.2, -5, 5)
    a5Vjet = RooRealVar('a5Vjet', 'width of the erf', -0.2, -5, 5)
    a6Vjet = RooRealVar('a6Vjet', 'width of the erf', -0.2, -5, 5)
    a7Vjet = RooRealVar('a7Vjet', 'width of the erf', -0.2, -5, 5)
    if channel=='XVZnnlp' and (('POL5'in fitFuncVjet) or ('POL4'in fitFuncVjet) or ('POL3'in fitFuncVjet)):#qua
        print "pol"
##try refitting twice
#        a0Vjet.setVal(-1.38130078)
#        a0Vjet.setConstant(True)
#        #a0Vjet.setMin(-1000.)
#        #a0Vjet.setMax(+1000.)
#        a1Vjet.setVal(0.)
#        #a1Vjet.setMin(-1000.)
#        #a1Vjet.setMax(+1000.)
#        a2Vjet.setVal(0.)
#        #a2Vjet.setMin(-1000.)
#        #a2Vjet.setMax(+1000.)
#        a3Vjet.setVal(0.)
#        a4Vjet.setVal(0.)
#converge pol5 out of the box:
##        a0Vjet.setVal(-1.)#.5870e+00)
##        a0Vjet.setConstant(True)
##        a0Vjet.setMin(-1.4)
##        a0Vjet.setMax(1)
#        a1Vjet.setVal(2.e-01)
##        a1Vjet.setConstant(True)
##        a1Vjet.setMin(0.1)
#        a1Vjet.setMax(1)
#        a2Vjet.setVal(-2.0759e-02)
##        a2Vjet.setConstant(True)
##        a2Vjet.setMin(-0.01)
#        a2Vjet.setMax(-0.1)
#        a3Vjet.setVal(1.4761e-02)
##        a3Vjet.setConstant(True)
#        a3Vjet.setMin(0.02)#(0.001)#(-0.001)#0.02 mi fa finalmente cambiare la concavita!
#        a3Vjet.setMax(0.1)
#        a4Vjet.setVal(-0.2)
##        a4Vjet.setConstant(True)
##        a4Vjet.setMin(-0.1)
#        a4Vjet.setMax(-0.001)#(-0.001)
#        a4Vjet.setVal(-0.20)

        #let's try something different!
        #a4Vjet.setVal(-0.12)
        #a4Vjet.setConstant(True)

    if channel=='XVZnnlp':
        constVjet   = RooRealVar('constVjet',   'slope of the exp',      -0.03, -1., 0.)#default
        const2Vjet   = RooRealVar('const2Vjet',   'slope of the 2nd exp',      -0.03, -1., 0.)#default
#        constVjet   = RooRealVar('constVjet',   'slope of the exp',      -0.040, -1., 0.)
        offsetVjet  = RooRealVar('offsetVjet',  'offset of the erf',     10., 0., 1000.)
#        offsetVjet  = RooRealVar('offsetVjet',  'offset of the erf',     70., 25., 200.)#default
#        widthVjet   = RooRealVar('widthVjet',   'width of the erf',      5., 0., 50.)#50 e 400
        widthVjet   = RooRealVar('widthVjet',   'width of the erf',      10., 0.01, 110.)#50 e 400#16 feb
#        widthVjet   = RooRealVar('widthVjet',   'width of the erf',      10., 0.1, 400.)#default
        constVjet2   = RooRealVar('constVjet2',   'slope of the exp',      -0.03, -1., 0.)#L#
        offsetVjet2  = RooRealVar('offsetVjet',  'offset of the erf',     10., 0., 1000.)#L#
        widthVjet2   = RooRealVar('widthVjet',   'width of the erf',      10., 0.01, 100.)#L#

    elif channel=='XVZnnhp':
        constVjet   = RooRealVar('constVjet',   'slope of the exp',      -0.0202, -1., 0.)
        const2Vjet   = RooRealVar('const2Vjet',   'slope of the 2nd exp',      -0.03, -1., 0.)
        offsetVjet  = RooRealVar('offsetVjet',  'offset of the erf',     84.8, 20., 200.)
        widthVjet   = RooRealVar('widthVjet',   'width of the erf',      43.7, 0.1, 400.)        
        constVjet2   = RooRealVar('constVjet2',   'slope of the exp',      -0.02, -1., 0.)#L#
        offsetVjet2  = RooRealVar('offsetVjet',  'offset of the erf',     84., 20., 200.)#L#
        widthVjet2   = RooRealVar('widthVjet',   'width of the erf',      45., 0.1, 400.)#L#
    

    g1meanVjet  = RooRealVar("g1meanVjet", "mean of gaus1", 150., 80., 220.)
    g1sigmaVjet = RooRealVar("g1sigmaVjet", "sigma of gaus1", 60, 10., 100.)
    gaus1Vjet   = RooGaussian('gaus1Vjet', 'gaus1Vjet2', J_mass, g1meanVjet, g1sigmaVjet)
#    g2meanVjet  = RooRealVar("g2meanVjet", "mean of gaus2", 54., 30., 120.)#default
#    g2sigmaVjet = RooRealVar("g2sigmaVjet", "sigma of gaus2", 40, 10., 100.)#default
    g2meanVjet  = RooRealVar("g2meanVjet", "mean of gaus2", 75., 70., 80.)#8feb
    g2meanVjet.setVal(75.)
    #g2meanVjet.setConstant(True)#8feb
    g2sigmaVjet = RooRealVar("g2sigmaVjet", "sigma of gaus2", 20, 0., 1000.)#8feb
    #g2sigmaVjet.setVal(10.)#9feb
    g2sigmaVjet.setConstant(True)
    g2meanVjetLess  = RooRealVar("g2meanVjetLess", "mean of gaus2", 42., 35., 55.)#8feb
    #g2meanVjetLess.setConstant(True)#8feb
    g2sigmaVjetLess = RooRealVar("g2sigmaVjetLess", "sigma of gaus2", 30, 0., 1000.)#8feb
    #g2sigmaVjetLess.setVal(5.)#9feb
    g2sigmaVjetLess.setConstant(True)
    gaus2Vjet   = RooGaussian('gaus2Vjet', 'gaus2Vjet2', J_mass, g2meanVjet, g2sigmaVjet)
    g12fracVjet = RooRealVar("g12fracVjet", "frac gaus2 gaus1", 0.3, 0., 1.)
    g2meanVjet2  = RooRealVar("g2meanVjet2", "mean of gaus2", 54., 30., 120.)#default
    g2sigmaVjet2 = RooRealVar("g2sigmaVjet2", "sigma of gaus2", 40, 10., 100.)#default

  
    expoVjet = RooExponential('expoVjet', 'expoVjet', J_mass, constVjet)
    gausVjet = RooGaussian('gausVjet', 'gausVjet', J_mass, g2meanVjet, g2sigmaVjet)
    gausVjetLess = RooGaussian('gausVjetLess', 'gausVjetLess', J_mass, g2meanVjetLess, g2sigmaVjetLess)#8feb
    fracVjet = RooRealVar('fracVjet',   'fraction of gaussian wrt exp', 0.3, 0., 1.)#default,se la cambio non cambia nulla
    #fracVjet.setVal(0.0)#9feb
    #fracVjet.setConstant(True)#9feb
    fracVjetLess = RooRealVar('fracVjetLess',   'fraction of gaussian wrt exp', 0.3, 0., 1.)#8feb
    expoVjet2 = RooExponential('expoVjet2', 'expoVjet2', J_mass, constVjet2)
    gausVjet2 = RooGaussian('gausVjet2', 'gausVjet2', J_mass, g2meanVjet2, g2sigmaVjet2)
    fracVjet2 = RooRealVar('fracVjet2',   'fraction of gaussian wrt exp', 0.3, 0., 1.)#default,se la cambio non cambia nulla

    Vjet_var1   = RooRealVar('Vjet_var1', 'Vjet_var1', 17, 5, 30)
    Vjet_var2   = RooRealVar('Vjet_var2', 'Vjet_var2', 0.034, 0.01, 0.05)
    Vjet_var3   = RooRealVar('Vjet_var3', 'Vjet_var3', 54., 40., 65.)
    Vjet_var4   = RooRealVar('Vjet_var4', 'Vjet_var4', 36., 25., 45.)

    # Define V+jets model
    if fitFuncVjet == 'ERFEXP': VjetMass = RooErfExpPdf('VjetMass', fitFuncVjet, J_mass, constVjet, offsetVjet, widthVjet)
    elif fitFuncVjet == 'EXP': VjetMass  = RooExponential('VjetMass',   fitFuncVjet, J_mass, constVjet)
    elif fitFuncVjet == 'EXPPOL':
        #VjetMass  = RooGenericPdf('VjetMass', fitFuncVjet, "exp(@0 * @1)*(1 + @0*@2 + @0*@0*@3 + @0*@0*@0*@4 + @0*@0*@0*@0*@5 + @0*@0*@0*@0*@0*@6)", RooArgList(J_mass, constVjet, a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet))
        VjetMass  = RooGenericPdf('VjetMass', fitFuncVjet, "exp(@0 * @1)*(1 + @0*@2 + @0*@0*@3 + @0*@0*@0*@4)", RooArgList(J_mass, constVjet, a0Vjet, a1Vjet, a2Vjet))
    elif fitFuncVjet == 'EXPGAUS': VjetMass  = RooAddPdf('VjetMass',   fitFuncVjet, RooArgList(gausVjet, expoVjet), RooArgList(fracVjet))
    elif fitFuncVjet == 'EXPGAUS2': VjetMass  = RooAddPdf('VjetMass',   fitFuncVjet, RooArgList(gausVjet, gausVjetLess, expoVjet), RooArgList(fracVjet, fracVjetLess))#8feb
    elif fitFuncVjet == 'POL3': VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
    elif fitFuncVjet == 'POL4':
        VjetMass3 = RooChebychev('VjetMass3', 'POL3', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
        VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
    elif fitFuncVjet == 'POL5':
        VjetMass3 = RooChebychev('VjetMass3', 'POL3', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
        VjetMass4 = RooChebychev('VjetMass4', 'POL4', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
        #this is the last fit:
        VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet))
    elif fitFuncVjet == 'POLT':
        VjetMass3 = RooChebychev('VjetMass3', 'POL3', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
        VjetMass4 = RooChebychev('VjetMass4', 'POL4', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
        #this is the last fit:
        VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet))
        VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet, a5Vjet))
    elif fitFuncVjet == 'POL6':
        VjetMass = RooChebychev('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet, a5Vjet))
    elif fitFuncVjet == 'POLY2':
        a0Vjet.setVal(0)
        a0Vjet.setMin(-1000)
        a0Vjet.setMax(1000)
        a1Vjet.setVal(0)
        a1Vjet.setMin(-1000)
        a1Vjet.setMax(1000)
        a2Vjet.setVal(0)
        a2Vjet.setMin(-1000)
        a2Vjet.setMax(1000)
        VjetMass = RooPolynomial('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet))
    elif fitFuncVjet == 'POLY3':
        a0Vjet.setVal(0)
        a1Vjet.setVal(0)
        a2Vjet.setVal(0)
        VjetMass = RooPolynomial('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
    elif fitFuncVjet == 'POLY4':
        VjetMass = RooPolynomial('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
    elif fitFuncVjet == 'POLY5':
        VjetMass = RooPolynomial('VjetMass', fitFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet, a4Vjet))
    elif fitFuncVjet == 'TEST':
        #p2Vjet = RooPolynomial('p2Vjet', 'p2Vjet', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
        #p3Vjet = RooChebychev('p3Vjet', 'p3Vjet', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
        #p3Vjet = RooPolynomial('p3Vjet', 'p3Vjet', J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet, a3Vjet))
        #VjetMass = RooAddPdf('VjetMass', fitFuncVjet, RooArgList(expoVjet, p2Vjet), RooArgList(fracVjet))
        #VjetMass = RooFFTConvPdf('VjetMass', fitFuncVjet, J_mass, expoVjet, p2Vjet)
        #a0Vjet.setVal(0)
        #a1Vjet.setVal(-20)
        #a2Vjet.setVal(1)
        #fracVjet.setVal(0.8)
        #erfexpVjet = RooErfExpPdf('erfexpVjet', 'erfexpVjet', J_mass, constVjet, offsetVjet, widthVjet)
        #VjetMass = RooAddPdf('VjetMass', fitFuncVjet, RooArgList(erfexpVjet, p2Vjet), RooArgList(fracVjet))
        VjetMass = RooErfPow2Pdf('VjetMass', fitFuncVjet, J_mass, a0Vjet, a1Vjet, offsetVjet, widthVjet)
        a0Vjet.setVal(10)
        a0Vjet.setMin(0)
        a0Vjet.setMax(20)
        a1Vjet.setVal(0.1)
    else:
        print '  ERROR! Pdf', fitFuncVjet, 'is not implemented for Vjets'
        exit()

    if fitAltFuncVjet == 'POL3': VjetMass2 = RooChebychev('VjetMass2', fitAltFuncVjet, J_mass, RooArgList(a0Vjet, a1Vjet, a2Vjet))
    elif fitAltFuncVjet == 'ERFEXP': VjetMass2 = RooErfExpPdf('VjetMass2', fitAltFuncVjet, J_mass, constVjet2, offsetVjet2, widthVjet2)#L#
    elif fitAltFuncVjet == 'EXPGAUS': VjetMass2  = RooAddPdf('VjetMass2',   fitAltFuncVjet, RooArgList(gausVjet2, expoVjet2), RooArgList(fracVjet2))#L#
    else:
        print '  ERROR! Pdf', fitAltFuncVjet, 'is not implemented for Vjets'
        exit()
   
    #Lisa: provo a simmetrizzare gli errori:
    if channel == 'XVZnnlp':
        #Funzione principale = EXPGAUS
        if fitFuncVjet == 'EXPGAUS':
            #g2meanVjet.setMin(0.)#8 feb
            #g2meanVjet.setMax(100.)#8 feb
            #g2meanVjet.setVal(50.)#8 feb
            #g2sigmaVjet.setMin(0.)#8 feb
            #g2sigmaVjet.setMax(80.)#8 feb
            #g2sigmaVjet.setVal(40.)#8 feb
            fracVjet.setMin(0.)
            fracVjet.setMax(1.)
            fracVjet.setVal(0.53)
            constVjet.setMin(-0.3)
            constVjet.setMax(-0.003)
            constVjet.setVal(-0.03)
            #setto costanti in cerca del colpevole di cotanto errore!
            #prima setto costanti le gaussiane
            g2meanVjet.setVal(49.913)
            g2meanVjet.setMin(38)
            g2meanVjet.setMax(64)
            #g2meanVjet.setConstant(True)
            #g2sigmaVjet.setVal(38.612)
            #g2sigmaVjet.setConstant(True)
        #Funzione principale = ERFEXP
        #if fitFuncVjet == 'EXPGAUS':
        #    constVjet.setVal(-0.04)
        #    constVjet.setMin(-1.)
        #    constVjet.setMax(0.)
        #    offsetVjet.setVal(50.)
        #    offsetVjet.setMin(40.)
        #    offsetVjet.setMax(100.)
        #    widthVjet.setVal(55.)
        #    widthVjet.setMin(0.)
        #    widthVjet.setMax(1000.)
        ######################
        #Funzione alternativa = EXPGAUS
        if fitAltFuncVjet == 'EXPGAUS':
            g2meanVjet2.setMin(0.)
            g2meanVjet2.setMax(100.)
            g2meanVjet2.setVal(50.)
            g2sigmaVjet2.setMin(0.)
            g2sigmaVjet2.setMax(80.)
            g2sigmaVjet2.setVal(40.)
            fracVjet2.setMin(0.)
            fracVjet2.setMax(1.)
            fracVjet2.setVal(0.53)
            constVjet2.setMin(-0.3)
            constVjet2.setMax(-0.003)
            constVjet2.setVal(-0.03)
            g2meanVjet2.setVal(49.913)
            g2meanVjet2.setMin(38)
            g2meanVjet2.setMax(64)
        #Funzione alternativa = ERFEXP
        #if fitAltFuncVjet == 'EXPGAUS':
        #    constVjet2.setVal(-0.04)
        #    constVjet2.setMin(-1.)
        #    constVjet2.setMax(0.)
        #    offsetVjet2.setVal(50.)
        #    offsetVjet2.setMin(40.)
        #    offsetVjet2.setMax(100.)
        #    widthVjet2.setVal(55.)
        #    widthVjet2.setMin(0.)
        #    widthVjet2.setMax(1000.)
    if channel == 'XVZnnhp':
        #Funzione principale = ERFEXP
        if fitFuncVjet == 'ERFEXP':
            constVjet.setVal(-0.0202)
            constVjet.setMin(-1.)
            constVjet.setMax(0.)
            offsetVjet.setVal(84.8)#(50.)
            offsetVjet.setMin(20.)#(40.)
            offsetVjet.setMax(200.)#(100.)
            widthVjet.setVal(43.7)#(55.)
            widthVjet.setMin(0.1)#(0.)
            widthVjet.setMax(400.)#(1000.)
        #Funzione principale = EXPGAUS
        if fitFuncVjet == 'EXPGAUS':
            #C#g2meanVjet.setMin(80.)
            #C#g2meanVjet.setMax(100.)
            #C#g2meanVjet.setVal(90.)
            g2sigmaVjet.setMin(10.)
            g2sigmaVjet.setMax(50.)
            g2sigmaVjet.setVal(30.)
            fracVjet.setMin(0.)
            fracVjet.setMax(1.)
            fracVjet.setVal(0.8)#0.5
            constVjet.setVal(-0.005)
            constVjet.setMin(-1.)
            constVjet.setMax(0.)
        ######################
        #Funzione alternativa = EXPGAUS
        if fitAltFuncVjet == 'EXPGAUS':#da uniformare!
            #C#g2meanVjet2.setMin(80.)
            #C#g2meanVjet2.setMax(100.)
            #C#g2meanVjet2.setVal(90.)
            g2sigmaVjet2.setMin(10.)
            g2sigmaVjet2.setMax(50.)
            g2sigmaVjet2.setVal(30.)
            fracVjet2.setMin(0.)
            fracVjet2.setMax(1.)
            fracVjet2.setVal(0.8)#0.5
            constVjet2.setVal(-0.005)
            constVjet2.setMin(-1.)
            constVjet2.setMax(0.)
        #Funzione alternativa = ERFEXP
        if fitAltFuncVjet == 'ERFEXP':
            constVjet2.setVal(-0.0202)
            constVjet2.setMin(-1.)
            constVjet2.setMax(0.)
            offsetVjet2.setVal(84.8)#(50.)
            offsetVjet2.setMin(20.)#(40.)
            offsetVjet2.setMax(200.)#(100.)
            widthVjet2.setVal(43.7)#(55.)
            widthVjet2.setMin(0.1)#(0.)
            widthVjet2.setMax(400.)#(1000.)

#    if fitFuncVjet == 'POL4':
#        frVjet3 = VjetMass3.fitTo(setVjet, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
#    if fitFuncVjet == 'POL5':
#        frVjet3 = VjetMass3.fitTo(setVjet, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
#        if VERBOSE or LISA: print '********** Fit result [POL3] *'+'*'*40, '\n', frVjet3.Print("v"), '\n', '*'*80
#        a2Vjet.setMax(0.5)
#        frVjet4 = VjetMass4.fitTo(setVjet, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
#        a3Vjet.setMin(0.)
#        #a4Vjet.setMax(0.)
#        a4Vjet.setVal(-0.0001)
#        #a4Vjet.setConstant(True)
#        if VERBOSE or LISA: print '********** Fit result [POL4] *'+'*'*40, '\n', frVjet4.Print("v"), '\n', '*'*80

    frVjet = VjetMass.fitTo(setVjet, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    if VERBOSE or LISA: print '********** Fit result [JET MASS Vjets] *'+'*'*40, '\n', frVjet.Print("v"), '\n', '*'*80


    frVjet2 = VjetMass2.fitTo(setVjet, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    #if VERBOSE or LISA: print '********** Fit result [JET MASS Vjets ALTERNATIVE] *'+'*'*40, '\n', frVjet2.Print("v"), '\n', '*'*80

    drawPlot('VjetMass', channel, J_mass, VjetMass, setVjet, [frVjet], -1, None, '', VjetMass2 if ALTERNATIVE else None, chi=1)

    iSBVjet_pre = VjetMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))
    iSRVjet_pre = VjetMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    iSBVjet2_pre = VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))
    iSRVjet2_pre = VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    print 'Debug: iSBVjet pre fit: ', iSBVjet_pre.getVal()
    print 'Debug: SB yield pre fit: ', iSBVjet_pre.getVal()*nMainEasy.getVal()
    print 'Debug: iSRVjet pre fit: ', iSRVjet_pre.getVal()
    print 'Debug: SR yield pre fit: ', iSRVjet_pre.getVal()*nMainEasy.getVal()
    print 'Debug: iSBVjet2 pre fit: ', iSBVjet2_pre.getVal()
    print 'Debug: SB yield2 pre fit: ', iSBVjet2_pre.getVal()*nMain2.getVal()    
    print 'Debug: iSRVjet2 pre fit: ', iSRVjet2_pre.getVal()
    print 'Debug: SR yield2 pre fit: ', iSRVjet2_pre.getVal()*nMain2.getVal()
    
    if SCAN:   
        likelihoodScan(VjetMass, setVjet, [constVjet, const2Vjet, offsetVjet, widthVjet],channel)
    #exit()

    #*******************************************************#
    #                                                       #
    #                 VV, VH normalization                  #
    #                                                       #
    #*******************************************************#
    
    # Variables for VV
    
    if channel=='XVZnnlp':
        offsetVV = RooRealVar('offsetVV', 'offset of the erf', 70.,     30., 300.)
        widthVV  = RooRealVar('widthVV',  'width of the erf',  50.,     1., 100.)
        constVV  = RooRealVar('constVV',  'slope of the exp',  -0.015, -0.1, 0.)
        meanVV   = RooRealVar('meanVV',   'mean of the gaussian', 98, 60., 100.)
        fracVV   = RooRealVar('fracVV',   'fraction of gaussian wrt erfexp', 0.65, 0.45, 1.)#0.002, 0., 1.)
        #sigmaVV  = RooRealVar('sigmaVV',  'sigma of the gaussian', 28., 6., 30.)
        #sigmaVV  = RooRealVar('sigmaVV',  'sigma of the gaussian', 8., 0., 30.)
        sigmaVV  = RooRealVar('sigmaVV',  'sigma of the gaussian', 8., 6., 11.)
        meanVH   = RooRealVar('meanVH',   'mean of the gaussian',    125.,   100., 150.)
        sigmaVH  = RooRealVar('sigmaVH',  'sigma of the gaussian',   10.,     5.,  50.)
        fracVH   = RooRealVar('fracVH',   'fraction of gaussian wrt erfexp',  1.5e-2, 0.,   1.)
    elif channel=='XVZnnhp':
        offsetVV = RooRealVar('offsetVV', 'offset of the erf', 70., 30., 300.)
        widthVV  = RooRealVar('widthVV',  'width of the erf',  50., 20., 100.)
        #constVV  = RooRealVar('constVV',  'slope of the exp',  -0.012, -0.1, 0.)#default
        constVV  = RooRealVar('constVV',  'slope of the exp',  -0.008, -0.1, 0.)
        meanVV   = RooRealVar('meanVV',   'mean of the gaussian', 90., 60., 100.)
        fracVV   = RooRealVar('fracVV',   'fraction of gaussian wrt erfexp', 0.6, 0., 1.)
        sigmaVV  = RooRealVar('sigmaVV',  'sigma of the gaussian', 27., 6., 30.)
        meanVH   = RooRealVar('meanVH',   'mean of the gaussian',           125.,   100., 150.)
        sigmaVH  = RooRealVar('sigmaVH',  'sigma of the gaussian',           10.,     5.,  50.)
        fracVH   = RooRealVar('fracVH',   'fraction of gaussian wrt erfexp',  1.5e-2, 0.,   1.)
        #devo prendere bene il trend vicino a 50-60 GeV altrimenti sovrastimo. innanzitutto fisso costanti e cerco di capire cosa alza e abbassa lo spettro all'inizio
        constVV.setVal(-8.2878e-03)#ed e' anche quellla con l'errore piu' alto!
        constVV.setMax(0.1)
        #constVV.setConstant(True)

    # Error function and exponential to model the bulk
    erfrVV   = RooErfExpPdf('baseVV', 'error function for VV jet mass', J_mass, constVV, offsetVV, widthVV)
    expoVV   = RooExponential('baseVV', 'error function for VV jet mass', J_mass, constVV)
    # gaussian for the V mass peak
    gausVV   = RooGaussian('gausVV',  'gaus for VV jet mass', J_mass, meanVV, sigmaVV)
    # gaussian for the H mass peak
    gausVH   = RooGaussian('gausVH',  'gaus for VH jet mass', J_mass, meanVH, sigmaVH)
        
    # Define VV model
    if fitFuncVV == 'ERFEXP': VVMass = RooErfExpPdf('VVMass', fitFuncVV, J_mass, constVV, offsetVV, widthVV)
    elif fitFuncVV == 'ERFEXPGAUS': VVMass  = RooAddPdf('VVMass',   fitFuncVV, RooArgList(gausVV, erfrVV), RooArgList(fracVV))
    elif fitFuncVV == 'ERFEXPGAUS2': VVMass  = RooAddPdf('VVMass',   fitFuncVV, RooArgList(gausVH, gausVV, erfrVV), RooArgList(fracVH, fracVV))
    elif fitFuncVV == 'EXPGAUS': VVMass  = RooAddPdf('VVMass',   fitFuncVV, RooArgList(gausVV, expoVV), RooArgList(fracVV))
    elif fitFuncVV == 'EXPGAUS2': VVMass  = RooAddPdf('VVMass',   fitFuncVV, RooArgList(gausVH, gausVV, expoVV), RooArgList(fracVH, fracVV))
    elif fitFuncVV == 'GAUS2': VVMass  = RooGaussian('VVMass',  fitFuncVV, J_mass, meanVV, sigmaVV)
    else:
        print '  ERROR! Pdf', fitFuncVV, 'is not implemented for VV'
        exit()

    # fit to secondary bkg in MC (whole range)
    frVV = VVMass.fitTo(setVV, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))

    drawPlot('VVMass', channel, J_mass, VVMass, setVV, [frVV],chi=1)
    
    if VERBOSE or LISA: print '********** Fit result [JET MASS VV] ****'+'*'*40, '\n', frVV.Print("v"), '\n', '*'*80
    
    #if SCAN:   
    #    likelihoodScan(VVMass, setVV, [offsetVV, widthVV],channel)    

    #*******************************************************#
    #                                                       #
    #                 Top, ST normalization                 #
    #                                                       #
    #*******************************************************#
        
    if channel=='XVZnnlp':
        #constTop  = RooRealVar('constTop',  'slope of the exp', -0.025, -1., 0.)#default
        constTop  = RooRealVar('constTop',  'slope of the exp', -0.0431, -1., 0.)#-0.05 pare bene;0.03 converge veloce;0.04 bene
        #fracT     = RooRealVar('fracT',     'fraction of gaussian wrt erfexp', 0.06, 0., 1.)#default
        #fracW     = RooRealVar('fracW',     'fraction of gaussian wrt erfexp', 0.05, 0., 1.)#default
        fracT     = RooRealVar('fracT',     'fraction of gaussian wrt erfexp', 0.018, 0., 1.)
        fracW     = RooRealVar('fracW',     'fraction of gaussian wrt erfexp', 0.15, 0., 1.)#0.06bene
        meanT     = RooRealVar('meanT',     'mean of the gaussian', 165., 150., 200.)#default
        #meanW     = RooRealVar('meanW',     'mean of the gaussian', 80., 60., 100.)#default
        meanW     = RooRealVar('meanW',     'mean of the gaussian', 80., 70., 90.)
#        offsetTop = RooRealVar('offsetTop', 'offset of the erf', 38.55, 10., 300.)#default
        offsetTop = RooRealVar('offsetTop', 'offset of the erf', 140, 10., 300.)#40 bene
        sigmaT    = RooRealVar('sigmaT',    'sigma of the gaussian', 15., 5., 40.)#default
#        sigmaW    = RooRealVar('sigmaW',    'sigma of the gaussian', 10., 6., 30.)#default
        sigmaW    = RooRealVar('sigmaW',    'sigma of the gaussian', 7., 4., 10.)#8 bene,10uguale
        widthTop  = RooRealVar('widthTop',  'width of the erf', 63., 10., 300.)#default 45
    elif channel=='XVZnnhp':
        constTop  = RooRealVar('constTop',  'slope of the exp', -0.025, -1., 0.)#default
        constTop  = RooRealVar('constTop',  'slope of the exp', -0.04, -1., 0.)
        fracT     = RooRealVar('fracT',     'fraction of gaussian wrt erfexp', 0.06, 0., 1.)
        fracW     = RooRealVar('fracW',     'fraction of gaussian wrt erfexp', 0.1, 0., 1.)
#        meanT     = RooRealVar('meanT',     'mean of the gaussian', 175., 150., 200.)#default
#        meanW     = RooRealVar('meanW',     'mean of the gaussian', 80., 60., 100.)#default
        #meanT     = RooRealVar('meanT',     'mean of the gaussian', 175., 170., 180.)#buono
        meanT     = RooRealVar('meanT',     'mean of the gaussian', 130., 120., 180.)#prova
        meanW     = RooRealVar('meanW',     'mean of the gaussian', 80., 70., 90.)#?
        offsetTop = RooRealVar('offsetTop', 'offset of the erf', 38.55, 10., 300.)
#        sigmaT    = RooRealVar('sigmaT',    'sigma of the gaussian', 40., 5., 50.)#default
#        sigmaW    = RooRealVar('sigmaW',    'sigma of the gaussian', 10., 6., 30.)#default
        sigmaT    = RooRealVar('sigmaT',    'sigma of the gaussian', 30., 5., 50.)#30-5-50
        sigmaW    = RooRealVar('sigmaW',    'sigma of the gaussian', 10., 1., 30.)
        widthTop  = RooRealVar('widthTop',  'width of the erf',  45., 10., 300.)
#Lisa: fisso un po' di roba e vediamo
        ###constTop.setVal(-4.e-02)#(-6.6233e-02)#compatibile con zero
        ###constTop.setMin(-1.e-01)
        ###constTop.setMax(-1.e-3)
        #constTop.setConstant(True)
        ###fracT.setVal(1.9278e-02)#errore maggiore del valore
        ###fracT.setMin(0.)
        ###fracT.setMax(0.1)
        #fracT.setConstant(True)
        ###meanT.setVal(130.)#(1.7043e+02)#molto asimmetrico
        ###meanT.setMax(180.)#(190.)#lo rende simmetrico ma peggiora la banda
        ###meanT.setMin(120.)#(150.)#lo rende simmetrico ma peggiora la banda
        #meanT.setConstant(True)
        ###sigmaW.setMin(1.)
        ###sigmaW.setMax(10.)
        ###sigmaT.setVal(1.2339e+01)#un po' asimmetrico
        #sigmaT.setConstant(True)

#    elif channel=='XVZnnhp':
#        constTop  = RooRealVar('constTop',  'slope of the exp', -0.015, -1., 0.)
#        fracT     = RooRealVar('fracT',     'fraction of gaussian wrt erfexp', 0.1, 0., 1.)
#        fracW     = RooRealVar('fracW',     'fraction of gaussian wrt erfexp', 0.3, 0., 1.)
#        meanT     = RooRealVar('meanT',     'mean of the gaussian', 175., 160., 190.)
#        meanW     = RooRealVar('meanW',     'mean of the gaussian', 90., 70., 110.)
#        offsetTop = RooRealVar('offsetTop', 'offset of the erf', 73., 50., 250.)
#        sigmaT    = RooRealVar('sigmaT',    'sigma of the gaussian', 20., 10., 50.)
#        sigmaW    = RooRealVar('sigmaW',    'sigma of the gaussian', 30., 5., 60.)
#        widthTop  = RooRealVar('widthTop',  'width of the erf',  40., 10., 300.)
    
    # Variables for Top
    # Error Function * Exponential to model the bulk
    gausTop   = RooGaussian('baseTop',  'gaus for Top jet mass', J_mass, offsetTop, widthTop)
    erfrTop   = RooErfExpPdf('baseTop', 'error function for Top jet mass', J_mass, constTop, offsetTop, widthTop)
    # gaussian for the W mass peak
    gausW     = RooGaussian('gausW',    'gaus for W jet mass', J_mass, meanW, sigmaW)
    # gaussian for the Top mass peak
    gausT     = RooGaussian('gausT',    'gaus for T jet mass', J_mass, meanT, sigmaT)
    
    # Define Top model
    if fitFuncTop == 'ERFEXP': TopMass = RooErfExpPdf('TopMass', fitFuncTop, J_mass, constTop, offsetTop, widthTop)
    elif fitFuncTop == 'ERFEXPGAUS2': TopMass = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausW, gausT, erfrTop), RooArgList(fracW, fracT))
    elif fitFuncTop == 'ERFEXPGAUS': 
        if channel=='XVZmmhp' or channel=='XVZeehp' or channel=='XVZnnhp':
            TopMass = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausW, erfrTop), RooArgList(fracW))
        else:
            TopMass = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausT, erfrTop), RooArgList(fracT))
    #elif fitFuncTop == 'EXPGAUS2': TopMass  = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausW, gausT,... expoVV), RooArgList(fracVH, fracVV))#L
    elif fitFuncTop == 'GAUS3': TopMass  = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausW, gausT, gausTop), RooArgList(fracW, fracT))
    elif fitFuncTop == 'GAUS2': TopMass  = RooAddPdf('TopMass',   fitFuncTop, RooArgList(gausT, gausTop), RooArgList(fracT))
    elif fitFuncTop == 'GAUS': TopMass  = RooGaussian('TopMass', fitFuncTop, J_mass, offsetTop, widthTop)
    else:
        print '  ERROR! Pdf', fitFuncTop, 'is not implemented for Top'
        exit()
    
    # fit to secondary bkg in MC (whole range)
    if channel =="XVZnnlp":
        constTop.setVal(-3.8399e-02)
        constTop.setConstant(True)#L + J#
    frTop = TopMass.fitTo(setTop, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.Minos(False if channel == 'XVZmmlp' else True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    ###frTop = TopMass.fitTo(setTop, RooFit.SumW2Error(True), RooFit.Range('V_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit'), RooFit.Minos(True), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    
    drawPlot('TopMass', channel, J_mass, TopMass, setTop, [frTop],chi=1)
    
    if VERBOSE or LISA: print '********** Fit result [JET MASS TOP] ***'+'*'*40, '\n', frTop.Print("v"), '\n', '*'*80
    
    #if SCAN:
    #    likelihoodScan(TopMass, setTop, [offsetTop, widthTop],channel)

    #*******************************************************#
    #                                                       #
    #                 All bkg normalization                 #
    #                                                       #
    #*******************************************************#
    
    if EXTRAPOLATE:
        if lowP:
            Vjet_var2.setConstant(True)
        if higP:
            Vjet_var2.setConstant(True)
        #constVjet.setConstant(False)
        #offsetVjet.setConstant(False)
        #widthVjet.setConstant(True)
        #a0Vjet.setConstant(False)
        #a1Vjet.setConstant(False)
        #a2Vjet.setConstant(False)

    if channel == 'XVZmmlp' : Vjet_var2.setConstant(True)

    constVV.setConstant(True)
    offsetVV.setConstant(True)
    widthVV.setConstant(True)
    meanVV.setConstant(True)
    sigmaVV.setConstant(True)
    fracVV.setConstant(True)
    meanVH.setConstant(True)
    sigmaVH.setConstant(True)
    fracVH.setConstant(True)
    
    constTop.setConstant(True)
    offsetTop.setConstant(True)
    widthTop.setConstant(True)
    meanW.setConstant(True)
    sigmaW.setConstant(True)
    fracW.setConstant(True)
    meanT.setConstant(True)
    sigmaT.setConstant(True)
    fracT.setConstant(True)
    
    nVV.setConstant(True)
    nTop.setConstant(True)
    nMain.setConstant(False)
    nMain2.setConstant(False)
    nMainEasy.setConstant(False)
    #cheating:
    #a0Vjet.setVal(-1.7)#alzato il val assoluto: tira su linearmente, non va bene
    #a0Vjet.setConstant(True)
    #a1Vjet.setVal(0.82)#se lo aumento in valore assoluto:
    #a1Vjet.setConstant(True)
    #a2Vjet.setConstant(True)
    #a3Vjet.setConstant(True)

    
    # Final background model by adding the main+secondary pdfs (using 'coef': ratio of the secondary/main, from MC)
    TopMass_ext  = RooExtendPdf('TopMass_ext',  'extended p.d.f', TopMass,  nTop)
    VVMass_ext   = RooExtendPdf('VVMass_ext',   'extended p.d.f', VVMass,   nVV)
    VjetMass_ext = RooExtendPdf('VjetMass_ext', 'extended p.d.f', VjetMass, nMain)
    VjetMassEasy_ext = RooExtendPdf('VjetMassEasy_ext', 'extended p.d.f', VjetMass, nMainEasy)
    VjetMass2_ext = RooExtendPdf('VjetMass2_ext', 'extended p.d.f', VjetMass2, nMain2)
    BkgMass = RooAddPdf('BkgMass', 'BkgMass', RooArgList(TopMass_ext, VVMass_ext, VjetMass_ext), RooArgList(nTop, nVV, nMain))
    BkgMass2 = RooAddPdf('BkgMass2', 'BkgMass2', RooArgList(TopMass_ext, VVMass_ext, VjetMass2_ext), RooArgList(nTop, nVV, nMain2))
    BkgMass.fixAddCoefRange('V_reasonable_range')
    BkgMass2.fixAddCoefRange('V_reasonable_range')
    #L
    BkgMassEasy = RooAddPdf('BkgMassEasy', 'BkgMassEasy', RooArgList(TopMass_ext, VVMass_ext, VjetMassEasy_ext), RooArgList(nTop, nVV, nMainEasy))#L
    #BkgMassEasy = RooAddPdf('BkgMassEasy', 'BkgMassEasy', RooArgList(VjetMassEasy_ext), RooArgList(nMainEasy))#L#questa funziona
    BkgMassEasy.fixAddCoefRange('V_reasonable_range')
    
    # Extended fit model to data in SB
    
    frMassEasy = BkgMassEasy.fitTo(setDataSB, RooFit.SumW2Error(True), RooFit.Extended(True), RooFit.Range('LSBrange,HSBrange'), RooFit.Strategy(2), RooFit.Minimizer('Minuit'), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    if VERBOSE or LISA: print '********** Fit result [JET MASS DATA] **'+'*'*40, '\n', frMassEasy.Print("v"), '\n', '*'*80

    frMass2 = BkgMass2.fitTo(setDataSB, RooFit.SumW2Error(True), RooFit.Extended(True), RooFit.Range('LSBrange,HSBrange'), RooFit.Strategy(2), RooFit.Minimizer('Minuit'), RooFit.Save(1), RooFit.PrintLevel(1 if VERBOSE else -1), RooFit.NumCPU(24))
    if VERBOSE or LISA: print '********** Fit result [Alt JET MASS DATA] **'+'*'*40, '\n', frMass2.Print("v"), '\n', '*'*80


    #L#frMass2.floatParsFinal().Print("s")#eliminare
    #L#frMassEasy.floatParsFinal().Print("s")#eliminare

    #if SCAN:
    #   likelihoodScan(VjetMass, setVjet, [constVjet, offsetVjet, widthVjet],channel)
    
    # Fix normalization and parameters of V+jets after the fit to data
    nMain.setConstant(True)
    nMainEasy.setConstant(True)
    nMain2.setConstant(True)
    print "Debug nMain2: ", nMain2.getVal() 
    
    constVjet.setConstant(True)
    offsetVjet.setConstant(True)
    widthVjet.setConstant(True)
    a0Vjet.setConstant(True)
    a1Vjet.setConstant(True)
    a2Vjet.setConstant(True)
    
    # integrals for global normalization
    # do not integrate the composte model: results have no sense
    
    # integral for normalization in the SB
    iSBVjet = VjetMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))
    iSBVjet2 = VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))
    iAllVjet = VjetMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg))
    iAllVjet2 = VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg))
    iSBVV = VVMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))
    iSBTop = TopMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('LSBrange,HSBrange'))

    ## integral for normalization in the SR
    iSRVjet = VjetMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    iSRVV = VVMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    iSRTop = TopMass.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    iSRVjet2 = VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))#try
    print "---------------"
    print "---------------"
    print "Debug: iSBVjet %.5f, iSRVjet %.5f, iAllVjet %.5f" % (iSBVjet.getVal(), iSRVjet.getVal(), iAllVjet.getVal())
    print "Debug: iSBVjet2 %.5f, iSRVjet2 %.5f, iAllVjet2 %.5f" % (iSBVjet2.getVal(), iSRVjet2.getVal(), iAllVjet2.getVal())
    print "---------------"    
    print "---------------"    
    #Dummy function just to check true integral
    #DummyexpoVjet2 = RooExponential('DummyexpoVjet2', 'expoVjet2', J_mass, constVjet2)
    #DummygausVjet2 = RooGaussian('DummygausVjet2', 'gausVjet2', J_mass, g2meanVjet2, g2sigmaVjet2)
    #DummyfracVjet2 = RooRealVar('DummyfracVjet2',   'fraction of gaussian wrt exp', fracVjet2.getVal(), 0., 1.)
    #DummyVjetMass2  = RooAddPdf('DummyVjetMass2', fitAltFuncVjet, RooArgList(DummygausVjet2, DummyexpoVjet2), RooArgList(DummyfracVjet2))
    #DummyiSRVjet2 = DummyVjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    #print "Debug DummyiSRVjet2: ", DummyiSRVjet2.getVal()
    #BooVjetMass2  = RooAddPdf('BooVjetMass2', fitAltFuncVjet, RooArgList(gausVjet2, expoVjet2), RooArgList(fracVjet2))
    #BooiSRVjet2 = BooVjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    #rif = VjetMass2.clone("rif")
    #rifaccio = rif.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))
    #print "Debug BooiSRVjet2: ", BooiSRVjet2.getVal()
    #print "========================================="
    #print "Debug iSRVjet2 print: ", iSRVjet2.Print("v")
    #print "Debug iSRVjet2 val: ", iSRVjet2.getVal()
    #print "Debug DummyiSRVjet2 print: ", DummyiSRVjet2.Print("v")
    #print "Debug DummyiSRVjet2 val: ", DummyiSRVjet2.getVal()
    #print "Debug BooiSRVjet2 print: ", BooiSRVjet2.Print("v")
    #print "Debug BooiSRVjet2 val: ", BooiSRVjet2.getVal()
    #print VjetMass2.expectedEvents(jetMassArg)
    #print BooVjetMass2.expectedEvents(jetMassArg)
    #print "rifaccio: ", rifaccio.getVal()
    #print "========================================="

    
    ## formual vars
    SByieldEasy = RooFormulaVar('SByieldEasy', 'extrapolation to SB', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSBVjet, nMainEasy, iSBVV, nVV, iSBTop, nTop))
    #SByieldEasy = RooFormulaVar('SByield', 'extrapolation to SB', '@0*@1', RooArgList(iSBVjet, nMainEasy))#questa funziona
    SByield2 = RooFormulaVar('SByield2', 'extrapolation to SB', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSBVjet2, nMain2, iSBVV, nVV, iSBTop, nTop))
    SByield = RooFormulaVar('SByield', 'extrapolation to SB', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSBVjet, nMain, iSBVV, nVV, iSBTop, nTop))
    #HRyield = RooFormulaVar('HRyield', 'extrapolation to HR', '@0*@1 + @2*@3 + @4*@5', RooArgList(iHRVjet, nMain, iHRVV, nVV, iHRTop, nTop))
    #SRyield = RooFormulaVar('SRyield', 'extrapolation to SR', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSRVjet, nMain, iSRVV, nVV, iSRTop, nTop))
    SRyieldEasy = RooFormulaVar('SRyieldEasy', 'extrapolation to SR', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSRVjet, nMainEasy, iSRVV, nVV, iSRTop, nTop))
    SRyield2 = RooFormulaVar('SRyield2', 'extrapolation to SR', '@0*@1 + @2*@3 + @4*@5', RooArgList(iSRVjet2, nMain2, iSRVV, nVV, iSRTop, nTop))
    #MainSRyield = RooFormulaVar('MainSRyield', 'extrapolation to SR', '@0*@1', RooArgList(iSRVjet, nMain))
    MainSRyieldEasy = RooFormulaVar('MainSRyieldEasy', 'extrapolation to SR', '@0*@1', RooArgList(iSRVjet, nMainEasy))#tried,works
    #print "Debug iSRVjet2: ", iSRVjet2.getVal() 
    MainSRyield2 = RooFormulaVar('MainSRyield2', 'extrapolation to SR', '@0*@1', RooArgList(iSRVjet2, nMain2))#try
    #print "Debug MainSRyield2: ", MainSRyield2.getVal() 
    VVSRyield = RooFormulaVar('VVSRyield', 'extrapolation to SR', '@0*@1', RooArgList(iSRVV, nVV))#tried,works
    TopSRyield = RooFormulaVar('TopSRyield', 'extrapolation to SR', '@0*@1', RooArgList(iSRTop, nTop))#tried,works
    #yield nelle SB
    MainSByieldEasy = RooFormulaVar('MainSByieldEasy', 'extrapolation to SB', '@0*@1', RooArgList(iSBVjet, nMainEasy))#tried,works
    MainSByield2 = RooFormulaVar('MainSByield2', 'extrapolation to SB', '@0*@1', RooArgList(iSBVjet2, nMain2))#try
    VVSByield = RooFormulaVar('VVSByield', 'extrapolation to SB', '@0*@1', RooArgList(iSBVV, nVV))#tried,works
    TopSByield = RooFormulaVar('TopSByield', 'extrapolation to SB', '@0*@1', RooArgList(iSBTop, nTop))#tried,works
    ## alternative method for main bakground estimation
    SByieldRatio = RooFormulaVar('SByieldRatio', 'extrapolation to SB','@0-@1*@3-@2*@4', RooArgList(entrySB, nVV, nTop, iSBVV, iSBTop))#tried,works
    SRyieldRatio = RooFormulaVar('SRyieldRatio', 'extrapolation to SR','@0*@1/@2', RooArgList(SByieldRatio, iSRVjet, iSBVjet))#tried,works
    BkgYieldRatio            = SRyieldRatio.getVal()#tried,works
    BkgYieldRatio_error      = SRyieldRatio.getPropagatedError(frVjet)#tried,works
    
    ## fractions
    fMainSB = RooRealVar('fMainSB', 'Fraction of Vjet events in SB', iSBVjet.getVal()*nMainEasy.getVal()/SByieldEasy.getVal(), 0., 1.)
    fVVSB = RooRealVar('fVVSB', 'Fraction of VV events in SB', iSBVV.getVal()*nVV.getVal()/SByieldEasy.getVal(), 0., 1.)
    fTopSB = RooRealVar('fTopSB', 'Fraction of Top events in SB', iSBTop.getVal()*nTop.getVal()/SByieldEasy.getVal(), 0., 1.)
    
    fMainSR = RooRealVar('fMainSR', 'Fraction of Vjet events in SR', iSRVjet.getVal()*nMainEasy.getVal()/SRyieldEasy.getVal(), 0., 1.)#tried,works
    fVVSR = RooRealVar('fVVSR', 'Fraction of VV events in SR', iSRVV.getVal()*nVV.getVal()/SRyieldEasy.getVal(), 0., 1.)
    fTopSR = RooRealVar('fTopSR', 'Fraction of Top events in SR', iSRTop.getVal()*nTop.getVal()/SRyieldEasy.getVal(), 0., 1.)
    
    ## final normalization values
    BkgYield            = SRyieldEasy.getVal()
    #BkgYield            = SRyield.getVal()
    #BkgYield2           = (VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))).getVal()*nMain2.getVal() + iSRVV.getVal()*nVV.getVal() + iSRTop.getVal()*nTop.getVal()#try
    BkgYield2            = SRyield2.getVal()
    #BkgYield2           = (VjetMass2.createIntegral(jetMassArg, RooFit.NormSet(jetMassArg), RooFit.Range('SRrange'))).getVal()*nMain2.getVal() + iSRVV.getVal()*nVV.getVal() + iSRTop.getVal()*nTop.getVal()
    #BkgYield_syst       = math.sqrt(SRyield.getPropagatedError(frVV)**2 + SRyield.getPropagatedError(frTop)**2)
    BkgYield_syst       = math.sqrt(SRyieldEasy.getPropagatedError(frVV)**2 + SRyieldEasy.getPropagatedError(frTop)**2)#tried,works
    BkgYield_systSB     = math.sqrt(SByieldEasy.getPropagatedError(frVV)**2 + SByieldEasy.getPropagatedError(frTop)**2)#tried,works
    #BkgYield_stat       = math.sqrt(SRyield.getPropagatedError(frMass)**2)
    BkgYield_stat       = math.sqrt(SRyieldEasy.getPropagatedError(frMassEasy)**2)#tried,works
    BkgYield_statSB       = math.sqrt(SByieldEasy.getPropagatedError(frMassEasy)**2)#tried,works
    #LISA#Warning! This is without alternative syst!
    BkgYield_alte       = abs(BkgYield - BkgYield2)#try
    BkgYield_alteSB     = abs(SByieldEasy.getVal() - SByield2.getVal())#try
    print 'Bkg yield ', BkgYield, ' ; Bk yield Alternative ', BkgYield2
    print 'Bkg yield alte abs diff ', BkgYield_alte
    #BkgYield_alte       = 0.#tried,works
    nBkgSR              = RooRealVar('nBkgSR', 'expected yield in SR', BkgYield, 0., 1.e7)
    #LISA#Warning! This is without alternative syst!
    nBkgSR.setError(math.sqrt( BkgYield_stat**2 + BkgYield_syst**2 + BkgYield_alte**2 + (BkgYield*fTopSR.getVal()*WSFErr[higP])**2 + (BkgYield*fVVSR.getVal()*WSFErr[higP])**2 ))

    ## backgrounds, separately
    #MainYield           = MainSRyield.getVal()
    #MainYieldErr        = MainSRyield.getPropagatedError(frMass)
    MainYield2          = MainSRyield2.getVal()
    MainYield2Err       = MainSRyield2.getPropagatedError(frMass2)
    MainYieldEasy       = MainSRyieldEasy.getVal()#tried,works
    MainYieldEasyErr    = MainSRyieldEasy.getPropagatedError(frMassEasy)#tried,works
    VVYield             = VVSRyield.getVal()#tried,works
    VVYieldErr          = VVSRyield.getPropagatedError(frVV)#tried,works
    TopYield            = TopSRyield.getVal()#tried,works
    TopYieldErr         = TopSRyield.getPropagatedError(frTop)#tried,works
    
    print '-----LISA-----'
    if channel=="XVZnnlp":
        print 'Nella SB: '
#    print '%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATI')
        print '%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATA')
        print '%.1f&\t%.1f&\t%.1f&\t%.1f&\t%.1f&\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.0f' % (MainSByieldEasy.getVal(), MainSByield2.getVal(),TopSByield.getVal(),VVSByield.getVal(), SByieldEasy.getVal(), BkgYield_statSB, BkgYield_systSB, BkgYield_alteSB, setDataSB.sumEntries())
        print 'Tot SB uncertainty: %.1f, in percent: %.1f' % (math.sqrt(BkgYield_statSB*BkgYield_statSB + BkgYield_systSB*BkgYield_systSB + BkgYield_alteSB*BkgYield_alteSB), math.sqrt(BkgYield_statSB*BkgYield_statSB + BkgYield_systSB*BkgYield_systSB + BkgYield_alteSB*BkgYield_alteSB)*100/SByieldEasy.getVal())
#    print '%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.0f' % (MainSByieldEasy.getVal(), MainSByield2.getVal(),TopSByield.getVal(),VVSByield.getVal(), SByieldEasy.getVal(), BkgYield_statSB, BkgYield_systSB, BkgYield_alteSB, setDataSB.sumEntries())
        print 'Nella SR: '
        print '%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATA')
        print '%.1f&\t%.1f&\t\t%.1f&\t%.1f&\t%.1f&\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.0f' % (MainSRyieldEasy.getVal(), MainSRyield2.getVal(),TopSRyield.getVal(),VVSRyield.getVal(), SRyieldEasy.getVal(), BkgYield_stat, BkgYield_syst, BkgYield_alte, setDataSR.sumEntries())
        print 'Tot SR uncertainty: %.1f, in percent: %.1f' % (math.sqrt(BkgYield_stat*BkgYield_stat + BkgYield_syst*BkgYield_syst + BkgYield_alte*BkgYield_alte), math.sqrt(BkgYield_stat*BkgYield_stat + BkgYield_syst*BkgYield_syst + BkgYield_alte*BkgYield_alte)*100/SRyieldEasy.getVal())
    if channel=="XVZnnhp":
        print 'Nella SB: '
#    print '%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s%-20s' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATI')
        print '%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATA')
        print '%.1f&\t%.1f&\t\t%.1f&\t%.1f&\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.0f' % (MainSByieldEasy.getVal(), MainSByield2.getVal(),TopSByield.getVal(),VVSByield.getVal(), SByieldEasy.getVal(), BkgYield_statSB, BkgYield_systSB, BkgYield_alteSB, setDataSB.sumEntries())
        print 'Tot SB uncertainty: %.1f, in percent: %.1f' % (math.sqrt(BkgYield_statSB*BkgYield_statSB + BkgYield_systSB*BkgYield_systSB + BkgYield_alteSB*BkgYield_alteSB), math.sqrt(BkgYield_statSB*BkgYield_statSB + BkgYield_systSB*BkgYield_systSB + BkgYield_alteSB*BkgYield_alteSB)*100/SByieldEasy.getVal())
#    print '%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.1f&%-20.0f' % (MainSByieldEasy.getVal(), MainSByield2.getVal(),TopSByield.getVal(),VVSByield.getVal(), SByieldEasy.getVal(), BkgYield_statSB, BkgYield_systSB, BkgYield_alteSB, setDataSB.sumEntries())
        print 'Nella SR: '
        print '%s\t%s\t%s\t%s\t%s\t\t%s\t%s\t%s\t%s\t' % ('v+jet', 'v+jet alt', 'top', 'VV', 'TOT', 'stat unc', 'syst unc.', 'syst alt.', 'DATA')
        print '%.1f&\t%.1f&\t\t%.1f&\t%.1f&\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.1f&\t\t%.0f' % (MainSRyieldEasy.getVal(), MainSRyield2.getVal(),TopSRyield.getVal(),VVSRyield.getVal(), SRyieldEasy.getVal(), BkgYield_stat, BkgYield_syst, BkgYield_alte, setDataSR.sumEntries())
        print 'Tot SR uncertainty: %.1f, in percent: %.1f' % (math.sqrt(BkgYield_stat*BkgYield_stat + BkgYield_syst*BkgYield_syst + BkgYield_alte*BkgYield_alte), math.sqrt(BkgYield_stat*BkgYield_stat + BkgYield_syst*BkgYield_syst + BkgYield_alte*BkgYield_alte)*100/SRyieldEasy.getVal())

    print '--- Channel', channel, '---'
    print '@', channel, ' bkg composition : V+jets %.3f+/-%.3f (%.1f%%),   VV %.3f+/-%.3f (%.1f%%),   Top %.3f+/-%.3f (%.1f%%)' % (MainYieldEasy, MainYieldEasyErr, fMainSR.getVal()*100, VVYield, VVYieldErr, fVVSR.getVal()*100, TopYield, TopYieldErr, fTopSR.getVal()*100)#tried,works
    print '@', channel, ' number of events: Integral = $%.1f$ & $\pm %.1f$ & $\pm %.1f$ & $\pm %.1f$ & $%.0f$' % (BkgYield, BkgYield_stat, BkgYield_syst, BkgYield_alte, setDataSR.sumEntries() if not False else -1)#tried,works
    #print '#', channel, ' main:', MainYield, '+-', MainYieldErr, ', alternate:', BkgYieldRatio, '+-', BkgYieldRatio_error, '(fit) +-', math.sqrt(SByieldRatio.getVal()), '(stat)'
    print '#', channel, ' main:', MainYieldEasy, '+-', MainYieldEasyErr, ', alternate:', BkgYieldRatio, '+-', BkgYieldRatio_error, '(fit) +-', math.sqrt(SByieldRatio.getVal()), '(stat)'#tried,works
    print '-'*11*2
    
    print 'Vjet in SB da MC: ', setVjetSB.sumEntries()
    print 'Data in SB: ', entrySB.getVal()

    drawPlot('JetMassEasy', channel, J_mass, BkgMassEasy, setDataSB if BLIND else setDataSRSB, [frMassEasy], SByieldEasy.getVal())
    drawPlot(channel+'_JetMass', channel, J_mass, BkgMassEasy, setDataSB if BLIND else setDataSRSB, [frMassEasy], SByieldEasy.getVal())

    #exit()
    # ====== CONTROL VALUE ======

    wx = RooWorkspace('Norm','workspace')
    getattr(wx, 'import')(fMainSB, RooFit.Rename(fMainSB.GetName()))
    getattr(wx, 'import')(fMainSR, RooFit.Rename(fMainSR.GetName()))
    getattr(wx, 'import')(fVVSB, RooFit.Rename(fVVSB.GetName()))
    getattr(wx, 'import')(fVVSR, RooFit.Rename(fVVSR.GetName()))
    getattr(wx, 'import')(fTopSB, RooFit.Rename(fTopSB.GetName()))
    getattr(wx, 'import')(fTopSR, RooFit.Rename(fTopSR.GetName()))
    BkgYieldVar = RooRealVar('BkgYieldVar','BkgYieldVar', BkgYield, 0., 1.e7)
    getattr(wx, 'import')(BkgYieldVar, RooFit.Rename(BkgYieldVar.GetName()))
    getattr(wx, 'import')(nBkgSR, RooFit.Rename(nBkgSR.GetName()))
    MainYieldEasyVar = RooRealVar('MainYieldEasyVar','MainYieldEasyVar',MainYieldEasy, 0., 1.e7)
    MainYieldEasyErrVar = RooRealVar('MainYieldEasyErrVar','MainYieldEasyErrVar',MainYieldEasyErr, 0., 1.e7)
    MainYield2Var = RooRealVar('MainYield2Var','MainYield2Var',MainYield2, 0., 1.e7)
    TopYieldVar = RooRealVar('TopYieldVar','TopYieldVar',TopYield, 0., 1.e7)
    TopYieldErrVar = RooRealVar('TopYieldErrVar','TopYieldErrVar',TopYieldErr, 0., 1.e7)
    VVYieldVar = RooRealVar('VVYieldVar','VVYieldVar',VVYield, 0., 1.e7)
    VVYieldErrVar = RooRealVar('VVYieldErrVar','VVYieldErrVar',VVYieldErr, 0., 1.e7)
    getattr(wx, 'import')(MainYieldEasyVar, RooFit.Rename(MainYieldEasyVar.GetName()))
    getattr(wx, 'import')(MainYieldEasyErrVar, RooFit.Rename(MainYieldEasyErrVar.GetName()))
    getattr(wx, 'import')(MainYield2Var, RooFit.Rename(MainYield2Var.GetName()))
    getattr(wx, 'import')(TopYieldVar, RooFit.Rename(TopYieldVar.GetName()))
    getattr(wx, 'import')(TopYieldErrVar, RooFit.Rename(TopYieldErrVar.GetName()))
    getattr(wx, 'import')(VVYieldVar, RooFit.Rename(VVYieldVar.GetName()))
    getattr(wx, 'import')(VVYieldErrVar, RooFit.Rename(VVYieldErrVar.GetName()))
    if not EXTRAPOLATE:
        wx.writeToFile('workspaces/%s_norm.root' % channel, True)
    else:
        wx.writeToFile('workspaces/%s_normExt.root' % channel, True)
    print "MainYieldEasy: ", MainYieldEasy
    print "MainYield2: ", MainYield2
    print "Norm workspace written: workspaces/%s_norm.root" % channel
    print "General debug:"
    print "fMainSB", fMainSB.getVal()
    print "fMainSR", fMainSR.getVal()
    print "fVVSB", fVVSB.getVal()
    print "fVVSR", fVVSR.getVal()
    print "fTopSB", fTopSB.getVal()
    print "fTopSR", fTopSR.getVal()
    print "BkgYieldVar", BkgYieldVar.getVal()
    print "nBkgSR", nBkgSR.getVal()
    print "MainYieldEasyVar", MainYieldEasyVar.getVal()
    print "MainYieldEasyErrVar", MainYieldEasyErrVar.getVal()
    print "MainYield2Var", MainYield2Var.getVal()
    print "TopYieldVar", TopYieldVar.getVal()
    print "TopYieldErrVar", TopYieldErrVar.getVal()
    print "VVYieldVar", VVYieldVar.getVal()
    print "VVYieldErrVar", VVYieldErrVar.getVal()
    

############################################33
#def funzioni
############################################33
 
def getColor(name, channel):
    if 'Top' in name: return 798
    elif 'VV' in name: return 602
    elif name in samples.keys(): return samples[name]['linecolor']
    elif 'XWh' in channel: return samples['WJetsToLNu']['linecolor']
    elif 'XZhnn' or 'XVZnn' in channel: return samples['DYJetsToNuNu_PtZ']['linecolor']
    else: return samples['DYJetsToLL']['linecolor']

def fixData(hist):
    for i in range(0, hist.GetN()):
        if (hist.GetX()[i]>=HMIN and hist.GetX()[i]<=HMAX and hist.GetY()[i]==0) or abs(hist.GetY()[i])<=1.e-6:
            hist.SetPoint(i, hist.GetX()[i], -1.e-4)
            hist.SetPointError(i, hist.GetErrorXlow(i), hist.GetErrorXhigh(i), 1.e-6, 1.e-6, )
        if hist.GetErrorXlow(i)<1.e-4:
            binwidth = hist.GetX()[1]-hist.GetX()[0]
            hist.SetPointEXlow(i, binwidth/2.)
            hist.SetPointEXhigh(i, binwidth/2.)


def drawPlot(name, channel, variable, model, dataset, fitRes=[], norm=-1, reg=None, cat='', alt=None, signal=None, snorm=-1, chi=0):
    isData = norm>0
    isMass = 'Mass' in name
    isSignal = 'X' in name and '_M' in name
    isCategory = reg is not None
    postfix = 'Mass' if isMass else ('SR' if 'SR' in name else ('SB' if 'SB' in name else ''))
    cut = 'reg==reg::'+cat if reg is not None else ''
    normRange = 'V_extended_reasonable_range' if isMass else 'X_reasonable_range'
    dataRange = 'LSBrange,HSBrange' if isMass and isData else normRange
    
    lumi = LUMI
    cmsLabel = 'Preliminary' if isData else 'Simulation Preliminary'
    pullRange = 5
    
    # ====== CONTROL PLOT ======
    c = TCanvas('c_'+name, 'Fitting '+name, 1000, 800)
    c.Divide(1, 2)
    c.cd(1)
    setTopPad(c.GetPad(1), RATIO)
    setBotPad(c.GetPad(2), RATIO)
    frame = variable.frame()
    setHistStyle(frame, 1.2)
    
    # Plot Data
    data, res = None, None
    if dataset is not None: data = dataset.plotOn(frame, RooFit.Cut(cut), RooFit.Binning(binsJmass if isMass else binsXmass), RooFit.DataError(RooAbsData.Poisson if isData else RooAbsData.SumW2), RooFit.MarkerStyle(21), RooFit.Range(dataRange), RooFit.Name('data_obs'))
    if data is not None and isData: fixData(data.getHist())
    
    # ---------- SIMPLE FIT ----------
    if isData:
        if isCategory:
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.DrawOption('F'), RooFit.LineColor(getColor(name, channel)), RooFit.FillColor(getColor(name, channel)), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('Vjet')) #RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), FIXME
            res = frame.pullHist()
            chi2 = frame.chiSquare()
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.Components('VV'+postfix+',Top'+postfix), RooFit.DrawOption('F'), RooFit.LineColor(798), RooFit.FillColor(798), RooFit.VLines(), RooFit.Name('Top'))
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.Components('VV'+postfix),  RooFit.DrawOption('F'), RooFit.LineColor(602), RooFit.FillColor(602), RooFit.VLines(), RooFit.Name('VV'))
            if type(fitRes) is list:
                for f in fitRes:
                    if f is not None: model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.VisualizeError(f, 1, False), RooFit.SumW2Error(True), RooFit.LineColor(0), RooFit.FillColor(1), RooFit.FillStyle(3005), RooFit.Name('Uncertainty'))
                    #model.paramOn(frame, RooFit.Label(model.GetTitle()), RooFit.Layout(0.5, 0.95, 0.94), RooFit.Format('NEAU')) #FIXME
            elif fitRes is not None: frame.addObject(fitRes, 'E3')
        else:
            if isMass:
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.LineColor(602), RooFit.DrawOption('F'), RooFit.FillColor(602), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('VV'))
                res = frame.pullHist()
                chi2 = frame.chiSquare()
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.Components('Vjet'+postfix+',Top'+postfix), RooFit.LineColor(798), RooFit.DrawOption('F'), RooFit.FillColor(798), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('Top'))
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.Components('Vjet'+postfix), RooFit.LineColor(getColor(name, channel)), RooFit.DrawOption('F'), RooFit.FillColor(getColor(name, channel)), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('Vjet'))
                #if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(921))#default
                if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(46))
                if type(fitRes) is list:
                    for f in fitRes:
                        if f is not None: model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.VisualizeError(f, 1, False), RooFit.SumW2Error(True), RooFit.LineColor(0), RooFit.FillColor(1), RooFit.FillStyle(3005), RooFit.Name('Uncertainty'))
                elif fitRes is not None: frame.addObject(fitRes, 'E3')
            else:
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.LineColor(getColor(name, channel)), RooFit.DrawOption('F'), RooFit.FillColor(getColor(name, channel)), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('Vjet'))
                res = frame.pullHist()
                chi2 = frame.chiSquare()
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.Components('VV'+postfix+',Top'+postfix), RooFit.LineColor(798), RooFit.DrawOption('F'), RooFit.FillColor(798), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('Top'))
                model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.Components('VV'+postfix), RooFit.LineColor(602), RooFit.DrawOption('F'), RooFit.FillColor(602), RooFit.FillStyle(1001), RooFit.VLines(), RooFit.Name('VV'))
                #if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(921))#default
                if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(46))
                if type(fitRes) is list:
                    for f in fitRes:
                        if f is not None: model.plotOn(frame, RooFit.Normalization(norm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.VisualizeError(f, 1, False), RooFit.SumW2Error(True), RooFit.LineColor(0), RooFit.FillColor(1), RooFit.FillStyle(3005), RooFit.Name('Uncertainty'))
                elif fitRes is not None: frame.addObject(fitRes, 'E3')
                if signal is not None: signal.plotOn(frame, RooFit.Normalization(snorm, RooAbsReal.NumEvent), RooFit.Range(normRange), RooFit.LineColor(629), RooFit.DrawOption('L'), RooFit.Name('Signal'))
                
    # Simple fit
    else:
        if isCategory:
            if type(fitRes) is list:
                for f in fitRes:
                    if f is not None: model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.VisualizeError(f, 1, False), RooFit.SumW2Error(True), RooFit.FillColor(1), RooFit.FillStyle(3005))
            elif fitRes is not None: frame.addObject(fitRes, 'E3')
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.LineColor(getColor(name, channel)))
            res = frame.pullHist()
            chi2 = frame.chiSquare()
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.LineColor(getColor(name, channel)), RooFit.LineStyle(2), RooFit.Components('baseTop'))
            model.plotOn(frame, RooFit.Slice(reg, cat), RooFit.ProjWData(RooArgSet(reg), dataset), RooFit.LineColor(getColor(name, channel)), RooFit.LineStyle(2), RooFit.Components('baseVV'))
        else:
            if type(fitRes) is list:
                for f in fitRes:
                    if f is not None: model.plotOn(frame, RooFit.VisualizeError(f, 1, False), RooFit.SumW2Error(True), RooFit.FillColor(1), RooFit.FillStyle(3005), RooFit.DrawOption('F'))
                    model.paramOn(frame, RooFit.Label(model.GetTitle()), RooFit.Layout(0.5, 0.95, 0.94), RooFit.Format('NEAU'))
            elif fitRes is not None: frame.addObject(fitRes, 'E3')
            model.plotOn(frame, RooFit.LineColor(getColor(name, channel)))
            res = frame.pullHist()
            chi2 = frame.chiSquare()
            model.plotOn(frame, RooFit.LineColor(getColor(name, channel)), RooFit.LineStyle(2), RooFit.Components('baseTop'))
            model.plotOn(frame, RooFit.LineColor(getColor(name, channel)), RooFit.LineStyle(2), RooFit.Components('baseVV'))
            #if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(921))#default
            if alt is not None: alt.plotOn(frame, RooFit.LineStyle(7), RooFit.LineColor(46))
    
    # Replot data
    if dataset is not None: data = dataset.plotOn(frame, RooFit.Cut(cut), RooFit.Binning(binsJmass if isMass else binsXmass), RooFit.DataError(RooAbsData.Poisson if isData else RooAbsData.SumW2), RooFit.Range(dataRange), RooFit.Name('data_obs'))
    if data is not None and isData: fixData(data.getHist())
    
    if not isMass and not isSignal: # Log scale
        frame.SetMaximum(frame.GetMaximum()*10)
        frame.SetMinimum(max(frame.SetMinimum(), 5.e-3 if isData else 1.e-4))
        c.GetPad(1).SetLogy()
    else:
        frame.GetYaxis().SetRangeUser(0, frame.GetMaximum())
        frame.SetMaximum(frame.GetMaximum()*1.25)
        frame.SetMinimum(0)
    if isSignal: frame.GetXaxis().SetRangeUser(0, 5000)
    frame.Draw()
    drawCMS(lumi, cmsLabel)
    drawAnalysis('VZinv')
    drawRegion(channel + ('' if isData and not isCategory else ('SR' if 'SR' in name else ('SB' if 'SB' in name else ''))), True)
    if isData and isMass:
        box = drawBox(HMIN, frame.GetMinimum(), HMAX, frame.GetMaximum()/1.30, '') #'(blind)'
        lineL = drawLine(LOWMAX, frame.GetMinimum(), LOWMAX, frame.GetMaximum()/1.30)
        lineSm = drawLine(SIGMIN, frame.GetMinimum(), SIGMIN, frame.GetMaximum()/1.30)
        lineSM = drawLine(SIGMAX, frame.GetMinimum(), SIGMAX, frame.GetMaximum()/1.30)
        lineHm = drawLine(HMIN, frame.GetMinimum(), HMIN, frame.GetMaximum()/1.30)
        lineHM = drawLine(HMAX, frame.GetMinimum(), HMAX, frame.GetMaximum()/1.30)
        lineU = drawLine(HIGMIN, frame.GetMinimum(), HIGMIN, frame.GetMaximum()/1.30)
        textL = drawText((LOWMIN+LOWMAX)/2, frame.GetMaximum()/1.35, 'LSB',2)
        textV = drawText((SIGMIN+SIGMAX)/2, frame.GetMaximum()/1.35, 'SR',2)
        textH = drawText((HMIN+HMAX)/2, frame.GetMaximum()/1.35, 'Higgs',2)
        textU = drawText(HIGMIN+20, frame.GetMaximum()/1.35, 'HSB',2)
    
    if isData:
        leg = TLegend(0.65, 0.6, 0.95, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0) #1001
        leg.SetFillColor(0)
        leg.AddEntry('data_obs', 'Data', 'PL')
        leg.AddEntry('Vjet', 'Z(#nu#nu) + jets, W(l#nu) + jets', 'F')
        leg.AddEntry('Top', 't#bar{t}, ST', 'F')
        leg.AddEntry('VV', 'VV', 'F')
        if type(fitRes) is list: leg.AddEntry('Uncertainty', 'Bkg. unc.', 'F')
        elif fitRes is not None: leg.AddEntry(fitRes, 'Bkg. unc.', 'F')
        if signal is not None: leg.AddEntry('Signal', signal.GetTitle(), 'L')
        leg.SetY1(0.9-leg.GetNRows()*0.06)
        leg.Draw()
        if signal is not None:
            latex = TLatex()
            latex.SetNDC()
            latex.SetTextSize(0.045)
            latex.SetTextFont(42)
            if SIGNALNAME == 'XZZ':
                latex.DrawLatex(0.725, leg.GetY1()-0.045, 'G_{Bulk} #rightarrow Z_{had}Z_{lep}')
            elif SIGNALNAME == 'XWZ': 
                latex.DrawLatex(0.725, leg.GetY1()-0.045, 'W\' #rightarrow W_{had}Z_{lep}') 
            else:
                latex.DrawLatex(0.725, leg.GetY1()-0.045, 'X #rightarrow V_{had}Z_{lep}') 
    
    c.cd(2)
    frame_res = variable.frame()
    setHistStyle(frame_res, 1.2)
    if dataset is not None: frame_res.addPlotable(res, 'P')
    setBotStyle(frame_res, RATIO, False)
    frame_res.GetYaxis().SetRangeUser(-pullRange, pullRange)
    frame_res.GetYaxis().SetTitle('Pulls')
    frame_res.Draw()
    
    if isData and isMass:
        #line_res = drawLine(frame_res.GetXaxis().GetXmin(), 0, frame_res.GetXaxis().GetXmax(), 0, 2)
        box_res = drawBox(HMIN, -pullRange, HMAX, pullRange)
        lineL_res = drawLine(LOWMAX, -pullRange, LOWMAX, pullRange)
        lineSm_res = drawLine(SIGMIN, -pullRange, SIGMIN, pullRange)
        lineSM_res = drawLine(SIGMAX, -pullRange, SIGMAX, pullRange)
        lineHm_res = drawLine(HMIN, -pullRange, HMIN, pullRange)
        lineHM_res = drawLine(HMAX, -pullRange, HMAX, pullRange)
        lineU_res = drawLine(HIGMIN, -pullRange, HIGMIN, pullRange)
    line_res = drawLine(frame_res.GetXaxis().GetXmin(), 0, frame_res.GetXaxis().GetXmax(), 0, 2)
    
    chilatex = TLatex()
    chilatex.SetNDC()
    chilatex.SetTextColor(getColor(name, channel))
    chilatex.SetTextFont(62)
    chilatex.SetTextSize(0.13)
    if(chi): chilatex.DrawLatex(0.6, 0.85, "Fit #chi^{2} = %.2f" % (chi2))

    c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.png')
    c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.pdf')
    if signal is not None:
        c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.root')
        c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.C')
    if name=='VjetMass':
        c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.root')
        c.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.C')        
    #if VERBOSE: raw_input('Press Enter to continue...')
    # ======   END PLOT   ======


def drawAlphaPlot(name, channel, variable, alpha, bkgSB, bkgSR, fitRes, alpha2=None, bkgSB2=None, bkgSR2=None, fitRes2=None):

    # ====== CONTROL PLOT ======
    c = TCanvas('c_'+name, 'Alpha function', 1000, 800)
    c.cd()
    frame = variable.frame()
    setHistStyle(frame, 1.1)
    #bkgSB.plotOn(frame, RooFit.LineColor(602)) #FIXME
    #bkgSR.plotOn(frame, RooFit.LineColor(2)) #FIXME
    alpha.plotOn(frame, RooFit.VisualizeError(fitRes, 2, False), RooFit.LineColor(400), RooFit.FillColor(400), RooFit.Name('2sigma'))
    alpha.plotOn(frame, RooFit.VisualizeError(fitRes, 1, False), RooFit.LineColor(416), RooFit.FillColor(416), RooFit.Name('1sigma'))
    alpha.plotOn(frame, RooFit.LineColor(1), RooFit.Name('alpha'))
    if ALTERNATIVE:
        #alpha2.plotOn(frame, RooFit.VisualizeError(fitRes2, 2, False), RooFit.DrawOption('L'), RooFit.LineColor(921), RooFit.LineStyle(8))
        #alpha2.plotOn(frame, RooFit.VisualizeError(fitRes2, 1, False), RooFit.DrawOption('L'), RooFit.LineColor(921), RooFit.LineStyle(7))
        alpha2.plotOn(frame, RooFit.LineColor(921), RooFit.LineStyle(7), RooFit.Name('alpha2'))
    frame.GetXaxis().SetRangeUser(XBINMIN, XBINMAX)
    frame.SetYTitle('')
    frame.Draw()
    #drawCMS(-1, 'Simulation')
    drawAnalysis('VZinv')
    drawRegion(channel, True)
    #drawCMS(LUMISILVER if 'XZhee' in channel or 'XZhmm' in channel or 'XZhll' in channel else LUMIGOLDEN, 'Preliminary')
    drawCMS(LUMI, 'Preliminary')
    
    leg = TLegend(0.5, 0.65, 0.95, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.AddEntry('alpha', alpha.GetTitle(), 'L')
    leg.AddEntry('1sigma', '#alpha function #pm 1#sigma', 'F')
    leg.AddEntry('2sigma', '#alpha function #pm 2#sigma', 'F')
    if ALTERNATIVE: leg.AddEntry('alpha2', alpha2.GetTitle(), 'L')
    leg.Draw()
    
    
    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaRatio.pdf')
    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaRatio.png')
    #c.SaveAs(PLOTDIR+'/'+channel+'/AlphaRatio.root')
#    frame.SetYTitle('')
#    #bkgSB.plotOn(frame, RooFit.LineColor(602))
#    #bkgSR.plotOn(frame, RooFit.LineColor(2))
#    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.pdf')
#    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.png')
#    #c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.root')
#    frame.SetMinimum(max(frame.SetMinimum(), 2.e-3))
#    c.GetPad(0).SetLogy()  
#    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.pdf')
#    c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.png')
#    #c.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.root')
    # ======   END PLOT   ======
    
    
    # ====== CONTROL PLOT ======
    c_alpha2 = TCanvas('c_alpha2', 'Alpha method', 1000, 800)
    c_alpha2.cd()
    frame_alpha2 = variable.frame()
    setHistStyle(frame, 1.1)
    #alpha.plotOn(frame_alpha2, RooFit.VisualizeError(fitRes, 2, False), RooFit.LineColor(400), RooFit.FillColor(400), RooFit.Name('2sigma'))
    #alpha.plotOn(frame_alpha2, RooFit.VisualizeError(fitRes, 1, False), RooFit.LineColor(416), RooFit.FillColor(416), RooFit.Name('1sigma'))
    alpha.plotOn(frame_alpha2, RooFit.LineColor(1), RooFit.Name('alpha'))
    if ALTERNATIVE:
        #alpha2.plotOn(frame_alpha2, RooFit.LineColor(921), RooFit.LineStyle(7), RooFit.Name('alpha2'))#default
        #bkgSB2.plotOn(frame_alpha2, RooFit.LineColor(602), RooFit.LineStyle(7), RooFit.Name('bkgSB2'))#default
        #bkgSR2.plotOn(frame_alpha2, RooFit.LineColor(2), RooFit.LineStyle(7), RooFit.Name('bkgSR2'))#default
        alpha2.plotOn(frame_alpha2, RooFit.LineColor(921), RooFit.LineStyle(7), RooFit.Name('alpha2'))
        bkgSB2.plotOn(frame_alpha2, RooFit.LineColor(66), RooFit.LineStyle(7), RooFit.Name('bkgSB2')) 
        bkgSR2.plotOn(frame_alpha2, RooFit.LineColor(8), RooFit.LineStyle(7), RooFit.Name('bkgSR2')) 
    bkgSB.plotOn(frame_alpha2, RooFit.LineColor(602), RooFit.Name('bkgSB')) 
    bkgSR.plotOn(frame_alpha2, RooFit.LineColor(2), RooFit.Name('bkgSR')) 
    frame_alpha2.GetXaxis().SetRangeUser(XBINMIN, XBINMAX)
    frame_alpha2.GetYaxis().SetRangeUser(0, frame_alpha2.GetMaximum())#L#
    frame_alpha2.Draw()
    frame_alpha2.SetYTitle('')
    drawAnalysis('VZinv')
    drawRegion(channel, True)
    #drawCMS(LUMISILVER if 'XZhee' in channel or 'XZhmm' in channel or 'XZhll' in channel else LUMIGOLDEN, 'Preliminary')
    drawCMS(LUMI, 'Preliminary')
    
    leg2 = TLegend(0.5, 0.55, 0.95, 0.9)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0) #1001
    leg2.SetFillColor(0)
    leg2.AddEntry('alpha', alpha.GetTitle(), 'L')
    leg2.AddEntry('bkgSB', 'bkg. fit in SB', 'L')
    leg2.AddEntry('bkgSR', 'bkg. pred. in SR', 'L')
    if ALTERNATIVE:
        leg2.AddEntry('alpha2', alpha2.GetTitle(), 'L')
        leg2.AddEntry('bkgSB2', 'alternative bkg. fit in SB', 'L')
        leg2.AddEntry('bkgSR2', 'alternative bkg. pred. in SR', 'L')
    leg2.Draw()
    
    c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.pdf')
    c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.png')
    #c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod.root')
    #frame_alpha2.SetMaximum(1.e+1)#default
    frame_alpha2.SetMaximum(2.e-2)#default
    frame_alpha2.SetMinimum(1.e-6)
    c_alpha2.GetPad(0).SetLogy()  
    c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.pdf')
    c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.png')
    #c_alpha2.SaveAs(PLOTDIR+'/'+channel+'/AlphaMethod_log.root')
    
    if VERBOSE: raw_input('Press Enter to continue...')
    # ======   END PLOT   ======
    


def drawMultiPlot(name, channel, variable, models):

    # ====== CONTROL PLOT ======
    c = TCanvas('c_'+name, 'Signal', 800, 600)
    c.cd()
    variable.setMin(0)
    variable.setMax(4000)
    frame = variable.frame()
    setHistStyle(frame, 1.1)
    for m, model in models.iteritems(): model.plotOn(frame, RooFit.Normalization(1, RooAbsReal.NumEvent), RooFit.LineColor(getColor('XZZ_M%d' % m, channel)), RooFit.LineStyle(1 if not getColor('XZZ_M%d' % m, channel)==1 else 3), RooFit.LineWidth(2), RooFit.DrawOption('L'), RooFit.Name(model.GetName()))
    frame.GetXaxis().SetRangeUser(0., 4000.)
    frame.SetYTitle('Signal normalization')
    frame.Draw()
    drawAnalysis('ZZ')
    drawRegion(channel, True)
    drawCMS(-1, 'Preliminary Simulation')
    
    leg = TLegend(0.65, 0.5, 0.99, 0.925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    #leg.SetNColumns(len(models) / 15)
    for m, model in sorted(models.iteritems()):
        if not getColor('XZZ_M%d' % m, channel)==1: leg.AddEntry(model.GetName(), model.GetTitle(), 'L')
    leg.Draw()
    c.SaveAs('plotsAlpha/'+channel+'/'+channel+'_'+name+'.pdf')
    c.SaveAs('plotsAlpha/'+channel+'/'+channel+'_'+name+'.png')
    
    if VERBOSE: raw_input('Press Enter to continue...')
    # ======   END PLOT   ======
    


def plotSys(name, channel, variable, model, sys):
    c_sys = TCanvas('c_sys_'+name, name, 800, 600)
    c_sys.cd()
    frame_sys = variable.frame()
    model.plotOn(frame_sys, RooFit.LineColor(1), RooFit.Name('Central'))
    sys.setVal(1)
    model.plotOn(frame_sys, RooFit.LineColor(634), RooFit.Name('Up'))
    sys.setVal(-1)
    model.plotOn(frame_sys, RooFit.LineColor(598), RooFit.Name('Down'))
    sys.setVal(0)
    frame_sys.GetXaxis().SetRangeUser(XBINMIN, XBINMAX)
    frame_sys.Draw()
    drawCMS(-1, 'Simulation')
    drawAnalysis('VZinv')
    drawRegion(channel)
    #drawCMS(LUMISILVER if 'XZhee' in channel or 'XZhmm' in channel or 'XZhll' in channel else LUMIGOLDEN, 'Preliminary')
    drawCMS(LUMI, 'Preliminary')
    frame_sys.SetMaximum(frame_sys.GetMaximum()*5)
    frame_sys.SetMinimum(1.e-4)
    c_sys.GetPad(0).SetLogy()
    leg = TLegend(0.65, 0.65, 0.95, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0) #1001
    leg.SetFillColor(0)
    leg.AddEntry('Up', '+1 #sigma', 'L')
    leg.AddEntry('Central', 'central', 'L')
    leg.AddEntry('Down', '-1 #sigma', 'L')
    leg.Draw()
    
    c_sys.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.pdf')
    c_sys.SaveAs(PLOTDIR+'/'+channel+'/'+name+'.png')




def likelihoodScan(model, dataset, par = [], chan = ''):

    nll = model.createNLL(dataset, RooFit.Range('X_reasonable_range'), RooFit.Strategy(2), RooFit.Minimizer('Minuit2'), RooFit.SumW2Error(False)) #RooFit.NumCPU(10)
    
    #gROOT.SetBatch(False)
    nv = (len(par)-1) / 2 + 1
    nh = len(par) / 2 + 1
    
    c_scan = TCanvas('c_scan', 'Likelihood scan', 600*nh, 600*nv)
    c_scan.Divide(nh, nv)
    frame = {}
    for i, p in enumerate(par):
        c_scan.cd(i+1)
        frame[i] = p.frame()
        nll.plotOn(frame[i], RooFit.ShiftToZero(), RooFit.PrintEvalErrors(-1), RooFit.EvalErrorValue(nll.getVal()+10))
        frame[i].GetXaxis().SetRangeUser(p.getMin(), p.getMax())
        frame[i].GetYaxis().SetRangeUser(0, 9)
        #c_scan.GetPad(i).SetLogx()
        #c_scan.GetPad(i).SetLogy()
        frame[i].Draw()
    c_scan.Print('plotsAlpha/'+chan+'/Scan_'+model.GetName()+'_'+dataset.GetName()+'.pdf')
    #raw_input('Press Enter to continue...')
    #if options.batch: gROOT.SetBatch(True)





##main function##
if __name__ == '__main__':
    if options.all:
        for c in channelList:
            p = multiprocessing.Process(target=alpha, args=(c,))
            jobs.append(p)
            p.start()
    else:
        if options.channel in channelList: alpha_norm(options.channel)
        else:
            print 'Channel not set or not recognized. Quitting...'
            exit()
