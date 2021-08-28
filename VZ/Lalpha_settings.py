
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
from ROOT import PdfDiagonalizer, RooAlphaExp, RooErfExpPdf, Roo2ExpPdf, RooAlpha42ExpPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooErfExpTailPdf, RooAlpha4ExpTailPdf, RooExpTail2Pdf, RooAlpha4ExpTail2Pdf, RooErfExpNPdf, RooExpN2Pdf, RooAlpha4ErfExpNPdf, RooAlpha4ExpN2Pdf, RooAlpha4ErfExpTailPdf, RooAlpha4ErfExpNPdf, RooAlpha, RooDoubleCrystalBall

from Analysis.ALPHA.LdrawUtils import *
from Analysis.ALPHA.variables import *
from Analysis.ALPHA.selectionsForAlpha import *
from Analysis.ALPHA.samples import samples

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
parser.add_option("-B", "--blind", action="store_true", default=False, dest="blind")
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

ALTERNATIVE = options.different
EXTRAPOLATE = options.extrapolate
SCAN        = options.scan
SIGNALNAME  = "XZZInv"
NTUPLEDIR   = '/lustre/cmswork/lbenato/VZAnalyses/CMSSW_8_0_25/src/Analysis/ALPHA/v3/Skim/'
#LUMI        = 12900 # in pb-1
LUMI        = 35867#36814 # in pb-1
PLOTDIR     = 'plotsAlpha' if not EXTRAPOLATE else 'plotsAlphaExt'
RATIO       = 4
SHOWERR     = True
BLIND       = False if (EXTRAPOLATE or options.unblind) else True
VERBOSE     = options.verbose
LISA        = True
CHAN        = options.channel


channelList = ['XVZnnlp', 'XVZnnhp']#['XVZmmlp', 'XVZmmhp', 'XVZeelp', 'XVZeehp']

signName  = 'XZZInv' # 'XZZ'

genPoints = [600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000 ]
#massPoints = [x for x in range(600, 4500+1, 100)] #if not HVTMODEL else genPoints
massPoints = [x for x in range(600, 3000+1, 100)] #if not HVTMODEL else genPoints
#massPoints = [600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000]


#gROOT.ProcessLine('.L %s/src/Analysis/ALPHA/plugins/Objects.h+' % os.environ['CMSSW_BASE'])
#from ROOT import LeptonType, JetType, FatJetType, MEtType, MEtFullType, CandidateType, LorentzType

LOWMIN = 30.
LOWMAX = 50. if EXTRAPOLATE else 65. 

SIGMIN = 50. if EXTRAPOLATE else 65.
SIGMAX = 65. if EXTRAPOLATE else 105.

HMIN = 105.
HMAX = 135.

HIGMIN = 135.
HIGMAX = 300.

VBINMIN= 30.
VBINMAX= 300.
VBINS  = 54

XBINMIN= 750.#550.
XBINMAX= 3150.
XBINS_PLOT  = 25
XBINS  = 1000

# Binning
binsJmass = RooBinning(VBINS, VBINMIN, VBINMAX)
binsXmass = RooBinning(XBINS_PLOT, XBINMIN, XBINMAX)
binsSignal = RooBinning(XBINS_PLOT*5, XBINMIN, XBINMAX)
binsSignal.addUniform(XBINS_PLOT*5, XBINMIN, XBINMAX)

# WSF[higP]76X with v2 JEC
#WSF = [1.154, 0.960, ]
#WSFErr = [0.305, 0.090, ]
# WSF[higP]80X with ICHEP JEC
WSF = [0.88, 1.03, ]
WSFErr = [0.24, 0.078,]
