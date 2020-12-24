import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint


signalAndBackground50KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")

signal0KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

signalMultiplicity = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root","READ")


hSeedMultiplicity_signalAndBackground50KeV = signalAndBackground50KeVFile.Get("hSeedMultiplicity")
hSeedMultiplicity_signal0KeV               = signal0KeVFile.Get("hSeedMultiplicity")
hSignalMultiplicity                        = signalMultiplicity.Get("hSignalMultiplicity")



FirstTH1 = [hSignalMultiplicity, hSeedMultiplicity_signalAndBackground50KeV, hSeedMultiplicity_signal0KeV]

PlotColor = [kGray, 2, 4]
LegendName = ['true signal','seeds from sig+bkg', 'seeds from sig']

xAxisName   = "Seed multiplicity"
yAxisName   = "BX"
xrange1down = 0
xrange1up   = 50
yrange1down = 0.001
yrange1up   = 1200
yline1low   = 1
yline1up    = 1

drawline    = False
logy        = False
latexName   = 'e+laser'
latexName2  = 'hics'
latexName3  = 'w0 3000 nm'
leftLegend  = False
doAtlas     = False
doLumi      = False
noRatio     = False
do80        = False
do59        = False
drawPattern = ""
logz        = False
logx        = False

DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "seedMultiplicityVariableECut", yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)
