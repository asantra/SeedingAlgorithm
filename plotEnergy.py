import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint


#signalAndBackground50KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfo_SignalAndBackground_50KeVCut.root", "READ")

signalAndBackground0KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_SignalAndBackground.root", "READ")

#signal0KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfo_Signal_0KeVCut.root", "READ")

signalMultiplicity = TFile("seedingInformationSignalOnly_hics_165gev_w0_3000nm_WIS_trackInfoClean.root","READ")


#hSeedMultiplicity_signalAndBackground50KeV = signalAndBackground50KeVFile.Get("hSeedMultiplicity")
hSeedMultiplicity_signalAndBackground0KeV  = signalAndBackground0KeVFile.Get("hSeedEnergy")
#hSeedMultiplicity_signal0KeV               = signal0KeVFile.Get("hSeedMultiplicity")
hSignalMultiplicity                        = signalMultiplicity.Get("hSigEnergy")


FirstTH1 = [hSignalMultiplicity, hSeedMultiplicity_signalAndBackground0KeV]

PlotColor = [kGray, 2, 4]
LegendName = ['energy from sig', 'seed energy from sig+bkg']

xAxisName   = "Energy"
yAxisName   = "entries"
xrange1down = 0
xrange1up   = 50
yrange1down = 0.001
yrange1up   = 8000
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

DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "seedMultiplicityNewEnergy", yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)
