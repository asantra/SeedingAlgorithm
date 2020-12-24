import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint



signalAndBackground0KeVFile = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")

signalFromBkgFile = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

signal0KeVFile = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root","READ")


#hSeedMultiplicity_signalAndBackground50KeV = signalAndBackground50KeVFile.Get("hSeedMultiplicity")
hSeedEnergy_signalAndBackground0KeV  = signalAndBackground0KeVFile.Get("hSeedEnergy")
#hSeedMultiplicity_signal0KeV               = signal0KeVFile.Get("hSeedMultiplicity")
hSignalEnergy                        = signal0KeVFile.Get("hSigEnergy")
hSeedEnergy_signal0KeV               = signalFromBkgFile.Get("hSeedEnergy")

plotSuffix = "BkgEBeam_SignalHics5000nmProvisional"
FirstTH1 = [hSignalEnergy, hSeedEnergy_signalAndBackground0KeV, hSeedEnergy_signal0KeV]

PlotColor = [kGray, 2, 4]
LegendName = ['energy (sig, from Geant4)', 'seed energy (reco, sig+bkg)', 'seed energy (reco, sig matched)']

xAxisName   = "E [GeV]"
yAxisName   = "Particles/4.0BX"
xrange1down = 0
xrange1up   = 20
yrange1down = 0.001
yrange1up   = 1400
yline1low   = 1
yline1up    = 1

drawline    = False
logy        = False
latexName   = 'sig from e+laser (updated)'
latexName2  = 'bkg from e-beam only'
latexName3  = ''
leftLegend  = False
doAtlas     = False
doLumi      = False
noRatio     = False
do80        = False
do59        = False
drawPattern = ""
logz        = False
logx        = False

h2 = hSeedEnergy_signal0KeV.Clone("h2")
h2.Reset()
h2.GetYaxis().SetTitle("#frac{seed energy (reco, sig+bkg)}{energy (sig, from Geant4)}")

h3 = hSeedEnergy_signal0KeV.Clone("h3")
h3.Reset()
h3.GetYaxis().SetTitle("#frac{seed energy (reco, sig matched)}{energy (sig, from Geant4)}")



DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "energyDistributionSeedAndSignal_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)

drawline    = True
DrawHistsRatioTwo(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, "energyDistributionSeedAndSignal_WithRatio_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3)
