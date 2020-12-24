import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint

def main():
    
    # give the signal sample you want to use, old with less signal tracks or new with more signal tracks
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-o', action="store_true", dest="needOldSignal", default=False)
    args = parser.parse_args()
    
    if(args.needOldSignal):
        #### signal and background file
        signalAndBackgroundBX1 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX2 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX3 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX4 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")

        #### only signal file
        signalFromBkgFileBX1 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX2 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX3 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX4 = TFile("seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

        plotSuffix = "BkgEBeam_SignalHics3000nmOld_PerBX"
        latexName   = 'sig from e+laser (old)'
        
    else:
        signalAndBackgroundBX1 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX2 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX3 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX4 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        
        
        signalFromBkgFileBX1 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX2 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX3 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX4 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        
        plotSuffix = "BkgEBeam_SignalHics5000nmProvisional_PerBX"
        latexName   = 'sig from e+laser (updated)'
    
    
    
    hSeedMultiplicitySignalAndBackgroundBX1    = signalAndBackgroundBX1.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalAndBackgroundBX2    = signalAndBackgroundBX2.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalAndBackgroundBX3    = signalAndBackgroundBX3.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalAndBackgroundBX4    = signalAndBackgroundBX4.Get("hSeedMultiplicity")


    hSeedMultiplicitySignalAndBackgroundBX1.Add(hSeedMultiplicitySignalAndBackgroundBX2)
    hSeedMultiplicitySignalAndBackgroundBX1.Add(hSeedMultiplicitySignalAndBackgroundBX3)
    hSeedMultiplicitySignalAndBackgroundBX1.Add(hSeedMultiplicitySignalAndBackgroundBX4)


    

    #signal0KeVFile = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

    #signalMultiplicity = TFile("seedingInformation_hics_165gev_w0_3000nm_WIS_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root","READ")

    hSeedMultiplicitySignalBX1               = signalFromBkgFileBX1.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalBX2               = signalFromBkgFileBX2.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalBX3               = signalFromBkgFileBX3.Get("hSeedMultiplicity")
    hSeedMultiplicitySignalBX4               = signalFromBkgFileBX4.Get("hSeedMultiplicity")

    hSeedMultiplicitySignalBX1.Add(hSeedMultiplicitySignalBX2)
    hSeedMultiplicitySignalBX1.Add(hSeedMultiplicitySignalBX3)
    hSeedMultiplicitySignalBX1.Add(hSeedMultiplicitySignalBX4)

    hSignalMultiplicityBX1                   = signalFromBkgFileBX1.Get("hSignalMultiplicity")
    hSignalMultiplicityBX2                   = signalFromBkgFileBX2.Get("hSignalMultiplicity")
    hSignalMultiplicityBX3                   = signalFromBkgFileBX3.Get("hSignalMultiplicity")
    hSignalMultiplicityBX4                   = signalFromBkgFileBX4.Get("hSignalMultiplicity")

    hSignalMultiplicityBX1.Add(hSignalMultiplicityBX2)
    hSignalMultiplicityBX1.Add(hSignalMultiplicityBX3)
    hSignalMultiplicityBX1.Add(hSignalMultiplicityBX4)


    FirstTH1 = [hSignalMultiplicityBX1, hSeedMultiplicitySignalAndBackgroundBX1, hSeedMultiplicitySignalBX1]

    PlotColor = [kGray, 2, 4]
    LegendName = ['true signal','seeds from sig+bkg', 'seeds from sig']

    xAxisName   = "Multiplicity"
    yAxisName   = "BX"
    xrange1down = 0
    xrange1up   = 2000
    yrange1down = 0.001
    yrange1up   = 4
    yline1low   = 1
    yline1up    = 1

    drawline    = False
    logy        = False
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

    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "seedMultiplicityDistributionSeedAndSignal_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)



if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
