import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint


def main():
    directory = "OldSeedPlots"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # give the signal sample you want to use, old with less signal tracks or new with more signal tracks
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-o', action="store_true", dest="needOldSignal", default=False)
    args = parser.parse_args()
    inDirectory = "NoRemovalOfDuplicateTracks"
    
    if(args.needOldSignal):
        #### signal and background file
        signalAndBackgroundBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")

        #### only signal file
        signalFromBkgFileBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

        yrange1up   = 30
        plotSuffix = "BkgEBeam0.826_SignalHics3000nmOld_PerBX_NoRemoveDuplicateTracks"
        latexName  = 'sig from e+laser (old)'
        
        
    else:
        #### signal and background file
        signalAndBackgroundBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")


        #### only signal file
        signalFromBkgFileBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

        yrange1up   = 450
        plotSuffix = "BkgEBeam0.826_SignalHics5000nmProvisional_PerBX_NoRemoveDuplicateTracks"
        latexName  = 'sig from e+laser (updated)'



    hSeedEnergy_signalAndBackgroundBX1  = signalAndBackgroundBX1.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX2  = signalAndBackgroundBX2.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX3  = signalAndBackgroundBX3.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX4  = signalAndBackgroundBX4.Get("hSeedEnergy")

    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX2)
    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX3)
    hSeedEnergy_signalAndBackgroundBX1.Add(hSeedEnergy_signalAndBackgroundBX4)
    
    hSeedEnergy_signalAndBackgroundBX1.Rebin(4)


    hSignalEnergyBX1                     = signalFromBkgFileBX1.Get("hSigEnergy")
    hSignalEnergyBX2                     = signalFromBkgFileBX2.Get("hSigEnergy")
    hSignalEnergyBX3                     = signalFromBkgFileBX3.Get("hSigEnergy")
    hSignalEnergyBX4                     = signalFromBkgFileBX4.Get("hSigEnergy")

    hSignalEnergyBX1.Add(hSignalEnergyBX2)
    hSignalEnergyBX1.Add(hSignalEnergyBX3)
    hSignalEnergyBX1.Add(hSignalEnergyBX4)
    
    hSignalEnergyBX1.Rebin(4)


    hSeedEnergyMatchedBX1                = signalFromBkgFileBX1.Get("hSeedEnergy")
    hSeedEnergyMatchedBX2                = signalFromBkgFileBX2.Get("hSeedEnergy")
    hSeedEnergyMatchedBX3                = signalFromBkgFileBX3.Get("hSeedEnergy")
    hSeedEnergyMatchedBX4                = signalFromBkgFileBX4.Get("hSeedEnergy")


    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX2)
    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX3)
    hSeedEnergyMatchedBX1.Add(hSeedEnergyMatchedBX4)
    
    hSeedEnergyMatchedBX1.Rebin(4)



    FirstTH1 = [hSignalEnergyBX1, hSeedEnergy_signalAndBackgroundBX1, hSeedEnergyMatchedBX1]

    PlotColor = [kGray, 2, 4]
    LegendName = ['sig, from Geant4', 'seed (sig+bkg)', 'seed (sig matched)']

    xAxisName   = "E [GeV]"
    yAxisName   = "Particles/4.0BX"
    xrange1down = 0
    xrange1up   = 17
    yrange1down = 0.001
    
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

    h2 = hSeedEnergy_signalAndBackgroundBX1.Clone("h2")
    h2.Reset()
    h2.GetYaxis().SetTitle("#frac{seed energy (reco, sig+bkg)}{energy (sig, from Geant4)}")

    h3 = hSeedEnergy_signalAndBackgroundBX1.Clone("h3")
    h3.Reset()
    h3.GetYaxis().SetTitle("#frac{seed energy (reco, sig matched)}{energy (sig, from Geant4)}")



    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignal_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)

    drawline    = True
    DrawHistsRatioTwo(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignal_WithRatio_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3)


if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
