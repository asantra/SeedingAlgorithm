import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint


def main():

    directory = "FitPlots3Hits"
    if not os.path.exists(directory):
        os.makedirs(directory)
    # give the signal sample you want to use, old with less signal tracks or new with more signal tracks
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-o', action="store_true", dest="needOldSignal", default=False)
    args = parser.parse_args()
    inDirectory = "."
    if(args.needOldSignal):
        #### signal and background file
        signalAndBackgroundBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")

        #### only signal file
        signalFromBkgFileBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        
        
        
        
        
        
        #### signal and background file
        signalAndBackgroundBX1WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX2WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX3WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX4WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")

        #### only signal file
        signalFromBkgFileBX1WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX2WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX3WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX4WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics3000nmOld_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        
        
        
        

        yrange1up   = 20
        plotSuffix = "BkgEBeam0.826_SignalHics3000nmOld_PerBX_SVDFit"
        latexName2  = 'Low signal multiplicity'
        
        
    else:
        #### signal and background file
        signalAndBackgroundBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalAndBackgroundBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithoutFit3HitsTracks.root", "READ")


        #### only signal file
        signalFromBkgFileBX1 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX2 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX3 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")
        signalFromBkgFileBX4 = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithoutFit3HitsTracks.root", "READ")


        #### signal and background file
        signalAndBackgroundBX1WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX2WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX3WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")
        signalAndBackgroundBX4WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_WithFit3HitsTracks.root", "READ")


        #### only signal file
        signalFromBkgFileBX1WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX2WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX3WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")
        signalFromBkgFileBX4WithFit = TFile(inDirectory+"/seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_WithFit3HitsTracks.root", "READ")

        
        
        
        
        yrange1up   = 300
        plotSuffix = "BkgEBeam0.826_SignalHics5000nmProvisional_PerBX_SVDFit"
        latexName2  = 'High signal multiplicity'


    ### without fit
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



    ##### with Fit
    
    hSeedEnergy_signalAndBackgroundBX1WithFit  = signalAndBackgroundBX1WithFit.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX2WithFit  = signalAndBackgroundBX2WithFit.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX3WithFit  = signalAndBackgroundBX3WithFit.Get("hSeedEnergy")
    hSeedEnergy_signalAndBackgroundBX4WithFit  = signalAndBackgroundBX4WithFit.Get("hSeedEnergy")

    hSeedEnergy_signalAndBackgroundBX1WithFit.Add(hSeedEnergy_signalAndBackgroundBX2WithFit)
    hSeedEnergy_signalAndBackgroundBX1WithFit.Add(hSeedEnergy_signalAndBackgroundBX3WithFit)
    hSeedEnergy_signalAndBackgroundBX1WithFit.Add(hSeedEnergy_signalAndBackgroundBX4WithFit)
    
    hSeedEnergy_signalAndBackgroundBX1WithFit.Rebin(4)


    hSignalEnergyBX1WithFit                     = signalFromBkgFileBX1WithFit.Get("hSigEnergy")
    hSignalEnergyBX2WithFit                     = signalFromBkgFileBX2WithFit.Get("hSigEnergy")
    hSignalEnergyBX3WithFit                     = signalFromBkgFileBX3WithFit.Get("hSigEnergy")
    hSignalEnergyBX4WithFit                     = signalFromBkgFileBX4WithFit.Get("hSigEnergy")

    hSignalEnergyBX1WithFit.Add(hSignalEnergyBX2WithFit)
    hSignalEnergyBX1WithFit.Add(hSignalEnergyBX3WithFit)
    hSignalEnergyBX1WithFit.Add(hSignalEnergyBX4WithFit)
    
    hSignalEnergyBX1WithFit.Rebin(4)


    hSeedEnergyMatchedBX1WithFit                = signalFromBkgFileBX1WithFit.Get("hSeedEnergy")
    hSeedEnergyMatchedBX2WithFit                = signalFromBkgFileBX2WithFit.Get("hSeedEnergy")
    hSeedEnergyMatchedBX3WithFit                = signalFromBkgFileBX3WithFit.Get("hSeedEnergy")
    hSeedEnergyMatchedBX4WithFit                = signalFromBkgFileBX4WithFit.Get("hSeedEnergy")


    hSeedEnergyMatchedBX1WithFit.Add(hSeedEnergyMatchedBX2WithFit)
    hSeedEnergyMatchedBX1WithFit.Add(hSeedEnergyMatchedBX3WithFit)
    hSeedEnergyMatchedBX1WithFit.Add(hSeedEnergyMatchedBX4WithFit)
    
    hSeedEnergyMatchedBX1WithFit.Rebin(4)

    
    
    
    
    
    
    
    
    


    FirstTH1 = [hSignalEnergyBX1, hSeedEnergy_signalAndBackgroundBX1, hSeedEnergyMatchedBX1, hSeedEnergy_signalAndBackgroundBX1WithFit, hSeedEnergyMatchedBX1WithFit]

    PlotColor = [kGray, 2, 4, 2, 4]
    LegendName = ['sig', 'seed without fit', 'seed without fit (sig matched)', 'seed with fit', 'seed with fit (sig matched)']

    xAxisName   = "E [GeV]"
    yAxisName   = "Tracks/BX"
    
    #### put the x axis label and y axis label
    for i in xrange(0,len(FirstTH1)):
        FirstTH1[i].Scale(1./4.0)
        FirstTH1[i].GetYaxis().SetTitle(yAxisName)
        FirstTH1[i].GetXaxis().SetTitle(xAxisName)

    
    xrange1down = 0
    xrange1up   = 17
    yrange1down = 0.001
    
    yline1low   = 1
    yline1up    = 1

    drawline    = False
    logy        = False
    
    latexName   = "Data equivalent to 4 BXs"
    latexName3  = "Requiring 3/4 hit tracks"
    latexName4  = "sig from e+laser"
    latexName5  = "bkg from beam-only"

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
    
    h4 = hSeedEnergy_signalAndBackgroundBX1.Clone("h4")
    h4.Reset()
    h4.GetYaxis().SetTitle("#frac{seed energy (reco, sig+bkg)}{energy (sig, from Geant4)}")

    h5 = hSeedEnergy_signalAndBackgroundBX1.Clone("h5")
    h5.Reset()
    h5.GetYaxis().SetTitle("#frac{seed energy (reco, sig matched)}{energy (sig, from Geant4)}")



    drawline    = True
    DrawHistsRatioFour(FirstTH1, LegendName, PlotColor, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/energyDistributionSeedAndSignal_3HitFit_WithRatio_"+plotSuffix, h2, h3, yline1low, yline1up, drawline, logy, False, False, latexName, latexName2, latexName3, latexName4, latexName5, h4, h5)


if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
