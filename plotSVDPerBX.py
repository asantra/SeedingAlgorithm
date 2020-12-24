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

        yrange1up   = 30
        plotSuffix = "BkgEBeam0.826_SignalHics3000nmOld_PerBX"
        latexName  = 'sig from e+laser (old)'
        
        
    else:
        #### signal and background file
        signalAndBackgroundBX1 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX2 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX3 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")
        signalAndBackgroundBX4 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide.root", "READ")


        #### only signal file
        signalFromBkgFileBX1 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX2 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX2_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX3 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX3_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")
        signalFromBkgFileBX4 = TFile("seedingInformation_BkgEBeam_SignalHics5000nmProvisional_BX4_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide.root", "READ")

        yrange1up   = 150
        plotSuffix = "BkgEBeam0.826_SignalHics5000nmProvisional_PerBX"
        latexName  = 'sig from e+laser (updated)'

    
    hSVD1_SignalAndBackgroundBX1        = signalAndBackgroundBX1.Get("hSVDValues1")
    hSVD1_SignalAndBackgroundBX2        = signalAndBackgroundBX2.Get("hSVDValues1")
    hSVD1_SignalAndBackgroundBX3        = signalAndBackgroundBX3.Get("hSVDValues1")
    hSVD1_SignalAndBackgroundBX4        = signalAndBackgroundBX4.Get("hSVDValues1")
    
    hSVD1_SignalAndBackgroundBX1.Add(hSVD1_SignalAndBackgroundBX2)
    hSVD1_SignalAndBackgroundBX1.Add(hSVD1_SignalAndBackgroundBX3)
    hSVD1_SignalAndBackgroundBX1.Add(hSVD1_SignalAndBackgroundBX4)
    
    
    hSVD2_SignalAndBackgroundBX1        = signalAndBackgroundBX1.Get("hSVDValues2")
    hSVD2_SignalAndBackgroundBX2        = signalAndBackgroundBX2.Get("hSVDValues2")
    hSVD2_SignalAndBackgroundBX3        = signalAndBackgroundBX3.Get("hSVDValues2")
    hSVD2_SignalAndBackgroundBX4        = signalAndBackgroundBX4.Get("hSVDValues2")
    
    hSVD2_SignalAndBackgroundBX1.Add(hSVD2_SignalAndBackgroundBX2)
    hSVD2_SignalAndBackgroundBX1.Add(hSVD2_SignalAndBackgroundBX3)
    hSVD2_SignalAndBackgroundBX1.Add(hSVD2_SignalAndBackgroundBX4)
    
    
    
    hSVD3_SignalAndBackgroundBX1        = signalAndBackgroundBX1.Get("hSVDValues3")
    hSVD3_SignalAndBackgroundBX2        = signalAndBackgroundBX2.Get("hSVDValues3")
    hSVD3_SignalAndBackgroundBX3        = signalAndBackgroundBX3.Get("hSVDValues3")
    hSVD3_SignalAndBackgroundBX4        = signalAndBackgroundBX4.Get("hSVDValues3")
    
    hSVD3_SignalAndBackgroundBX1.Add(hSVD3_SignalAndBackgroundBX2)
    hSVD3_SignalAndBackgroundBX1.Add(hSVD3_SignalAndBackgroundBX3)
    hSVD3_SignalAndBackgroundBX1.Add(hSVD3_SignalAndBackgroundBX4)





    hSVD1_signalFromBkgFileBX1        = signalFromBkgFileBX1.Get("hSVDValues1")
    hSVD1_signalFromBkgFileBX2        = signalFromBkgFileBX2.Get("hSVDValues1")
    hSVD1_signalFromBkgFileBX3        = signalFromBkgFileBX3.Get("hSVDValues1")
    hSVD1_signalFromBkgFileBX4        = signalFromBkgFileBX4.Get("hSVDValues1")
    
    hSVD1_signalFromBkgFileBX1.Add(hSVD1_signalFromBkgFileBX2)
    hSVD1_signalFromBkgFileBX1.Add(hSVD1_signalFromBkgFileBX3)
    hSVD1_signalFromBkgFileBX1.Add(hSVD1_signalFromBkgFileBX4)
    
    
    hSVD2_signalFromBkgFileBX1        = signalFromBkgFileBX1.Get("hSVDValues2")
    hSVD2_signalFromBkgFileBX2        = signalFromBkgFileBX2.Get("hSVDValues2")
    hSVD2_signalFromBkgFileBX3        = signalFromBkgFileBX3.Get("hSVDValues2")
    hSVD2_signalFromBkgFileBX4        = signalFromBkgFileBX4.Get("hSVDValues2")
    
    hSVD2_signalFromBkgFileBX1.Add(hSVD2_signalFromBkgFileBX2)
    hSVD2_signalFromBkgFileBX1.Add(hSVD2_signalFromBkgFileBX3)
    hSVD2_signalFromBkgFileBX1.Add(hSVD2_signalFromBkgFileBX4)
    
    
    
    hSVD3_signalFromBkgFileBX1        = signalFromBkgFileBX1.Get("hSVDValues3")
    hSVD3_signalFromBkgFileBX2        = signalFromBkgFileBX2.Get("hSVDValues3")
    hSVD3_signalFromBkgFileBX3        = signalFromBkgFileBX3.Get("hSVDValues3")
    hSVD3_signalFromBkgFileBX4        = signalFromBkgFileBX4.Get("hSVDValues3")
    
    hSVD3_signalFromBkgFileBX1.Add(hSVD3_signalFromBkgFileBX2)
    hSVD3_signalFromBkgFileBX1.Add(hSVD3_signalFromBkgFileBX3)
    hSVD3_signalFromBkgFileBX1.Add(hSVD3_signalFromBkgFileBX4)




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
    
    
    
    FirstTH1 = [hSVD1_SignalAndBackgroundBX1, hSVD1_signalFromBkgFileBX1]

    PlotColor = [kGray, 2]
    LegendName = ['SVD[0] sig+bkg', 'SVD[0] bkg']

    xAxisName   = "SVD[0]"
    yAxisName   = "Entries"
    xrange1down = 150
    xrange1up   = 250
    yrange1down = 0.001

    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "SVD0_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)
    
    
    
    xAxisName   = "SVD[1]"
    yAxisName   = "Entries"
    xrange1down = 0
    xrange1up   = 0.6
    yrange1down = 0.001

    FirstTH1 = [hSVD2_SignalAndBackgroundBX1, hSVD2_signalFromBkgFileBX1]
    LegendName = ['SVD[1] sig+bkg', 'SVD[1] bkg']
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "SVD2_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)
    
    
    
    xAxisName   = "SVD[2]"
    yAxisName   = "Entries"
    xrange1down = 0
    xrange1up   = 0.6
    yrange1down = 0.001

    FirstTH1 = [hSVD3_SignalAndBackgroundBX1, hSVD3_signalFromBkgFileBX1]
    LegendName = ['SVD[2] sig+bkg', 'SVD[2] bkg']
    
    DrawHists(FirstTH1, LegendName, PlotColor,xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, "SVD3_"+plotSuffix, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx)
    


if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
