import os, sys, glob, time
from ROOT import *
import argparse
from copy import copy, deepcopy
sys.path.insert(0, '/Users/arkasantra/arka/include')
from Functions import *
import pprint




### the seed pY width
PseedMin = -0.005 # GeV
PseedMax = 0.005 # GeV

def getCutValue(name):
    cutValue = {}
    cutValue['hNminus1d_High_Multiplicity_ELt4_TrackInclusive'] = {'d':0.005, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_High_Multiplicity_EGt4_TrackInclusive'] = {'d':0.003, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Medium_Multiplicity_ELt4_TightTrack'] = {'d':0.006, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Medium_Multiplicity_EGt4_TightTrack'] = {'d':0.004, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Medium_Multiplicity_ELt4_LooseTrack'] = {'d': 0.005, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Medium_Multiplicity_EGt4_LooseTrack'] = {'d':0.003, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Low_Multiplicity_ELt4_TightTrack'] = {'d':0.006, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Low_Multiplicity_EGt4_TightTrack'] = {'d':0.004, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Low_Multiplicity_ELt4_LooseTrack'] = {'d': 0.005, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1d_Low_Multiplicity_EGt4_LooseTrack'] = {'d': 0.005, 'xlow':0.0, 'xhigh':0.1}
    cutValue['hNminus1Py_High_Medium_Multiplicity_EInclusive_TrackInclusive'] = {'d':PseedMax, 'dmin': -PseedMax,'xlow':-0.02, 'xhigh':0.02}
    cutValue['hNminus1Py_Low_Multiplicity_ELt4_TightTracks'] = {'d':PseedMax*2, 'dmin': -PseedMax*2,'xlow':-0.02, 'xhigh':0.02}
    cutValue['hNminus1Py_Low_Multiplicity_ELt4_LooseTracks'] = {'d':PseedMax*1.9, 'dmin':-PseedMax*1.9,'xlow':-0.02, 'xhigh':0.02}
    cutValue['hNminus1Py_Low_Multiplicity_EGt4_TightTracks'] = {'d':PseedMax*7, 'dmin':-PseedMax*7,'xlow':-0.06, 'xhigh':0.06}
    cutValue['hNminus1Py_Low_Multiplicity_EGt4_LooseTracks'] = {'d':PseedMax*6.2, 'dmin':-PseedMax*6.2,'xlow':-0.07, 'xhigh':0.07}
    cutValue['hNminus1FitPar1_High_Multiplicity_ELt4_TrackInclusive'] = {'d':0.08,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_High_Multiplicity_ELt4_TrackInclusive'] = {'d':0.02,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_High_Multiplicity_EGt4_TrackInclusive'] = {'d':0.05,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_High_Multiplicity_EGt4_TrackInclusive'] = {'d':0.01,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Medium_Multiplicity_ELt4_TrackInclusive'] = {'d':0.06,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Medium_Multiplicity_ELt4_TrackInclusive'] = {'d':0.05,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Medium_Multiplicity_EGt4_TrackInclusive'] = {'d':0.06,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Medium_Multiplicity_EGt4_TrackInclusive'] = {'d':0.06,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Low_Multiplicity_ELt4_TightTrack'] = {'d':0.1,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Low_Multiplicity_ELt4_TightTrack'] = {'d':0.1,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Low_Multiplicity_EGt4_TightTrack'] = {'d':0.1,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Low_Multiplicity_EGt4_TightTrack'] = {'d':0.1,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Low_Multiplicity_ELt4_LooseTrack'] = {'d':0.2,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Low_Multiplicity_ELt4_LooseTrack'] = {'d':0.05,'xlow':0.0, 'xhigh':0.2}
    cutValue['hNminus1FitPar1_Low_Multiplicity_EGt4_LooseTrack'] = {'d':0.065,'xlow':0.0, 'xhigh':0.3}
    cutValue['hNminus1FitPar2_Low_Multiplicity_EGt4_LooseTrack'] = {'d':0.065,'xlow':0.0, 'xhigh':0.2}
    
    return cutValue[name]


def getLatexNames(name):
    name1 = " "
    name2 = " "
    name3 = " "
    name4 = " "
    name5 = " "
    
    if('High_Multiplicity' in name):
        name1 = "high track multiplicity"
    elif('Medium_Multiplicity' in name):
        name1 = "medium track multiplicity"
    elif('Low_Multiplicity' in name):
        name1 = "low track multiplicity"
        
    if("EGt4" in name):
        name2 = "E > 4 GeV"
    else:
        name2 = "E < 4 GeV"
        
    if("TightTrack" in name):
        name3 = "tight tracks"
    elif("LooseTrack" in name):
        name3 = "loose tracks"
    else:
        name3 = "tight+loose tracks"
        
    return name1, name2, name3, name4, name5
    

def main():
    gROOT.SetBatch()
    directory = "CutFlowPlots"  ## _LooseSVDanddCuts_December29_2020
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    inDirectory = "."
    
    
    #### background file
    ePlusLaserBkgDistanceD            = TFile(inDirectory+"/seedingInformation_EBeamOnlyNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_distancedNminus1.root", "READ")
    ePlusLaserBkgTrackPy              = TFile(inDirectory+"/seedingInformation_EBeamOnlyNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_trackPyNminus1.root", "READ")
    ePlusLaserBkgFitParameter         = TFile(inDirectory+"/seedingInformation_EBeamOnlyNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_SVDParNminus1.root", "READ")
    
    ### signal
    ePlusLaserSigDistanceD            = TFile(inDirectory+"/seedingInformation_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_distancedNminus1.root", "READ")
    ePlusLaserSigTrackPy              = TFile(inDirectory+"/seedingInformation_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_trackPyNminus1.root", "READ")
    ePlusLaserSigFitParameter         = TFile(inDirectory+"/seedingInformation_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_SVDParNminus1.root", "READ")
    
    
    ### g+laser electron
    gPlusLaserBkgElectronDistanceD    = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_ElectronSide_distancedNminus1.root", "READ")
    gPlusLaserBkgElectronTrackPy      = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_ElectronSide_trackPyNminus1.root", "READ")
    gPlusLaserBkgElectronFitParameter = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_ElectronSide_SVDParNminus1.root", "READ")
    
    ### g+laser positron
    gPlusLaserBkgPositronDistanceD    = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_distancedNminus1.root", "READ")
    gPlusLaserBkgPositronTrackPy      = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_trackPyNminus1.root", "READ")
    gPlusLaserBkgPositronFitParameter = TFile(inDirectory+"/seedingInformation_gPlusLaserBkgNewSamples_AllBX_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_SVDParNminus1.root", "READ")
    
    #### get the distanceD cutflow first
    ### d parameter N-1 plot
    hCutFlow_ePlusLaserBkg_DistanceD = {}
    for keys in ePlusLaserBkgDistanceD.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1d' not in checkKeys: continue
        hist = ePlusLaserBkgDistanceD.Get(keys.GetName())
        hist.Scale(1./160.0)
        hCutFlow_ePlusLaserBkg_DistanceD[keys.GetName()] = hist
        
    hCutFlow_ePlusLaserSig_DistanceD = {}
    for keys in ePlusLaserSigDistanceD.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1d' not in checkKeys: continue
        hist = ePlusLaserSigDistanceD.Get(keys.GetName())
        hist.Scale(1./494.)
        hCutFlow_ePlusLaserSig_DistanceD[keys.GetName()] = hist
        
    hCutFlow_gPlusLaserBkgPositron_DistanceD = {}
    for keys in gPlusLaserBkgPositronDistanceD.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1d' not in checkKeys: continue
        hist = gPlusLaserBkgPositronDistanceD.Get(keys.GetName())
        hist.Scale(1./30.0)
        hCutFlow_gPlusLaserBkgPositron_DistanceD[keys.GetName()] = hist
        
    hCutFlow_gPlusLaserBkgElectron_DistanceD = {}
    for keys in gPlusLaserBkgElectronDistanceD.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1d' not in checkKeys: continue
        hist = gPlusLaserBkgElectronDistanceD.Get(keys.GetName())
        hist.Scale(1./30.0)
        hCutFlow_gPlusLaserBkgElectron_DistanceD[keys.GetName()] = hist
        
        
        
        
        
        
    #### get the TrackPy cutflow first
    ### pY parameter N-1 plot
    hCutFlow_ePlusLaserBkg_TrackPy = {}
    for keys in ePlusLaserBkgTrackPy.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1Py' not in checkKeys: continue
        hist = ePlusLaserBkgTrackPy.Get(keys.GetName())
        hist.Scale(1./160.0)
        hCutFlow_ePlusLaserBkg_TrackPy[keys.GetName()] = hist
        
    hCutFlow_ePlusLaserSig_TrackPy = {}
    for keys in ePlusLaserSigTrackPy.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1Py' not in checkKeys: continue
        hist = ePlusLaserSigTrackPy.Get(keys.GetName())
        hist.Scale(1./494.)
        hCutFlow_ePlusLaserSig_TrackPy[keys.GetName()] = hist
        
    hCutFlow_gPlusLaserBkgPositron_TrackPy = {}
    for keys in gPlusLaserBkgPositronTrackPy.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1Py' not in checkKeys: continue
        hist = gPlusLaserBkgPositronTrackPy.Get(keys.GetName())
        hist.Scale(1./30.)
        hCutFlow_gPlusLaserBkgPositron_TrackPy[keys.GetName()] = hist 
        
    hCutFlow_gPlusLaserBkgElectron_TrackPy = {}
    for keys in gPlusLaserBkgElectronTrackPy.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1Py' not in checkKeys: continue
        hist = gPlusLaserBkgElectronTrackPy.Get(keys.GetName())
        hist.Scale(1./30.)
        hCutFlow_gPlusLaserBkgElectron_TrackPy[keys.GetName()] = hist
    
    
    
    
    #### get the FitParameter cutflow first
    ### SVD parameter N-1 plot
    hCutFlow_ePlusLaserBkg_FitParameter = {}
    for keys in ePlusLaserBkgFitParameter.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1FitPar' not in checkKeys: continue
        hist = ePlusLaserBkgFitParameter.Get(keys.GetName())
        hist.Scale(1./160.)
        hCutFlow_ePlusLaserBkg_FitParameter[keys.GetName()] = hist
        
    hCutFlow_ePlusLaserSig_FitParameter = {}
    for keys in ePlusLaserSigFitParameter.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1FitPar' not in checkKeys: continue
        hist = ePlusLaserSigFitParameter.Get(keys.GetName())
        hist.Scale(1./494.)
        hCutFlow_ePlusLaserSig_FitParameter[keys.GetName()] = hist
        
    hCutFlow_gPlusLaserBkgPositron_FitParameter = {}
    for keys in gPlusLaserBkgPositronFitParameter.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1FitPar' not in checkKeys: continue
        hist = gPlusLaserBkgPositronFitParameter.Get(keys.GetName())
        hist.Scale(1./30.0)
        hCutFlow_gPlusLaserBkgPositron_FitParameter[keys.GetName()] = hist
        
    hCutFlow_gPlusLaserBkgElectron_FitParameter = {}
    for keys in gPlusLaserBkgElectronFitParameter.GetListOfKeys():
        checkKeys = str(keys)
        if 'hNminus1FitPar' not in checkKeys: continue
        hist = gPlusLaserBkgElectronFitParameter.Get(keys.GetName())
        hist.Scale(1./30.0)
        hCutFlow_gPlusLaserBkgElectron_FitParameter[keys.GetName()] = hist
    
    
    
    drawline    = False
    logy        = True
    leftLegend  = False
    doAtlas     = False
    doLumi      = False
    noRatio     = False
    do80        = False
    do59        = False
    drawPattern = ""
    logz        = False
    logx        = False
    drawline    = True
    energyPlots = True
    
    
    
    xAxisName   = "d [m]"
    yAxisName   = "Entries/BX"
    
    yrange1down = 0.001
    
    
    
    
    
    
    
    PlotColor   = [2, 4, kGreen+3, kAzure+10]
    LegendName  = ["e+laser (Sig)", "e+laser (Bkg)", "g+laser (Bkg) e^{+} side", "g+laser (Bkg) e^{-} side"]
    
    
    hNameList = ['hNminus1d_High_Multiplicity_ELt4_TrackInclusive', 'hNminus1d_High_Multiplicity_EGt4_TrackInclusive', 'hNminus1d_Medium_Multiplicity_ELt4_TightTrack', 'hNminus1d_Medium_Multiplicity_EGt4_TightTrack', 'hNminus1d_Medium_Multiplicity_ELt4_LooseTrack', 'hNminus1d_Medium_Multiplicity_EGt4_LooseTrack', 'hNminus1d_Low_Multiplicity_ELt4_TightTrack', 'hNminus1d_Low_Multiplicity_EGt4_TightTrack', 'hNminus1d_Low_Multiplicity_ELt4_LooseTrack', 'hNminus1d_Low_Multiplicity_EGt4_LooseTrack']
    
    ### for arrow and lines
    for hName in hNameList:
        if 'Medium_Multiplicity_EGt4_TightTrack' in hName:
            yrange1up = 1e3
            yline1up  = 1e1
        else:
            yrange1up   = 1e2
            yline1up    = 1e0
            
        latexName, latexName2, latexName3, latexName4, latexName5   = getLatexNames(hName)
        
        yline1low   = 0.0
        xline1low   = getCutValue(hName)['d']
        xline1high  = getCutValue(hName)['d']
        xrange1down = getCutValue(hName)['xlow']
        xrange1up   = getCutValue(hName)['xhigh']
        
        FirstTH1    = [hCutFlow_ePlusLaserSig_DistanceD[hName], hCutFlow_ePlusLaserBkg_DistanceD[hName], hCutFlow_gPlusLaserBkgPositron_DistanceD[hName], hCutFlow_gPlusLaserBkgElectron_DistanceD[hName]]
        
        DrawHistsCutValueOneArrow(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/"+hName, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4, xline1low, xline1high)
        
        
    xAxisName   = "p'_{Y} [GeV]"
    yAxisName   = "Entries/BX"
    hNameList = ['hNminus1Py_High_Medium_Multiplicity_EInclusive_TrackInclusive', 'hNminus1Py_Low_Multiplicity_ELt4_TightTracks', 'hNminus1Py_Low_Multiplicity_ELt4_LooseTracks', 'hNminus1Py_Low_Multiplicity_EGt4_TightTracks', 'hNminus1Py_Low_Multiplicity_EGt4_LooseTracks']
    
    ### for arrow and lines
    for hName in hNameList:
        latexName, latexName2, latexName3, latexName4, latexName5   = getLatexNames(hName)
        if "High_Medium" in hName:
            yrange1up = 1e5
        else:
            yrange1up   = 1e3
        yline1low   = 0.0
        yline1up    = 1e0
        xline1low   = getCutValue(hName)['d']
        xline1high  = getCutValue(hName)['d']
        negativexline1high  = getCutValue(hName)['dmin']
        negativexline1low   = getCutValue(hName)['dmin']
        xrange1down = getCutValue(hName)['xlow']
        xrange1up   = getCutValue(hName)['xhigh']
        
        FirstTH1    = [hCutFlow_ePlusLaserSig_TrackPy[hName], hCutFlow_ePlusLaserBkg_TrackPy[hName], hCutFlow_gPlusLaserBkgPositron_TrackPy[hName], hCutFlow_gPlusLaserBkgElectron_TrackPy[hName]]
        
        DrawHistsCutValueTwoArrows(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/"+hName, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4, xline1low, xline1high, negativexline1low, negativexline1high)
        
        
        
    xAxisName   = "Second SVD Fit Parameter"
    yAxisName   = "Entries/BX"
    hNameList = ['hNminus1FitPar1_High_Multiplicity_ELt4_TrackInclusive', 'hNminus1FitPar1_High_Multiplicity_EGt4_TrackInclusive', 'hNminus1FitPar1_Medium_Multiplicity_ELt4_TrackInclusive', 'hNminus1FitPar1_Medium_Multiplicity_EGt4_TrackInclusive', 'hNminus1FitPar1_Low_Multiplicity_ELt4_TightTrack', 'hNminus1FitPar1_Low_Multiplicity_EGt4_TightTrack', 'hNminus1FitPar1_Low_Multiplicity_ELt4_LooseTrack', 'hNminus1FitPar1_Low_Multiplicity_EGt4_LooseTrack']
    
    ### for arrow and lines
    for hName in hNameList:
        latexName, latexName2, latexName3, latexName4, latexName5   = getLatexNames(hName)
        yrange1up   = 1e3
        yrange1down = 0.0001
        #print latexName, " ", latexName2, " ", latexName3, " ", latexName4, " ", latexName5
        yline1low   = 0.0
        yline1up    = 1e0
        xline1low   = getCutValue(hName)['d']
        xline1high  = getCutValue(hName)['d']
        xrange1down = getCutValue(hName)['xlow']
        xrange1up   = getCutValue(hName)['xhigh']
        
        FirstTH1    = [hCutFlow_ePlusLaserSig_FitParameter[hName], hCutFlow_ePlusLaserBkg_FitParameter[hName], hCutFlow_gPlusLaserBkgPositron_FitParameter[hName], hCutFlow_gPlusLaserBkgElectron_FitParameter[hName]]
        
        DrawHistsCutValueOneArrow(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/"+hName, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4, xline1low, xline1high)
        
        
        
    xAxisName   = "Third SVD Fit Parameter"
    yAxisName   = "Entries/BX"
    hNameList = ['hNminus1FitPar2_High_Multiplicity_ELt4_TrackInclusive', 'hNminus1FitPar2_High_Multiplicity_EGt4_TrackInclusive', 'hNminus1FitPar2_Medium_Multiplicity_ELt4_TrackInclusive', 'hNminus1FitPar2_Medium_Multiplicity_EGt4_TrackInclusive', 'hNminus1FitPar2_Low_Multiplicity_ELt4_TightTrack', 'hNminus1FitPar2_Low_Multiplicity_EGt4_TightTrack', 'hNminus1FitPar2_Low_Multiplicity_ELt4_LooseTrack', 'hNminus1FitPar2_Low_Multiplicity_EGt4_LooseTrack']
    
    ### for arrow and lines
    for hName in hNameList:
        latexName, latexName2, latexName3, latexName4, latexName5   = getLatexNames(hName)
        yrange1up   = 1e3
        yrange1down = 0.0001
        yline1low   = 0.0
        if 'High_Multiplicity_EGt4_TrackInclusive' in hName:
            yline1up    = 1e1
        else:
            yline1up    = 1e0
        xline1low   = getCutValue(hName)['d']
        xline1high  = getCutValue(hName)['d']
        xrange1down = getCutValue(hName)['xlow']
        xrange1up   = getCutValue(hName)['xhigh']
        
        FirstTH1    = [hCutFlow_ePlusLaserSig_FitParameter[hName], hCutFlow_ePlusLaserBkg_FitParameter[hName], hCutFlow_gPlusLaserBkgPositron_FitParameter[hName], hCutFlow_gPlusLaserBkgElectron_FitParameter[hName]]
        
        DrawHistsCutValueOneArrow(FirstTH1, LegendName, PlotColor, xAxisName, yAxisName, xrange1down, xrange1up, yrange1down, yrange1up, directory+"/"+hName, yline1low, yline1up, drawline, logy, latexName, latexName2, latexName3, leftLegend, doAtlas, doLumi, noRatio, do80, do59, drawPattern, logz, logx, latexName4, xline1low, xline1high)
    
    
    
    
    

if __name__=="__main__":
    start = time.time()
    main()
    print("The time taken: ", time.time() - start, " s")
