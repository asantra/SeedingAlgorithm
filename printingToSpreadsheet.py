### code to divide the background in separate BXs
### run: python getBXForEBeamOnly.py <list of background file names>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse


def main():
    listOfTracks = [1, 5, 10, 20, 30, 50, 80, 100, 130, 150, 170, 185, 200, 220]
    bxList       = [1,2,3,4]
    inDirectory  = "." #"LooseSVDanddCuts_December29_2020"
    
    outCSVFile   = open("hybridSVDCutsNumbersUpdated.csv", "w")
    suffixToRootFile = "WithFit3or4HitsTracksAndDistanceCut"
    
    trackList = {}
    for tracks in listOfTracks:
        bxValueList = {}
        seedMultiplicityPrelimSignal = 0
        seedMultiplicityPrelimBkg    = 0
        seedMultiplicityPrelimCombi  = 0
        for bx in bxList:
            rootFileName = "seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
            
            rootFileNameOnlySignal = "seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+suffixToRootFile+".root"
            
            rootFileOnlyBackground = "seedingInformation_EBeamOnlyWIS_DividedByBX"+str(bx)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
            
            inRootFile                       = TFile(inDirectory+"/"+rootFileName)
            inRootFileSignal                 = TFile(inDirectory+"/"+rootFileNameOnlySignal)
            inRootFileBkg                    = TFile(inDirectory+"/"+rootFileOnlyBackground)
            
            ### prepare the signal only scenario
            hSignalMultiplicity              = inRootFileSignal.Get("hSignalMultiplicity")
            hSeedEnergySignal                = inRootFileSignal.Get("hSeedEnergy")
            hSeedEnergyLooseSignal           = inRootFileSignal.Get("hSeedEnergyLoose")
            hSeedEnergyTightSignal           = inRootFileSignal.Get("hSeedEnergyTight")
            hSeedMultiplicityPrelimSignal    = inRootFileSignal.Get("hSeedMultiplicityPrelim")
            
            hSeedEnergyBkg                   = inRootFileBkg.Get("hSeedEnergy")
            hSeedEnergyLooseBkg              = inRootFileBkg.Get("hSeedEnergyLoose")
            hSeedEnergyTightBkg              = inRootFileBkg.Get("hSeedEnergyTight")
            hSeedMultiplicityPrelimBkg       = inRootFileBkg.Get("hSeedMultiplicityPrelim")
            
            
            hSeedEnergy                      = inRootFile.Get("hSeedEnergy")
            hSeedEnergyLoose                 = inRootFile.Get("hSeedEnergyLoose")
            hSeedEnergyTight                 = inRootFile.Get("hSeedEnergyTight")
            hSeedMultiplicityPrelim          = inRootFile.Get("hSeedMultiplicityPrelim")
            
            
            ### irrespective of signal+background combination
            trueSignal = 0
            for nxBins in xrange(0, hSignalMultiplicity.GetNbinsX()+1):
                if(hSignalMultiplicity.GetBinContent(nxBins)!=0):
                    trueSignal = (nxBins - 1)
                    break
              
            #### find the seeds before the fit
            ### signal case
            prelimSigMultiSignal = 0
            for nxBins in xrange(0, hSeedMultiplicityPrelimSignal.GetNbinsX()+1):
                if(hSeedMultiplicityPrelimSignal.GetBinContent(nxBins)!=0):
                    prelimSigMultiSignal = (nxBins - 1)
                    break
            
            seedMultiplicityPrelimSignal += prelimSigMultiSignal
            ### background case
            prelimSigMultiBkg = 0
            for nxBins in xrange(0, hSeedMultiplicityPrelimBkg.GetNbinsX()+1):
                if(hSeedMultiplicityPrelimBkg.GetBinContent(nxBins)!=0):
                    prelimSigMultiBkg = (nxBins - 1)
                    break
                
            seedMultiplicityPrelimBkg += prelimSigMultiBkg
            ### signal and background case
            prelimSigMultiCombi = 0
            for nxBins in xrange(0, hSeedMultiplicityPrelim.GetNbinsX()+1):
                if(hSeedMultiplicityPrelim.GetBinContent(nxBins)!=0):
                    prelimSigMultiCombi = (nxBins - 1)
                    break
                
            seedMultiplicityPrelimCombi += prelimSigMultiCombi
                
            bxValueList[bx] = {"trueSignal": str(trueSignal), "signalOnlyInclusive":  str(hSeedEnergySignal.Integral()), "signalOnlyTight": str(hSeedEnergyTightSignal.Integral()), "signalOnlyLoose": str(hSeedEnergyLooseSignal.Integral()), "bkgOnlyInclusive": str(hSeedEnergyBkg.Integral()), "bkgOnlyTight": str(hSeedEnergyTightBkg.Integral()), "bkgOnlyLoose": str(hSeedEnergyLooseBkg.Integral()), "combinedInclusive": str(hSeedEnergy.Integral()), "combinedTight": str(hSeedEnergyTight.Integral()), "combinedLoose": str(hSeedEnergyLoose.Integral())}
            
        avgPrelimSeedMultiSignal   = seedMultiplicityPrelimSignal/4.0
        avgPrelimSeedMultiBkg      = seedMultiplicityPrelimBkg/4.0
        avgPrelimSeedMultiCombined = seedMultiplicityPrelimCombi/4.0
        
        ### storing the values in the ditionary
        trackList[tracks] = bxValueList
        trackList[tracks].update({"prelimSeedSignal": str(avgPrelimSeedMultiSignal), "prelimSeedBkg": str(avgPrelimSeedMultiBkg), "prelimSeedCombined": str(avgPrelimSeedMultiCombined)})
           
    
    ### first print the background
    outCSVFile.write(",,,,"+trackList[1][1]["bkgOnlyInclusive"]+","+trackList[1][1]["bkgOnlyTight"]+","+trackList[1][1]["bkgOnlyLoose"]+",,,,"+trackList[1][2]["bkgOnlyInclusive"]+","+trackList[1][2]["bkgOnlyTight"]+","+trackList[1][2]["bkgOnlyLoose"]+",,,,"+trackList[1][3]["bkgOnlyInclusive"]+","+trackList[1][3]["bkgOnlyTight"]+","+trackList[1][3]["bkgOnlyLoose"]+",,,,"+trackList[1][4]["bkgOnlyInclusive"]+","+trackList[1][4]["bkgOnlyTight"]+","+trackList[1][4]["bkgOnlyLoose"]+"\n")
    
    for tracks in listOfTracks:
        #print trackList[tracks][1][0], " ", trackList[tracks][1][1], " ", trackList[tracks][1][2], " ", trackList[tracks][1][3]
        outCSVFile.write(trackList[tracks][1]["trueSignal"]+","+trackList[tracks][1]["signalOnlyInclusive"]+","+trackList[tracks][1]["signalOnlyTight"]+","+trackList[tracks][1]["signalOnlyLoose"]+","+trackList[tracks][1]["combinedInclusive"]+","+trackList[tracks][1]["combinedTight"]+","+trackList[tracks][1]["combinedLoose"]+","+trackList[tracks][2]["trueSignal"]+","+trackList[tracks][2]["signalOnlyInclusive"]+","+trackList[tracks][2]["signalOnlyTight"]+","+trackList[tracks][2]["signalOnlyLoose"]+","+trackList[tracks][2]["combinedInclusive"]+","+trackList[tracks][2]["combinedTight"]+","+trackList[tracks][2]["combinedLoose"]+","+trackList[tracks][3]["trueSignal"]+","+trackList[tracks][3]["signalOnlyInclusive"]+","+trackList[tracks][3]["signalOnlyTight"]+","+trackList[tracks][3]["signalOnlyLoose"]+","+trackList[tracks][3]["combinedInclusive"]+","+trackList[tracks][3]["combinedTight"]+","+trackList[tracks][3]["combinedLoose"]+","+trackList[tracks][4]["trueSignal"]+","+trackList[tracks][4]["signalOnlyInclusive"]+","+trackList[tracks][4]["signalOnlyTight"]+","+trackList[tracks][4]["signalOnlyLoose"]+","+trackList[tracks][4]["combinedInclusive"]+","+trackList[tracks][4]["combinedTight"]+","+trackList[tracks][4]["combinedLoose"]+",,,"+trackList[tracks]["prelimSeedBkg"]+",,,"+trackList[tracks]["prelimSeedCombined"]+",,,"+trackList[tracks]["prelimSeedSignal"]+"\n")
    
    
    outCSVFile.close()
    
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
