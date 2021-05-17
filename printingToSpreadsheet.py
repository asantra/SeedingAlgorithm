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
    
    parser = argparse.ArgumentParser(description='Code to find seed tracks')
    parser.add_argument('-p', action="store", dest="Particle", type=str, default="Positron")
    parser.add_argument('-process', action="store", dest="Process", type=str, default="hics")
    args = parser.parse_args()
    
    particle     = args.Particle
    process      = args.Process
    
    if(process == "hics"):
        ### low signal multiplicity tracks
        #listOfTracks = [1, 5, 10, 20, 30, 50, 80, 100, 130, 150, 170, 185, 200, 220]
        #outFileName  = "hybridSVDCutsNumbersUpdated.csv"
        ### high signal multiplicity tracks
        listOfTracks = [250, 300, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500]
        outFileName  = "hybridSVDCutsNumbersUpdatedELaserHighMultiplicity.csv"
    else:
        listOfTracks = [1, 4, 7, 10, 13, 16]
        outFileName  = "hybridSVDCutsNumbersgPlusLaser"+particle+".csv"
        
    bxList       = [1,2,3,4]
    inDirectory  = "." #"LooseSVDanddCuts_December29_2020"
    
    outCSVFile   = open(outFileName, "w")
    suffixToRootFile = "WithFit3or4HitsTracksAndDistanceCut"
    
    trackList = {}
    for tracks in listOfTracks:
        bxValueList = {}
        seedMultiplicityPrelimSignal = 0
        seedMultiplicityPrelimBkg    = 0
        seedMultiplicityPrelimCombi  = 0
        for bx in bxList:
            if(process == "hics"):
                ##### this is for electron beam plus laser case
                ##rootFileName = "seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
                
                ##rootFileNameOnlySignal = "seedingInformation_BkgEBeam_SignalHics3000nmOldForSignalMultiplicityLessThan20Or5000nmForSignalMultiplicityMoreThan20_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+suffixToRootFile+".root"
                
                ##rootFileOnlyBackground = "seedingInformation_EBeamOnlyWIS_DividedByBX"+str(bx)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
                
                
                rootFileName = "seedingInformationFiles/seedingInformation_BkgEBeam_SignalPositronhics3000nm_jeti40_122020_9550dac4_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
                
                rootFileNameOnlySignal = "seedingInformationFiles/seedingInformation_BkgEBeam_SignalPositronhics3000nm_jeti40_122020_9550dac4_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_OnlySignal_PositronSide_"+suffixToRootFile+".root"
                
                rootFileOnlyBackground = "seedingInformationFiles/seedingInformation_ePlusLaserBkgKaptonWindowNewSamplesMarch62021_DividedByBX"+str(bx)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_PositronSide_"+suffixToRootFile+".root"
                
            else:
                #### this is for g+laser beam case
                
                rootFileName = "seedingInformation_BkgGBeam_Signal"+particle+"bppp3000nmOr5000nm_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_"+particle+"Side_"+suffixToRootFile+".root"
                
                rootFileNameOnlySignal = "seedingInformation_BkgGBeam_Signal"+particle+"bppp3000nmOr5000nm_BX"+str(bx)+"_SignalTracks"+str(tracks)+"_trackInfoClean_VariableEnergyCut_OnlySignal_"+particle+"Side_"+suffixToRootFile+".root"
                
                rootFileOnlyBackground = "seedingInformation_gPlusLaserBkgNewSamplesJan262021_DividedByBX"+str(bx)+"_trackInfoClean_VariableEnergyCut_SignalAndBackground_"+particle+"Side_"+suffixToRootFile+".root"
            
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
            for nxBins in range(hSignalMultiplicity.GetNbinsX()+1):
                if(hSignalMultiplicity.GetBinContent(nxBins)!=0):
                    trueSignal = (nxBins - 1)
                    break
              
            #### find the seeds before the fit
            ### signal case
            prelimSigMultiSignal = 0
            for nxBins in range(hSeedMultiplicityPrelimSignal.GetNbinsX()+1):
                if(hSeedMultiplicityPrelimSignal.GetBinContent(nxBins)!=0):
                    prelimSigMultiSignal = (nxBins - 1)
                    break
            
            seedMultiplicityPrelimSignal += prelimSigMultiSignal
            ### background case
            prelimSigMultiBkg = 0
            for nxBins in range(hSeedMultiplicityPrelimBkg.GetNbinsX()+1):
                if(hSeedMultiplicityPrelimBkg.GetBinContent(nxBins)!=0):
                    prelimSigMultiBkg = (nxBins - 1)
                    break
                
            seedMultiplicityPrelimBkg += prelimSigMultiBkg
            ### signal and background case
            prelimSigMultiCombi = 0
            for nxBins in range(hSeedMultiplicityPrelim.GetNbinsX()+1):
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
        
        #pprint.pprint(trackList)
           
    

    firstTrack = listOfTracks[0]
    outCSVFile.write(",,,,"+trackList[firstTrack][1]["bkgOnlyInclusive"]+","+trackList[firstTrack][1]["bkgOnlyTight"]+","+trackList[firstTrack][1]["bkgOnlyLoose"]+",,,,"+trackList[firstTrack][2]["bkgOnlyInclusive"]+","+trackList[firstTrack][2]["bkgOnlyTight"]+","+trackList[firstTrack][2]["bkgOnlyLoose"]+",,,,"+trackList[firstTrack][3]["bkgOnlyInclusive"]+","+trackList[firstTrack][3]["bkgOnlyTight"]+","+trackList[firstTrack][3]["bkgOnlyLoose"]+",,,,"+trackList[firstTrack][4]["bkgOnlyInclusive"]+","+trackList[firstTrack][4]["bkgOnlyTight"]+","+trackList[firstTrack][4]["bkgOnlyLoose"]+"\n")
    
    for tracks in listOfTracks:
        #print trackList[tracks][1][0], " ", trackList[tracks][1][1], " ", trackList[tracks][1][2], " ", trackList[tracks][1][3]
        outCSVFile.write(trackList[tracks][1]["trueSignal"]+","+trackList[tracks][1]["signalOnlyInclusive"]+","+trackList[tracks][1]["signalOnlyTight"]+","+trackList[tracks][1]["signalOnlyLoose"]+","+trackList[tracks][1]["combinedInclusive"]+","+trackList[tracks][1]["combinedTight"]+","+trackList[tracks][1]["combinedLoose"]+","+trackList[tracks][2]["trueSignal"]+","+trackList[tracks][2]["signalOnlyInclusive"]+","+trackList[tracks][2]["signalOnlyTight"]+","+trackList[tracks][2]["signalOnlyLoose"]+","+trackList[tracks][2]["combinedInclusive"]+","+trackList[tracks][2]["combinedTight"]+","+trackList[tracks][2]["combinedLoose"]+","+trackList[tracks][3]["trueSignal"]+","+trackList[tracks][3]["signalOnlyInclusive"]+","+trackList[tracks][3]["signalOnlyTight"]+","+trackList[tracks][3]["signalOnlyLoose"]+","+trackList[tracks][3]["combinedInclusive"]+","+trackList[tracks][3]["combinedTight"]+","+trackList[tracks][3]["combinedLoose"]+","+trackList[tracks][4]["trueSignal"]+","+trackList[tracks][4]["signalOnlyInclusive"]+","+trackList[tracks][4]["signalOnlyTight"]+","+trackList[tracks][4]["signalOnlyLoose"]+","+trackList[tracks][4]["combinedInclusive"]+","+trackList[tracks][4]["combinedTight"]+","+trackList[tracks][4]["combinedLoose"]+",,,"+trackList[tracks]["prelimSeedBkg"]+",,,"+trackList[tracks]["prelimSeedCombined"]+",,,"+trackList[tracks]["prelimSeedSignal"]+"\n")
    
    
    outCSVFile.close()
    
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
