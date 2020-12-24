#### this code mixes the background tracks and signal tracks per bunch crossing
#### run: python3 getSigBkgTextFile.py <bxNumberWanted>
import os
import sys
import time
import pprint
import math
from ROOT import *
from collections import OrderedDict
import argparse


def main():
    numberTxtFile = int(sys.argv[1])
    bkgFileName   = open("EBeamOnlyWIS_DividedByBX"+str(numberTxtFile)+"_trackInfoClean.txt")
    #sigFileName   = open("list_root_hics_165gev_w0_5000nm_provisional_10xi3_6cee466a_trackInfoClean.txt")
    #sigFileName   = open("list_root_hics_165gev_w0_8000nm_provisional_10xi3_6cee466a_trackInfoClean.txt")
    sigFileName   = open("list_root_hics_165gev_w0_3000nm_WIS_trackInfoClean.txt")
    
    #outFile       = open("BkgEBeam_SignalHics5000nmProvisional_BX"+str(numberTxtFile)+"_trackInfoClean.txt","w")
    #outFile       = open("BkgEBeam_SignalHics8000nmProvisional_BX"+str(numberTxtFile)+"_trackInfoClean.txt","w")
    outFile       = open("BkgEBeam_SignalHics3000nmOldReverse_BX"+str(numberTxtFile)+"_trackInfoClean.txt","w")
    
    ### write the bkg as it is
    for lines in bkgFileName.readlines():
        outFile.write(lines)
    print("Printed bkg to mix")
    
    ### work on the signal, first collect all the signal lines in a dictionary
    position = []
    particleInEachBX = {}
    for lines2 in sigFileName.readlines():
        lines2   = lines2.rstrip()
        if "#" in lines2: continue
        eachWord = lines2.split()
        
        bxNumber  = int(eachWord[0])
        trackId   = int(eachWord[2])
        pdgId     = int(eachWord[1])
        staveId   = int(eachWord[3])
        
        position.append([bxNumber, pdgId, trackId, staveId])
        particleInEachBX.setdefault(bxNumber, []).append(lines2) ### this ensures that the new line for the same bxnumber does not erase the previous line, the new line will be appended as a list 
        
        
    ### count the number of signal positrons in each BX and then put BX: number of signals in a dictionary
    print("Now counting signals in each BX")
    outRootFile = TFile("signalTrackMultiplicity_w0_5000nm_provisional_BX"+str(numberTxtFile)+".root", "RECREATE")
    outRootFile.cd()
    signalTrackMultiplicity = TH1F("signalTrackMultiplicity", "track multiplicity; number of tracks; BX",80, 0, 400)
    signalNumber     = {}
    for bx in range(0,494):
        counter = 0
        print("BX: ", bx)
        for tracks in position:
            if tracks[0] == bx:
                if(tracks[1] == -11 and tracks[2] == 1 and (tracks[3]==1000 or tracks[3]==1001)):
                    counter += 1
        signalNumber[bx] = counter ### store the number of signal for each BX
                
    ### sort the list according to the number of signals, descending order
    #sortedSignalNumber = {k: v for k, v in sorted(signalNumber.items(), key=lambda item: item[1], reverse=True)}
    ### this is with ascending order
    sortedSignalNumber = {k: v for k, v in sorted(signalNumber.items(), key=lambda item: item[1])}
    
    print(sortedSignalNumber)
    ### write one bx of signal to the output file
    writeNumber = 0
    for key in sortedSignalNumber:
        print(key,sortedSignalNumber[key])
        writeNumber += 1
        signalTrackMultiplicity.Fill(sortedSignalNumber[key])
        if(writeNumber!=numberTxtFile): continue
        print("writing only the ",numberTxtFile," th signal BX")
        for eachLine in particleInEachBX[key]:
            findPDG     = int(eachLine.split()[1])
            findTrackId = int(eachLine.split()[2])
            if findPDG == -11 and findTrackId == 1:
                outFile.write(eachLine+"\n")
    
    outFile.close()
    outRootFile.Write()
    outRootFile.Close()
                
            
    
if __name__=="__main__":
    start = time.time()
    main()
    print("--- The time taken: ", time.time() - start, " s")
